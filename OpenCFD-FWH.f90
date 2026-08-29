!-----------------------------------------------------------
! OpenCFD-FWH (MPI+OpenMP)
! Copyright by ZhangKeli
! Reference: https://arxiv.org/abs/2312.16263
! Permeable surface FWH for wind tunnel cases using Garrick Triangle
! Openmp paralle support
! Version 2: 边Read FWH data 边计算Qn Lr Lm以降低内存需求
! Version 3: 以Observer为单位进行循环计算，降低内存压力
! Version mpi: 支持mpi+openmp混合并行加速（同时避免在超算上单个节点爆内存
! 2023-12-22: a bug remove(Mi=-U0i/C0); support formatted FWH data file
!------------------------------------------------------------
  module Const_Variables
  implicit none
  include "mpif.h"
  ! integer,parameter:: PRE_EC=4           ! Single precision
  integer,parameter:: PRE_EC=8           ! Double Precision
  ! integer,parameter:: FWH_DATA_TYPE=MPI_REAL             ! Single precision
  integer,parameter:: FWH_DATA_TYPE=MPI_DOUBLE_PRECISION ! Double precision
  real(PRE_EC),parameter:: PI=3.14159265358979d0
end module Const_Variables

  module Global_Variables
  use Const_Variables
  implicit none
 ! global parameter (for all Faces) 
  integer,save:: Kstep_start,Kstep_end,NUM_Obs,delta_step,FWH_data_Format
  real(PRE_EC),save:: Ma,AoA,delta_t
  integer,save:: NUM_Face,NUM_Face_ALL                  ! 本mpi进程的网格面数，总网格面数  
  integer,save:: NUM_THREADS,NUM_Frame,Total_Faces      ! OpenMP线程数，总实例数，总子面数
  real(PRE_EC),pointer,dimension(:):: Obs_x,Obs_y,Obs_z
  real(PRE_EC),pointer,dimension(:):: t_interp,pp

!-----------------------------------------------------------------------------------------

!  核心变量——每面网格存储的信息 （全局变量）
   TYPE Face_TYPE           !  variables for each Face 
     integer :: nx,ny,nz
     real(PRE_EC),pointer,dimension(:,:,:):: x,y,z,dS,n1,n2,n3  ! coordinates of cell center, 网格中心坐标 
     real(PRE_EC),pointer,dimension(:,:,:):: R,r1,r2,r3         ! R & r1,r2,r3 
	 real(PRE_EC),pointer,dimension(:,:,:,:):: d,u,v,w,p,Lm,Qn,Lr
   real(PRE_EC),pointer,dimension(:,:,:,:):: pp_ret,t_obs,pp_interp
   End TYPE Face_TYPE  
   TYPE (Face_TYPE), save,dimension(:),allocatable,target:: Face
!-----------mpi data ----------------------------------------------------------- 
   integer:: my_id,NUM_proc                   ! my_id (本进程号), NUM_proc 总进程数
   integer,pointer,dimension(:):: F_Proc, F_n    ! F_proc(m) m面所在的进程号; F_n(m) m面所在进程中的内部编号
   integer,pointer,dimension(:):: my_Faces       ! 本进程包含的面号
  end module Global_Variables
!------------------------------------------------------------------


!-------------------------------------------
  program compute_farfield_noise
     use Global_Variables
     implicit none
     integer:: ierr
     integer :: obs,frame
     real(8) :: t_start,t_end
     character(len=50):: filename
!---------------init mpi----------------
     call mpi_init(ierr)                                     ! 初始化MPI
     call mpi_comm_rank(MPI_COMM_WORLD,my_id,ierr)           ! 获取本进程编号
     call mpi_comm_size(MPI_COMM_WORLD,NUM_proc,ierr)   
!---------------------------------------
     if(my_id .eq. 0) t_start = MPI_Wtime()
     call init  

     if(my_id .eq. 0) then
      t_end = MPI_Wtime()
      print*, "Cpu wall time for initialization: ",t_end-t_start
      call system('mkdir FWH_result-mpi')
     endif

      do obs = 1,NUM_Obs
      
        if(my_id .eq. 0) print*, "Processing Observer :",obs

        call compute_R(obs)
        if(my_id .eq. 0) print*, "Compute R OK "


        call compute_Noise
        if(my_id .eq. 0) print*, "Compute pp_Ret OK "


        call compute_ObserverTime
        if(my_id .eq. 0) print*, "Compute Observer Time OK "


        call Interp_PressureSignal
        if(my_id .eq. 0) print*, "Interp Pressure Signal OK "


        call Integrate_Sources
        if(my_id .eq. 0) then
          print*, "Integrate Sources OK, write pressure signals ... "

          write(filename,"('./FWH_result-mpi/p_observer-'I3.3'.log')") obs
          open(222,file=trim(filename))
          do frame =1,NUM_Frame
          write(222,"(2E30.16)") t_interp(frame),pp(frame)
          enddo
          close(222)
          print *,"Done writing",filename
          print *
        endif

      enddo

    if(my_id .eq. 0) then
      t_end = MPI_Wtime()
      print*, "Cpu wall time: ",t_end-t_start
    endif
    
    call MPI_Finalize(ierr)
  end

!=====================================================================================
!------------------------------------------------------------------------------     
! Read the message of the mesh and the initial flow;
! mpi partition，读取网格，初始流场信息; 
!------------------------------------------------------------------------------
   subroutine init
   use  Global_Variables
   implicit none
   integer,allocatable,dimension(:,:):: F_grid    ! number of cells and Face id
   integer,dimension(:),pointer:: Fproc           ! proc id
   integer:: m,G0,mg,n,t1,t2,mf,m0,mp
   integer:: i,j,k,nx,ny,nz,m1,m2
   real(PRE_EC):: grid_av
   integer,allocatable,dimension(:):: NI,NJ,NK,Pgrid ! nx,ny,nz, number of cells in proc
   Type (Face_TYPE),pointer:: F
   real(PRE_EC),dimension(:,:,:,:),pointer:: U
   integer:: ierr, status(MPI_status_size)
   character(len=50):: filename

!------------------------------------------------------------------
  call read_parameter_FWH
!$ call omp_set_num_threads(NUM_THREADS)   ! 设置OpenMP的运行线程数 （并行数目）， 本语句对openmp编译器不是注释!
!测试一下运行的进程 （openmp编译时，不是注释）
!$ if(my_id ==0) then 
!$OMP Parallel
!$  print*, "omp run ..."
!$OMP END parallel
!$ endif 
  if(my_id .eq. 0) then
    print*,  "-------------------------- OpenCFD-FWH (MPI-OpenMP version) ------------------------------"
    print*,  "                      Copyright by Zhang Keli, zkl70@buaa.edu.cn                          "
    print*,  "                         Programming by Zhang Keli  2023-12                               "
    print*,  " Ma,AoA,delta_t=", Ma,AoA,delta_t
    print*,  " NUM_Obs,NUM_Frame,NUM_THREADS=",NUM_Obs,NUM_Frame,NUM_THREADS
    print*,  "------------------------------------------------------------------------------------------" 
    print*
  endif
!---------------------------- mpi partition --------------------------
  if(my_id .eq. 0) then
    print*, "partition ......"
    open(100,file="FWH_Surface_Geo.dat")
    read(100,*)
    read(100,*) NUM_Face_ALL         ! 总面数
    print*, "NUM_Face=", NUM_Face_ALL
    allocate(NI(NUM_Face_ALL),NJ(NUM_Face_ALL),NK(NUM_Face_ALL) )   ! 每面的大小
    do m=1,NUM_Face_ALL      
      read(100,*) NI(m),NJ(m),NK(m) 
        do k=1,NK(m)
        do j=1,NJ(m)
        do i=1,NI(m)
           read(100,*) 
        enddo
        enddo
        enddo    
    enddo
    close(100)
!    读入每面网格数，按从多到少次序排列
    allocate(F_grid(NUM_Face_ALL,2))
    allocate(Pgrid(NUM_proc),Fproc(NUM_Face_ALL))
    Total_Faces=0
    do m=1,NUM_Face_ALL
      F_grid(m,1)=NI(m)*NJ(m)*NK(m) 
      F_grid(m,2)=m
      Total_Faces=Total_Faces + F_grid(m,1)
    enddo
    print*, "Total_Faces =",Total_Faces
!    按网格点从多到少的次序排序
    do m=1,NUM_Face_ALL
      G0=F_grid(m,1)     ! cell points
       mg=m
       do n=m+1,NUM_Face_ALL   ! Find maximum
          if(F_grid(n,1) .gt. G0 ) then
          G0=F_grid(n,1)
          mg=n
          endif
       enddo
 !    mg面与m面交换
      t1=F_grid(mg,1)
      t2=F_grid(mg,2)
      F_grid(mg,1)=F_grid(m,1)
      F_grid(mg,2)=F_grid(m,2)
      F_grid(m,1)=t1
      F_grid(m,2)=t2
    enddo

!----- partition algorithm ------
    Pgrid(:)=0
    do m=1,NUM_Face_ALL
     mf=F_grid(m,2) 
     !寻找网格数目最小的进程
     mg=Pgrid(1)
     m0=1
     do mp=1,NUM_Proc
      if(Pgrid(mp) .lt. mg) then
        mg=Pgrid(mp)
        m0=mp
      endif
     enddo
     !将该面放入网格数目最小的进程中
     Pgrid(m0)=Pgrid(m0)+F_grid(m,1)
     Fproc(mf)=m0-1
    enddo
!--------------------------------
    grid_av=1.d0*Total_Faces/NUM_Proc
    print*, "Max faces number is ", maxval(Pgrid),  " ... rato to mean faces is", maxval(Pgrid)/grid_av
    print*, "Min faces number is ", minval(Pgrid),  " ... rato to mean faces is", minval(Pgrid)/grid_av
    print*, "rato max to min is ", 1.0*maxval(Pgrid)/minval(Pgrid)  
  endif
  
  call MPI_Bcast(NUM_Face_ALL,1,MPI_Integer,0, MPI_COMM_WORLD,ierr)
  if(my_id .ne. 0) allocate(NI(NUM_Face_ALL),NJ(NUM_Face_ALL),NK(NUM_Face_ALL) )   ! 每面的大小
  call MPI_Bcast(NI,NUM_Face_ALL,MPI_Integer,0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(NJ,NUM_Face_ALL,MPI_Integer,0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(NK,NUM_Face_ALL,MPI_Integer,0, MPI_COMM_WORLD,ierr)

  allocate (F_Proc(NUM_Face_ALL),F_n(NUM_Face_ALL))   

  if(my_id .eq. 0) then
    F_Proc(1:NUM_Face_ALL)=Fproc(1:NUM_Face_ALL)
  endif
  call MPI_Bcast(F_Proc(1),NUM_Face_ALL,MPI_Integer,0,  MPI_COMM_WORLD,ierr)

! 计算F_n(m), 第m面在该进程的内部编号
  do m=1, NUM_Face_ALL
    F_n(m)=0
    do k=1,m  ! m面前面有多数个面在F_proc(m)进程
     if(F_Proc(k) .eq. F_Proc(m))  F_n(m)=F_n(m)+1
    enddo
  enddo 

!  统计每个进程包含的面数
  NUM_Face=0
  do m=1, NUM_Face_ALL
    if(F_proc(m) .eq. my_id )   NUM_Face=NUM_Face+1    ! my_id进程包含的面数
  enddo
  
  allocate (my_Faces(NUM_Face))       ! 本进程包含的面号（数组）
  k=1
  do m=1, NUM_Face_ALL
     if(F_proc(m) .eq. my_id )  then
      my_Faces(k)=m
      k=k+1
     endif
  enddo

!---------------------------------------------------------------------------------
! 分配内存
  allocate(Face(NUM_Face))                      
  do m=1,NUM_Face        
    F => Face(m)
    mf = my_Faces(m)
    F%nx=NI(mf); F%ny=NJ(mf) ; F%nz=NK(mf)   ! nx,ny,nz 每面的大小
    nx=F%nx ; ny= F%ny ; nz=F%nz
    ! print*, "my_id, mf, nx, ny, nz:", my_id,mf,nx,ny,nz  
   allocate(F%x(nx,ny,nz),F%y(nx,ny,nz),F%z(nx,ny,nz),F%dS(nx,ny,nz))
   allocate(F%n1(nx,ny,nz),F%n2(nx,ny,nz),F%n3(nx,ny,nz),F%R(nx,ny,nz))
   allocate(F%r1(nx,ny,nz),F%r2(nx,ny,nz),F%r3(nx,ny,nz))
   allocate(F%d(nx,ny,nz,NUM_Frame),F%u(nx,ny,nz,NUM_Frame), &
   F%v(nx,ny,nz,NUM_Frame),F%w(nx,ny,nz,NUM_Frame),F%p(nx,ny,nz,NUM_Frame))
  enddo

! 表面几何信息读取
if(my_id .eq. 0) then
   print*, "read FWH_Surface_Geo.dat... "
    open(100,file="FWH_Surface_Geo.dat")
    read(100,*)
    read(100,*)
! ---------------------读取每个面的几何信息---------------------------   
  do m=1,NUM_Face_ALL        
    read(100,*) 
    nx=NI(m); ny=NJ(m) ; nz=NK(m)   ! nx,ny,nz 每面的大小
    allocate(U(7,nx,ny,nz))
    print*, "m=", m
    print*, "nx,ny,nz=", nx,ny,nz
    do k=1,nz
    do j=1,ny
    do i=1,nx
       read(100,*) U(1,i,j,k),U(2,i,j,k),U(3,i,j,k),U(4,i,j,k),U(5,i,j,k),U(6,i,j,k),U(7,i,j,k)
    enddo
    enddo
    enddo

    if(F_proc(m) .eq. 0) then            ! 这些面属于根进程
      mf=F_n(m)                          ! 该面在进程内部的编号
	    F=>Face(mf)
	    do k=1,nz
	    do j=1,ny
	    do i=1,nx
	      F%x(i,j,k) =U(1,i,j,k)
	      F%y(i,j,k) =U(2,i,j,k)
	      F%z(i,j,k) =U(3,i,j,k)
        F%n1(i,j,k)=U(4,i,j,k)
	      F%n2(i,j,k)=U(5,i,j,k)
	      F%n3(i,j,k)=U(6,i,j,k)
        F%dS(i,j,k)=U(7,i,j,k)
	    enddo
	    enddo
 	    enddo
    else                        ! 将该面数据发送出
	    call MPI_send(U,7*nx*ny*nz,FWH_DATA_TYPE, F_proc(m), F_n(m), MPI_COMM_WORLD,ierr )
    endif
    deallocate(U)
  enddo
  close(100)
else ! 非根节点
  do m=1,NUM_Face
    F=>Face(m)
    nx=F%nx; ny=F%ny; nz=F%nz
    allocate(U(7,nx,ny,nz))
    call MPI_Recv(U,7*nx*ny*nz,FWH_DATA_TYPE, 0, m, MPI_COMM_WORLD,Status,ierr )
    do k=1,nz
	  do j=1,ny
	  do i=1,nx
	    F%x(i,j,k) =U(1,i,j,k)
	    F%y(i,j,k) =U(2,i,j,k)
	    F%z(i,j,k) =U(3,i,j,k)
      F%n1(i,j,k)=U(4,i,j,k)
	    F%n2(i,j,k)=U(5,i,j,k)
	    F%n3(i,j,k)=U(6,i,j,k)
      F%dS(i,j,k)=U(7,i,j,k)
	  enddo
	  enddo
 	  enddo
    deallocate(U)
  enddo
endif
call MPI_Barrier(MPI_COMM_WORLD,ierr)

  allocate(Obs_x(NUM_Obs),Obs_y(NUM_Obs),Obs_z(NUM_Obs))
  if(my_id .eq. 0) then
    print*, "read Observers.dat... "
    open(88,file="Observers.dat")
    do m=1,NUM_Obs
      read(88,*) Obs_x(m),Obs_y(m),Obs_z(m)
      ! print*, " Observers No. =",m
      ! print*, " Observers x,y,z =",Obs_x(m),Obs_y(m),Obs_z(m)
    enddo
    close(88)
  endif
  call MPI_Bcast(Obs_x,NUM_Obs,FWH_DATA_TYPE,0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Obs_y,NUM_Obs,FWH_DATA_TYPE,0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Obs_z,NUM_Obs,FWH_DATA_TYPE,0, MPI_COMM_WORLD,ierr)
  if(my_id .eq. 0) print*, "read Observers.dat ok "
  
  do m=Kstep_start,Kstep_end,delta_step
    m2=(m-Kstep_start)/delta_step+1
    if(my_id .eq. 0) then
      write(filename, "('FWH-'I8.8'.dat')") m
      print*, "read ",filename
      if( FWH_data_Format .eq. 0) then   ! 二进制文件
        open(77,file=trim(filename),form="unformatted", access='stream', status='old')
      else
        open(77,file=trim(filename))
      endif
      do m1=1,NUM_Face_ALL
        nx=NI(m1); ny=NJ(m1) ; nz=NK(m1)
        allocate(U(5,nx,ny,nz))
        if( FWH_data_Format .eq. 0) then   ! 二进制文件
          read(77) (((U(1,i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                   (((U(2,i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                   (((U(3,i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                   (((U(4,i,j,k),i=1,nx),j=1,ny),k=1,nz) , &
                   (((U(5,i,j,k),i=1,nx),j=1,ny),k=1,nz)
        else
          do k=1,nz
          do j=1,ny
          do i=1,nx
            read(77,*) U(1,i,j,k),U(2,i,j,k),U(3,i,j,k),U(4,i,j,k),U(5,i,j,k)
            ! print *,"frame i,j,k rou u v d p",m2,i,j,k,U(1,i,j,k),U(2,i,j,k),U(3,i,j,k),U(4,i,j,k),U(5,i,j,k)
          enddo
          enddo
          enddo
        endif
        
        
        if(F_proc(m1) .eq. 0) then            ! 这些面属于根进程
          mf=F_n(m1)                          ! 该面在进程内部的编号
          F=>Face(mf)
          do k=1,nz
          do j=1,ny
          do i=1,nx
            F%d(i,j,k,m2)=U(1,i,j,k)
            F%u(i,j,k,m2)=U(2,i,j,k)
            F%v(i,j,k,m2)=U(3,i,j,k)
            F%w(i,j,k,m2)=U(4,i,j,k)
            F%p(i,j,k,m2)=U(5,i,j,k)
          enddo
          enddo
          enddo
        else                        ! 将该面数据发送出
          ! print*, "send to id,mf,nx, ny, nz:", F_proc(m1),F_n(m1),nx,ny,nz  
          call MPI_send(U,5*nx*ny*nz,FWH_DATA_TYPE, F_proc(m1), F_n(m1), MPI_COMM_WORLD,ierr )
        endif
        deallocate(U)
      enddo
      close(77)
    else ! 非根节点
      do m1=1,NUM_Face
        F=>Face(m1)
        nx=F%nx; ny=F%ny; nz=F%nz
        allocate(U(5,nx,ny,nz))
        call MPI_Recv(U,5*nx*ny*nz,FWH_DATA_TYPE, 0, m1, MPI_COMM_WORLD,Status,ierr)
        ! print*, "my_id,mf,nx, ny, nz, m, m2:", my_id,F_n(m1),nx,ny,nz,m,m2  
        do k=1,nz
        do j=1,ny
        do i=1,nx
          F%d(i,j,k,m2)=U(1,i,j,k)
          F%u(i,j,k,m2)=U(2,i,j,k)
          F%v(i,j,k,m2)=U(3,i,j,k)
          F%w(i,j,k,m2)=U(4,i,j,k)
          F%p(i,j,k,m2)=U(5,i,j,k)
        enddo
        enddo
        enddo
        deallocate(U)
      enddo
    endif
    call MPI_Barrier(MPI_COMM_WORLD,ierr)  
  enddo
    
  allocate(t_interp(NUM_Frame))
  if(my_id .eq. 0) allocate(pp(NUM_Frame))

end subroutine init  

!------------read parameter (Namelist type)---------------- 
  subroutine read_parameter_FWH 
   use Global_Variables
   implicit none
   integer:: Ipara(10),ierr
   real(PRE_EC):: rpara(10)
   namelist /control_FWH/ Ma, AoA,Kstep_start, Kstep_end, &
              NUM_Obs, delta_t, delta_step,NUM_THREADS,FWH_data_Format
	

!---------default--------
        Ma=0.2d0
        AoA=0.d0
        delta_t=5.8d-1
        NUM_Obs=1
        Kstep_start=1
        Kstep_end=100
        delta_step=1
        NUM_THREADS=5
        FWH_data_Format=0 ! 默认FWH data为二进制文件
!---------------------------------
  if(my_id .eq. 0) then
	  open(99,file="control.fwh")
	  read(99,nml=control_FWH)
    close(99)
    AoA=AoA*PI/180.d0   ! Angle of attack
    NUM_Frame=(Kstep_end-Kstep_start)/delta_step+1

    rpara(1)=Ma 
    rpara(2)=AoA
    rpara(3)=delta_t
    
    Ipara(1)=NUM_Obs
    Ipara(2)=Kstep_start
    Ipara(3)=Kstep_end
    Ipara(4)=delta_step
    Ipara(5)=NUM_THREADS
    Ipara(6)=NUM_Frame
    Ipara(7)=FWH_data_Format
  endif

  call MPI_Bcast(rpara,10,FWH_DATA_TYPE,0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Ipara,10,MPI_Integer,0, MPI_COMM_WORLD,ierr)
    Ma=rpara(1)    
    AoA=rpara(2)
    delta_t=rpara(3)
    
    NUM_Obs=Ipara(1)
    Kstep_start=Ipara(2)
    Kstep_end=Ipara(3)
    delta_step=Ipara(4)
    NUM_THREADS=Ipara(5)
    NUM_Frame=Ipara(6)
    FWH_data_Format=Ipara(7)

 end  subroutine read_parameter_FWH 
 
 !------compute effective acoustic distance for 2D motion cases---------------- 
 subroutine compute_R(o)
  use Global_Variables
  implicit none
  integer :: i,j,k,o,m,nx,ny,nz
  real(PRE_EC):: d1_o,d2_o,d3_o,R_s,beta2
  real(PRE_EC):: d1,d2,d3,r1,r2,r3
  Type (Face_TYPE),pointer:: F
  if(my_id .eq. 0) print*,"Obs_x Obs_y Obs_z",Obs_x(o),Obs_y(o),Obs_z(o)

  beta2=1-Ma**2
  do m=1,NUM_Face        
    F => Face(m)
    nx=F%nx ; ny= F%ny ; nz=F%nz
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,d1_o,d2_o,d3_o,d1,d2,d3,R_s,r1,r2,r3)	 	 
    do k=1,nz
    do j=1,ny
    do i=1,nx
        d1_o=Obs_x(o)-F%x(i,j,k)
        d2_o=Obs_y(o)-F%y(i,j,k)
        d3_o=Obs_z(o)-F%z(i,j,k)
        ! d1_o=3
        ! d2_o=2
        ! d3_o=1
        ! 坐标变换
        d1=d1_o*cos(AoA)+d2_o*sin(AoA)
        d2=-d1_o*sin(AoA)+d2_o*cos(AoA)
        d3=d3_o
        R_s=sqrt(d1**2+beta2*(d2**2+d3**2))  ! R*
        F%R(i,j,k)=(-Ma*d1+R_s)/beta2
        r1=(-Ma*R_s+d1)/(beta2*F%R(i,j,k))
        r2=d2/F%R(i,j,k)
        r3=d3/F%R(i,j,k)
        ! 坐标变换回去
        F%r1(i,j,k)=r1*cos(-AoA)+r2*sin(-AoA)
        F%r2(i,j,k)=-r1*sin(-AoA)+r2*cos(-AoA)
        F%r3(i,j,k)=r3
        ! print*, "face observe r1 r2 r3 R=",m,o,F%r1(i,j,k,o),F%r2(i,j,k,o),F%r3(i,j,k,o),F%R(i,j,k,o)
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO 
 enddo
end  subroutine compute_R 

!---------- compute pp_ret ----------------
! 不储存Qn和Lm以减小内存消耗
subroutine compute_Noise
  use Global_Variables
  implicit none
  integer :: i,j,k,frame,m,nx,ny,nz
  real(PRE_EC):: M1,M2,vel_n,L1,L2,L3,MR,Qn_dot,Lr_dot,P_T,P_L
  Type (Face_TYPE),pointer:: F

  M1=-Ma*cos(AoA)
  M2=-Ma*sin(AoA)

  do m=1,NUM_Face        
    F => Face(m)
    nx=F%nx ; ny= F%ny ; nz=F%nz

!------------ compute Qn Lm Lr ------------------------------------
    allocate(F%Qn(nx,ny,nz,NUM_Frame),F%Lm(nx,ny,nz,NUM_Frame),F%Lr(nx,ny,nz,NUM_Frame))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(frame,i,j,k,vel_n,L1,L2,L3)
    do frame=1,NUM_Frame
      do k=1,nz
      do j=1,ny
      do i=1,nx
        ! print*, "Face d u v w p=",m,d(i,j,k),u(i,j,k),v(i,j,k),w(i,j,k),p(i,j,k)
        vel_n=F%u(i,j,k,frame)*F%n1(i,j,k)+F%v(i,j,k,frame)*F%n2(i,j,k)+F%w(i,j,k,frame)*F%n3(i,j,k)
        L1=F%p(i,j,k,frame)*F%n1(i,j,k)+F%d(i,j,k,frame)*(F%u(i,j,k,frame)-cos(AoA))*vel_n
        L2=F%p(i,j,k,frame)*F%n2(i,j,k)+F%d(i,j,k,frame)*(F%v(i,j,k,frame)-sin(AoA))*vel_n
        L3=F%p(i,j,k,frame)*F%n3(i,j,k)+F%d(i,j,k,frame)*F%w(i,j,k,frame)*vel_n
        F%Lm(i,j,k,frame)=L1*M1+L2*M2
        ! print*, "Face Lm=",m,F%Lm(i,j,k,frame1)
        F%Qn(i,j,k,frame)=-(cos(AoA)*F%n1(i,j,k)+sin(AoA)*F%n2(i,j,k))+F%d(i,j,k,frame)*vel_n
        ! print*, "Face Qn=",m,F%Qn(i,j,k,frame1)
        F%Lr(i,j,k,frame)=L1*F%r1(i,j,k)+L2*F%r2(i,j,k)+L3*F%r3(i,j,k)
      enddo
      enddo
      enddo
    enddo
!$OMP END PARALLEL DO 

!----------------------- Compute pp_ret ---------------------------------
    allocate(F%pp_ret(nx,ny,nz,NUM_Frame))
! compute Lr_dot & Qn_dot
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(frame,i,j,k,Qn_dot,Lr_dot,MR,P_T,P_L)
    do frame=1,NUM_Frame
      do k=1,nz
      do j=1,ny
      do i=1,nx
        if (frame .eq. 1) then
          Qn_dot = (-F%Qn(i,j,k,frame+2)+4.d0*F%Qn(i,j,k,frame+1)-3.d0*F%Qn(i,j,k,frame))/(2.d0*delta_t)
          Lr_dot = (-F%Lr(i,j,k,frame+2)+4.d0*F%Lr(i,j,k,frame+1)-3.d0*F%Lr(i,j,k,frame))/(2.d0*delta_t)
        elseif (frame .eq. NUM_Frame) then
          Qn_dot = (3.d0*F%Qn(i,j,k,frame)-4.d0*F%Qn(i,j,k,frame-1)+F%Qn(i,j,k,frame-2))/(2.d0*delta_t)
          Lr_dot = (3.d0*F%Lr(i,j,k,frame)-4.d0*F%Lr(i,j,k,frame-1)+F%Lr(i,j,k,frame-2))/(2.d0*delta_t)
        else
          Qn_dot = (F%Qn(i,j,k,frame+1)-F%Qn(i,j,k,frame-1))/(2.d0*delta_t)
          Lr_dot = (F%Lr(i,j,k,frame+1)-F%Lr(i,j,k,frame-1))/(2.d0*delta_t)
        endif

        MR=M1*F%r1(i,j,k)+M2*F%r2(i,j,k)
        ! c_0=1/Ma
        P_T=Qn_dot/(F%R(i,j,k)*(1.d0-MR)**2) + F%Qn(i,j,k,frame)*(1.d0/Ma)*(MR-Ma**2)/(F%R(i,j,k)**2*(1.d0-MR)**3)
        P_L=Lr_dot*Ma/(F%R(i,j,k)*(1.d0-MR)**2) + & 
            (F%Lr(i,j,k,frame)-F%Lm(i,j,k,frame))/(F%R(i,j,k)**2*(1.d0-MR)**2) + &
            (F%Lr(i,j,k,frame)*(Mr-Ma**2))/(F%R(i,j,k)**2*(1.d0-MR)**3)
        F%pp_ret(i,j,k,frame)=(P_T+P_L)*F%dS(i,j,k)/(4.d0*PI)
      enddo
      enddo
      enddo
    enddo
!$OMP END PARALLEL DO
    deallocate(F%Qn,F%Lm,F%Lr)
    ! print*, "Done computing pp_ret at Face:", m
  enddo

 

end  subroutine compute_Noise 

!------compute observer time---------------- 
subroutine compute_ObserverTime 
  use Global_Variables
  implicit none
  integer :: frame,i,j,k,m,nx,ny,nz,ierr
  real(PRE_EC):: r_min_F(NUM_Face),r_max_F(NUM_Face),r_min,r_max
  real(PRE_EC):: t_interp_start,t_interp_end,dt_interp
  Type (Face_TYPE),pointer:: F

! compute r_min & r_max
!$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(NUM_Face,r_min_F,r_max_F)
  do m=1,NUM_Face        
      F => Face(m)
      r_min_F(m)=minval(F%R(:,:,:))
      r_max_F(m)=maxval(F%R(:,:,:))
      ! print*, "Face Observer r_min r_max",m,o,r_min_F(m),r_max_F(m)
  enddo
!$OMP END PARALLEL DO 

  call MPI_ALLREDUCE(minval(r_min_F),r_min,1,FWH_DATA_TYPE,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(maxval(r_max_F),r_max,1,FWH_DATA_TYPE,MPI_MAX,MPI_COMM_WORLD,ierr)
  ! print*, "r_min r_max",r_min,r_max

! compute t_obs
do m=1,NUM_Face        
    F => Face(m)
    nx=F%nx ; ny= F%ny ; nz=F%nz
    allocate(F%t_obs(nx,ny,nz,NUM_Frame))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(frame,i,j,k)	 	 
  do frame=1,NUM_Frame
    do k=1,nz
    do j=1,ny
    do i=1,nx
				F%t_obs(i,j,k,frame) = delta_t*(frame-1) + F%R(i,j,k)*Ma  ! c_0=1/Ma
    enddo
    enddo
    enddo
  enddo
!$OMP END PARALLEL DO 
enddo

! compute t_interp
  t_interp_start = r_max*Ma  ! c_0=1/Ma
  t_interp_end = delta_t*(NUM_Frame-1) + r_min*Ma ! c_0=1/Ma
  dt_interp = (t_interp_end - t_interp_start)/(NUM_Frame-1)
  ! print*, "t_interp_start t_interp_end dt_interp",t_interp_start,t_interp_end,dt_interp
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(frame)	 	 
  do frame=1,NUM_Frame
    t_interp(frame) = t_interp_start + dt_interp*(frame-1)
    ! print*, "frame t_interp",frame,t_interp(frame)
  enddo
!$OMP END PARALLEL DO 

end  subroutine compute_ObserverTime 

!-------- Interp Pressure Signal ------------- 
subroutine Interp_PressureSignal 
  use Global_Variables
  implicit none
  integer :: i,j,k,m,nx,ny,nz
  Type (Face_TYPE),pointer:: F

  do m=1,NUM_Face        
    F => Face(m)
    nx=F%nx ; ny= F%ny ; nz=F%nz
    allocate(F%pp_interp(nx,ny,nz,NUM_Frame))
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)	 	 
    do k=1,nz
    do j=1,ny
    do i=1,nx
        call interp1d_Cubic(F%t_obs(i,j,k,:),F%pp_ret(i,j,k,:),NUM_Frame,t_interp(:),F%pp_interp(i,j,k,:))
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO 
    deallocate(F%t_obs,F%pp_ret)
    ! print*, "Done Interp Pressure Signal at Face:", m
 enddo
end  subroutine Interp_PressureSignal 

!-------- Integrate Sources ------------- 
subroutine Integrate_Sources 
  use Global_Variables
  implicit none
  integer :: frame,m,ierr
  real(PRE_EC) :: pp_f(NUM_Face),pp0(NUM_Frame)
  Type (Face_TYPE),pointer:: F

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(frame,m,F,pp_f)	 	 
    do frame=1,NUM_Frame
        do m=1,NUM_Face        
          F => Face(m)
          pp_f(m) = sum(F%pp_interp(:,:,:,frame))
        enddo
        pp0(frame)=sum(pp_f)
    enddo
!$OMP END PARALLEL DO 
    call MPI_Reduce(pp0,pp,NUM_Frame,FWH_DATA_TYPE,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    do m=1,NUM_Face        
      F => Face(m)
      deallocate(F%pp_interp)
    enddo

end  subroutine Integrate_Sources 

!-------------- Cubic Spline Interpolation ----------------------
! Reference :http://phys.uri.edu/nigh/NumRec/bookfpdf/f3-3.pdf
subroutine spline(x,y,n,y2)
  use Global_Variables
  implicit none
  integer, intent(in) :: n
  real(PRE_EC), intent(in) :: x(n),y(n)
  real(PRE_EC), intent(out) :: y2(n)
  integer :: i,k
  real(PRE_EC) :: p,qn,sig,un,u(n)
    y2(1)=0.d0
    u(1)=0.d0
    do i=2,n-1      ! This is the decomposition loop of the tridiagonal algorithm.
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
         /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    qn=0.d0
    un=0.d0
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
    do k=n-1,1,-1 
        y2(k)=y2(k)*y2(k+1)+u(k) ! This is the backsubstitution loop of the tridiagonal algorithm.
    enddo 
end
    
subroutine interp1d_Cubic(xa,ya,n,x,y)
  use Global_Variables
  implicit none
    integer, intent(in) :: n
    real(PRE_EC), intent(in) :: x(n),xa(n),ya(n)
    real(PRE_EC), intent(out) :: y(n)
    real(PRE_EC) :: y2a(n)
    integer :: i,k,khi,klo
    real(PRE_EC) :: a,b,h
    call spline(xa,ya,n,y2a)
    klo=1    
    khi=n-1
    do i = 1,n
        if (khi .lt. n) then
        khi = khi + 1
        else
        khi = n
        endif
        ! print *, khi, khi-klo
1   if (khi-klo .gt. 1) then
        k=(khi+klo)/2
        if(xa(k) .gt. x(i))then
            khi=k
        else
            klo=k
        endif
        ! print *,klo, khi
    goto 1
    endif
    ! print *,klo, khi, khi-klo
    
    h=xa(khi)-xa(klo)
    a=(xa(khi)-x(i))/h
    b=(x(i)-xa(klo))/h
    y(i)=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
    enddo
end subroutine interp1d_Cubic
