!-----------------------------------------------------------
! OpenCFD-FWH (OpenMP only)
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
  integer,parameter:: PRE_EC=8
end module Const_Variables

  module Global_Variables
  use Const_Variables
  implicit none
 ! global parameter (for all Faces) 
  real(PRE_EC),parameter:: PI=3.14159265358979d0
  integer,save:: Kstep_start,Kstep_end,NUM_Obs,delta_step,FWH_data_Format
  real(PRE_EC),save:: Ma,AoA,delta_t
  integer,save:: NUM_Face,NUM_THREADS,NUM_Frame,Total_Faces
  real(PRE_EC),pointer,dimension(:):: Obs_x,Obs_y,Obs_z
  real(PRE_EC),pointer,dimension(:):: t_interp,pp

!-----------------------------------------------------------------------------------------

!  核心变量——每块网格存储的信息 （全局变量）
   TYPE Face_TYPE           !  variables for each Face 
     integer :: nx,ny,nz
     real(PRE_EC),pointer,dimension(:,:,:):: x,y,z,dS,n1,n2,n3  ! coordinates of cell center, 网格中心坐标 
     real(PRE_EC),pointer,dimension(:,:,:):: R,r1,r2,r3         ! R & r1,r2,r3 
	 real(PRE_EC),pointer,dimension(:,:,:,:):: d,u,v,w,p,Lm,Qn,Lr
   real(PRE_EC),pointer,dimension(:,:,:,:):: pp_ret,t_obs,pp_interp
   End TYPE Face_TYPE  
   TYPE (Face_TYPE), save,dimension(:),allocatable,target:: Face
  
  end module Global_Variables
!------------------------------------------------------------------


!-------------------------------------------
  program compute_farfield_noise   
     use Global_Variables
     implicit none
     integer :: obs,frame,cr,cm,c1,c2
     real(PRE_EC) :: rate
     character(len=50):: filename
      
     ! initialize the system_clock
      call system_clock(count_rate=cr)
      call system_clock(count_max=cm)
      rate = REAL(cr)
      call system_clock(c1)
	    
      call init
      call system('mkdir FWH_result')

    do obs = 1,NUM_Obs
      
      print*, "Processing Observer :",obs

      call compute_R(obs)
      print*, "Compute R OK "

      call compute_Noise
      print*, "Compute pp_Ret OK "

      call compute_ObserverTime
      print*, "Compute Observer Time OK "

      call Interp_PressureSignal
      print*, "Interp Pressure Signal OK "

      call Integrate_Sources
      print*, "Integrate Sources OK, write pressure signals ... "
      
      write(filename,"('./FWH_result/p_observer-'I3.3'.log')") obs
      open(222,file=trim(filename))
      do frame =1,NUM_Frame
      write(222,"(2E30.16)") t_interp(frame),pp(frame)
      enddo
      close(222)
      print *,"Done writing",filename
      print *

    enddo
      
      call system_clock(c2)
      print*, "run time : ",(c2 - c1)/rate

  end

!=====================================================================================
!------------------------------------------------------------------------------     
! Read the message of the mesh and the initial flow;
! 读取网格，初始流场信息; 
! 分配内存变量；
! 计算几何量；
!------------------------------------------------------------------------------
   subroutine init
   use  Global_Variables
   implicit none
   integer :: i,j,k,m,nx,ny,nz,m1,m2
   integer,allocatable,dimension(:):: NI,NJ,NK
   Type (Face_TYPE),pointer:: F
   character(len=50):: filename

!------------------------------------------------------------------
    call read_parameter_FWH
!$ call omp_set_num_threads(NUM_THREADS)   ! 设置OpenMP的运行线程数 （并行数目）， 本语句对openmp编译器不是注释!
!测试一下运行的进程 （openmp编译时，不是注释）
!$OMP Parallel
!$  print*, "omp run ..."
!$OMP END parallel
    print*,  "---------------------------- OpenCFD-FWH (OpenMP version) --------------------------------"
    print*,  "                      Copyright by Zhang Keli, zkl70@buaa.edu.cn                          "
    print*,  "                         Programming by Zhang Keli  2023-12                               "
    print*,  " Ma,AoA,delta_t=", Ma,AoA,delta_t
    print*,  " NUM_Obs,NUM_Frame,NUM_THREADS=",NUM_Obs,NUM_Frame,NUM_THREADS
    print*,  "------------------------------------------------------------------------------------------" 
    print*
!  表面几何信息文件   
   print*, "read FWH_Surface_Geo.dat... "
    open(100,file="FWH_Surface_Geo.dat")
    read(100,*)
    read(100,*) NUM_Face         ! 总块数
   print*, "NUM_Face=", NUM_Face
   
    allocate(Face(NUM_Face))             
    allocate(NI(NUM_Face),NJ(NUM_Face),NK(NUM_Face) )   ! 每面的大小
! 读取每块信息----------------------------------------   
    do m=1,NUM_Face        
     F => Face(m)
     read(100,*) NI(m),NJ(m),NK(m) 
     F%nx=NI(m); F%ny=NJ(m) ; F%nz=NK(m)   ! nx,ny,nz 每面的大小
     nx=F%nx ; ny= F%ny ; nz=F%nz
! ----------  几何量 -----------------------------------------------
    allocate(F%x(nx,ny,nz),F%y(nx,ny,nz),F%z(nx,ny,nz),F%dS(nx,ny,nz))
    allocate(F%n1(nx,ny,nz),F%n2(nx,ny,nz),F%n3(nx,ny,nz),F%R(nx,ny,nz))
    allocate(F%r1(nx,ny,nz),F%r2(nx,ny,nz),F%r3(nx,ny,nz))
    allocate(F%d(nx,ny,nz,NUM_Frame),F%u(nx,ny,nz,NUM_Frame), &
    F%v(nx,ny,nz,NUM_Frame),F%w(nx,ny,nz,NUM_Frame),F%p(nx,ny,nz,NUM_Frame))

   print*, "m=", m
   print*, "nx,ny,nz=", nx,ny,nz
    do k=1,nz
    do j=1,ny
    do i=1,nx
       read(100,*) F%x(i,j,k),F%y(i,j,k),F%z(i,j,k),F%n1(i,j,k),F%n2(i,j,k),F%n3(i,j,k),F%dS(i,j,k)
    enddo
    enddo
    enddo

  enddo
  close(100)

  Total_Faces=0
  do m=1,NUM_Face
    Total_Faces=Total_Faces + NI(m)*NJ(m)*NK(m) 
  enddo
  print*, "Total_Faces =",Total_Faces

  print*, "read Observers.dat... "
  allocate(Obs_x(NUM_Obs),Obs_y(NUM_Obs),Obs_z(NUM_Obs))
  open(88,file="Observers.dat")
  do m=1,NUM_Obs
    read(88,*) Obs_x(m),Obs_y(m),Obs_z(m)
  enddo
  close(88)
  print*, "read Observers.dat ok "

  do m=Kstep_start,Kstep_end,delta_step
    write(filename, "('FWH-'I8.8'.dat')") m
    m2=(m-Kstep_start)/delta_step+1
    print*, "read ",filename
    if( FWH_data_Format .eq. 0) then   ! 二进制文件
      open(77,file=trim(filename),form="unformatted")
    else
      open(77,file=trim(filename))
    endif
      do m1=1,NUM_Face
        F => Face(m1)
        nx=F%nx ; ny= F%ny ; nz=F%nz
        if( FWH_data_Format .eq. 0) then   ! 二进制文件
          read(77)   (((F%d(i,j,k,m2),i=1,nx),j=1,ny),k=1,nz) , &
                     (((F%u(i,j,k,m2),i=1,nx),j=1,ny),k=1,nz) , &
                     (((F%v(i,j,k,m2),i=1,nx),j=1,ny),k=1,nz) , &
                     (((F%w(i,j,k,m2),i=1,nx),j=1,ny),k=1,nz) , &
                     (((F%p(i,j,k,m2),i=1,nx),j=1,ny),k=1,nz)
        else
          do k=1,nz
          do j=1,ny
          do i=1,nx
            read(77,*) F%d(i,j,k,m2),F%u(i,j,k,m2),F%v(i,j,k,m2),F%w(i,j,k,m2),F%p(i,j,k,m2)
          enddo
          enddo
          enddo
        endif
      enddo
      close(77)  
    enddo

    allocate(t_interp(NUM_Frame),pp(NUM_Frame))

end subroutine init  

!------read parameter (Namelist type)---------------- 
  subroutine read_parameter_FWH 
   use Global_Variables
   implicit none
   namelist /control_FWH/ Ma, AoA,Kstep_start, Kstep_end, &
              NUM_Obs, delta_t, delta_step,NUM_THREADS,FWH_data_Format
	

!---- default--------
        Ma=0.2d0
        AoA=0.d0
        delta_t=5.8d-1
        NUM_Obs=1
        Kstep_start=1
        Kstep_end=100
        delta_step=1
        NUM_THREADS=10
        FWH_data_Format=0 ! 默认FWH data为二进制文件
!---------------------------------
	open(99,file="control.fwh")
	read(99,nml=control_FWH)
  close(99)
  AoA=AoA*PI/180.d0   ! Angle of attack
  NUM_Frame=(Kstep_end-Kstep_start)/delta_step+1
 end  subroutine read_parameter_FWH 
 
 !------compute effective acoustic distance for 2D motion cases---------------- 
 subroutine compute_R(o)
  use Global_Variables
  implicit none
  integer :: i,j,k,o,m,nx,ny,nz
  real(PRE_EC):: d1_o,d2_o,d3_o,R_s,beta2
  real(PRE_EC):: d1,d2,d3,r1,r2,r3
  Type (Face_TYPE),pointer:: F
  print*,"Obs_x Obs_y Obs_z",Obs_x(o),Obs_y(o),Obs_z(o)

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

  do m=1,Num_Face        
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
  integer :: frame,i,j,k,m,nx,ny,nz
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

  r_min=minval(r_min_F)
  r_max=maxval(r_max_F)
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
  integer :: frame,m
  real(PRE_EC) :: pp_f(NUM_Face)
  Type (Face_TYPE),pointer:: F

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(frame,m,F,pp_f)	 	 
    do frame=1,NUM_Frame
        do m=1,NUM_Face        
          F => Face(m)
          pp_f(m) = sum(F%pp_interp(:,:,:,frame))
        enddo
        pp(frame)=sum(pp_f)
    enddo
!$OMP END PARALLEL DO 
    
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