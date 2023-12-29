# Welcome to OpenCFD-FWH

OpenCFD-FWH is a far-file noise prediction code based on a permeable surface nondimensional FW-H (Ffowcs Williams-Hawkings) acoustics analogy with convective effect and AoA corrections. It is developed to use as a post-processing code for our finite volume CFD solver [OpenCFD-EC](https://opencfd.cn/downloads/)， but **it can be used with other solvers** with the right data format, or by modifying the data reading procedure of the code.

More detail information about the methodology, framework,  validation, and input data format of the code pleased refer to our preprint paper: [An MPI-OpenMP mixing parallel open source FW-H code for aeroacoustics calculation](https://arxiv.org/abs/2312.16263). 

### Features:
- Nondimensional FW-H acoustics analogy for permeable data surface
- Convective and AoA effect taken into considertion
- MPI-OpenMP mixing parallel
- Only required an MPI library and a Fortran 90 compilation environment (Intel MPI library is recommended)
- Only OpenMP parallel version is provided for single workstation users

# Tutorials

In the *Tutorials* folder Matlab programs for monopole and dipole validation cases are provided for tutorial purpose.

## For Linux users
1. Make sure you install a MPI library. You can run```sudo apt install openmpi``` to install OpenMPI library, or  install [Intel MPI library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html), which is more recommended for better performance.

2. Compile the OpenCFD-FWH.f90 by:
```bash
	mpif90 -o OpenCFD-FWH.o OpenCFD-FWH.f90 -fopenmp
```
3. Move the OpenCFD-FWH. o file to the FWH-monpole or *FWH-dipole* folder within the *monopole* or *dipole* folder.
4. Enter the *monopole* or *dipole* folder, there are two .m file in each folder.  Create a new folder *FWH-monopole* or *FWH-dipole*, according to which folder you are in.
5.  Run the _samlpe.m files by Matlab, it will generate the control.fwh, FWH_Surface_Geo.dat, Observers.dat, and sampling dataset on the permable FW-H surface in the *FWH-monpole/dipole* folder for OpenCFD-FWH. Also a matlab.control file is created for the other .m file to read in the parameters.
6.  After the completion of the execution of the _sample.m file (about 30s with the default setting on a lapton with i7-13620H CPU), run the following command if you are in the *monpole* folder:
```bash
	cd FWH-monpole
	mpirun -n 1 ./OpenCFD-FWH.o
```
7.  Only one MPI processor is used due to the permeable FW-H data surface is generated in one Face by the _samlpe.m files. The running time of the code is just in  2.5 senconds with 12 OpenMP threads. You can change the number of OpenMP threads through the value of NUM_THREADS in the control.fwh file. And the pressure signal result of each observer is in the *FWH_result-mpi* folder.
8. Go back to the *monpole* or *dipole* folder, run the monpole.m/dipole .m file by Matlab, it will compare the results of the code with the exact analytical solution.
9. You can change the parameters in the  _samlpe.m files for different monpole/dipole conditions and FW-H data sampling strategies.
 
## For Windows users
- You can follow the same procedure of [For Linux users](https://github.com/Z-K-L/OpenCFD-FWH?tab=readme-ov-file#for-linux-users), expect compile the .f90 file to a .exe file and use ```mpiexe -n 1 .\OpenCFD-FWH.exe``` to excute the code.
- Or you can place the OpenCFD-FWH_OpenMP-only.exe file in the *FWH-monopole* or *FWH-dipole* folder, and excute it (*even you don't have OpenMP install in your system*) after the finish of the _sample.m file's execution 
 in the *monopole/dipole* folder. Note the the only OpenMP version will output the results in the *FWH_result* folder. You need to rename it as *FWH_result-mpi* or change the *FWH_result-mpi* to *FWH_result* in the FWH results reading part of the monopole.m/dipole.m file. Then folowing the 8. and 9. step in [For Linux users](https://github.com/Z-K-L/OpenCFD-FWH?tab=readme-ov-file#for-linux-users).

# Hybrid parallel acceleration

Initialization, computation, and overall execution time of  OpenCFD-FWH runs in different MPI processors and OpenMP threads for 18252 subfaces, 6535 sampling frames, and 40 observers of a 30P30N case on the CAS SunRising platform with Intel MPI library being deployed：
MPI processors|OMP threads|init time/s|computing time/s|total time/s|computing acceleration ratio|total acceleration ratio
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
1 |  ×  | 753.7 | 18120.3 | 18873.9 |  \  |  \ 
31 |  ×  | 887.2 | 811.3 | 1698.5 | 22.3 | 11.1 
1 | 32 | 823.9 | 774.3 | 1598.1 | 23.4 | 11.8 
31 | 32 | 717.5 | 33.7 | 751.2 | 538.5 | 25.1 

# Citation of OpenCFD-FWH
If you use the OpenCFD-FWH code for academic research, please cite the following paper:
(*under review*)

# For problems
Any problems encounter with the code, you can publish them on Issues.
