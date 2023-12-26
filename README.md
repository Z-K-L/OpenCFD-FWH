# Welcome to OpenCFD-FWH

OpenCFD-FWH is a far-file noise prediction code based on a permeable surface nondimensional FW-H (Ffowcs Williams-Hawkings) acoustics analogy with convective effect and AoA corrections. It is developed to use as a post-processing code for our finite volume CFD solver [OpenCFD-EC](https://opencfd.cn/downloads/)， but **it can be used with other solvers** with the right data format, or by modifying the data reading procedure of the code.

More detail information about the methodology, framework,  validation, and input data format of the code pleased refer to [our preprint paper](222) (will be released soon). 

### Features:
- Nondimensional FW-H acoustics analogy for permeable data surface
- Convective and AoA effect taken into considertion
- MPI-OpenMP mixing parallel
- Only required an MPI library and a Fortran 90 compilation environment (Intel MPI library is recommended)

# Tutorials

In the Tutorials folder Matlab programs for monopole and dipole validation cases are provided for tutorial purpose.
1. Make sure you install a MPI library. If you are in a Linux system, you can run
```bash
	sudo apt install openmpi
```
&emsp;&emsp; or you can install [Intel MPI library](https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html), which is more recommended for better performance.

2. Compile the OpenCFD-FWH.f90 by :
```bash
	mpif90 -o OpenCFD-FWH.o OpenCFD-FWH.f90 -fopenmp
```
3. Move the OpenCFD-FWH. o file to the FWH-monpole or FWH-dipole folder within the monopole or diople folder.
4. Enter the monopole or diople folder, there are two .m file in each folder.  Create a new folder :FWH-monopole or FWH-dipole, according to which folder you are in.
7.  Run the _samlpe.m files by Matlab will generate the control.fwh, FWH_Surface_Geo.dat, Observers.dat, and the sampling dataset on the permable FW-H surface in the FWH-monpole/dipole folder for OpenCFD-FWH. Also a matlab.control file is generated for the other .m file to read in the parameters.
8.  After the completion of the execution of the _sample.m file (about 30s with the default setting on a lapton with i7-13620H CPU), run the following command if you are in the monpole folder:
```bash
	cd FWH-monpole
	mpirun -n 1 ./OpenCFD-FWH.o
```
9.  Only one MPI processor is used due to the permeable FW-H data surface is generated in one Face by the _samlpe.m files. The running time of the code is just in  2.5 senconds with 12 OpenMP threads. You can change the number of OpenMP threads through the value of NUM_THREADS in the control.fwh file. And the pressure signal result of each observer is in the FWH_result-mpi folder.
10. Go back to the monpole or dipole folder, run the monpole.m/dipole .m file by Matlab, it will compare the results of the code with the exact analytical solution.
11. You can change the parameters in the  _samlpe.m files for different monpole/dipole conditions and FW-H data sampling strategies.
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
Any problems encounter with the code, you can publish themon Issues。
