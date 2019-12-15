# XFLAT
XFLAT, is developed to study neutrino oscillations in supernovae. 
XFLAT is designed to utilize multiple levels of parallelism through MPI, OpenMP, and SIMD instructions (vectorization). 
It can run on both the CPU and the Xeon Phi co-processor, the latter of which is based on the Intel Many Integrated Core Architecture (MIC).

## Compilation and Build Instructions
XFLAT is a command line application intended to study neutrino flavor oscillations in supernovae environments. The code is C++ implementation with the hybrid architecture that exploits SIMD, OpenMP and MPI for performance acceleration. It is capable to be run on heterogeneous supercomputers and can utilize both traditional CPUs and the newer Intel Many Integrated Core Architecture or Intel MIC (Xeon Phi).

The code contains several modules that can be swapped in or out from the build using the provided switches via the Makefile. In addition to the modules, other features can be switched on or off from the Makefile as well. Features such as the usage of SIMD, OpenMP, and MPI are controllable from the Makefile. Turning off optimizations can help the debugging process. 

The following code features can be switched on or off from the second line of the provided Makefile starting with CXXFLAGS, and using -D switch:
```
SIMD # when defined the code will use SIMD instructions
OMP # when defined the code will use OpenMP threads
MMPI # when defined the code will use MPI
```

In addition there are multiple modules that can be changed from the same line of the Makefile. These modules can be categorized as follow: 

Geometry related modules:
```
SA # Single Angle supernova module
MA # Multi Angle supernova module
MAA # Multi Azimutal Angle extended supernova module
CLN # Cylindrical module
PNT # Plane module
LIN # Multi-source line module IO
```

related modules:
```
IOF # performs file IO for each node
IOFI # performs indirect IO in which MIC sends its data to CPU first
```

An example of the line with the usage of IOF and MAA modules which utilizes SIMD,
OpenMP and MPI is shown as below:
```
CXXFLAGS = -O3 -openmp -DIOF -DMAA -DSIMD -DOMP -DMMPI
```

The code employs the NetCDF(either version 3 or 4) library for its current IO modules, however, if the NetCDF4 is used, HDF5 library is required as well.
In order to build binary for CPU the following commands must be issued from console:
```
$ cp Makefile.cpu Makefile
$ make all
```
OR
```
$ make all -f Makefile.cpu
```

Likewise, in order to build the code for the Xeon Phi the following commands must be issued:
```
$ cp Makefile.mic Makefile
$ make all
```
OR
```
$ make all -f Makefile.mic
```
Consequently, the CPU binary will be called XFLAT.cpu and the Xeon Phi binary is XFLAT.mic. Please note, if the OpenMP feature is switched on, the -openmp flag also must be added to the compiler flags and if the MPI feature is switched on inside the Makefile, a few MPI scripts are required so as to run it on multiple nodes.

In order to optimize the code for any Intel CPU, one has to add -xHOST to the compiler flags set in the Makefile. Likewise, in order to build the code for first generation MIC the -mmic flag should be added to the compiler flags set. For the compilation on the second generation of MIC the -xMIC-AVX512 flag shall be used.