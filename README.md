## About this repository 

This repository implements a variant on the Force Coupling Method for Stokes suspensions, based on the exponential of a semicircle (ES) kernel, see:
[CITE THE ARTICLE HERE]
This repo also contains code and scripts to reproduce the data for the figures in the applications section of the article (see the README inside each folder for information about each simulation).

We provide a simple python interface but the key performance-sensitive pieces (namely FFTs, BVP solvers, and spreading and interpolation) are implemented in C++/CUDA using FFTW/cuFFT and LAPACK. One can use either the GPU and CPU versions as needed without changing calls.

The doubly periodic Stokes solver for GPUs is provided by [UAMMD](https://github.com/RaulPPelaez/uammd), which is included as a submodule inside the `source/gpu` folder. **Make sure to clone this repo recursively** by doing:
```shell
git clone --recursive https://github.com/stochasticHydroTools/DPStokesTests
```
At the moment, only a python interface to use the module is provided.

The CPU version of the solver is included in `source/cpu`. In particular this implements an OpenMP-based C++ spreading and interpolation library in 3D that supports also non-uniform grids in the z direction. One can use the C++ library directly if desired, but there is a python interface as well.

## Installation

To be able to use either the GPU or CPU versions on demand you will need to have reasonably new versions of CUDA, GNU or Intel C++ compilers, and python 3. For example, on Courant machines you can use something like:
```shell
module load cuda-10.2
module load intel-2019
module load gcc-6.4.0
module load python-3.8
```
From DoublyPeriodicStokes, Running 
```shell
make 
```
will compile and set up both the CPU and GPU python interfaces for the solver.

The top level `Makefile` in DoublyPeriodicStokes contains a section where a user
can specify the dependency library names/paths, install paths and the like.
[Donev: Moved this here] Users should source the bash script `cpuconfig.sh` before using either 
the GPU or CPU Python interface in a new shell, and can edit the thread environment 
settings therin as needed. The PYTHONPATH and LD_LIBRARY_PATH environment variables
are appeneded to so that the modules can be used anywhere within the filesystem.
By default, the script will exist in the $INSTALL_DIR specified in the top Makefile.

If you want to use the Intel compiler for the CPU code, prefix the call to make as
```shell
CPU=Intel make
``` 
Note, even if using the Intel compiler, you must load the module for gcc-6.4.0 or higher, 
as the compiler relies on GNU headers. Also note that by default with Intel compilers the (MKL library)[https://en.wikipedia.org/wiki/Math_Kernel_Library] is used to provide LAPACK/BLAS functionality.
 
You can compile both CPU and GPU libraries in debug node through
```shell
DEBUG=True make
```
Both CPU and DEBUG can also be set from within the Makefile, though the 
command line setting will override the ones in the Makefile.

## Python Interface

A common python interface is provided that allows to compute the hydrodynamic displacements for a group of positions with forces and/or torques acting on them in different geometries, mainly:  

	* Triply periodic, using a standard FFT-based spectral Stokes solver.
	* Doubly periodic (either with no walls, a bottom wall or a slit channel), see paper for details.

The modules will use the ES kernel (called BM in UAMMD). Other kernels (e.g., Gaussian, as in the traditional FCM method) can be used, but this simple interface does not allow them for the moment.  

Hydrodynamic displacements coming from forces and torques can be computed. 
For the GPU interface, if the torque-related arguments are ommited, the computations related to them are skipped entirely.
For the CPU interface, the user must specify whether torques are involved with a boolean swith, like `has_torque=True`.
        
The file `python_interface/common_interface_wrapper.py` is the joint CPU-GPU interface. 
The file `python_interface/dpstokes_common.py` is an example of using the joint interface. 

### GPU Python interface

The GPU interface requires [pybind11](https://github.com/pybind/pybind11) to compile, [Donev: Raul, aren't we changing this?] which is included as a submodule and will be automatically downloaded if this repo is cloned recursively (see "About this repo" above).  
In order to use it you must compile the python wrappers using make inside python_interface.  
A file called uammd.*.so will be created and then "import uammd" can be used inside python. 
See `python_interface/dpstokesGPU.py` for a usage example. Once compiled and imported you can use "help(uammd)" in python for additional usage information.  

### CPU Python interface

See the `source/cpu/README.md` for details. Note, the build instructions contained therein are for using cmake3 as the build system. 
The section can be ignored, or followed analogously through the provided top level Makefile. The file `python_interface/dpstokesCPU.py` contains an example. 
OpenMP is used for parallelization and users should control the number of threads via the bash script `cpuconfig.sh`.
