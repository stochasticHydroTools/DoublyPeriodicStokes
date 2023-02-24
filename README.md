## About this repository 

This repository implements a variant on the Force Coupling Method for Stokes suspensions, based on the exponential of a semicircle (ES) kernel, see:	
https://doi.org/10.48550/arXiv.2210.01837  
This repo also contains scripts and raw data to regenerate the majority of the figures (see figures/). The [examples on microroller suspensions](https://github.com/stochasticHydroTools/RigidMultiblobsWall/tree/New_FCM_Lubrication/Lubrication/Lubrication_Examples) are in the [RigidMultiblobWall branch New_FCM_Lubrication](https://github.com/stochasticHydroTools/RigidMultiblobsWall/tree/New_FCM_Lubrication/Lubrication/Lubrication_Examples) repo. We are working on fully integrating the FCM doubly periodic Stokes solver into that repo.

We provide a simple python interface, but the key performance-sensitive pieces (namely FFTs, BVP solvers, and spreading and interpolation) are implemented in C++/CUDA using FFTW/cuFFT and LAPACK. One can use either the GPU and CPU versions as needed without changing calls.

Using these codes requires carefully selecting the parameters of the ES kernel and the number of cells to achieve a required hydrodynamic radius for the particles (only monodisperse blobs are supported) with sufficient translational invariance, and to optimize the FFT performance. We provide [this python routine](https://github.com/stochasticHydroTools/DoublyPeriodicStokes/blob/main/source/cpu/python/GridAndKernelConfig.py) for this purpose. For GPUs, the choice of FFT-friendly grid sizes is GPU-specific and made in [this CUDA routine](https://github.com/RaulPPelaez/UAMMD/blob/v2.x/src/utils/Grid.cuh#L125-L176).

The doubly periodic Stokes solver for GPUs is provided by [UAMMD](https://github.com/RaulPPelaez/uammd), which is included as a submodule inside the `source/gpu` folder. **Make sure to clone this repo recursively** by doing:
```shell
git clone --recursive https://github.com/stochasticHydroTools/DoublyPeriodicStokes
```
This repository provides a joint CPU/GPU python interface to use the module. 
A CUDA/C++ version is available in UAMMD, see [here](https://uammd.readthedocs.io/en/latest/Integrator/BrownianHydrodynamics.html#doubly-periodic-stokes-dpstokes). The class DPStokesUAMMD in source/gpu/python_wrapper/uammd_wrapper.cu can also be used as an example on how to call the CUDA/C++ interface.

The CPU version of the solver is included in `source/cpu`. In particular, this implements an OpenMP-based C++ spreading and interpolation library in 3D that supports also non-uniform grids in the z direction. One can use the C++ library directly if desired, but there is a python interface as well.

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
Users should source the bash script `cpuconfig.sh` before using either 
the GPU or CPU Python interface in a new shell, and can edit the thread environment 
settings therin as needed. The PYTHONPATH and LD_LIBRARY_PATH environment variables
are appeneded to so that the modules can be used anywhere within the filesystem.
By default, the script will exist in the $INSTALL_DIR specified in the top Makefile.

If you want to use the Intel compiler for the CPU code, prefix the call to make as
```shell
CPU=Intel make
``` 
Note, even if using the Intel compiler, you must load the module for gcc-6.4.0 or higher, 
as the compiler relies on GNU headers. Also, note that by default with Intel compilers, the [MKL library](https://en.wikipedia.org/wiki/Math_Kernel_Library) is used to provide LAPACK/BLAS functionality. MKL (or other LAPACK/BLAS implementations) can be used with GNU compilers as well by specifying the relevant paths in the top level `Makefile`.
 
You can compile both CPU and GPU libraries in debug mode through
```shell
DEBUG=True make
```
Both CPU and DEBUG can also be set from within the `Makefile`, though the 
command line setting will override the ones in the `Makefile`.

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

The GPU interface requires [pybind11](https://github.com/pybind/pybind11) to compile, the default Makefile expects pip to be available and pybind11 being installed via "pip install pybind11". There is an option in the Makefile if neither pip nor pybind11 are available in the system.  
In order to use it you must compile the python wrappers using make in the root directory.  
A file called uammd.*.so will be created and then "import uammd" can be used inside python. The module name can be changed in the Makefile.
See `python_interface/dpstokesGPU.py` for a usage example. Once compiled and imported you can use "help(uammd)" in python for additional usage information.  

### CPU Python interface

See the `source/cpu/README.md` for details. Note, the build instructions contained therein are for using cmake3 as the build system. 
The section can be ignored, or followed analogously through the provided top level Makefile. The file `python_interface/dpstokesCPU.py` contains an example.
One can specify particles with differing radii in the CPU Python interface, though we only expose single radius setting in the joint interface. 
OpenMP is used for parallelization and users should control the number of threads via the bash script `cpuconfig.sh`.

