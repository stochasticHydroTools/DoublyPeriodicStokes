## About this repository 

This repo contains code and scripts to reproduce the data for the figures in the applications section of the article [CITE THE ARTICLE HERE].  

See the README inside each folder for information about each simulation  

The doubly periodic Stokes solver for GPUs is provided by [UAMMD](https://github.com/RaulPPelaez/uammd), which is included as a submodule inside the `source/gpu` folder. **Make sure to clone this repo recursively.**  
This means that you must clone using  
```shell
git clone --recursive https://github.com/stochasticHydroTools/DPStokesTests
```
At the moment, only a python interface to use the module is provided.

The CPU version of the solver is included in `source/cpu`.

## USAGE:  

Running make will compile and set up both the CPU and GPU python interfaces for the solver.
You may have to adapt `compile_flags.txt`, `python_interface/compile_flags.txt`, and the `python_interface/Makefile` to build
the GPU solver. The CPU solver should compile as is, but there is more information for customizing the build
in the source/cpu/README.md. 

If you want to use the Intel compiler for the CPU code, prefix the call to make as
```shell
CPU=Intel make
```  

NOTE: Not useful yet, since reproduce scripts are not included.
  
### GPU Python interface

A python interface is provided that allows to compute the hydrodynamic displacements for a group of positions with forces and/or torques acting on them in different geometries, mainly:  

	* Doubly periodic (either with no walls, a bottom wall or a slit channel)  
	* Triply periodic (using force coupling method)  

The GPU interface requires [pybind11](https://github.com/pybind/pybind11) to compile, which is included as a submodule and will be automatically downloaded if this repo is cloned recursively (see "About this repo" above).  
In order to use it you must compile the python wrappers using make inside python_interface.  
A file called uammd.*.so will be created and then "import uammd" can be used inside python. 
See `python_interface/dpstokesGPU.py` for a usage example. Once compiled and imported you can use "help(uammd)" in python for additional usage information.  

### CPU Python interface

See the `source/cpu/README.md` for details. Importantly, a file `cpuconfig.sh` will be generated in the `python_interface` folder. 
Users must source this bash script before using the CPU Python interface in a new shell, and can edit the thread environment 
settings therein as needed.

## About the spreading/interpolation kernels in the python interface

The modules will use the ES kernel (called BM in UAMMD).  

Other kernels can be used, but the interface does not allow them for the moment.  

Hydrodynamic displacements coming from forces and torques can be computed. 
For the GPU interface, if the torque-related arguments are ommited, the computations related to them are skipped entirely.
For the CPU interface, the user must specify whether torques are involved with a boolean swith, like `has_torque=True`.

The files python_interface/dpstokesGPU.py and python_interface/dpstokesCPU.py contain more info.  
