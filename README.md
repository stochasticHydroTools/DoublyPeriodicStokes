## About this repository 

This repo contains code and scripts to reproduce the data for the figures in the applications section of the article [CITE THE ARTICLE HERE].  

See README inside each folder for information about each simulation  

The doubly periodic Stokes solver is provided by [UAMMD](https://github.com/RaulPPelaez/uammd), which is included as a submodule inside the source folder. **Make sure to clone this repo recursively.**  
This means that you must clone using  
```shell
git clone --recursive https://github.com/stochasticHydroTools/DPStokesTests
```
At the moment, only a python interface to use the module is provided.

## USAGE:  

Run make to compile, you might have to adapt it to your particular system before.  
  
### Python interface

A python interface is provided that allows to compute the hydrodynamic displacements for a group of positions with forces and torques acting on them in different geometries, mainly:  

	* Doubly periodic (either with no walls, a bottom wall or a slit channel)  
	* Triply periodic (using force coupling method)  

Note that for now, FCM cannot deal with torques, just forces.

This interface requires [pybind11](https://github.com/pybind/pybind11) to compile, which is included as a submodule and will be automatically downloaded if this repo is cloned recursively (see "About this repo" above).  
In order to use it you must compile the python wrappers using make (doing ```make python``` or ```make all``` here will also compile the python library).  
A file called uammd.*.so will be created and then "import uammd" can be used inside python. Notice that you might have to customize python\_interface/Makefile for your particular system.  
See python_interface/dpstokes.py for a usage example.  
Once compiled you can use "help(uammd)" for additional usage information.  

## About the spreading/interpolation kernels in the python interface

The doubly periodic options will use the ES kernel (called BM in UAMMD), in this case the tolerance parameter will be ignored and the support must be specified.  
The gw parameter will always be ignored (it comes from when a Gaussian was used).  

The triply periodic mode (FCM) will use a Gaussian, the tolerance parameter can be specified instead of the support (which will be autocomputed from it).  

Other kernels can be used, but the interface does not allow them for the moment, since the specific kernel parameters are still a work in progress.  

The file dpstokes.py contains more info about this.
