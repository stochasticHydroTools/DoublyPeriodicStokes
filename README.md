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

A python interface is provided that allows to compute the product between the mobility tensor (from a group of Gaussian sources) and a vector (i.e the force).  
This interface requires [pybind11](https://github.com/pybind/pybind11) to compile, which is included as a submodule and will be automatically downloaded if this repo is cloned recursively (see "About this repo" above).  
In order to use it you must compile the python wrappers using make (doing ```make python``` or ```make all``` here will also compile the python library).  
A file called uammd.*.so will be created and then "import uammd" can be used inside python. Notice that you might have to customize python\_interface/Makefile for your particular system.  
See python_interface/dpstokes.py for a usage example.  
Once compiled you can use "help(uammd)" for additional usage information.  

