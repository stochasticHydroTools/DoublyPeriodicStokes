################ USAGE #########################

# [CPU=Intel] [DEBUG=True] make           

# i.e) make, CPU=Intel make, DEBUG=True, CPU=Intel DEBUG=True make
# CPU=Intel will toggle the Intel compiler (icpc) for cpu code
# DEBUG=True will compile in debug mode for both gpu and cpu code

# Users shoud edit the USER EDIT section as needed

########################## BEGIN USER EDIT ##############################
# specify location of shared libraries
export INSTALL_DIR      = $(shell pwd)/python_interface/lib
# root of uammd
export UAMMD_ROOT       = source/gpu/uammd/
# This variable can be commented if the system provides pybind11
export PYBIND_ROOT      = python_interface/pybind11/
# Location/name of the python3 executable
export PYTHON3          = python3
# cuda compiler and bin 
export NVCC             = nvcc
export CUDA_ROOT        = "$(shell dirname `which nvcc`)"/../

# specify where lapacke.h and lapacke.so 
#   openblas
export LAPACKE_FLAGS    = -I/usr/include/openblas -L/usr/lib64
export LAPACKE_LIBS     = -lopenblas -llapacke
#   mkl
#export LAPACKE_FLAGS   = -I/opt/intel/mkl/include -L/opt/intel/mkl/lib/intel64 -DUSE_MKL
#export LAPACKE_LIBS    = -lmkl_rt -lm -ldl

# specify where fftw*.h/.so are and select whether to use fftw wisdom (see README)
# change USE_FFTW_MEASURE to USE_FFTW_PATIENT for more optimal fftw plans
ifneq ($(cpu), Intel)
  export FFTW_FLAGS     = -I/usr/include -L/usr/lib64 -DENABLE_WISDOM -DUSE_FFTW_MEASURE
endif
# use stack instead of heap for part of cpu spreading algorithm
export SPREAD_FLAGS     = -DUSE_STACK

# uncomment for double precision - UAMMD is compiled in single by default
#export DOUBLEPRECISION = -DDOUBLE_PRECISION 

# UAMMD can be quite verbose, 5 shows only some messages at initialization/exit
# 0 will only print critical errors, 1 will print non crashing errors, 
# 2 will print messages for recoverable exceptions, 
# 3 will print warnings, 15 will print A LOT.
export VERBOSITY        = 3
# name of gpu module
export GPU_MODULE_NAME  = uammd

# if you want to use the Intel compiler, set the following
ifeq ($(cpu),Intel)
  # set fftw_install dir
  export FFTW_INSTALL   = $(PWD)/source/cpu/fftw_install
  # set mkl install dir
  export MKLROOT        = /opt/intel/mkl/include
  # change USE_FFTW_MEASURE to USE_FFTW_PATIENT for more optimal fftw plans
  export FFTW_FLAGS     = -I$(FFTW_INSTALL)/include -L$(FFTW_INSTALL)/lib -DENABLE_WISDOM -DUSE_FFTW_MEASURE
  export LAPACKE_FLAGS  = -I$(MKLROOT)/include -L$(MKLROOT)/intel64 -DUSE_MKL
  export LAPACKE_LIBS   = -lmkl_rt -lpthread -ldl 
endif
################################ END USER EDIT ##################################

# 'CPU=Intel make' will use intel compiler for cpu code
export cpu      = $(CPU)
# 'DEBUG=True make' will compile in debug mode for cpu and gpu code
export debug    = $(DEBUG)
export clonedir = $(PWD)

all: python

python: python_cpu python_gpu

python_cpu:
  make dpstokesCPU -C python_interface

python_gpu:
  make dpstokesGPU -C python_interface

clean:
  make -C python_interface clean

clean_cpu:
  make clean_cpu -C python_interface

clean_gpu:
  make clean_gpu -C python_interface
