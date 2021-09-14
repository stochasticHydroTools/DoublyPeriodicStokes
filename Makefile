################ USAGE #########################

# [CPU=Intel] [DEBUG=True] make           

# i.e) make, CPU=Intel make, DEBUG=True make, CPU=Intel DEBUG=True make
# CPU=Intel will toggle the Intel compiler (icpc) for cpu code
# DEBUG=True will compile in debug mode for both gpu and cpu code
# These vars can be specified below as well.

# Users shoud edit the USER EDIT section as needed

########################## BEGIN USER EDIT ##############################
# specify location of DoublyPeriodicStokes directory
export DPSTOKES_ROOT    = $(PWD)
# specify desired location of shared libraries
export DPSTOKES_INSTALL = $(DPSTOKES_ROOT)/python_interface/lib
# name of the python3 executable
export PYTHON3          = python3
# if True, compilation will use debug mode for cpu and gpu code
export DEBUG           ?= False

# Donev split this into two sections for improved organization (but you can use #### if you want):

#---------------------------------------------------
# CPU/OpenMP settings
#---------------------------------------------------
# specify compiler type - GNU|Intel (only used for cpu build)
export CPU             ?= Intel

# specify where lapacke.h and liblapacke.so 
# For openblas (see CPU=Intel below for Intel's MKL)
export LAPACKE_FLAGS    = -I/usr/include/openblas -L/usr/lib64
export LAPACKE_LIBS     = -lopenblas -llapacke
# Donev removed MKL stuff here since it is below (I don't like commented out sections in general)

# specify where fftw*.h/.so are and select whether to use fftw wisdom (see README)
ifeq ($(CPU),Intel)
  # set fftw_install dir (custom install needed for intel cc)
  export FFTW_INSTALL   = $(DPSTOKES_ROOT)/source/cpu/fftw_install
  # set mkl install dir
  export MKLROOT       ?= /opt/intel/mkl
  export FFTW_FLAGS     = -I$(FFTW_INSTALL)/include -L$(FFTW_INSTALL)/lib -DENABLE_WISDOM -DUSE_FFTW_MEASURE
  # override lapack flags/libs to always use MKL for best performance
  export LAPACKE_FLAGS  = -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -DUSE_MKL
  export LAPACKE_LIBS   = -lmkl_rt -lpthread -ldl
else # Donev changed this to else
  # change USE_FFTW_MEASURE to USE_FFTW_PATIENT for more optimal fftw plans
  export FFTW_FLAGS     = -I/usr/include -L/usr/lib64 -DENABLE_WISDOM -DUSE_FFTW_MEASURE  
endif

# use stack instead of heap for part of cpu spreading algorithm
export SPREAD_FLAGS     = -DUSE_STACK

#---------------------------------------------------
# GPU settings
#---------------------------------------------------
# cuda compiler and bin 
export NVCC             = nvcc
export CUDA_ROOT        = "$(shell dirname `which $(NVCC)`)"/..

# uncomment for double precision - UAMMD is compiled in single by default
# Donev for Raul: This doesn't really work at present unless you change self.precision = np.float32 in the common interface
#export DOUBLEPRECISION = -DDOUBLE_PRECISION 

# UAMMD can be quite verbose, 5 shows only some messages at initialization/exit
# 0 will only print critical errors, 1 will print non crashing errors, 
# 2 will print messages for recoverable exceptions,
# 3 will print warnings,
# 5 will print messages (initialization, parameters, ...)
# 15 will print A LOT.
export VERBOSITY        = 5

# name of gpu module
export GPU_MODULE_NAME  = uammd

# root of uammd (change if compiled separately)
export UAMMD_ROOT       = $(DPSTOKES_ROOT)/source/gpu/uammd

# This variable can be commented if the system provides pybind11
export PYBIND_ROOT      = $(DPSTOKES_ROOT)/python_interface/pybind11

################################ END USER EDIT ##################################

export calling_from_parent = True
all: python

python: python_cpu python_gpu

python_cpu:
	make dpstokesCPU -C $(DPSTOKES_ROOT)/python_interface
	@sed -i "/DPSTOKES_ROOT=/c DPSTOKES_ROOT=$(DPSTOKES_ROOT)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh
	@sed -i "/DPSTOKES_INSTALL=/c DPSTOKES_INSTALL=$(DPSTOKES_INSTALL)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh
	@sed -i "/CPU=/c CPU=$(CPU)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh

python_gpu:
	make dpstokesGPU -C $(DPSTOKES_ROOT)/python_interface
	@sed -i "/DPSTOKES_ROOT=/c DPSTOKES_ROOT=$(DPSTOKES_ROOT)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh
	@sed -i "/DPSTOKES_INSTALL=/c DPSTOKES_INSTALL=$(DPSTOKES_INSTALL)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh
	@sed -i "/CPU=/c CPU=$(CPU)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh

clean:
	make clean -C $(DPSTOKES_ROOT)/python_interface

clean_cpu:
	make clean_cpu -C $(DPSTOKES_ROOT)/python_interface

clean_gpu:
	make clean_gpu -C $(DPSTOKES_ROOT)/python_interface
