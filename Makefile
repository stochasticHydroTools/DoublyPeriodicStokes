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

#---------------------------------------------------
# CPU/OpenMP settings
#---------------------------------------------------
# specify compiler type - GNU|Intel (only used for cpu build)
export CPU             ?= GNU

# specify where lapacke.h and liblapacke.so 
# For openblas (see CPU=Intel below for Intel's MKL)
export LAPACKE_FLAGS    = -I/usr/include/openblas -L/usr/lib64
export LAPACKE_LIBS     = -lopenblas -llapacke

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
else 
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

# Uncomment for double precision - UAMMD is compiled in single by default
#export DOUBLEPRECISION = -DDOUBLE_PRECISION 

# UAMMD can be quite verbose, 5 shows only some messages at initialization/exit
# 0 will only print critical errors, 1 will print non crashing errors, 
# 2 will print messages for recoverable exceptions,
# 3 will print warnings,
# 5 will print messages (initialization, parameters, ...)
# 15 will print A LOT.
export VERBOSITY        = 5

# name of gpu module (changing this will require to change the common_wrapper*py code accordingly)
export GPU_MODULE_NAME  = uammd

# root of uammd (uammd should be automatically downloaded, if you want to use a local version, change this variable accordingly). Note that UAMMD does not have to be compiled on its own (it is a header only library).
export UAMMD_ROOT       = $(DPSTOKES_ROOT)/source/gpu/uammd

#Location of pybind11 root directory.
#If pybind11 was installed via "pip install pybind11" this snippet will probably work out of the box.
#export PYBIND_ROOT      := $(shell pip list -v 2>&1 | grep '^pybind11[[:space:]]' | awk '{print $$3"/pybind11"}')
#If you do not pip available you can comment the above line and uncomment the two lines below to auto download it:
DOWNLOAD_PYBIND11 := 1
export PYBIND_ROOT = $(DPSTOKES_ROOT)/source/gpu/python_wrapper/pybind11
################################ END USER EDIT ##################################

export calling_from_parent = True
all: python

python: python_cpu python_gpu

python_cpu: envconfig
ifeq ($(CPU), GNU)
	make -f Makefile.GNU -C $(DPSTOKES_ROOT)/source/cpu; 
endif
ifeq ($(CPU), Intel)
	make -f Makefile.Intel -C $(DPSTOKES_ROOT)/source/cpu;
endif

python_gpu: envconfig $(DOWNLOAD_PYBIND11)
	make dpstokesGPU -C $(DPSTOKES_ROOT)/source/gpu/python_wrapper

$(DOWNLOAD_PYBIND11): $(PYBIND_ROOT)

$(PYBIND_ROOT):
	git clone --depth=1 https://github.com/pybind/pybind11 $@

envconfig:
	@sed -i "/DPSTOKES_ROOT=/c DPSTOKES_ROOT=$(DPSTOKES_ROOT)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh
	@sed -i "/DPSTOKES_INSTALL=/c DPSTOKES_INSTALL=$(DPSTOKES_INSTALL)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh
	@sed -i "/CPU=/c CPU=$(CPU)" $(DPSTOKES_ROOT)/python_interface/cpuconfig.sh

clean: clean_cpu clean_gpu

clean_cpu:
	rm -rf $(DPSTOKES_INSTALL)/lib*.so

clean_gpu:
	rm -rf $(DPSTOKES_ROOT)/source/gpu/python_wrapper/*.o $(DPSTOKES_INSTALL)/uammd*.so
