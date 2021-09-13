################ USAGE #########################

# [CPU=Intel] [DEBUG=True] make           

# i.e) make, CPU=Intel make, DEBUG=True make, CPU=Intel DEBUG=True make
# CPU=Intel will toggle the Intel compiler (icpc) for cpu code
# DEBUG=True will compile in debug mode for both gpu and cpu code
# These vars can be specified below as well.

# Users shoud edit the USER EDIT section as needed

########################## BEGIN USER EDIT ##############################
# specify compiler type - GNU|Intel (only used for cpu build)
export CPU             ?= Intel
# specify location of DoublyPeriodicStokes directory
export SRC_DIR          = /home/srn324/DoublyPeriodicStokes
# specify desired location of shared libraries
export INSTALL_DIR      = $(SRC_DIR)/python_interface/lib
# name of the python3 executable
export PYTHON3          = python3
# cuda compiler and bin 
export NVCC             = nvcc
export CUDA_ROOT        = /usr/local/stow/cuda-10.2/bin
# specify where lapacke.h and liblapacke.so 
#   openblas
export LAPACKE_FLAGS    = -I/usr/include/openblas -L/usr/lib64
export LAPACKE_LIBS     = -lopenblas -llapacke
#   mkl
#export LAPACKE_FLAGS   = -I/opt/intel/mkl/include -L/opt/intel/mkl/lib/intel64 -DUSE_MKL
#export LAPACKE_LIBS    = -lmkl_rt -lm -ldl

# specify where fftw*.h/.so are and select whether to use fftw wisdom (see README)
ifneq ($(CPU), Intel)
	# change USE_FFTW_MEASURE to USE_FFTW_PATIENT for more optimal fftw plans
  export FFTW_FLAGS     = -I/usr/include -L/usr/lib64 -DENABLE_WISDOM -DUSE_FFTW_MEASURE
endif
ifeq ($(CPU),Intel)
  # set fftw_install dir (custom install needed for intel cc)
  export FFTW_INSTALL   = $(SRC_DIR)/source/cpu/fftw_install
  # set mkl install dir
  export MKLROOT       ?= /opt/intel/mkl
  export FFTW_FLAGS     = -I$(FFTW_INSTALL)/include -L$(FFTW_INSTALL)/lib -DENABLE_WISDOM -DUSE_FFTW_MEASURE
	# override lapack flags/libs to always use mkl
  export LAPACKE_FLAGS  = -I$(MKLROOT)/include -L$(MKLROOT)/lib/intel64 -DUSE_MKL
  export LAPACKE_LIBS   = -lmkl_rt -lpthread -ldl
endif

# use stack instead of heap for part of cpu spreading algorithm
export SPREAD_FLAGS     = -DUSE_STACK

# uncomment for double precision - UAMMD is compiled in single by default
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

# if True, compilation will use debug mode for cpu and gpu code
export DEBUG           ?= False

# root of uammd
export UAMMD_ROOT       = $(SRC_DIR)/source/gpu/uammd

# This variable can be commented if the system provides pybind11
export PYBIND_ROOT      = $(SRC_DIR)/python_interface/pybind11

################################ END USER EDIT ##################################

export calling_from_parent = True
all: python

python: python_cpu python_gpu

python_cpu:
	make dpstokesCPU -C $(SRC_DIR)/python_interface
	@sed -i "/SRC_DIR=/c SRC_DIR=$(SRC_DIR)" $(SRC_DIR)/python_interface/cpuconfig.sh
	@sed -i "/INSTALL_DIR=/c INSTALL_DIR=$(INSTALL_DIR)" $(SRC_DIR)/python_interface/cpuconfig.sh
	@sed -i "/CPU=/c CPU=$(CPU)" $(SRC_DIR)/python_interface/cpuconfig.sh

python_gpu:
	make dpstokesGPU -C $(SRC_DIR)/python_interface
	@sed -i "/SRC_DIR=/c SRC_DIR=$(SRC_DIR)" $(SRC_DIR)/python_interface/cpuconfig.sh
	@sed -i "/INSTALL_DIR=/c INSTALL_DIR=$(INSTALL_DIR)" $(SRC_DIR)/python_interface/cpuconfig.sh
	@sed -i "/CPU=/c CPU=$(CPU)" $(SRC_DIR)/python_interface/cpuconfig.sh

clean:
	make clean -C $(SRC_DIR)/python_interface

clean_cpu:
	make clean_cpu -C $(SRC_DIR)/python_interface

clean_gpu:
	make clean_gpu -C $(SRC_DIR)/python_interface
