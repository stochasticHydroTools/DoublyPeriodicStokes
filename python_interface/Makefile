# Donev: What is compile_flags.txt about?

UAMMD_ROOT=../source/gpu/uammd/
PYBIND_ROOT=pybind11/
# Donev: Is this PYTHON settings stuff shared between the CPU and GPU? Surely one has to use the same python for both to use the common interface?
PYTHON=python3
# Donev: Why do you need CUDA_ROOT? It doesn't work with "module load". I suggest instead
# NVCC=nvcc
# It is my understanding that when you use nvcc, -I$(CUDA_ROOT)/include is automatic, i.e., it knows where to link from?
# This is similar to how the mpi compiler scripts do the right thing for -I and -L
# It seems to me there is nothing here to control which MPI is used which suggests this is done just be making sure the right mpicc is in the PATH?
CUDA_ROOT=/usr/local/cuda
# Donev: I am confused here how the cpu version finds lapack? It seems to me LAPACKE_FLAGS is only used in make GPU
# I guess the two codes can use different versions of lapack?
LAPACKE_FLAGS=-llapack -llapacke -I/usr/include/lapacke
# Donev: What does commenting out versus not commenting out the next line do?
#MKL_FLAGS=-DUSE_MKL -DMKL_ILP64 -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64  -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_avx512 -lmkl_def -lpthread -lm -ldl

# Donev: It appears to me the rest of these are specific to UAMMD, right? If so, say this clearly and separate by empty lines or otherwise delineate
# Stuff like "DEBUG" could apply to either/both CPU and GPU versions
NVCC=$(CUDA_ROOT)/bin/nvcc
#UAMMD can be quite verbose, 5 shows only some messages at initialization/exit, 0 will only print critical errors, 15 will print A LOT.
# Donev: What does 1 do?
VERBOSITY=1
#Uncomment for double precision, UAMMD is compiled in single by default
#DOUBLEPRECISION=-DDOUBLE_PRECISION
#In case you prefer to import with other name
MODULE_NAME=uammd
INCLUDE_FLAGS= -I$(CUDA_ROOT)/include -I$(UAMMD_ROOT)/src -I$(UAMMD_ROOT)/src/third_party `$(PYTHON)-config --includes` -I $(PYBIND_ROOT)/include/

DEBUG= -O3 #-DUSE_NVTX #-DUAMMD_DEBUG #  -g  -G -O0

LIBRARY_NAME=$(MODULE_NAME)`$(PYTHON)-config --extension-suffix`
FILE=uammd_wrapper.cu

all: dpstokesGPU dpstokesCPU

# Donev: I don't like or don't get either of these Make rules. What you have written is a script not a Makefile
# A makefile should have dependencies and should only recompile stuff that needs to be recompiled if one changes a file
# Here, it looks like this will rm build directory for CPU every time you recompile, maybe even every time you run make?
# How does it decide when it needs to re-run the command without any dependence after the colon?
# Maybe the way to do this is to make part of this Makefile a new Makefile inside the GPU directory and then write a build.sh or some such SCRIPT
# instead of writing this as a "pretend" Makefile?

# Why is this not written as
# dpstokesGPU: $(FILE)
dpstokesGPU:
	$(NVCC) -w -shared -std=c++14 -DMAXLOGLEVEL=$(VERBOSITY) $(DOUBLEPRECISION) $(INCLUDE_FLAGS) $(DEBUG) -Xcompiler "-fPIC -w"  $(FILE) -o $(LIBRARY_NAME) $(LAPACKE_FLAGS) $(MKL_FLAGS) -lcufft -lnvToolsExt

# Somehow flags like -j6 to make seem wrong (hard-wired). I think there is some mechanism for "recursively calling make" that do this better
# Or, define a new variable "MAKE_FLAGS" inside the build.sh

dpstokesCPU:
# Donev: Who defines the variable "cpu"? What if someone has Intel but doesn't want to deal with MKL?
# Do you really mean "USE_MKL" here?
ifneq ($(cpu), Intel)
	cd ../source/cpu; rm -rf build; mkdir build; cd build; cmake3 -Dfftw_wisdom=on -Duse_stack=on ..; make -j6; make install; cd ../../../python_interface; cp ../source/cpu/examples/cpuconfig.sh .
else
# Sachin: I am a bit confused by fftw here. Is it required to recompile it with icc to link it? (should not be the case)
# Or is it about performance? I thought that Intel provides its own version of FFT libraries or FFTW interface
	cd ../source/cpu; rm -rf build; mkdir build; cd build;\
	if [ ! -d "../fftw_install" ]; then\
		wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.9.tar.gz;\
		tar xzf fftw-3.3.9.tar.gz;\
		cd fftw-3.3.9;\
		sed -i 's/fopenmp/qopenmp/g' configure;\
		CC=icc F77=ifort ./configure --prefix=$(clonedir)/source/cpu/fftw_install --enable-shared --enable-openmp --enable-sse2 --enable-avx --enable-avx2;\
		make -j6; make install; cd ..;\
	fi;\
	cmake3 -Dfftw_wisdom=on -Duse_stack=on -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc ..; make -j6; make install; cd ../../../python_interface; cp ../source/cpu/examples/cpuconfig.sh .
endif

clean:
	rm $(LIBRARY_NAME); cd ../source/cpu; rm -rf build lib
