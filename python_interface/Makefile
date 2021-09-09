#Raul: Donev, keep in mind that this Makefile was written just to compile uammd_wrapper.cu, Sachin then injected a call to his CMake file as an additional rule (not using any of the defined variables)

# Donev: What is compile_flags.txt about?
# A, Raul: an spureous file that got accidentally commited, removed.
UAMMD_ROOT=../source/gpu/uammd/
PYBIND_ROOT=pybind11/
#Location/name of the python3 executable
PYTHON3=python3
# Donev: Why do you need CUDA_ROOT? It doesn't work with "module load". I suggest instead
NVCC=nvcc
# It is my understanding that when you use nvcc, -I$(CUDA_ROOT)/include is automatic, i.e., it knows where to link from?
# A, Raul: This is true for the libraries, but I need to put -I$(CUDA_ROOT)/include, otherwise there are conflicts stemming from CUDA including a certain library (cub).
#          Some versions ago, CUDA did not include this library, so uammd has a local copy. If I remove the -I, I have to give up old versions of CUDA (8.0 and before AFAIK).
#          I now use $(which nvcc) to get CUDA_ROOT, this works for any installation directory hierarchy made by the official CUDA installer (including module load in principle). Additionally, if CUDA_ROOT discovery fails because of a non-standard installation, the Makefile should still just work.

# This is similar to how the mpi compiler scripts do the right thing for -I and -L
# It seems to me there is nothing here to control which MPI is used which suggests this is done just be making sure the right mpicc is in the PATH?
# A, Raul: I do not use MPI, I believe Sachin does it in his CMake file.
CUDA_ROOT="$(shell dirname `which nvcc`)"/../
# Donev: I am confused here how the cpu version finds lapack? It seems to me LAPACKE_FLAGS is only used in make GPU
# I guess the two codes can use different versions of lapack?
LAPACKE_FLAGS=-llapack -llapacke -I/usr/include/lapacke
# Donev: What does commenting out versus not commenting out the next line do?
#MKL_FLAGS=-DUSE_MKL -DMKL_ILP64 -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64  -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_avx512 -lmkl_def -lpthread -lm -ldl
# A, Raul: These two vars are proper to UAMMD, commenting LAPACKE_FLAGS and uncommenting MKL_FLAGS makes uammd use mkl instead of lapacke/openblas/etc.
#          The CPU version discovers lapacke on its own via Sachin's CMake, AFAIK Sachin's will look openblas.
#          I guess, given that CMake is a more sophisticated build system, I should port all uammd compilation into Sachin's CMakeFiles.txt.

# Donev: It appears to me the rest of these are specific to UAMMD, right? If so, say this clearly and separate by empty lines or otherwise delineate
# Stuff like "DEBUG" could apply to either/both CPU and GPU versions
# A, Raul: All defined variables in this Makefile are just for UAMMD, Sachin does not use them.
#	   Regarding debug, I separated it between GPU and CPU variables (the flags are different)

#UAMMD can be quite verbose, 5 shows only some messages at initialization/exit
# 0 will only print critical errors, 1 will print non crashing errors, 2 will print messages for recoverable exceptions, 3 will print warnings, 15 will print A LOT.
VERBOSITY=3 #Raul: It would be useful that Sachin shares this variable
#Uncomment for double precision, UAMMD is compiled in single by default
#Raul: I assume it is meaningless for the CPU impl to be compiled in single, should this variable be uammd only then?
#DOUBLEPRECISION=-DDOUBLE_PRECISION 
#In case you prefer to import with other name
GPU_MODULE_NAME=uammd
GPU_INCLUDE_FLAGS= -I$(CUDA_ROOT)/include -I$(UAMMD_ROOT)/src -I$(UAMMD_ROOT)/src/third_party `$(PYTHON3)-config --includes` -I $(PYBIND_ROOT)/include/

#DEBUG= -g -O0
ifdef $(DEBUG)
GPU_DEBUG= $(DEBUG) -G -DUSE_NVTX -DUAMMD_DEBUG
else
GPU_OPTIMIZATION= -O3
endif
#Raul: Sachin, can you make your CMake aware of the DEBUG variable?

GPU_LIBRARY_NAME=$(GPU_MODULE_NAME)`$(PYTHON)-config --extension-suffix`

all: dpstokesGPU dpstokesCPU

# Donev: I don't like or don't get either of these Make rules. What you have written is a script not a Makefile
# A makefile should have dependencies and should only recompile stuff that needs to be recompiled if one changes a file
# Here, it looks like this will rm build directory for CPU every time you recompile, maybe even every time you run make?
# How does it decide when it needs to re-run the command without any dependence after the colon?
# A, Raul: It recompiles everytime

# Maybe the way to do this is to make part of this Makefile a new Makefile inside the GPU directory and then write a build.sh or some such SCRIPT
# instead of writing this as a "pretend" Makefile?
# A, Raul: Sachin uses CMake, so for him using a Makefile is awkward, the only way he can do it is like he is doing.
#	   If you are not cool with a pretend Makefile, my suggestion is to put uammd compilation inside a CMakeList.txt here.
#	   This CMakeList.txt would compile uammd and also include Sachin's cmake via "add_subdirectory". My CMake foo is limited, though.
#	   Alternatively, this Makefile could simply have two lines, calling Makefiles under source/gpu and sources/cpu.

# Why is this not written as
# dpstokesGPU: $(FILE)
#A, Raul: It is now.
dpstokesGPU: uammd_wrapper.cu
	$(NVCC) -w -shared -std=c++14 -DMAXLOGLEVEL=$(VERBOSITY) $(GPU_OPTIMIZATION) $(DOUBLEPRECISION) $(GPU_INCLUDE_FLAGS) $(GPU_DEBUG) -Xcompiler "-fPIC -w"  $< -o $(GPU_LIBRARY_NAME) $(LAPACKE_FLAGS) $(MKL_FLAGS) -lcufft -lnvToolsExt

# Somehow flags like -j6 to make seem wrong (hard-wired). I think there is some mechanism for "recursively calling make" that do this better
# Or, define a new variable "MAKE_FLAGS" inside the build.sh

dpstokesCPU:
# Donev: Who defines the variable "cpu"? What if someone has Intel but doesn't want to deal with MKL?
# Do you really mean "USE_MKL" here?
# A, Raul: USE_MKL is uammd proper, Sachin does not use mkl AFAIK. I Sachin intended "cpu" to be defined at calling: "make -Dcpu=Intel"

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
# Raul: Sachin, I think it would make sense to put the "build" and "lib" folders here, in the folder where "make" was run. What do you think?
endif

clean:
	rm $(GPU_LIBRARY_NAME); cd ../source/cpu; rm -rf build lib
