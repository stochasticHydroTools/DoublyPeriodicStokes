UAMMD_ROOT=../source/gpu/uammd/
#This variable can be commented if the system provides pybind11
PYBIND_ROOT=pybind11/
#Location/name of the python3 executable
PYTHON3=python3

NVCC=nvcc
CUDA_ROOT="$(shell dirname `which nvcc`)"/../

LAPACKE_LIBS=-llapack -llapacke 
LAPACKE_FLAGS=-I/usr/include/lapacke 

#These two vars are proper to UAMMD, commenting LAPACKE_FLAGS and uncommenting MKL_FLAGS makes uammd use mkl instead of lapacke/openblas/etc.
#MKL_FLAGS=-DUSE_MKL -DMKL_ILP64 -m64 -I${MKLROOT}/include
#MKL_LIBS=-L${MKLROOT}/lib/intel64  -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_avx512 -lmkl_def -lpthread -lm -ldl

#UAMMD can be quite verbose, 5 shows only some messages at initialization/exit
# 0 will only print critical errors, 1 will print non crashing errors, 2 will print messages for recoverable exceptions, 3 will print warnings, 15 will print A LOT.
VERBOSITY=3
#Uncomment for double precision, UAMMD is compiled in single by default
#DOUBLEPRECISION=-DDOUBLE_PRECISION 
#In case you prefer to import with other name
GPU_MODULE_NAME=uammd
GPU_INCLUDE_FLAGS= -I$(CUDA_ROOT)/include -I$(UAMMD_ROOT)/src -I$(UAMMD_ROOT)/src/third_party

#DEBUG= -g -O0
ifdef $(DEBUG)
GPU_DEBUG= $(DEBUG) -G -DUSE_NVTX -DUAMMD_DEBUG
else
GPU_OPTIMIZATION= -O3
endif
#Raul: Sachin, can you make your CMake aware of the DEBUG variable?

GPU_LIBRARY_NAME=$(GPU_MODULE_NAME)$(shell $(PYTHON3)-config --extension-suffix)

all: dpstokesGPU dpstokesCPU

dpstokesGPU: $(GPU_LIBRARY_NAME)

python_gpu/uammd_gpu.o: python_gpu/uammd_wrapper.cu
	$(NVCC) -w -std=c++14 -DMAXLOGLEVEL=$(VERBOSITY) $(GPU_OPTIMIZATION) $(DOUBLEPRECISION) $(GPU_INCLUDE_FLAGS) $(GPU_DEBUG) -Xcompiler "-fPIC -w" -c $< -o $@  $(LAPACKE_FLAGS) $(MKL_FLAGS)

python_gpu/uammd_python.o: python_gpu/uammd_python.cpp
	$(CXX) -std=c++14 -O3 $(DEBUG) -fPIC -w `$(PYTHON3)-config --includes` -I $(PYBIND_ROOT)/include/ -c $< -o $@

$(GPU_LIBRARY_NAME): python_gpu/uammd_gpu.o python_gpu/uammd_python.o
	$(NVCC) $(GPU_OPTIMIZATION) $(GPU_DEBUG) -w -shared $^ -o $@ -lcufft -lnvToolsExt $(LAPACKE_LIBS) $(MKL_LIBS)

dpstokesCPU:

ifneq ($(cpu), Intel)
	cd ../source/cpu; rm -rf build; mkdir build; cd build; cmake3 -Dfftw_wisdom=on -Duse_stack=on ..; make -j6; make install; cd ../../../python_interface; cp ../source/cpu/examples/cpuconfig.sh .
else
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
	rm -f  python_gpu/uammd_gpu.o python_gpu/uammd_python.o $(GPU_LIBRARY_NAME); (cd ../source/cpu; rm -rf build lib) 
