UAMMD_ROOT=../source/gpu/uammd/
PYBIND_ROOT=pybind11/
PYTHON=python3
CUDA_ROOT=/usr/local/cuda
LAPACKE_FLAGS=-llapack -llapacke -I/usr/include/lapacke
#MKL_FLAGS=-DUSE_MKL -DMKL_ILP64 -m64 -I${MKLROOT}/include -L${MKLROOT}/lib/intel64  -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lmkl_avx512 -lmkl_def -lpthread -lm -ldl
NVCC=$(CUDA_ROOT)/bin/nvcc
#UAMMD can be quite verbose, 5 shows only some messages at initialization/exit, 0 will only print critical errors, 15will print A LOT.
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

dpstokesGPU:
	$(NVCC) -w -shared -std=c++14 -DMAXLOGLEVEL=$(VERBOSITY) $(DOUBLEPRECISION) $(INCLUDE_FLAGS) $(DEBUG) -Xcompiler "-fPIC -w"  $(FILE) -o $(LIBRARY_NAME) $(LAPACKE_FLAGS) $(MKL_FLAGS) -lcufft -lnvToolsExt

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
	rm $(LIBRARY_NAME); cd ../source/cpu; rm -rf build lib
