#Print an error if this Makefile is called directly
ifndef calling_from_parent
  errormsg:=$(shell printf '\033[31m\033[1mThis Makefile must be called by the Makefile in the root directory. Do not try to run it directly\033[0m\n')
  $(error $(errormsg))
endif
GPU_INCLUDE_FLAGS   = -I$(CUDA_ROOT)/include -I$(UAMMD_ROOT)/src -I$(UAMMD_ROOT)/src/third_party
DEBUG_FLAGS         = -g -O0
ifeq ($(DEBUG), True)
  GPU_DEBUG = $(DEBUG_FLAGS) -G -DUSE_NVTX -DUAMMD_DEBUG
  export CPU_DEBUG  = $(DEBUG_FLAGS) -DDEBUG 
else
GPU_OPTIMIZATION    = -O3
CPU_OPTIMIZATION    = -O3
endif

GPU_LIBRARY_NAME    = $(INSTALL_DIR)/$(GPU_MODULE_NAME)$(shell $(PYTHON3)-config --extension-suffix)

all: dpstokesCPU dpstokesGPU

dpstokesGPU: $(GPU_LIBRARY_NAME)

python_gpu/uammd_gpu.o: python_gpu/uammd_wrapper.cu
	$(NVCC) -w -std=c++14 -DMAXLOGLEVEL=$(VERBOSITY) $(GPU_OPTIMIZATION) $(DOUBLEPRECISION) $(GPU_INCLUDE_FLAGS) $(GPU_DEBUG) -Xcompiler "-fPIC -w -fno-lto" -c $< -o $@  $(LAPACKE_FLAGS)

python_gpu/uammd_python.o: python_gpu/uammd_python.cpp
	$(CXX) -std=c++14 -O3 $(CPU_DEBUG) $(CPU_OPTIMIZATION) -flto -fPIC -w `$(PYTHON3)-config --includes` -I $(PYBIND_ROOT)/include/ -c $< -o $@

$(GPU_LIBRARY_NAME): python_gpu/uammd_gpu.o python_gpu/uammd_python.o
	mkdir -p $(INSTALL_DIR)
	$(NVCC) $(GPU_OPTIMIZATION) $(GPU_DEBUG) -w -shared $^ -o $@ -lcufft -lnvToolsExt $(LAPACKE_LIBS) $(MKL_LIBS)

dpstokesCPU:

ifeq ($(CPU), GNU)
	make -f Makefile.GNU -C ../source/cpu; 
endif
ifeq ($(CPU), Intel)
	make -f Makefile.Intel -C ../source/cpu;
endif

clean: clean_cpu clean_gpu

clean_cpu:
	rm -rf $(INSTALL_DIR)/lib*.so

clean_gpu:
	rm -rf python_gpu/uammd_gpu.o python_gpu/uammd_python.o $(INSTALL_DIR)/uammd*.so
