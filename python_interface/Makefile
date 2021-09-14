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

GPU_LIBRARY_NAME    = $(DPSTOKES_INSTALL)/$(GPU_MODULE_NAME)$(shell $(PYTHON3)-config --extension-suffix)

PYTHON_WRAPPER_DIR :=$(UAMMD_ROOT)/../python_wrapper
all: dpstokesCPU dpstokesGPU

dpstokesGPU: $(GPU_LIBRARY_NAME)

uammd_gpu.o: $(PYTHON_WRAPPER_DIR)/uammd_wrapper.cu
	$(NVCC) -w -std=c++14 -DMAXLOGLEVEL=$(VERBOSITY) $(GPU_OPTIMIZATION) $(DOUBLEPRECISION) $(GPU_INCLUDE_FLAGS) $(GPU_DEBUG) -Xcompiler "-fPIC -w" -c $< -o $@  $(LAPACKE_FLAGS)

uammd_python.o: $(PYTHON_WRAPPER_DIR)/uammd_python.cpp
	$(CXX) -std=c++14 -O3 $(CPU_DEBUG) $(CPU_OPTIMIZATION) -fPIC -w `$(PYTHON3)-config --includes` -I $(PYBIND_ROOT)/include/ -c $< -o $@

$(GPU_LIBRARY_NAME): uammd_gpu.o uammd_python.o
	mkdir -p $(DPSTOKES_INSTALL)
	$(NVCC) $(GPU_OPTIMIZATION) $(GPU_DEBUG) -w -shared $^ -o $@ -lcufft -lnvToolsExt $(LAPACKE_LIBS) $(MKL_LIBS)

dpstokesCPU:

ifeq ($(CPU), GNU)
	make -f Makefile.GNU -C $(DPSTOKES_ROOT)/source/cpu; 
endif
ifeq ($(CPU), Intel)
	make -f Makefile.Intel -C $(DPSTOKES_ROOT)/source/cpu;
endif

clean: clean_gpu

clean_cpu:
	rm -rf $(DPSTOKES_INSTALL)/lib*.so

clean_gpu:
	rm -rf uammd_gpu.o uammd_python.o $(DPSTOKES_INSTALL)/uammd*.so
