ifndef calling_from_parent
  errormsg:=$(shell printf '\033[31m\033[1mThis Makefile must be called by the Makefile in the root directory. Do not try to run it directly\033[0m\n')
  $(error $(errormsg))
endif

CURDIR          = $(shell pwd)
CXX             = g++ 

ifneq ($(DEBUG), True)
  CXXFLAGS      = -I$(CURDIR)/include -w -O3 -march=native -shared -L$(DPSTOKES_INSTALL)
else
  CXXFLAGS      = -I$(CURDIR)/include -w -shared -g -O0 -DDEBUG -L$(DPSTOKES_INSTALL)
endif

chebSRC         = src/Quadrature.cpp
chebINC         = include/Quadrature.h
gridSRC         = src/Grid.cpp wrapper/GridWrapper.cpp
gridINC         = include/Grid.h include/common.h include/exceptions.h include/Quadrature.h 
particlesSRC    = src/ParticleList.cpp wrapper/ParticleListWrapper.cpp
particlesINC    = include/ParticleList.h include/ESKernels.h include/Grid.h include/Quadrature.h include/exceptions.h
spreadInterpSRC = src/SpreadInterp.cpp wrapper/SpreadInterpWrapper.cpp
spreadInterpINC = include/SpreadInterp.h include/Grid.h include/ParticleList.h include/exceptions.h
transformSRC    = src/Transform.cpp wrapper/TransformWrapper.cpp
transformINC    = include/Transform.h include/exceptions.h include/common.h
linSolveSRC     = src/LinearSolvers.cpp
linSolveINC     = include/LinearSolvers.h
dpToolsSRC      = src/DPTools.cpp
dpToolsINC      = include/DPTools.h include/exceptions.h
bcSRC           = wrapper/BCWrapper.cpp
bcINC           = include/BoundaryConditions.h include/common.h include/Grid.h include/ParticleList.h
LIBS_           = libcheb.so libgrid.so libparticles.so libspreadInterp.so libspreadInterpDerivX.so\
                  libspreadInterpDerivY.so libspreadInterpDerivZ.so libtransform.so liblinSolve.so\
                  libdpTools.so libBC.so

LIBS            = $(patsubst %,$(DPSTOKES_INSTALL)/%,$(LIBS_))

all: $(LIBS)

$(DPSTOKES_INSTALL)/libcheb.so: $(chebSRC) $(chebINC)
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(chebSRC) $(CXXFLAGS) -fPIC -fopenmp

$(DPSTOKES_INSTALL)/libgrid.so: $(gridSRC) $(gridINC) 
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(gridSRC) $(CXXFLAGS) $(FFTW_FLAGS) -lfftw3 -lm -fPIC -fopenmp

$(DPSTOKES_INSTALL)/libparticles.so: $(particlesSRC) $(particlesINC) $(DPSTOKES_INSTALL)/libgrid.so $(DPSTOKES_INSTALL)/libcheb.so
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(particlesSRC) $(CXXFLAGS) -lfftw3 -fPIC -fopenmp -lgrid -lcheb

$(DPSTOKES_INSTALL)/libspreadInterp.so: $(spreadInterpSRC) $(spreadInterpINC) $(DPSTOKES_INSTALL)/libparticles.so
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(spreadInterpSRC) $(CXXFLAGS) $(SPREAD_FLAGS) -lfftw3 -fPIC -fopenmp -lparticles

$(DPSTOKES_INSTALL)/libspreadInterpDerivX.so: $(spreadInterpSRC) $(spreadInterpINC) $(DPSTOKES_INSTALL)/libparticles.so
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(spreadInterpSRC) $(CXXFLAGS) $(SPREAD_FLAGS) -lfftw3 -DTORQUE_X -fPIC -fopenmp -lparticles
  
$(DPSTOKES_INSTALL)/libspreadInterpDerivY.so: $(spreadInterpSRC) $(spreadInterpINC) $(DPSTOKES_INSTALL)/libparticles.so
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(spreadInterpSRC) $(CXXFLAGS) $(SPREAD_FLAGS) -lfftw3 -DTORQUE_Y -fPIC -fopenmp -lparticles

$(DPSTOKES_INSTALL)/libspreadInterpDerivZ.so: $(spreadInterpSRC) $(spreadInterpINC) $(DPSTOKES_INSTALL)/libparticles.so
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(spreadInterpSRC) $(CXXFLAGS) $(SPREAD_FLAGS) -lfftw3 -DTORQUE_Z -fPIC -fopenmp -lparticles

$(DPSTOKES_INSTALL)/libtransform.so: $(transformSRC) $(transformINC)
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(transformSRC) $(CXXFLAGS) $(FFTW_FLAGS) -lfftw3_omp -lfftw3 -lm -fopenmp -fPIC

$(DPSTOKES_INSTALL)/liblinSolve.so: $(linSolveSRC) $(linSolveINC)
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(linSolveSRC) $(CXXFLAGS) $(FFTW_FLAGS) $(LAPACKE_FLAGS) -lm $(LAPACKE_LIBS) -lfftw3 -fopenmp -fPIC

$(DPSTOKES_INSTALL)/libdpTools.so: $(dpToolsSRC) $(dpToolsINC) $(DPSTOKES_INSTALL)/liblinSolve.so
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(dpToolsSRC) $(CXXFLAGS) $(FFTW_FLAGS) $(LAPACKE_FLAGS) -lm $(LAPACKE_LIBS) \
	-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE -fopenmp -lfftw3_omp -lfftw3 -fPIC -llinSolve

$(DPSTOKES_INSTALL)/libBC.so: $(bcSRC) $(bcINC) $(DPSTOKES_INSTALL)/libgrid.so $(DPSTOKES_INSTALL)/libparticles.so
	@mkdir -p $(DPSTOKES_INSTALL)
	$(CXX) -o $@ $(bcSRC) $(CXXFLAGS) -fPIC -fopenmp -lgrid -lparticles 
