CURDIR          = $(shell pwd)
CXX             = g++ 

ifneq ($(debug), True)
  CXXFLAGS      = -I$(CURDIR)/include -w -O3 -march=native -shared -L$(INSTALL_DIR)
else
  CXXFLAGS      = -I$(CURDIR)/include -w -shared $(CPU_DEBUG) -L$(INSTALL_DIR)
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

LIBS            = $(patsubst %,$(INSTALL_DIR)/%,$(LIBS_))

all: $(LIBS)

$(INSTALL_DIR)/libcheb.so: $(chebSRC) $(chebINC)
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(chebSRC) $(CXXFLAGS) -fPIC -fopenmp

$(INSTALL_DIR)/libgrid.so: $(gridSRC) $(gridINC) 
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(gridSRC) $(CXXFLAGS) $(FFTW_FLAGS) -lfftw3 -lm -fPIC -fopenmp

$(INSTALL_DIR)/libparticles.so: $(particlesSRC) $(particlesINC) $(INSTALL_DIR)/libgrid.so $(INSTALL_DIR)/libcheb.so
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(particlesSRC) $(CXXFLAGS) -fPIC -fopenmp -lgrid -lcheb

$(INSTALL_DIR)/libspreadInterp.so: $(spreadInterpSRC) $(spreadInterpINC) $(INSTALL_DIR)/libparticles.so
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(spreadInterpSRC) $(CXXFLAGS) $(SPREAD_FLAGS) -fPIC -fopenmp -lparticles

$(INSTALL_DIR)/libspreadInterpDerivX.so: $(spreadInterpSRC) $(spreadInterpINC) $(INSTALL_DIR)/libparticles.so
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(spreadInterpSRC) $(CXXFLAGS) $(SPREAD_FLAGS) -DTORQUE_X -fPIC -fopenmp -lparticles
  
$(INSTALL_DIR)/libspreadInterpDerivY.so: $(spreadInterpSRC) $(spreadInterpINC) $(INSTALL_DIR)/libparticles.so
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(spreadInterpSRC) $(CXXFLAGS) $(SPREAD_FLAGS) -DTORQUE_Y -fPIC -fopenmp -lparticles

$(INSTALL_DIR)/libspreadInterpDerivZ.so: $(spreadInterpSRC) $(spreadInterpINC) $(INSTALL_DIR)/libparticles.so
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(spreadInterpSRC) $(CXXFLAGS) $(SPREAD_FLAGS) -DTORQUE_Z -fPIC -fopenmp -lparticles

$(INSTALL_DIR)/libtransform.so: $(transformSRC) $(transformINC)
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(transformSRC) $(CXXFLAGS) $(FFTW_FLAGS) -lfftw3_omp -lfftw3 -lm -fopenmp -fPIC

$(INSTALL_DIR)/liblinSolve.so: $(linSolveSRC) $(linSolveINC)
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(linSolveSRC) $(CXXFLAGS) $(FFTW_FLAGS) $(LAPACKE_FLAGS) -lm $(LAPACKE_LIBS) -lfftw3 -fopenmp -fPIC

$(INSTALL_DIR)/libdpTools.so: $(dpToolsSRC) $(dpToolsINC) $(INSTALL_DIR)/liblinSolve.so
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(dpToolsSRC) $(CXXFLAGS) $(FFTW_FLAGS) $(LAPACKE_FLAGS) -lm $(LAPACKE_LIBS) \
	-DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_STRUCTURE -fopenmp -lfftw3_omp -lfftw3 -fPIC -llinSolve

$(INSTALL_DIR)/libBC.so: $(bcSRC) $(bcINC) $(INSTALL_DIR)/libgrid.so $(INSTALL_DIR)/libparticles.so
	@mkdir -p $(INSTALL_DIR)
	$(CXX) -o $@ $(bcSRC) $(CXXFLAGS) -fPIC -fopenmp -lgrid -lparticles 

clean:
	rm -rf $(INSTALL_DIR)
