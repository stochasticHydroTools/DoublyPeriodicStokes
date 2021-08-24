/* Raul P. Pelaez 2021. Doubly Periodic Stokes python bindings
   Allows to call the DPStokes module from python to compute the product between the mobility tensor and a list forces and torques acting on a group of positions.
   For additional info use:
   import uammd
   help(uammd)

*/
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include <uammd.cuh>
//Doubly Periodic FCM implementation (currently without noise)
#include <Integrator/BDHI/DoublyPeriodic/DPStokesSlab.cuh>
//Triply Periodic FCM implementation
#include <Integrator/BDHI/BDHI_FCM.cuh>

//Some convenient aliases
namespace py = pybind11;
using FCM_BM = uammd::BDHI::FCM_ns::Kernels::BarnettMagland;
using FCM = uammd::BDHI::FCM_impl<FCM_BM, FCM_BM>;
using DPStokesSlab = uammd::DPStokesSlab_ns::DPStokes;
using uammd::DPStokesSlab_ns::WallMode;
using uammd::System;
using real = uammd::real;
struct PyParameters{
  //The number of cells in each direction
  //If -1, they will be autocomputed from the tolerance if possible (DP cannot do it, FCM can)
  int nx = -1;
  int ny = -1;
  int nz = -1;
  real dt;
  real viscosity;
  real Lx;
  real Ly;
  real zmin, zmax;
  //Tolerance will be ignored in DP mode, TP will use only tolerance and nxy/nz
  real tolerance = 1e-5;
  real w, w_d;
  real hydrodynamicRadius = -1;
  real beta = -1;
  real beta_d = -1;
  real alpha = -1;
  real alpha_d = -1;
  //Can be either none, bottom, slit or periodic
  std::string mode;
};

//Helper functions and objects
struct Real3ToReal4{
  __host__ __device__ uammd::real4 operator()(uammd::real3 i){
    auto pr4 = uammd::make_real4(i);
    return pr4;
  }
};
struct Real4ToReal3{
  __host__ __device__ uammd::real3 operator()(uammd::real4 i){
    auto pr3 = uammd::make_real3(i);
    return pr3;
  }
};

struct Real3ToReal4SubstractOriginZ{
  real origin;
  Real3ToReal4SubstractOriginZ(real origin):origin(origin){}
  __host__ __device__ uammd::real4 operator()(uammd::real3 i){
    auto pr4 = uammd::make_real4(i);
    pr4.z -= origin;
    return pr4;
  }
};

FCM::Parameters createFCMParameters(PyParameters pypar){
  FCM::Parameters par;
  par.temperature = 0; //FCM can compute fluctuations, but they are turned off here
  par.viscosity = pypar.viscosity;
  par.tolerance = pypar.tolerance;
  par.box = uammd::Box({pypar.Lx, pypar.Ly, pypar.zmax- pypar.zmin});
  par.cells = {pypar.nx, pypar.ny, pypar.nz};
  par.kernel = std::make_shared<FCM_BM>(pypar.w, pypar.alpha, pypar.beta, pypar.Lx/pypar.nx);
  par.kernelTorque = std::make_shared<FCM_BM>(pypar.w_d, pypar.alpha_d, pypar.beta_d, pypar.Lx/pypar.nx);
  return par;
}

WallMode stringToWallMode(std::string str){
  if(str.compare("nowall") == 0){
    return WallMode::none;
  }
  else if(str.compare("slit") == 0){
    return WallMode::slit;
  }
  else if(str.compare("bottom") == 0){
    return WallMode::bottom;
  }
  else return WallMode::none;
}

DPStokesSlab::Parameters createDPStokesParameters(PyParameters pypar){
  DPStokesSlab::Parameters par;
  par.nx         = pypar.nx;
  par.ny         = pypar.ny;
  par.nz	  = pypar.nz;
  par.dt	  = pypar.dt;
  par.viscosity	  = pypar.viscosity;
  par.Lx	  = pypar.Lx;
  par.Ly	  = pypar.Ly;
  par.H		  = pypar.zmax-pypar.zmin;
  par.w = pypar.w;
  par.w_d = pypar.w_d;
  par.hydrodynamicRadius = pypar.hydrodynamicRadius;
  par.beta = pypar.beta;
  par.beta_d = pypar.beta_d;
  par.alpha = pypar.alpha;
  par.alpha_d = pypar.alpha_d;
  par.mode = stringToWallMode(pypar.mode);
  return par;
}

//Wrapper to UAMMD's TP and DP hydrodynamic modules, python interface is below
struct DPStokesUAMMD {
private:
  auto computeHydrodynamicDisplacements(bool useTorque){
    auto force = pd->getForce(uammd::access::gpu, uammd::access::read);
    auto pos = pd->getPos(uammd::access::gpu, uammd::access::read);
    auto torque = pd->getTorqueIfAllocated(uammd::access::gpu, uammd::access::read);
    auto d_torques_ptr = useTorque?torque.raw():nullptr;
    if(fcm){
      return fcm->computeHydrodynamicDisplacements(pos.raw(), force.raw(),
						   d_torques_ptr, numberParticles, st);
    }
    else if(dpstokes){
      return dpstokes->Mdot(pos.raw(), force.raw(),
			    d_torques_ptr, numberParticles, st);
    }
  }
public:
  std::shared_ptr<DPStokesSlab> dpstokes;
  std::shared_ptr<FCM> fcm;
  std::shared_ptr<uammd::System> sys;
  std::shared_ptr<uammd::ParticleData> pd;
  int numberParticles;
  cudaStream_t st;
  thrust::device_vector<uammd::real3> tmp;
  real zOrigin;

  DPStokesUAMMD(PyParameters pypar, int numberParticles): numberParticles(numberParticles){
    this->sys = std::make_shared<uammd::System>();
    this->pd = std::make_shared<uammd::ParticleData>(numberParticles, sys);
    if(pypar.mode.compare("periodic")==0){
      auto par = createFCMParameters(pypar);
      this->fcm = std::make_shared<FCM>(par);
      zOrigin = 0;
    }
    else{
      auto par = createDPStokesParameters(pypar);
      this->dpstokes = std::make_shared<DPStokesSlab>(par);
      zOrigin = pypar.zmin + par.H*0.5;
    }
    CudaSafeCall(cudaStreamCreate(&st));
  }

  //Copy positions to UAMMD's ParticleData
  void setPositions(py::array_t<real> h_pos){
    tmp.resize(numberParticles);
    auto pos = pd->getPos(uammd::access::gpu, uammd::access::write);
    thrust::copy((uammd::real3*)h_pos.data(), (uammd::real3*)h_pos.data() + numberParticles,
		 tmp.begin());
    thrust::transform(thrust::cuda::par.on(st), tmp.begin(), tmp.end(),
		      pos.begin(), Real3ToReal4SubstractOriginZ(zOrigin));
  }

  //Compute the hydrodynamic displacements due to a series of forces and/or torques acting on the particles
  void Mdot(py::array_t<real> h_forces, py::array_t<real> h_torques,
	    py::array_t<real> h_MF,
	    py::array_t<real> h_MT){
    // static int uses = 0;
    // uses++;
    //if(uses>=10) isNVTXEnabled = true;
    tmp.resize(numberParticles);
    bool useTorque = h_torques.size() != 0;
    {
      auto force = pd->getForce(uammd::access::gpu, uammd::access::write);
      thrust::copy((uammd::real3*)h_forces.data(), (uammd::real3*)h_forces.data() + numberParticles, tmp.begin());
      thrust::transform(thrust::cuda::par.on(st),
			tmp.begin(), tmp.end(), force.begin(), Real3ToReal4());
    }
    if(useTorque){
      auto torque = pd->getTorque(uammd::access::gpu, uammd::access::write);
      thrust::copy((uammd::real3*)h_torques.data(), (uammd::real3*)h_torques.data() + numberParticles, tmp.begin());
      thrust::transform(thrust::cuda::par, tmp.begin(), tmp.end(), torque.begin(), Real3ToReal4());
    }
    auto mob = this->computeHydrodynamicDisplacements(useTorque);
    thrust::copy(mob.first.begin(), mob.first.end(), (uammd::real3*)h_MF.mutable_data());   
    if(mob.second.size()){
      thrust::copy(mob.second.begin(), mob.second.end(), (uammd::real3*)h_MT.mutable_data());
    }    
  }
  
  ~DPStokesUAMMD(){
    cudaDeviceSynchronize();
    cudaStreamDestroy(st);
  }

};


//Python interface for the DPStokes module, see the accompanying example for more information
/*Usage:
  1- Call initialize with a set of parameters
  2- Call setPositions (the format must be [x0 y0 z0 x1 y1 z1,...])
  3- Call Mdot
  4- Call clear to free any memory allocated by the module and ensure a gracious finish

initialize can be called again in order to change the parameters.
Calling initialize twice is cheaper than calling initialize, then clear, then initialize again.

 */
class DPStokesPython{
  std::shared_ptr<DPStokesUAMMD> dpstokes;
public:

  //Initialize the modules with a certain set of parameters
  //Reinitializes if the module was already initialized
  void initialize(PyParameters pypar, int numberParticles){
    dpstokes = std::make_shared<DPStokesUAMMD>(pypar, numberParticles);
  }

  //Clears all memory allocated by the module.
  //This leaves the module in an unusable state until initialize is called again.
  void clear(){
    dpstokes->sys->finish();
    dpstokes.reset();
  }

  //Set positions to compute mobility matrix
  void setPositions(py::array_t<real> h_pos){
    throwIfInvalid();
    dpstokes->setPositions(h_pos);
  }

  //Compute the dot product of the mobility matrix with the forces and/or torques acting on the previously provided positions
  void Mdot(py::array_t<real> h_forces, py::array_t<real> h_torques,
	    py::array_t<real> h_MF,
	    py::array_t<real> h_MT){
    throwIfInvalid();
    dpstokes->Mdot(h_forces, h_torques, h_MF, h_MT);
  }

private:
  void throwIfInvalid(){
    if(not dpstokes){
      throw std::runtime_error("DPStokes is not initialized. Call Initialize first");
    }
  }
};

using namespace pybind11::literals;

//Transform between the enumerator for selecting a mode and a string
std::string wallModeToString(WallMode mode){
  switch(mode){
  case WallMode::none:
    return "no wall";
  case WallMode::slit:
    return "slit channel";
  case WallMode::bottom:
    return "bottom wall";
  };
}


//Pybind bindings
PYBIND11_MODULE(uammd, m) {
  m.doc() = "UAMMD DPStokes Python interface";
  py::class_<DPStokesPython>(m, "DPStokes").
    def(py::init()).
    def("initialize", &DPStokesPython::initialize,
	"Initialize the DPStokes module, can be called on an already initialize module to change the parameters.",
	"Parameters"_a, "numberParticles"_a).
    def("clear", &DPStokesPython::clear, "Release all memory allocated by the module").
    def("setPositions", &DPStokesPython::setPositions, "Set the positions to compute the mobility matrix",
	"positions"_a).
    def("Mdot", &DPStokesPython::Mdot, "Computes the product of the Mobility tensor with the provided forces and torques. If torques are not present, they are assumed to be zero and angular displacements will not be computed",
	"forces"_a, "torques"_a = py::array_t<real>(),
	"velocities"_a, "angularVelocities"_a = py::array_t<real>());
  
  py::class_<PyParameters>(m, "StokesParameters").
    def(py::init([](uammd::real viscosity,
		    uammd::real  Lx, uammd::real Ly, uammd::real zmin, uammd::real zmax,
		    uammd::real w, uammd::real w_d,
		    uammd::real alpha, uammd::real alpha_d,
		    uammd::real beta, uammd::real beta_d,
		    int Nx, int Ny, int nz, std::string mode) {
      auto tmp = std::unique_ptr<PyParameters>(new PyParameters);
      tmp->viscosity = viscosity;
      tmp->Lx = Lx;
      tmp->Ly = Ly;
      tmp->zmin = zmin;
      tmp->zmax = zmax;
      tmp->nx = Nx;
      tmp->ny = Ny;
      tmp->nz = nz;
      tmp->mode = mode;
      tmp->w = w;
      tmp->w_d = w_d;
      tmp->beta =beta;
      tmp->beta_d = beta_d;
      tmp->alpha = alpha;
      tmp->alpha_d = alpha_d;
      return tmp;	
    }),"viscosity"_a  = 1.0,"Lx"_a = 0.0, "Ly"_a = 0.0, "zmin"_a = 0.0,"zmax"_a = 0.0,
	"w"_a=1.0, "w_d"_a=1.0,
	"alpha"_a = -1.0, "alpha_d"_a=-1.0,
	"beta"_a = -1.0, "beta_d"_a=-1.0,
	"nx"_a = -1,"ny"_a = -1, "nz"_a = -1, "mode"_a="none").
    def_readwrite("viscosity", &PyParameters::viscosity, "Viscosity").
    def_readwrite("Lx", &PyParameters::Lx, "Domain size in the plane").
    def_readwrite("Ly", &PyParameters::Ly, "Domain size in the plane").
    def_readwrite("zmin", &PyParameters::zmin, "Minimum height of a particle (or bottom wall location)").
    def_readwrite("zmax", &PyParameters::zmax, "Maximum height of a particle (or top wall location)").
    def_readwrite("mode", &PyParameters::mode, "Domain walls mode, can be any of: none (no walls), bottom (wall at the bottom), slit (two walls) or periodic (uses force coupling method).").
    def_readwrite("nz", &PyParameters::nz, "Number of cells in Z").
    def_readwrite("nx", &PyParameters::nx, "Number of cells in X").
    def_readwrite("ny", &PyParameters::ny, "Number of cells in Y").
    def_readwrite("alpha", &PyParameters::alpha, "ES kernel monopole alpha").
    def_readwrite("alpha_d", &PyParameters::alpha_d, "ES kernel dipole alpha").
    def_readwrite("beta", &PyParameters::beta, "ES kernel monopole beta").
    def_readwrite("beta_d", &PyParameters::beta_d, "ES kernel dipole beta").
    def_readwrite("w", &PyParameters::w, "ES kernel monopole width").
    def_readwrite("w_d", &PyParameters::w_d, "ES kernel dipole width").
    def("__str__", [](const PyParameters &p){
      return"viscosity = " + std::to_string(p.viscosity) +"\n"+
	"box (L = " + std::to_string(p.Lx) +
	"," + std::to_string(p.Ly) + "," +
	std::to_string(p.zmin) + ":" + std::to_string(p.zmax) +" )\n"+
	"Nx = " + std::to_string(p.nx) + "\n" +
	"Ny = " + std::to_string(p.ny) + "\n" +
	"nz = " + std::to_string(p. nz) + "\n" +
	"mode = " + p.mode + "\n";
    });
    
}
