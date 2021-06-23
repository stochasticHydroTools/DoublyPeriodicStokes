/* Raul P. Pelaez 2021. Doubly Periodic Stokes python bindings
   Allows to call the DPStokes module from python to compute the product between the mobility tensor and a list forces and torques acting on a group of positions.
   For additional info use:
   import uammd
   help(uammd)

*/
#include "Integrator/BDHI/DoublyPeriodic/StokesSlab/utils.cuh"
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include <uammd.cuh>
#include <Integrator/BDHI/DoublyPeriodic/DPStokesSlab.cuh>
#include <Integrator/BDHI/BDHI_FCM.cuh>


namespace py = pybind11;
using uammd::BDHI::FCM;
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
  real tolerance = 1e-7;
  real w, w_d;
  real hydrodynamicRadius;
  real beta = -1;
  real beta_d = -1;
  real alpha = -1;
  real alpha_d = -1;
  //Can be either none, bottom, slit or periodic
  std::string mode;
};


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

struct UAMMD {
  std::shared_ptr<DPStokesSlab> dpstokes;
  std::shared_ptr<FCM> fcm;
  std::shared_ptr<uammd::System> sys;
  std::shared_ptr<uammd::ParticleData> pd;
  int numberParticles;
  cudaStream_t st;
  thrust::device_vector<uammd::real3> tmp;
  real zOrigin;
  UAMMD(PyParameters pypar, int numberParticles): numberParticles(numberParticles){
    this->sys = std::make_shared<uammd::System>();
    this->pd = std::make_shared<uammd::ParticleData>(numberParticles, sys);
    if(pypar.mode.compare("periodic")==0){
      auto par = createFCMParameters(pypar);
      this->fcm = std::make_shared<FCM>(pd, sys, par);
      zOrigin = 0;
    }
    else{
      auto par = createDPStokesParameters(pypar);
      this->dpstokes = std::make_shared<DPStokesSlab>(par);
      zOrigin = pypar.zmin + par.H*0.5;
    }
    CudaSafeCall(cudaStreamCreate(&st));
  }

  void Mdot(py::array_t<real> h_pos, py::array_t<real> h_forces, py::array_t<real> h_torques,
	    py::array_t<real> h_MF,
	    py::array_t<real> h_MT){
    tmp.resize(numberParticles);
    bool useTorque = h_torques.size() != 0;
    {
      auto pos = pd->getPos(uammd::access::gpu, uammd::access::write);
      auto force = pd->getForce(uammd::access::gpu, uammd::access::write);
      thrust::copy((uammd::real3*)h_pos.data(), (uammd::real3*)h_pos.data() + numberParticles, tmp.begin());
      thrust::transform(thrust::cuda::par, tmp.begin(), tmp.end(), pos.begin(), Real3ToReal4SubstractOriginZ(zOrigin));
      thrust::copy((uammd::real3*)h_forces.data(), (uammd::real3*)h_forces.data() + numberParticles, tmp.begin());
      thrust::transform(thrust::cuda::par, tmp.begin(), tmp.end(), force.begin(), Real3ToReal4());
      if(useTorque){
	auto torque = pd->getTorque(uammd::access::gpu, uammd::access::write);
        thrust::copy((uammd::real3*)h_torques.data(), (uammd::real3*)h_torques.data() + numberParticles, tmp.begin());
	thrust::transform(thrust::cuda::par, tmp.begin(), tmp.end(), torque.begin(), Real3ToReal4());
      }
    }
    auto force = pd->getForce(uammd::access::gpu, uammd::access::read);
    auto pos = pd->getPos(uammd::access::gpu, uammd::access::read);
    if(fcm){
      if(h_torques.size() != 0){
    	System::log<System::EXCEPTION>("Cannot process torques in triply periodic mode");
    	throw std::runtime_error("Invalid mode");
      }
      auto tmp_ptr = thrust::raw_pointer_cast(tmp.data());
      fcm->Mdot(tmp_ptr, force.raw(), 0);
      thrust::copy(tmp.begin(), tmp.end(), (uammd::real3*)h_MF.mutable_data());
    }
    if(dpstokes){
      auto torque = pd->getTorqueIfAllocated(uammd::access::gpu, uammd::access::read);
      auto d_torques_ptr = useTorque?torque.raw():nullptr;
      //mob is a tuple containing MF and MT. The mobilities for translational and rotational contributions
      auto mob = dpstokes->Mdot(pos.raw(), force.raw(), d_torques_ptr, numberParticles, st);
      if(mob.second.size()){
	auto MT_real3 = thrust::make_transform_iterator(mob.second.begin(), Real4ToReal3());
	thrust::copy(MT_real3, MT_real3 + numberParticles, (uammd::real3*)h_MT.mutable_data());
      }
      auto MF_real3 = thrust::make_transform_iterator(mob.first.begin(), Real4ToReal3());
      thrust::copy(MF_real3, MF_real3 + numberParticles, (uammd::real3*)h_MF.mutable_data());
    }
  }
  
  ~UAMMD(){
    cudaDeviceSynchronize();
    cudaStreamDestroy(st);
  }
};

using namespace pybind11::literals;


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

PYBIND11_MODULE(uammd, m) {
  m.doc() = "UAMMD DPStokes Python interface";
  py::class_<UAMMD>(m, "DPStokes").
    def(py::init<PyParameters, int>(),"Parameters"_a, "numberParticles"_a).
    def("Mdot", &UAMMD::Mdot, "Computes the product of the Mobility tensor with the provided forces and torques. If torques are not present, they are assumed to be zero and angular displacements will not be computed",
	"positions"_a,"forces"_a, "torques"_a = py::array_t<real>(),
	"velocities"_a, "angularVelocities"_a = py::array_t<real>());
  
  py::class_<PyParameters>(m, "StokesParameters").
    def(py::init([](uammd::real viscosity,
		    uammd::real  Lx, uammd::real Ly, uammd::real zmin, uammd::real zmax,
		    uammd::real w, uammd::real w_d,
		    uammd::real alpha, uammd::real alpha_d,
		    uammd::real beta, uammd::real beta_d,
		    uammd::real hydrodynamicRadius,
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
      tmp->hydrodynamicRadius = hydrodynamicRadius;
      tmp->beta =beta;
      tmp->beta_d = beta_d;
      tmp->alpha = alpha;
      tmp->alpha_d = alpha_d;
      return tmp;	
    }),"viscosity"_a  = 1.0,"Lx"_a = 0.0, "Ly"_a = 0.0, "zmin"_a = 0.0,"zmax"_a = 0.0,
	"w"_a=1.0, "w_d"_a=1.0,
	"alpha"_a = -1.0, "alpha_d"_a=-1.0,
	"beta"_a = -1.0, "beta_d"_a=-1.0,
	"hydrodynamicRadius"_a = 1.0,
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
    def_readwrite("hydrodynamicRadius", &PyParameters::hydrodynamicRadius, "Hydrodynamic radius").
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
