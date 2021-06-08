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
  int nxy = -1;
  int nz = -1;
  real dt;
  real viscosity;
  real Lxy;
  real H;
  //Tolerance will be ignored in DP mode
  real tolerance = 1e-7;
  real gw; //Gaussian width, unused with the BM kernel and in TP
  int support = -1; //-1 means auto compute from tolerance if possible (FCM can do this)
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

FCM::Parameters createFCMParameters(PyParameters pypar){
  FCM::Parameters par;
  par.temperature = 0; //FCM can compute fluctuations, but they are turned off here
  par.viscosity = pypar.viscosity;
  par.tolerance = pypar.tolerance;
  par.box = uammd::Box({pypar.Lxy, pypar.Lxy, pypar.H});
  par.cells = {pypar.nxy, pypar.nxy, pypar.nz};
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
  par.nxy         = pypar.nxy;
  par.nz	  = pypar.nz;
  par.dt	  = pypar.dt;
  par.viscosity	  = pypar.viscosity;
  par.Lxy	  = pypar.Lxy;
  par.H		  = pypar.H;
  par.gw	  = pypar.gw;
  par.support 	  = pypar.support;
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
  UAMMD(PyParameters pypar, int numberParticles): numberParticles(numberParticles){
    this->sys = std::make_shared<uammd::System>();
    this->pd = std::make_shared<uammd::ParticleData>(numberParticles, sys);
    if(pypar.mode.compare("periodic")==0){
      auto par = createFCMParameters(pypar);
      this->fcm = std::make_shared<FCM>(pd, sys, par);
    }
    else{
      auto par = createDPStokesParameters(pypar);
      this->dpstokes = std::make_shared<DPStokesSlab>(par);
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
      thrust::transform(thrust::cuda::par, tmp.begin(), tmp.end(), pos.begin(), Real3ToReal4());
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
		    uammd::real  Lxy, uammd::real H,
		    uammd::real gw,
		    int support, int Nxy, int nz, std::string mode) {
      auto tmp = std::unique_ptr<PyParameters>(new PyParameters);
      tmp->viscosity = viscosity;
      tmp->Lxy = Lxy;
      tmp->H = H;      
      tmp->gw = gw;
      tmp->support = support;
      tmp->nxy = Nxy;
      tmp->nz = nz;
      tmp->mode = mode;
      return tmp;	
    }),"viscosity"_a  = 1.0,"Lxy"_a = 0.0,"H"_a = 0.0,"gw"_a=1.0, "support"_a = -1, "Nxy"_a=-1, "nz"_a = -1, "mode"_a="none").
    def_readwrite("viscosity", &PyParameters::viscosity, "Viscosity").
    def_readwrite("gw", &PyParameters::gw, "Gaussian width of the sources").
    def_readwrite("Lxy", &PyParameters::Lxy, "Domain size in the plane").
    def_readwrite("H", &PyParameters::H, "Domain width").
    def_readwrite("support", &PyParameters::support, "Number of support cells for spreading/interpolation").
    def_readwrite("mode", &PyParameters::mode, "Domain walls mode, can be any of: none (no walls), bottom (wall at the bottom), slit (two walls) or periodic (uses force coupling method).").
    def_readwrite("nz", &PyParameters::nz, "Number of cells in Z").
    def_readwrite("nxy", &PyParameters::nxy, "Number of cells in XY").
    def("__str__", [](const PyParameters &p){
      return"viscosity = " + std::to_string(p.viscosity) +"\n"+
	"gw = " + std::to_string(p. gw)+ "\n" +
	"box (L = " + std::to_string(p.Lxy) +
	"," + std::to_string(p.Lxy) + "," + std::to_string(p.H) + ")\n"+
	"support = " + std::to_string(p. support)+ "\n" + 
	"Nxy = " + std::to_string(p. nxy) + "\n" +
	"nz = " + std::to_string(p. nz) + "\n" +
	"mode = " + p.mode + "\n";
    });
    
}
