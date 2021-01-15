/* Raul P. Pelaez 2021. Doubly Periodic Stokes python bindings
   Allows to call the DPStokes module from python to compute the product between the mobility tensor and a list forces acting on a group of Gaussian sources.
   For additional info use:
   import uammd
   help(uammd)

*/
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
#include <uammd.cuh>
#include <Integrator/BDHI/DoublyPeriodic/DPStokesSlab.cuh>

namespace py = pybind11;
using DPStokesSlab = uammd::DPStokesSlab_ns::DPStokes;
using Parameters = DPStokesSlab::Parameters;

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

struct UAMMD {
  using real = uammd::real;
  std::shared_ptr<DPStokesSlab> dpstokes;
  std::shared_ptr<uammd::System> sys;
  int numberParticles;
  cudaStream_t st;
  thrust::device_vector<uammd::real4> pos, force;
  thrust::device_vector<uammd::real3> tmp;
  UAMMD(Parameters par, int numberParticles): numberParticles(numberParticles){
    this->sys = std::make_shared<uammd::System>();
    this->dpstokes = std::make_shared<DPStokesSlab>(par);
    CudaSafeCall(cudaStreamCreate(&st));
  }

  void Mdot(py::array_t<real> h_pos, py::array_t<real> h_forces, py::array_t<real> h_MF){
    pos.resize(numberParticles);
    force.resize(numberParticles);
    tmp.resize(numberParticles);
    thrust::copy((uammd::real3*)h_pos.data(), (uammd::real3*)h_pos.data() + numberParticles, tmp.begin());
    thrust::transform(tmp.begin(), tmp.end(), pos.begin(), Real3ToReal4());
    thrust::copy((uammd::real3*)h_forces.data(), (uammd::real3*)h_forces.data() + numberParticles, tmp.begin());
    thrust::transform(tmp.begin(), tmp.end(), force.begin(), Real3ToReal4());
    auto d_pos = thrust::raw_pointer_cast(pos.data());
    auto d_force = thrust::raw_pointer_cast(force.data());
    auto MF = dpstokes->Mdot(d_pos, d_force, numberParticles, st);
    auto MF_real3 = thrust::make_transform_iterator(MF.begin(), Real4ToReal3());
    thrust::copy(MF_real3, MF_real3 + numberParticles, (uammd::real3*)h_MF.mutable_data());
  }
  
  ~UAMMD(){
    cudaDeviceSynchronize();
    cudaStreamDestroy(st);
  }
};



using namespace pybind11::literals;

PYBIND11_MODULE(uammd, m) {
  m.doc() = "UAMMD DPStokes Python interface";
  py::class_<UAMMD>(m, "DPStokes").
    def(py::init<Parameters, int>(),"Parameters"_a, "numberParticles"_a).
    def("Mdot", &UAMMD::Mdot, "Computes the product of the Mobility tensor with a provided array",
	"positions"_a,"forces"_a,"result"_a);
  
  py::class_<uammd::Box>(m, "Box").
    def(py::init<uammd::real>()).
    def(py::init([](uammd::real x, uammd::real y, uammd::real z) {
      return std::unique_ptr<uammd::Box>(new uammd::Box(uammd::make_real3(x,y,z)));
    }));

  py::class_<Parameters>(m, "DPStokesParameters").
    def(py::init([](uammd::real viscosity,
		    uammd::real  Lxy, uammd::real H,
		    uammd::real gw,
		    int support, int Nxy, int nz) {             
      auto tmp = std::unique_ptr<Parameters>(new Parameters);
      tmp->viscosity = viscosity;
      tmp->box = uammd::Box(uammd::make_real3(Lxy, Lxy, H));
      tmp->box.setPeriodicity(1,1,0);
      tmp->gw = gw;
      tmp->support = support;
      tmp->cells = make_int3(Nxy, Nxy, nz);
      return tmp;	
    }),"viscosity"_a  = 1.0,"Lxy"_a = 0.0,"H"_a = 0.0,"gw"_a=1.0, "support"_a = -1, "Nxy"_a=-1, "nz"_a = -1).
    def_readwrite("viscosity", &Parameters::viscosity).
    def_readwrite("gw", &Parameters::gw).
    def_readwrite("box", &Parameters::box).
    def_readwrite("support", &Parameters::support).
    def("__str__", [](const Parameters &p){
      return"viscosity = " + std::to_string(p.viscosity) +"\n"+
	"gw = " + std::to_string(p. gw)+ "\n" +
	"box (L = " + std::to_string(p.box.boxSize.x) +
	"," + std::to_string(p.box.boxSize.y) + "," + std::to_string(p.box.boxSize.z) + ")\n"+
	"support = " + std::to_string(p. support)+ "\n" + 
	"Nxy = " + std::to_string(p. cells.x) + "\n" +
	"nz = " + std::to_string(p. cells.z) + "\n";
    });
    
}
