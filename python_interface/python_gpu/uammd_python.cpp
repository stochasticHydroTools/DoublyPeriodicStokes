#include"uammd_interface.h"
#include<pybind11/pybind11.h>
#include<pybind11/numpy.h>
namespace py = pybind11;
using namespace pybind11::literals;

class DPStokesPython2{
  std::shared_ptr<DPStokesPython> dpstokes;
public:

  //Initialize the modules with a certain set of parameters
  //Reinitializes if the module was already initialized
  void initialize(PyParameters pypar, int numberParticles){
    dpstokes = std::make_shared<DPStokesPython>();
    dpstokes->initialize(pypar, numberParticles);
  }

  //Clears all memory allocated by the module.
  //This leaves the module in an unusable state until initialize is called again.
  void clear(){dpstokes->clear();}
  //Set positions to compute mobility matrix
  void setPositions(py::array_t<real> h_pos){
    dpstokes->setPositions(h_pos.data());
  }

  //Compute the dot product of the mobility matrix with the forces and/or torques acting on the previously provided positions
  void Mdot(py::array_t<real> h_forces, py::array_t<real> h_torques,
	    py::array_t<real> h_MF,
	    py::array_t<real> h_MT){
    dpstokes->Mdot(h_forces.data(), h_torques.data(),
		   h_MF.mutable_data(), h_MT.mutable_data());
  }
};

//Pybind bindings
PYBIND11_MODULE(uammd, m) {
  m.doc() = "UAMMD DPStokes Python interface";
  py::class_<DPStokesPython2>(m, "DPStokes").
    def(py::init()).
    def("initialize", &DPStokesPython2::initialize,
	"Initialize the DPStokes module, can be called on an already initialize module to change the parameters.",
	"Parameters"_a, "numberParticles"_a).
    def("clear", &DPStokesPython2::clear, "Release all memory allocated by the module").
    def("setPositions", &DPStokesPython2::setPositions, "Set the positions to compute the mobility matrix",
	"positions"_a).
    def("Mdot", &DPStokesPython2::Mdot, "Computes the product of the Mobility tensor with the provided forces and torques. If torques are not present, they are assumed to be zero and angular displacements will not be computed",
	"forces"_a, "torques"_a = py::array_t<real>(),
	"velocities"_a, "angularVelocities"_a = py::array_t<real>());
  
  py::class_<PyParameters>(m, "StokesParameters").
    def(py::init([](real viscosity,
		    real  Lx, real Ly, real zmin, real zmax,
		    real w, real w_d,
		    real alpha, real alpha_d,
		    real beta, real beta_d,
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
