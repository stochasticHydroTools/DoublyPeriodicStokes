
#include<string>
#include<memory>
using real = float;
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

class DPStokesUAMMD;
class DPStokesGlue{
  std::shared_ptr<DPStokesUAMMD> dpstokes;
public:

  //Initialize the modules with a certain set of parameters
  //Reinitializes if the module was already initialized
  void initialize(PyParameters pypar, int numberParticles);

  //Clears all memory allocated by the module.
  //This leaves the module in an unusable state until initialize is called again.
  void clear();
  //Set positions to compute mobility matrix
  void setPositions(const real* h_pos);

  //Compute the dot product of the mobility matrix with the forces and/or torques acting on the previously provided positions
  void Mdot(const real* h_forces, const real* h_torques,
	    real* h_MF,
	    real* h_MT);

private:
  void throwIfInvalid();
};

namespace uammd_wrapper{
  std::string getPrecision();
}
