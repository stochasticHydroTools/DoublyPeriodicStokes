#include"Grid.h"
#include"exceptions.h"
#include<fftw3.h>
/* C wrapper for calling from Python. Any functions
   defined here should also have their prototypes 
   and wrappers defined in Grid.py */
extern "C"
{
  Grid* MakeGrid()
  {
    Grid* grid = new Grid();
    return grid;
  }
  
  void SetupGrid(Grid* grid) {grid->setup();}
  /* set interior grid extent */
  void SetL(Grid* grid, const double Lx, const double Ly, const double Lz,
            const double minX, const double minY, const double minZ) 
  {
    grid->Lx = Lx; grid->Ly = Ly; grid->Lz = Lz;
    grid->minX = minX; grid->minY = minY; grid->minZ = minZ;
  }
  /* set interior grid size */
  void SetN(Grid* grid, const unsigned int Nx, const unsigned int Ny, 
            const unsigned int Nz) 
  {
    grid->Nx = Nx; grid->Ny = Ny; grid->Nz = Nz;
  }
  /* set interior grid spacing, with hz=0 indicating a non-uniform z grid */
  void Seth(Grid* grid, const double hx, const double hy, const double hz) 
  { 
    grid->hx = hx; grid->hy = hy; grid->hz = hz;
  }
  /* set z grid and weights, if hz is set to 0 */
  void SetZ(Grid* grid, const double* zpts, const double* zwts) 
  {
    grid->zG = (double*) fftw_malloc(grid->Nz * sizeof(double));
    grid->zG_wts = (double*) fftw_malloc(grid->Nz * sizeof(double));  
    std::copy(zpts, zpts + grid->Nz, grid->zG);
    std::copy(zwts, zwts + grid->Nz, grid->zG_wts); 
  } 
  /* set periodicity of each axis */ 
  void SetPeriodicity(Grid* grid, bool x, bool y, bool z)
  {
    grid->isperiodic[0] = x; grid->isperiodic[1] = y;
    grid->isperiodic[2] = z; grid->has_periodicity = true;
  }
  /* set boundary conditions for end of each axis */
  void SetBCs(Grid* grid, unsigned int* BCs)
  {
    //if (grid->BCs) {free(grid->BCs); grid->BCs = 0;}
    if (!grid->BCs) {grid->BCs = (BC*) malloc(grid->dof * 6 * sizeof(BC));}
    #pragma omp simd
    for (unsigned int i = 0; i < 6 * grid->dof; ++i)
    {
      grid->BCs[i] = reinterpret_cast<BC*>(BCs)[i];
    }
    grid->has_bc = true; 
  } 
  void Setdof(Grid* grid, const unsigned int dof) {grid->dof = dof;}
  void ZeroExtGrid(Grid* grid){grid->zeroExtGrid();}
  void ZeroExtGrid_ghost(Grid* grid){grid->zeroExtGrid_ghost();}
  void ZeroIntGrid(Grid* grid){grid->zeroIntGrid();}
  void CleanGrid(Grid* g) {g->cleanup();}
  void DeleteGrid(Grid* g) {if(g) {delete g; g = 0;}} 
  double* GetData(Grid* g) {return g->fG;}
  void CopyData(Grid* g, double* fG) {std::copy(g->fG, g->fG + g->dof * g->Nx * g->Ny * g->Nz, fG);}
  void SetData(Grid* g, double* f){g->setData(f);} 

  void SetNumThreads(Grid* g, int num_threads) {g->n_threads =  num_threads;}
  void WriteGrid(Grid* g, const char* fname) {g->writeGrid(fname);}
  void WriteCoords(Grid* g, const char* fname) {g->writeCoords(fname);}  
}

