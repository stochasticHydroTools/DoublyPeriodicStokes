#include"BoundaryConditions.h"
#include"Grid.h"
#include"ParticleList.h"

/* C wrapper for calling BoundaryConditions methods from Python. Any functions
   defined here should also have their prototypes 
   and wrappers defined in ParticleList.py */
extern "C"
{
  
  // fold spread data from ghost region of extended grid into the interior
  // according to periodicity or boundary condition for each data component
  void DeGhostify(Grid* grid, ParticleList* particles)
  {
    omp_set_num_threads(grid->n_threads);
    fold(grid->fG_unwrap, grid->fG, grid->Nx, grid->Ny, 
         grid->Nz, grid->Nxeff, grid->Nyeff, grid->Nzeff, 
         grid->dof, grid->isperiodic, grid->BCs);
  }

  // copy spread data from interior grid to ghost region of extended grid
  // according to periodicity or boundary condition for each data component
  void Ghostify(Grid* grid, ParticleList* particles)
  {
    omp_set_num_threads(grid->n_threads);
    copy(grid->fG_unwrap, grid->fG, grid->Nx, grid->Ny, 
         grid->Nz, grid->Nxeff, grid->Nyeff, grid->Nzeff, 
         grid->dof, grid->isperiodic, grid->BCs);
  }
}
