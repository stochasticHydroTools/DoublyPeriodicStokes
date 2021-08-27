#include"SpreadInterp.h"
#include"ParticleList.h"
#include"Grid.h"
#include"exceptions.h"
#include<fftw3.h>
#include<omp.h>


/* C wrapper for calling from Python. Any functions
   defined here should also have their prototypes 
   and wrappers defined in SpreadInterp.py */

extern "C"
{
  // Spread and get pointer to data
  void Spread(ParticleList* s, Grid* g) 
  { 
    omp_set_num_threads(g->n_threads);
    // spread from particles onto grid
    if (g->unifZ) {spreadUnifZ(*s, *g);}
    else {spreadNonUnifZ(*s, *g);} 
  }

  // Spread and get pointer to data
  void Interpolate(ParticleList* s, Grid* g) 
  {
    omp_set_num_threads(g->n_threads);
    // interpolate data from grid onto particles
    if (g->unifZ) {interpUnifZ(*s, *g);}
    else {interpNonUnifZ(*s, *g);}
  }
  
  /* set the dof to a new one so we can re-use search structures
     for spread and interp for any dof.

     the grid must already be setup with a non-zero dof. 
     the particles must already be setup with a non-zero dof.
     if the new and old dof are the same, we do nothing.
  
     the grid must also already have a locator */
  void Resetdof(Grid* g, ParticleList* s, const unsigned int dof)
  { 
    if (s->dof != dof)
    {
      if (s->fP) {fftw_free(s->fP); s->fP = 0;}
      s->fP = (double*) fftw_malloc(s->nP * dof * sizeof(double));
      s->dof = dof;
    }
    if (g->dof && g->has_locator)
    {
      if (g->dof != dof)
      {
        if (g->fG) {fftw_free(g->fG); g->fG = 0;}
        if (g->fG_unwrap) {fftw_free(g->fG_unwrap); g->fG_unwrap = 0;}
        g->fG = (double*) fftw_malloc(g->Nx * g->Ny * g->Nz * dof * sizeof(double));
        g->fG_unwrap = (double*) fftw_malloc(g->Nxeff * g->Nyeff * g->Nzeff * dof * sizeof(double));
        g->dof = dof;
      }
    }
    else
    {
      exitErr("Grid has not been setup");
    }
  }

  void addDx(Grid* gm, Grid* gd)
  {
    unsigned int N = gm->Nz * gm ->Ny * gm->Nx;
    double* fmy = gm->fG + 1;
    double* fmz = gm->fG + 2;
    double* fdy = gd->fG + 1;
    double* fdz = gd->fG + 2;
    #pragma omp parallel for num_threads(gm->n_threads)
    for (unsigned int i = 0; i < 3 * N; i += 3)
    {
      fmy[i] -= 0.5 * fdz[i];    
      fmz[i] += 0.5 * fdy[i];    
    }
  }
  void addDy(Grid* gm, Grid* gd)
  {
    unsigned int N = gm->Nz * gm ->Ny * gm->Nx;
    double* fmx = gm->fG;
    double* fmz = gm->fG + 2;
    double* fdx = gd->fG;
    double* fdz = gd->fG + 2;
    #pragma omp parallel for num_threads(gm->n_threads) 
    for (unsigned int i = 0; i < 3 * N; i += 3)
    {
      fmx[i] += 0.5 * fdz[i];
      fmz[i] -= 0.5 * fdx[i];
    }
  }
  void addDz(Grid* gm, Grid* gd)
  {
    unsigned int N = gm->Nz * gm ->Ny * gm->Nx;
    double* fmx = gm->fG;
    double* fmy = gm->fG + 1;
    double* fdx = gd->fG;
    double* fdy = gd->fG + 1;
    #pragma omp parallel for num_threads(gm->n_threads)
    for (unsigned int i = 0; i < 3 * N; i += 3)
    {
      fmx[i] -= 0.5 * fdy[i];
      fmy[i] += 0.5 * fdx[i];
    }
  }
}
