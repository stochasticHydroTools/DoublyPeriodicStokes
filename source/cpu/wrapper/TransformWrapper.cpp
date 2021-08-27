#include<omp.h>
#include "Transform.h"

/* C wrapper for calling from Python. Any functions
   defined here should also have their prototypes 
   and wrappers defined in Transform.py */
extern "C"
{

  /*
    This function initializes fftw with threads. Specifically,
      - single buffer of size (1 + 2*(Nx*Ny*Nz*dof)) is allocated
      - the first Nx*Ny*Nz*dof is the real part, the next is the
        complex part of the input/output (which are aliased for in-place transform)
      - forward and backward plans are created
      - If wisdom exists on disk (see init def for details on naming), then the 
        plans use the wisdom (much faster).
      - If wisdom doesn't exist, a plan is created using the FFTW_MEASURE flag and the
        wisdom is saved to disk as 'fftw_wisdom/fftw_wisdom_forward|backward_Nx[Nx]_Ny[Ny]_Nz[Nz]_dof[dof]'
      - We can change FFTW_MEASURE to FFTW_PATIENT or FFTW_EXHAUSTIVE to compute more 
        optimal plans, but with longer precompute cost.

      This function only has to be called once. After the call, 
        setFData, setBData, fTransform, bTransform can be called any number of times
      And, at the end of the program, one must call cleanup()     

    Parameters:
      Nx, Ny, Nz, dof - dimensions of the problem (grid and degrees of freedom)
      n_threads - number of threads for fftw
      type - type of transform
           - 0 = 3D Fourier transform on each dof
           - 1 = 3D Fourier-Chebyshev transfom (Fourier in x,y Cheb in z) on each dof
    Side Effects: 
      Buffers (in_real/in_imag) are allocated, forward/backward plans are determined (from disk if possible)   

    Returns:
      pointer to initialized transform struct
  */
  Transform* InitTransforms(const unsigned int Nx, const unsigned int Ny,
                            const unsigned int Nz, const unsigned int dof,
                            const unsigned int n_threads, transformType type)
  {
    if (not fftw_init_threads())
    {
      exitErr("Could not initialize threads for FFTW");
    }
    omp_set_num_threads(n_threads);
    Transform* t = new Transform();
    t->init(Nx, Ny, Nz, dof, n_threads, type);
    return t; 
  }
  
  /*
    This function sets the input data for fftw forward, based on transformType
    We assume forward transforms are always on real data
    
    If transformType is 0, the data is just copied into the internal buffers of t (regular fft)
    If transformType is 1, the data in x,y is copied, and the data in z reflected and copied (Fourier-Chebyshev)
  */
  void SetFData(Transform* t, double* in_real){t->setFData(in_real);}
  /*
    This function sets the input data for fftw backward, based on transformType
    
    If transformType is 0, the data is just copied into the internal buffers of t (regular fft)
    If transformType is 1, the data in x,y is copied, and the data in z reflected and copied (Fourier-Chebyshev)
  */
  void SetBData(Transform* t, double* in_real, double* in_imag){t->setBData(in_real, in_imag);}
  /* 
    Execute forward or backward transforms, with correct mode dispatched based on transformType
    
    If transformType = 0, t->out_real and t->out_imag have the output
    If transformType = 1, t->out_real_tr and t->out_imag_tr have the (truncated in z) output
      
  */
  void Ftransform(Transform* t) {t->fTransform();}
  void Btransform(Transform* t){t->bTransform();}
  
  /* get real output given transformType*/ 
  double* getRealOut(Transform* t)
  {
    if (t->type == fourier) {return t->out_real;}
    else if (t->type == fourierCheb) {return t->out_real_tr;}
    else {exitErr("Invalid transform type");}
  }
  
  /* get complex output given transformType*/ 
  double* getComplexOut(Transform* t) 
  {
    if (t->type == fourier) {return t->out_imag;}
    else if (t->type == fourierCheb) {return t->out_imag_tr;}
    else {exitErr("Invalid transform type");}
  }

  /* Delete all memory owned by t */
  void CleanTransform(Transform* t) {t->cleanup();}
  /* Delete t itself */
  void DeleteTransform(Transform* t) {if(t) {delete t; t = 0;}}
}
