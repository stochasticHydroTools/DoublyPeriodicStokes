#ifndef TRANSFORM_H
#define TRANSFORM_H
#include<fftw3.h>
#include"exceptions.h"
#include"common.h"

/* Forward and Backward Fourier transform object
   which wraps desired features of fftw */

enum transformType {fourier, fourierCheb, unset};
struct Transform
{
  // real and imag input
  double *in_real, *in_imag;
  // real and imag output (these are aliased to input ptrs)
  double *out_real, *out_imag;
  // truncated output for Fourier-Chebyshev transforms
  double *out_real_tr, *out_imag_tr;
  // forward and backward plans
  fftw_plan pF, pB;
  // structs for configuring mem layout
  fftw_iodim *dims, *howmany_dims; 
  // mem layout
  unsigned int Nx, Ny, Nz;
  // degrees of freedom in the input, dimension of the problem
  // eg. if 4-component vector field in 3D, dof = 4 and rank = 3
  unsigned int dof, rank;
  // internal flag indicating whether we do a forward or back transform 
  int mode;
  // num threads for fftw
  int n_threads;
  transformType type;

  // null constructor
  Transform();
  
  /*
    This function initializes fftw with threads. Specifically,
      - single buffer of size (1 + 2*(Nx*Ny*Nz*dof)) is allocated
      - the first Nx*Ny*Nz*dof is the real part, the next is the
        imag part of the input/output (which are aliased for in-place transform)
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
    
    Note: mem layout is (Ny, Nx, Nz, dof) with inner dim the fastest

  */
  void init(const unsigned int Nx, const unsigned int Ny,
            const unsigned int Nz, const unsigned int dof,
            const unsigned int n_threads, transformType type);
  /*
    This function sets the input data for fftw forward, based on transformType
    We assume forward transforms are always on real data
  */
  void setFData(const double* in_real); 
  /*
    This function sets the input data for fftw backward, based on transformType
  */
  void setBData(const double* in_real, const double* in_imag);
  /* execute forward or backward transforms, with correct mode dispatched based on transformType */
  void fTransform();
  void bTransform();
  // configure memory layout (one for fftw, the other for internal use)
  void configDims();
  void setMemLayout();

  void cleanup();

};

#endif 
