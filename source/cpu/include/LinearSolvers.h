#ifndef _LINEAR_SOLVERS_H
#define _LINEAR_SOLVERS_H
#include<complex.h>
#include<fftw3.h>
extern "C"
{
  typedef double _Complex Complex;  
  /*************************************************************************
    Given a block linear system of the form 
    
    |A B||a |   |f|
    |C D||c0| = |alpha| 
         |d0|   |beta|,  
                
    
    where A is Nz x Nz and banded, B is Nz x 2, C is 2 x Nz and D is 2x2 (all real), 
    this function computes the inverse of the schur complement of A,
    
        i.e. Ginv = (C A^{-1} B - D)^{-1},
    
    as well as A^{-1}B, and info for later applications of A^{-1}. It uses the
    KBPENTA solver for pentadiagonal systems, specialized for those with 
    only 3 non-zero diags
      
      http://dx.doi.org/10.1080/00207160802326507
    
    The computations are performed for every k \in [0, Ny * Nx]

    Parameters:
      A - I - k^2 * Phi2, where Phi2 is the Chebyshev second integral matrix
        - this only stores the 3 non-zero diagonals at each k 
      B - Nz x 2 x Nyx tensor (stored in Fortran order)
      C - 2 x Nz x Nyx tensor (stored in Fortran order)
      D - 2 x 2 x Nyx tensor (stored in Fortran order)
      G - 2 x 2 x Nyx tensor of zeros (stored in Fortran order)
      Ginv - 2 x 2 x Nyx tensor of zeros (stored in Fortran order)
      Nyx - Nx * Ny (total points in x-y plane)
      Nz -  num points in z
    
    Side Effects:
      A_kbp is overwritten with info for future solves with A
      B is overwritten with A^{-1}B (from KBPENTA solve)
      G is overwritten with (C A^{-1} B - D)
      Ginv is overwritten with G^{-1} explicitly computed
  */
  void  precomputeLinOps(const double* A_diags, double* storage_kbp, double* B, 
                         const double* C, const double* D, double* G, 
                         double* G_inv, int Nyx, int Nz, int n_threads);


  // precompute storage required for kbpenta solve
  inline void precompute_KBPENTA(const double* diagonal, const double* diagonal_p2, 
                                 const double* diagonal_m2, double* storage, int nz)
  {	
    double* beta = (double*) fftw_malloc((nz + 1) * sizeof(double));
    beta[0] = 0;
    beta[1] = diagonal[0];
    beta[2] = diagonal[1];
    #pragma omp simd
    for(unsigned int i = 3; i <= nz; ++i)
    {
      beta[i] = diagonal[i - 1] - diagonal_m2[i - 3] * diagonal_p2[i - 3] / beta[i - 2];
    }
    for(unsigned int i = 0; i <= nz; ++i)
    {
      storage[i] = beta[i];
      if(i < nz - 2)
      {
	      storage[i + nz + 1] = diagonal_p2[i];
	      storage[i + 2 * nz - 1] = diagonal_m2[i];
      }
    }
    fftw_free(beta); beta = 0;
  }
  // solve for real system
  inline void solve_KBPENTA_real(double* xpenta, double* rhs, const double* storage, int nz, int nrhs)
  {

    const double* beta = storage;
    const double* diagonal_p2 = storage + nz + 1;
    const double* diagonal_m2 = storage + 2 * nz - 1;
    for (unsigned int irhs = 0; irhs < nrhs; ++irhs)
    {
      double* b = &(rhs[irhs * nz]);
      double* x = &(xpenta[irhs * nz]);
      x[0] = b[0];
      x[1] = b[1];
      #pragma omp simd
      for(unsigned int i = 2; i < nz; ++i)
      {
        x[i] = b[i] - diagonal_m2[i - 2] * x[i - 2] / beta[i - 1];
      }
      x[nz - 1] /= beta[nz];
      x[nz - 2] /= beta[nz - 1];
      b[nz - 1] = x[nz - 1]; 
      b[nz - 2] = x[nz - 2];
      for(int i = nz - 3; i >= 0; --i)
      {
          x[i] = (x[i] - diagonal_p2[i] * x[i + 2]) / beta[i + 1];
          b[i] = x[i];
      }
    }
  }

  // solve for complex system
  inline void solve_penta(Complex* xpenta, const Complex* rhs, const double* storage, int nz, int nrhs)
  {

    const double* beta = storage;
    const double* diagonal_p2 = storage + nz + 1;
    const double* diagonal_m2 = storage + 2 * nz - 1;
    Complex *x, *xmp, d2, bmp;
    const Complex* b;
    for (unsigned int i = 0; i < 2; ++i)
    {
      b = &(rhs[i * nrhs]); x = &(xpenta[i * nrhs]);
      #pragma omp simd
      for (unsigned int j = 0; j < nrhs; ++j) {x[j] = b[j];} 
    }
    for (unsigned int i = 2; i < nz; ++i)
    {
      b = &(rhs[i * nrhs]); x = &(xpenta[i * nrhs]); xmp = &(xpenta[(i - 2) * nrhs]);
      d2 = diagonal_m2[i - 2]; bmp = beta[i - 1]; 
      #pragma omp simd
      for (unsigned int j = 0; j < nrhs; ++j){x[j] = b[j] - d2 * xmp[j] / bmp;}
    }
    for (unsigned int j = 0; j < 2; ++j)
    {
      x = &(xpenta[(nz - j - 1) * nrhs]); 
      bmp = beta[nz -j];
      #pragma omp simd
      for (unsigned int i = 0; i < nrhs; ++i) {x[i] /= bmp;}
    }
    for(int i = nz - 3; i >= 0; --i)
    {
      x = &(xpenta[i * nrhs]); bmp = beta[i + 1];
      d2 = diagonal_p2[i]; xmp = &(xpenta[(i + 2) * nrhs]);
      #pragma omp simd
      for (unsigned int j = 0; j < nrhs; ++j) {x[j] = (x[j] - d2 * xmp[j]) / bmp;}
    }
  }

} // end extern

#endif
