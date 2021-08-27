#include<omp.h>
#include"LinearSolvers.h"
#ifdef INTEL
#include<mkl_blas.h>
#else
#include<cblas.h>
#endif

extern "C"
{
  void  precomputeLinOps(const double* A_diags, double* storage_kbp, double* B, 
                         const double* C, const double* D, double* G, 
                         double* G_inv, int Nyx, int Nz, int n_threads)
  {
    const double one = 1;
    const double zero = 0;
    const char notrans = 'N';
    int nrhs = 2; int Na = 3 * Nz - 4, Ns = Na + 1;
    double* Xpenta = (double*) fftw_malloc(Nyx * nrhs * Nz * sizeof(Complex));
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int i = 1; i < Nyx; ++i)
    {
      unsigned int offset_a = i * Na;
      unsigned int offset_stor = i * Ns;
      const double* d0 = &(A_diags[offset_a]);
      const double* dp2 = d0 + Nz;
      const double* dm2 = dp2 + Nz - 2;
      double* stor = &(storage_kbp[offset_stor]);
      // precompute for kbpenta solver
      precompute_KBPENTA(d0, dp2, dm2, stor, Nz);
      unsigned int offset_bc = nrhs * Nz * i;
      unsigned int offset_dg = nrhs * nrhs * i;
      double* b = &(B[offset_bc]);
      double* x = &(Xpenta[offset_bc]);
      const double* c = &(C[offset_bc]);
      const double* d = &(D[offset_dg]);
      double* g = &(G[offset_dg]);
      double* g_inv = &(G_inv[offset_dg]);
      // solve Ax = b, b is overwritten with A^-1 b
      solve_KBPENTA_real(x, b, stor, Nz, nrhs);
      #ifdef INTEL
      dgemm(&notrans, &notrans, &nrhs, &nrhs, &Nz, &one, 
            c, &nrhs, b, &Nz, &zero, g, &nrhs);
      #else
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                  nrhs, nrhs, Nz, 1, c, nrhs, b, Nz, 0, g, nrhs);
      #endif    
  
      g[0] -= d[0]; g[1] -= d[1]; g[2] -= d[2]; g[3] -= d[3];
      double det = 1.0 / (g[0] * g[3] - g[2] * g[1]);
      g_inv[0] = g[3] * det;
      g_inv[1] = -g[1] * det;
      g_inv[2] = -g[2] * det;
      g_inv[3] = g[0] * det;
    }
    fftw_free(Xpenta);
  }
}
