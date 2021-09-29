#include<complex.h> 

typedef double _Complex Complex;
extern "C"
{
  /* 
  Solve triply periodic Stokes eq in Fourier domain given the Fourier
  coefficients of the forcing.
  
  Parameters:
    U_hat_r[i], P_hat_r[i] - real and complex part of vel and pressure Fourier
                             coeffs (output)
    fG_hat_r, fG_hat_i - real and complex part of Fourier coefficients
                         of spread forces on the grid. These are both
                         arrays of doubles (not complex).
    eta - viscocity
    Ntotal - total number of points in x,y,z
    Kx, Ky, Kz (Dx..) - wave numbers and Fourier derivative operators 
  
  Returns: None - U_hat_r[i], P_hat_r[i] are written into

  Note: We assume the net force on the unit cell is 0 by *ignoring* 
        the k = 0 mode. That is, the k=0 mode of the output solution
        will be 0.
  */
  void TriplyPeriodicStokes(double* U_hat_r, double* U_hat_i, double* P_hat_r, double* P_hat_i, 
                            const double* fG_hat_r, const double* fG_hat_i, const Complex* Dx, 
                            const Complex* Dy, const Complex* Dz, const double* Ksq, double eta, 
                            unsigned int Ntotal, int n_threads);
}
