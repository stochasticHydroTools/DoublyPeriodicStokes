#include"TPTools.h"
#include<omp.h>

extern "C"
{
  void TriplyPeriodicStokes(double* U_hat_r, double* U_hat_i, double* P_hat_r, double* P_hat_i, 
                            const double* fG_hat_r, const double* fG_hat_i, const Complex* Dx, 
                            const Complex* Dy, const Complex* Dz, const double* Ksq, double eta, 
                            unsigned int Ntotal, int n_threads)
  {
    const double* f_hat_r = &(fG_hat_r[0]);
    const double* g_hat_r = &(fG_hat_r[1]);
    const double* h_hat_r = &(fG_hat_r[2]);
    double* u_hat_r = &(U_hat_r[0]);
    double* v_hat_r = &(U_hat_r[1]);
    double* w_hat_r = &(U_hat_r[2]);
    const double* f_hat_i = &(fG_hat_i[0]);
    const double* g_hat_i = &(fG_hat_i[1]);
    const double* h_hat_i = &(fG_hat_i[2]);
    double* u_hat_i = &(U_hat_i[0]);
    double* v_hat_i = &(U_hat_i[1]);
    double* w_hat_i = &(U_hat_i[2]);
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int i = 3; i < Ntotal; i+=3)
    {
      Complex dx, dy, dz; double ksq;
      dx = Dx[i / 3]; dy = Dy[i / 3]; dz = Dz[i / 3]; ksq = Ksq[i / 3];
      double val = 1.0 / (eta * ksq);
      Complex f_hat = f_hat_r[i] + I * f_hat_i[i];
      Complex g_hat = g_hat_r[i] + I * g_hat_i[i];
      Complex h_hat = h_hat_r[i] + I * h_hat_i[i];
      Complex p_hat = -(dx * f_hat + dy * g_hat + dz * h_hat) / ksq;
      Complex u_hat = val * (f_hat - (dx * p_hat)); 
      Complex v_hat = val * (g_hat - (dy * p_hat)); 
      Complex w_hat = val * (h_hat - (dz * p_hat));
    
      P_hat_r[i / 3] = creal(p_hat);
      u_hat_r[i] = creal(u_hat);      
      v_hat_r[i] = creal(v_hat);      
      w_hat_r[i] = creal(w_hat);      
      P_hat_i[i / 3] = cimag(p_hat);
      u_hat_i[i] = cimag(u_hat);      
      v_hat_i[i] = cimag(v_hat);      
      w_hat_i[i] = cimag(w_hat);      
    }   
  }
}
