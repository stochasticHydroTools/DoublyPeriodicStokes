#include"DPTools.h"
#ifdef INTEL
#include<mkl_blas.h>
#include<mkl_lapack.h>
#else
#include<lapacke.h>
#endif
#include<math.h>
#include"exceptions.h"
#include<sstream>
#include<omp.h>

const int _one = 1;
const int _two = 2;
const char notrans = 'N';
const char trans = 'T';

#ifndef INTEL
typedef lapack_complex_double lpComplex;
#endif

DoublyPeriodicStokes::DoublyPeriodicStokes(const Complex* _Dx, const Complex* _Dy,
                                           const double* _storage_kbp, const Complex* _C,
                                           const Complex* _Ginv, const Complex* _Ainv_B,
                                           const Complex* _FIMat, const Complex* _SIMat,
                                           const double* _K, const double _H, const double _eta,
                                           const unsigned int _Nyx, const unsigned int _Nz,
                                           const unsigned int _dof, const unsigned int _k0, 
                                           const int _n_threads)
  : Cf(0), Cf_r(0), Cf_i(0), Dx(_Dx), Dy(_Dy), storage_kbp(_storage_kbp), C(_C), 
    Ginv(_Ginv), Ainv_B(_Ainv_B), FIMat(_FIMat), SIMat(_SIMat), 
    K(_K), H(_H), eta(_eta), Nyx(_Nyx), Nz(_Nz), dof(_dof), k0(_k0), 
    p_rhs_k0_r(0), p_rhs_k0_i(0), u_rhs_k0_r(0), u_rhs_k0_i(0), fplan(0), 
    n_threads(_n_threads), corr_sol_hat(0), phi_out_i(0), phi_out_r(0)
{ 
  stokes_out = (Complex**) fftw_malloc(Nyx * sizeof(Complex*));
  X = (Complex**) fftw_malloc(Nyx * sizeof(Complex*));
  Y = (Complex**) fftw_malloc(Nyx * sizeof(Complex*));
  Xpenta = (Complex**) fftw_malloc(Nyx * sizeof(Complex*));
  U_BC_RHS = (Complex**) fftw_malloc(Nyx * sizeof(Complex*)); 
  str_df = str_cu = str_urhs = 3 * Nz;
  str_dp = str_cp = str_prhs = Nz;
  str_x = 3 * (Nz + 2);
  str_y = str_ubcrhs = 6;
  str_xp = 3 * Nz;
  
  offset_df = 0;
  offset_prhs = offset_df + str_df;
  offset_cp = offset_prhs + str_prhs;
  offset_dp = offset_cp + str_cp;
  offset_urhs = offset_dp + str_dp;
  offset_cu = offset_urhs + str_urhs;
  Nk = str_df + str_prhs + str_cp + str_dp + str_urhs + str_cu; 
  #pragma omp parallel for num_threads(n_threads)
  for (unsigned int i = 0; i < Nyx; ++i)
  {
    stokes_out[i] = (Complex*) fftw_malloc(Nk * sizeof(Complex));
    X[i] = (Complex*) fftw_malloc(str_x * sizeof(Complex));
    Y[i] = (Complex*) fftw_malloc(str_y * sizeof(Complex));
    Xpenta[i] = (Complex*) fftw_malloc(str_xp * sizeof(Complex));
    U_BC_RHS[i] = (Complex*) fftw_malloc(str_ubcrhs * sizeof(Complex));
  }
  U_hat_r = (double*) fftw_malloc(Nyx * Nz * dof * sizeof(double));
  U_hat_i = (double*) fftw_malloc(Nyx * Nz * dof * sizeof(double));
  P_hat_r = (double*) fftw_malloc(Nyx * Nz * sizeof(double));
  P_hat_i = (double*) fftw_malloc(Nyx * Nz * sizeof(double));
}

inline void DoublyPeriodicStokes::zero_init()
{
  #pragma omp parallel for num_threads(n_threads)
  for (unsigned int i = 0; i < Nyx; ++i)
  {
    memset(stokes_out[i], 0, Nk * sizeof(Complex));
    memset(X[i], 0, str_x * sizeof(Complex));
    memset(Y[i], 0, str_y * sizeof(Complex));
  }
  memset(U_hat_r, 0, sizeof(double) * Nz * dof); 
  memset(U_hat_i, 0, sizeof(double) * Nz * dof); 
  memset(P_hat_r, 0, sizeof(double) * Nz); 
  memset(P_hat_i, 0, sizeof(double) * Nz); 
}


inline void DoublyPeriodicStokes::doublyPeriodicSolve_no_wall()
{
  memset(U_hat_r, 0, sizeof(double) * Nz * dof); 
  memset(U_hat_i, 0, sizeof(double) * Nz * dof); 
  memset(P_hat_r, 0, sizeof(double) * Nz); 
  memset(P_hat_i, 0, sizeof(double) * Nz); 
  #pragma omp parallel for num_threads(n_threads)
  for (unsigned int i = 1; i < Nyx; ++i)
  {
    memset(stokes_out[i], 0, Nk * sizeof(Complex));
    memset(X[i], 0, str_x * sizeof(Complex));
    memset(Y[i], 0, str_y * sizeof(Complex));
    // compute cheb coefficients of z derivatives of each component of f
    chebCoeffDiff(i, 3);            
    // assemble rhs for pressure (Dx * Cf + Dy * Cg + Dh)
    assemble_p_rhs(i);
    // solve for coeffs of pressure and its derivative
    schurSolve_p(i);
    // assemble rhs for velocity
    assemble_u_rhs(i);
    // assemble bc rows for vel solve
    assemble_u_bc_rhs(i);
    // solve for vel coeffs
    schurSolve_u(i);  
    // save split output
    const Complex* Cu = &(stokes_out[i][offset_cu]);
    const Complex* Cp = &(stokes_out[i][offset_cp]);
    double* u_hat_r = &(U_hat_r[dof * Nz * i]);
    double* u_hat_i = &(U_hat_i[dof * Nz * i]);
    double* p_hat_r = &(P_hat_r[Nz * i]);
    double* p_hat_i = &(P_hat_i[Nz * i]);
    #pragma omp simd
    for (unsigned int j = 0; j < dof * Nz; ++j)
    {
      u_hat_r[j] = creal(Cu[j]); u_hat_i[j] = cimag(Cu[j]);
    } 
    #pragma omp simd
    for (unsigned int j = 0; j < Nz; ++j)
    {
      p_hat_r[j] = creal(Cp[j]); p_hat_i[j] = cimag(Cp[j]);
    }
  }
  // assemble p_rhs for k = 0 if requested
  if (k0)
  {
    chebCoeffDiff(0, 3);
    assemble_p_rhs(0);
  }
}

inline void DoublyPeriodicStokes::chebCoeffDiff(unsigned int i, unsigned int _dof)
{
  Complex* _df = &(stokes_out[i][offset_df]);
  //const Complex* _cf = &(Cf[i * str_df]);
  const double* _cf_r = &(Cf_r[i * str_df]);
  const double* _cf_i = &(Cf_i[i * str_df]);
  double fac = (Nz - 1) * 2 / H, fac1;
  //const Complex* cf;
  const double* cf_r;
  const double* cf_i;
  Complex *df1, *df2;
  
  //cf = &(_cf[_dof * (Nz - 1)]); 
  cf_r = &(_cf_r[_dof * (Nz - 1)]); 
  cf_i = &(_cf_i[_dof * (Nz - 1)]); 

  df1 = &(_df[_dof * (Nz - 2)]); 

  //#pragma omp simd
  //for (unsigned int l = 0; l < _dof; ++l){df1[l] = fac * cf[l];}
  #pragma omp simd
  for (unsigned int l = 0; l < _dof; ++l){df1[l] = fac * (cf_r[l] + I * cf_i[l]);}
  for (unsigned int k = 2; k < Nz; ++k)
  {
    df1 = &(_df[_dof * (Nz - k - 1)]);
    df2 = &(_df[_dof * (Nz - k + 1)]);
    //cf =  &(_cf[_dof * (Nz - k)]);
    cf_r =  &(_cf_r[_dof * (Nz - k)]);
    cf_i =  &(_cf_i[_dof * (Nz - k)]);
    fac1 = (Nz - k) * 2 / H;
    //#pragma omp simd
    //for (unsigned int l = 0; l < _dof; ++l){df1[l] = df2[l] + fac1 * cf[l];}
    #pragma omp simd
    for (unsigned int l = 0; l < _dof; ++l){df1[l] = df2[l] + fac1 * (cf_r[l] + I * cf_i[l]);}
  }
  #pragma omp simd
  for (unsigned int l = 0; l < _dof; ++l)
  {
    _df[l] /= 2;
  }
}

//p_RHS = Dx * Cf + Dy * Cg + Dh
inline void DoublyPeriodicStokes::assemble_p_rhs(unsigned int i)
{
  const Complex* df = &(stokes_out[i][offset_df]);
  //const Complex* cf = &(Cf[i * str_df]);
  const double* cf_r = &(Cf_r[i * str_df]);
  const double* cf_i = &(Cf_i[i * str_df]);
  Complex* p_rhs = &(stokes_out[i][offset_prhs]);
  Complex dx = Dx[i], dy = Dy[i];
  unsigned int j = 0;
  #pragma omp simd
  for (unsigned int k = 0; k < Nz; ++k)
  {
    //p_rhs[k] = dx * cf[j] + dy * cf[j + 1] + df[j + 2]; j += 3;
    p_rhs[k] = dx * (cf_r[j] + I * cf_i[j]) + 
               dy * (cf_r[j + 1] + I * cf_i[j + 1]) + 
               df[j + 2]; 
    j += 3;
  }
}

inline void DoublyPeriodicStokes::assemble_u_rhs(unsigned int i)
{
  const Complex* cp = &(stokes_out[i][offset_cp]);
  const Complex* dp = &(stokes_out[i][offset_dp]);
  //const Complex* cf = &(Cf[i * str_df]);
  const double* cf_r = &(Cf_r[i * str_df]);
  const double* cf_i = &(Cf_i[i * str_df]);
  const Complex dx = Dx[i], dy = Dy[i];
  Complex* u_rhs = &(stokes_out[i][offset_urhs]);
  #pragma omp simd
  for (unsigned int j = 0; j < Nz; ++j)
  {
    //u_rhs[dof * j] = (dx * cp[j] - cf[dof * j]) / eta;  
    u_rhs[dof * j] = (dx * cp[j] - (cf_r[dof * j] + I * cf_i[dof * j])) / eta;  
    //u_rhs[1 + dof * j] = (dy * cp[j] - cf[1 + dof * j]) / eta;  
    u_rhs[1 + dof * j] = (dy * cp[j] - (cf_r[1 + dof * j] + I * cf_i[1 + dof * j])) / eta;  
    //u_rhs[2 + dof * j] = (dp[j] - cf[2 + dof * j]) / eta;  
    u_rhs[2 + dof * j] = (dp[j] - (cf_r[2 + dof * j] + I * cf_i[2 + dof * j])) / eta;  
  }
}

inline void DoublyPeriodicStokes::assemble_u_bc_rhs(unsigned int i)
{
  const Complex* cp = &(stokes_out[i][offset_cp]);
  Complex sum, alt_sum; sum = alt_sum = 0;
  #pragma omp simd reduction(+:sum) reduction(+:alt_sum)
  for (unsigned int j = 0; j < Nz; ++j)
  {
    sum += cp[j];
    alt_sum += cp[j] * pow(-1.0,j);
  }
  double fac1 = 2 * eta;
  double fac = fac1 * K[i];
  Complex* bc_rhs = U_BC_RHS[i];
  bc_rhs[0] = -Dx[i] * sum / fac; bc_rhs[1] = Dx[i] * alt_sum / fac;
  bc_rhs[2] = -Dy[i] * sum / fac; bc_rhs[3] = Dy[i] * alt_sum / fac;
  bc_rhs[4] = sum / fac1; bc_rhs[5] = alt_sum / fac1;
}

inline void DoublyPeriodicStokes::schurSolve_p(unsigned int i)
{
  const int nrhs = 1; int Ns = 3 * Nz - 3;
  int Nzp2 = Nz+2; 
 
  const double* stor = &(storage_kbp[i * Ns]);
  const Complex* rhs = &(stokes_out[i][offset_prhs]);
  const Complex* c = &(C[2 * Nz * i]);
  const Complex* ainvb = &(Ainv_B[2 * Nz * i]);
  const Complex* ginv = &(Ginv[4 * i]);

  Complex* xpenta = Xpenta[i];
  Complex* x      = X[i]; 
  Complex* y = Y[i]; 
  Complex* cp = &(stokes_out[i][offset_cp]);    
  Complex* dp = &(stokes_out[i][offset_dp]);    
  // compute xpenta = A^{-1}*rhs
  solve_penta(xpenta, rhs, stor, Nz, 1);
  // compute  y = c * xpenta (2x1)
  Complex tot1, tot2; tot1 = tot2 = 0;
  #pragma omp simd reduction(+:tot1,tot2)
  for (unsigned int j = 0; j < Nz; ++j) 
  { 
    tot1 += xpenta[j] * c[2 * j]; tot2 += xpenta[j] * c[1 + 2 * j]; 
  }
  y[0] = tot1; y[1] = tot2; 
  // compute y = g^{-1} * y  (2 x 1)
  Complex y0 = ginv[0] * y[0] + ginv[2] * y[1];
  Complex y1 = ginv[1] * y[0] + ginv[3] * y[1];
  y[0] = y0; y[1] = y1;
  // compute xpenta - lu^{-1}*b*y and save into x
  #pragma omp simd 
  for (unsigned int j = 0; j < Nz; ++j) 
  {
    x[j] = xpenta[j] - (ainvb[j] * y[0] + ainvb[j + Nz] * y[1]);
  } 
  x[Nz] = y[0]; x[Nz + 1] = y[1];
  // compute sol and its derivative
  for (unsigned int m = 0; m < Nz; ++m)
  {
    tot1 = tot2 = 0;
    #pragma omp simd reduction(+:tot1,tot2)
    for (unsigned int j = 0; j < Nz + 2; ++j)
    {
      tot1 += SIMat[m + j * Nz] * x[j];
      tot2 += FIMat[m + j * Nz] * x[j];
    }
    cp[m] = tot1; dp[m] = tot2;
  }
}

inline void DoublyPeriodicStokes::schurSolve_u(unsigned int i)
{

  int Nzp2 = Nz+2; 
  const int nrhs = 3, Ns = 3 * Nz - 3;
  const double* stor = &(storage_kbp[i * Ns]);
  const Complex* rhs    = &(stokes_out[i][offset_urhs]);
  const Complex* c      = &(C[2 * Nz * i]);
  const Complex* ainvb  = &(Ainv_B[2 * Nz * i]);
  const Complex* ginv   = &(Ginv[4 * i]);
  const Complex* bc_rhs = U_BC_RHS[i];
  Complex* xpenta = Xpenta[i];
  Complex* x = X[i];    
  Complex* y = Y[i]; 
  Complex* cp = &(stokes_out[i][offset_cu]);
  // compute xpenta = A^{-1}*rhs (3 x Nz)
  solve_penta(xpenta, rhs, stor, Nz, 3);
  // compute y = c * xpenta (3 x 2)
  Complex tot1, tot2, tot3; 
  for (unsigned int m = 0; m < 3; ++m)
  {
    tot1 = tot2 = 0;
    #pragma omp simd reduction(+:tot1,tot2)
    for (unsigned int j = 0; j < Nz; ++j) 
    {
      tot1 += xpenta[m + 3 * j] * c[2 * j];
      tot2 += xpenta[m + 3 * j] * c[1 + 2 * j];
    }
    y[m] = tot1; y[m + 3] = tot2; 
  }
  // compute y = g^{-1}*(y - (alpha,beta))  (3 x 2)
  Ginv_y(ginv, y, bc_rhs, 3);
  Ginv_y(ginv, &y[1], &bc_rhs[2], 3);
  Ginv_y(ginv, &y[2], &bc_rhs[4], 3);  
  // compute xpenta - A^{-1}*b*y and save into x (3 x Nz)
  #pragma omp simd
  for (unsigned int j = 0; j < Nz; ++j)
  {
    x[j]              = xpenta[3 * j] - (ainvb[j] * y[0] + ainvb[j + Nz] * y[3]);
    x[j + Nz + 2]     = xpenta[1 + 3 * j] - (ainvb[j] * y[1] + ainvb[j + Nz] * y[4]);
    x[j + 2 * Nz + 4] = xpenta[2 + 3 * j] - (ainvb[j] * y[2] + ainvb[j + Nz] * y[5]);
  }
  x[Nz]         = y[0]; x[Nz + 1]     = y[3];
  x[2 * Nz + 2] = y[1]; x[2 * Nz + 3] = y[4];
  x[3 * Nz + 4] = y[2]; x[3 * Nz + 5] = y[5];
  // compute sol
  for (unsigned int m = 0; m < Nz; ++m)
  {
    tot1 = tot2 = tot3 = 0;
    #pragma omp simd reduction(+:tot1,tot2,tot3)
    for (unsigned int j = 0; j < Nz + 2; ++j)
    {
      tot1 += SIMat[m + j * Nz] * x[j];
      tot2 += SIMat[m + j * Nz] * x[j + Nz + 2];
      tot3 += SIMat[m + j * Nz] * x[j + 2 * Nz + 4];
    }
    cp[m * 3] = tot1; cp[1 + m * 3] = tot2; cp[2 + m * 3] = tot3;
  }
}


inline void DoublyPeriodicStokes::evalTheta(double theta, const int ind)
{
  if (!phi_out_r)
  {
    phi_out_r = (double*) fftw_malloc(2 * Nyx * dof * sizeof(double)); 
    phi_out_i = (double*) fftw_malloc(2 * Nyx * dof * sizeof(double)); 
  }
  memset(phi_out_r + ind * Nyx * dof, 0, Nyx * dof * sizeof(double)); 
  memset(phi_out_i + ind * Nyx * dof, 0, Nyx * dof * sizeof(double)); 
  const double* in_r = U_hat_r;
  const double* in_i = U_hat_i;
  #pragma omp parallel for num_threads(n_threads)
  for (unsigned int i = 0; i < Nyx; ++i)
  {
    double* out_xy_r = &(phi_out_r[dof * i + ind * Nyx * dof]);
    double* out_xy_i = &(phi_out_i[dof * i + ind * Nyx * dof]);
    for (unsigned int j = 0; j < Nz; ++j)
    {
      const double* in_xyz_r = &(in_r[dof * (j + Nz * i)]);
      const double* in_xyz_i = &(in_i[dof * (j + Nz * i)]);
      double alpha = cos(j * theta);
      #pragma omp simd
      for (unsigned int l = 0; l < dof; ++l)
      {
        //out_xy[l] += (in_xyz_r[l] + I * in_xyz_i[l]) * alpha; 
        out_xy_r[l] += in_xyz_r[l] * alpha;
        out_xy_i[l] += in_xyz_i[l] * alpha; 
      }
    }
  }
}

inline void DoublyPeriodicStokes::plan_fft()
{
  if (!this->fplan)
  {
    int n = 2 * Nz - 2;
    fftw_complex* in = (fftw_complex*) fftw_malloc(n * sizeof(fftw_complex));
    fftw_complex* out = in;
    fftw_plan_with_nthreads(1);
    fplan = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_PATIENT);
    if (!fplan) {exitErr("FFTW planning failed");}
    fftw_free(in);
  }
}

inline void DoublyPeriodicStokes::clean()
{
  if (stokes_out)
  {
    for (unsigned int i = 0; i < Nyx; ++i) {fftw_free(stokes_out[i]); stokes_out[i] = 0;}
    fftw_free(stokes_out); stokes_out = 0;
  }
  if (X)
  {
    for (unsigned int i = 0; i < Nyx; ++i) {fftw_free(X[i]); X[i] = 0;}
    fftw_free(X); X = 0;
  }
  if (Y)
  {
    for (unsigned int i = 0; i < Nyx; ++i) {fftw_free(Y[i]); Y[i] = 0;}
    fftw_free(Y); Y = 0;
  }
  if (Xpenta)
  {
    for(unsigned int i = 0; i < Nyx; ++i) {fftw_free(Xpenta[i]); Xpenta[i] = 0;}
    fftw_free(Xpenta); Xpenta = 0;
  }
  if (U_BC_RHS)
  {
    for(unsigned int i = 0; i < Nyx; ++i) {fftw_free(U_BC_RHS[i]); U_BC_RHS[i] = 0;}
    fftw_free(U_BC_RHS); U_BC_RHS = 0;
  }

  if (U_hat_r) {fftw_free(U_hat_r); U_hat_r = 0;}
  if (U_hat_i) {fftw_free(U_hat_i); U_hat_i = 0;}
  if (P_hat_r) {fftw_free(P_hat_r); P_hat_r = 0;}
  if (P_hat_i) {fftw_free(P_hat_i); P_hat_i = 0;}
  if (k0)
  {
    if (p_rhs_k0_r) {fftw_free(p_rhs_k0_r); p_rhs_k0_r = 0;}
    if (p_rhs_k0_i) {fftw_free(p_rhs_k0_i); p_rhs_k0_i = 0;}
    if (u_rhs_k0_r) {fftw_free(u_rhs_k0_r); u_rhs_k0_r = 0;}
    if (u_rhs_k0_i) {fftw_free(u_rhs_k0_i); u_rhs_k0_i = 0;}
  }
  if (fplan){fftw_destroy_plan(fplan);}
  if (corr_sol_hat)
  {
    for (unsigned int i = 0; i < Nyx-1; ++i) {fftw_free(corr_sol_hat[i]); corr_sol_hat[i] = 0;}
    fftw_free(corr_sol_hat); corr_sol_hat = 0;
  }
  if (phi_out_r)
  {
    fftw_free(phi_out_r); phi_out_r = 0;
    fftw_free(phi_out_i); phi_out_i = 0;
  }
}

extern "C"
{
  DoublyPeriodicStokes* Solver(const Complex* _Dx, const Complex* _Dy,
                               const double* _storage_kbp, const Complex* _C,
                               const Complex* _Ginv, const Complex* _Ainv_B,
                               const Complex* _FIMat, const Complex* _SIMat,
                               const double* _K, const double _H, const double _eta,
                               const unsigned int _Nyx, const unsigned int _Nz,
                               const unsigned int _dof, const unsigned int _k0,
                               const int _n_threads)
  {
    DoublyPeriodicStokes* solver = new DoublyPeriodicStokes(_Dx, _Dy, _storage_kbp, _C, 
                                                            _Ginv, _Ainv_B, _FIMat, _SIMat, 
                                                            _K, _H, _eta, _Nyx, _Nz,  
                                                            _dof, _k0, _n_threads);
    return solver;
  }

  void DoublyPeriodicSolve_no_wall(DoublyPeriodicStokes* solver) {solver->doublyPeriodicSolve_no_wall();}
  void DoublyPeriodicSolve_bottom_wall(DoublyPeriodicStokes* solver, 
                                      const double* Kx, const double* Ky,
                                      const double* z, double eta, unsigned int Nyx,
                                      unsigned int Nz, unsigned int dof, int n_threads)
  {
    solver->doublyPeriodicSolve_no_wall();
    evalCorrectionSol_bottomWall(solver, Kx, Ky, z, eta, Nyx, Nz, dof, n_threads);
  }
  void DoublyPeriodicSolve_slit_channel(DoublyPeriodicStokes* solver, 
                                      const double* Kx, const double* Ky,
                                      const double* z, double H, double eta, unsigned int Nyx,
                                      unsigned int Nz, unsigned int dof, int n_threads)
  {
    solver->doublyPeriodicSolve_no_wall();
    evalCorrectionSol_slitChannel(solver, Kx, Ky, z, H, eta, Nyx, Nz, dof, n_threads);
  }



  void ZeroInit(DoublyPeriodicStokes* solver) {solver->zero_init();}
  void SetRHS(DoublyPeriodicStokes* solver, const Complex* _Cf){solver->Cf = _Cf;}
  void SetRHS_split(DoublyPeriodicStokes* solver, const double* _Cf_r, const double* _Cf_i)
  {
    solver->Cf_r = _Cf_r; solver->Cf_i = _Cf_i;
  }
  void Clean(DoublyPeriodicStokes* solver) {if (solver){solver->clean(); delete solver; solver = 0;}}
  double* GetU_real(DoublyPeriodicStokes* solver) {return solver->U_hat_r;}
  double* GetU_imag(DoublyPeriodicStokes* solver) {return solver->U_hat_i;}
  double* GetP_real(DoublyPeriodicStokes* solver) {return solver->P_hat_r;}
  double* GetP_imag(DoublyPeriodicStokes* solver) {return solver->P_hat_i;}
  double* GetPhi_out_real(DoublyPeriodicStokes* solver) {return solver->phi_out_r;}
  double* GetPhi_out_imag(DoublyPeriodicStokes* solver) {return solver->phi_out_i;}
  
  double* GetP_RHS_real(DoublyPeriodicStokes* solver, unsigned int k)
  {
    if(not solver->p_rhs_k0_r)
    {
      solver->p_rhs_k0_r = (double*) fftw_malloc(solver->str_prhs * sizeof(double));
    }
    const Complex* p_rhs_k = &(solver->stokes_out[k][solver->offset_prhs]);
    #pragma omp simd
    for (unsigned int j = 0; j < solver->str_prhs; ++j)
    {
      solver->p_rhs_k0_r[j] = creal(p_rhs_k[j]);
    }
    return solver->p_rhs_k0_r;
  }   
  double* GetP_RHS_imag(DoublyPeriodicStokes* solver, unsigned int k)
  {
    if(not solver->p_rhs_k0_i)
    {
      solver->p_rhs_k0_i = (double*) fftw_malloc(solver->str_prhs * sizeof(double));
    }
    const Complex* p_rhs_k = &(solver->stokes_out[k][solver->offset_prhs]);
    #pragma omp simd
    for (unsigned int j = 0; j < solver->str_prhs; ++j)
    {
      solver->p_rhs_k0_i[j] = cimag(p_rhs_k[j]);
    }
    return solver->p_rhs_k0_i;
  }
  double* GetU_RHS_real(DoublyPeriodicStokes* solver, unsigned int k, bool assemble)
  {
    if(not solver->u_rhs_k0_r)
    {
      solver->u_rhs_k0_r = (double*) fftw_malloc(solver->str_urhs * sizeof(double));
    }
    if (assemble) {solver->assemble_u_rhs(k);}
    const Complex* u_rhs_k = &(solver->stokes_out[k][solver->offset_urhs]);
    #pragma omp simd
    for (unsigned int j = 0; j < solver->str_urhs; ++j)
    {
      solver->u_rhs_k0_r[j] = creal(u_rhs_k[j]);
    }
    return solver->u_rhs_k0_r;
  }   
  double* GetU_RHS_imag(DoublyPeriodicStokes* solver, unsigned int k, bool assemble)
  {
    if(not solver->u_rhs_k0_i)
    {
      solver->u_rhs_k0_i = (double*) fftw_malloc(solver->str_urhs * sizeof(double));
    }
    if (assemble) {solver->assemble_u_rhs(k);}
    const Complex* u_rhs_k = &(solver->stokes_out[k][solver->offset_urhs]);
    #pragma omp simd
    for (unsigned int j = 0; j < solver->str_urhs; ++j)
    {
      solver->u_rhs_k0_i[j] = cimag(u_rhs_k[j]);
    }
    return solver->u_rhs_k0_i;
  }
  
  void SetP(DoublyPeriodicStokes* solver, const Complex* cp,  unsigned int k)    
  {
    Complex* cp_k = &(solver->stokes_out[k][solver->offset_cp]);
    #pragma omp simd
    for (unsigned int j = 0; j < solver->str_cp; ++j)
    {
      cp_k[j] = cp[j];
    }
  } 
  void SetdP(DoublyPeriodicStokes* solver, const Complex* dp, unsigned int k)
  {
    Complex* dp_k = &(solver->stokes_out[k][solver->offset_dp]);
    #pragma omp simd
    for (unsigned int j = 0; j < solver->str_dp; ++j)
    {
      dp_k[j] = dp[j];
    }
  }
  void SetU(DoublyPeriodicStokes* solver, const Complex* cu,  unsigned int k)    
  {
    Complex* cu_k = &(solver->stokes_out[k][solver->offset_cu]);
    #pragma omp simd
    for (unsigned int j = 0; j < solver->str_cu; ++j)
    {
      cu_k[j] = cu[j];
    }
  } 

  void evalTheta(DoublyPeriodicStokes* solver, double theta)
  {
    solver->evalTheta(theta, 0);
  }

  void chebTransform(const Complex* in, Complex* out, const fftw_plan plan, unsigned int N)
  {
    unsigned int ext = 2 * N - 2;
    Complex* in_ext = (Complex*) fftw_malloc(ext * sizeof(Complex));
    #pragma omp simd
    for (unsigned int i = 0; i < N; ++i)
    {
      in_ext[i] = in[i];
    }
    #pragma omp simd
    for (unsigned int i = N; i < ext; ++i)
    {
      in_ext[i] = in[ext - i];
    }
    fftw_execute_dft(plan, reinterpret_cast<fftw_complex*>(in_ext), 
                     reinterpret_cast<fftw_complex*>(in_ext));
    out[0] = in_ext[0] / ext;
    #pragma omp simd
    for (unsigned int i = 1; i < N - 1; ++i)
    {
      out[i] = (in_ext[i] + in_ext[ext - i]) / ext;
    }
    out[N-1] = in_ext[N-1] / ext; 
    fftw_free(in_ext);
  }

  void evalCorrectionSol_bottomWall(DoublyPeriodicStokes* solver, 
                                    const double* Kx, const double* Ky,
                                    const double* z, double eta, unsigned int Nyx,
                                    unsigned int Nz, unsigned int dof, int n_threads)
  {
    // get velocities at bottom wall for BCs of correction problem
    solver->evalTheta(M_PI, 0); 
    // plan nyx-1 ffts of length 2nz-2, or use wisdom (no op if plan exists)
    solver->plan_fft(); 
    int n = 2 * Nz - 2;
    if (!solver->corr_sol_hat)
    {
      solver->corr_sol_hat = (Complex**) fftw_malloc((Nyx - 1) * sizeof(Complex*));
      #pragma omp parallel for num_threads(n_threads)
      for (unsigned int i = 0; i < Nyx-1; ++i)
      {
        solver->corr_sol_hat[i] = (Complex*) fftw_malloc(n * 4 * sizeof(Complex));     
      }
    }
    #pragma omp parallel for num_threads(n_threads) 
    for (unsigned int i = 1; i < Nyx; ++i)
    {
      // real/imaginary components of correction to pressure and vel
      Complex* Cp = &(solver->corr_sol_hat[i-1][0]); 
      Complex* Cu = &(solver->corr_sol_hat[i-1][n]);
      Complex* Cv = &(solver->corr_sol_hat[i-1][2 * n]); 
      Complex* Cw = &(solver->corr_sol_hat[i-1][3 * n]);
      
      unsigned int ix = dof * i, iy = ix + 1, iz = iy + 1;
      double kx = Kx[i], ky = Ky[i], k = solver->K[i];
      double enkz, zenkz;
      Complex fhat_x, fhat_y, fhat_z;
      fhat_x = -(solver->phi_out_r[ix] + I * solver->phi_out_i[ix]);
      fhat_y = -(solver->phi_out_r[iy] + I * solver->phi_out_i[iy]);
      fhat_z = -(solver->phi_out_r[iz] + I * solver->phi_out_i[iz]);
      Complex alpha = kx * fhat_x + ky * fhat_y + I * k * fhat_z;
      // compute corrections
      #pragma omp simd
      for (unsigned int j = 0; j < Nz; ++j)
      {
        enkz = exp(-k * z[j]); zenkz = z[j] * enkz;
        // pressure correction
        Cp[j] = -I * 2 * eta * alpha * enkz;
        // x vel correction
        Cu[j] = -kx * alpha * zenkz / k + fhat_x * enkz;
        // y vel correction
        Cv[j] = -ky * alpha * zenkz / k + fhat_y * enkz;
        // z vel correction
        Cw[j] = -I * alpha * zenkz + fhat_z * enkz;
      }
      // reflection for cheb transform
      #pragma omp simd
      for (unsigned int j = Nz; j < n; ++j)
      {
        Cp[j] = Cp[n - j];
        Cu[j] = Cu[n - j];
        Cv[j] = Cv[n - j];
        Cw[j] = Cw[n - j];
      }
      // forward transform in z to get cheb coeffs
      fftw_execute_dft(solver->fplan, reinterpret_cast<fftw_complex*>(Cp), 
                       reinterpret_cast<fftw_complex*>(Cp));
      fftw_execute_dft(solver->fplan, reinterpret_cast<fftw_complex*>(Cu), 
                       reinterpret_cast<fftw_complex*>(Cu));
      fftw_execute_dft(solver->fplan, reinterpret_cast<fftw_complex*>(Cv), 
                       reinterpret_cast<fftw_complex*>(Cv));
      fftw_execute_dft(solver->fplan, reinterpret_cast<fftw_complex*>(Cw), 
                     reinterpret_cast<fftw_complex*>(Cw));
      // truncate and copy output
      unsigned int offset = i * Nz;
      double* Cp_out_r = &(solver->P_hat_r[offset]);
      double* Cu_out_r = &(solver->U_hat_r[dof * offset]);
      double* Cv_out_r = &(solver->U_hat_r[1 + dof * offset]);
      double* Cw_out_r = &(solver->U_hat_r[2 + dof * offset]);
      double* Cp_out_i = &(solver->P_hat_i[offset]);
      double* Cu_out_i = &(solver->U_hat_i[dof * offset]);
      double* Cv_out_i = &(solver->U_hat_i[1 + dof * offset]);
      double* Cw_out_i = &(solver->U_hat_i[2 + dof * offset]);
      Cp_out_r[0] += creal(Cp[0]) / n;
      Cu_out_r[0] += creal(Cu[0]) / n;
      Cv_out_r[0] += creal(Cv[0]) / n;
      Cw_out_r[0] += creal(Cw[0]) / n;
      Cp_out_i[0] += cimag(Cp[0]) / n;
      Cu_out_i[0] += cimag(Cu[0]) / n;
      Cv_out_i[0] += cimag(Cv[0]) / n;
      Cw_out_i[0] += cimag(Cw[0]) / n;
      #pragma omp simd
      for (unsigned int j = 1; j < Nz - 1; ++j)
      {
        Cp_out_r[j]       += creal(Cp[j] + Cp[n - j]) / n;
        Cu_out_r[dof * j] += creal(Cu[j] + Cu[n - j]) / n;
        Cv_out_r[dof * j] += creal(Cv[j] + Cv[n - j]) / n;
        Cw_out_r[dof * j] += creal(Cw[j] + Cw[n - j]) / n;
        Cp_out_i[j]       += cimag(Cp[j] + Cp[n - j]) / n;
        Cu_out_i[dof * j] += cimag(Cu[j] + Cu[n - j]) / n;
        Cv_out_i[dof * j] += cimag(Cv[j] + Cv[n - j]) / n;
        Cw_out_i[dof * j] += cimag(Cw[j] + Cw[n - j]) / n;
      }
      Cp_out_r[Nz-1]         += creal(Cp[Nz-1]) / n; 
      Cu_out_r[dof * (Nz-1)] += creal(Cu[Nz-1]) / n; 
      Cv_out_r[dof * (Nz-1)] += creal(Cv[Nz-1]) / n; 
      Cw_out_r[dof * (Nz-1)] += creal(Cw[Nz-1]) / n; 
      Cp_out_i[Nz-1]         += cimag(Cp[Nz-1]) / n; 
      Cu_out_i[dof * (Nz-1)] += cimag(Cu[Nz-1]) / n; 
      Cv_out_i[dof * (Nz-1)] += cimag(Cv[Nz-1]) / n; 
      Cw_out_i[dof * (Nz-1)] += cimag(Cw[Nz-1]) / n; 
    }
  }

  inline unsigned int at(unsigned int i, unsigned int j){return i + 8 * j;}

  void evalCorrectionSol_slitChannel(DoublyPeriodicStokes* solver,
                                     const double* Kx, const double* Ky, 
                                     const double* z, double H, double eta, unsigned int Nyx, 
                                     unsigned int Nz, unsigned int dof, int n_threads)
  {
    // get velocities at walls for BCs of correction problem
    solver->evalTheta(M_PI, 0); 
    solver->evalTheta(0, 1); 
    // plan nyx-1 ffts of length z, or use wisdom (no op if plan exists)
    solver->plan_fft(); 
    double fac = 1.0 / 2.0 / eta;
    int one, eight; one = 1; eight = 8; int info; 
    int n = 2 * Nz - 2;
    if (!solver->corr_sol_hat)
    {
      solver->corr_sol_hat = (Complex**) fftw_malloc((Nyx - 1) * sizeof(Complex*));
      #pragma omp parallel for num_threads(n_threads)
      for (unsigned int i = 0; i < Nyx-1; ++i)
      {
        solver->corr_sol_hat[i] = (Complex*) fftw_malloc(n * 4 * sizeof(Complex));     
      }
    }
    #pragma omp parallel for num_threads(n_threads)
    for (unsigned int i = 1; i < Nyx; ++i)
    {
      unsigned int ix = dof * i, iy = ix + 1, iz = iy + 1;
      double kx = Kx[i], ky = Ky[i], k = sqrt(kx * kx + ky * ky);
      // pivot storage for lapack
      int* piv = (int*) fftw_malloc(8 * sizeof(int));
      // matrix to solve for coeffs of exponentials in correction sol
      Complex* coeffA = (Complex*) fftw_malloc(64 * sizeof(Complex));
      // right hand side of coeffA_r c_r = x_r (real and complex parts)
      Complex* x = (Complex*) fftw_malloc(8 * sizeof(Complex));
      #pragma omp simd
      for (unsigned int j = 0; j < 64; ++j) {coeffA[j] = 0;} 

      double enkh = exp(-k * H);
      coeffA[at(0,0)] = fac;
      coeffA[at(0,1)] = enkh * fac;
      coeffA[at(0,6)] = -k;
      coeffA[at(0,7)] = k * enkh;
      coeffA[at(1,0)] = (1 - k * H) * enkh * fac;
      coeffA[at(1,1)] = (1 + k * H) * fac;
      coeffA[at(1,6)] = -k * enkh;
      coeffA[at(1,7)] = k;
      coeffA[at(2,2)] = 1;
      coeffA[at(2,3)] = enkh;
      coeffA[at(3,4)] = 1;
      coeffA[at(3,5)] = enkh;
      coeffA[at(4,6)] = 1;
      coeffA[at(4,7)] = enkh;
      coeffA[at(5,0)] = -I * kx * H * enkh * fac / k;
      coeffA[at(5,1)] = I * kx * H * fac / k;
      coeffA[at(5,2)] = enkh;
      coeffA[at(5,3)] = 1;
      coeffA[at(6,0)] = -I * ky * H * enkh * fac / k;
      coeffA[at(6,1)] = I * ky * H * fac / k;
      coeffA[at(6,4)] = enkh;
      coeffA[at(6,5)] = 1;
      coeffA[at(7,0)] = H * enkh * fac;
      coeffA[at(7,1)] = H * fac;
      coeffA[at(7,6)] = enkh;
      coeffA[at(7,7)] = 1; 

      Complex fbhat_x = -(solver->phi_out_r[ix] + I * solver->phi_out_i[ix]);
      Complex fbhat_y = -(solver->phi_out_r[iy] + I * solver->phi_out_i[iy]);
      Complex fbhat_z = -(solver->phi_out_r[iz] + I * solver->phi_out_i[iz]);
      Complex fthat_x = -(solver->phi_out_r[ix + Nyx * dof] + I * solver->phi_out_i[ix + Nyx * dof]);
      Complex fthat_y = -(solver->phi_out_r[iy + Nyx * dof] + I * solver->phi_out_i[iy + Nyx * dof]);
      Complex fthat_z = -(solver->phi_out_r[iz + Nyx * dof] + I * solver->phi_out_i[iz + Nyx * dof]);
  
      x[0] = -I * kx * fbhat_x - I * ky * fbhat_y;
      x[1] = -I * kx * fthat_x - I * ky * fthat_y;
      x[2] = fbhat_x; x[3] = fbhat_y; 
      x[4] = fbhat_z; x[5] = fthat_x; 
      x[6] = fthat_y; x[7] = fthat_z; 
      // solve for coefficients of exponentials   
      #ifdef INTEL
      zgesv(&eight, &one, coeffA, &eight, piv, x, &eight, &info);      
      #else
      LAPACKE_zgesv(LAPACK_COL_MAJOR, 8, 1, reinterpret_cast<lapack_complex_double*>(coeffA), 
                    8, piv, reinterpret_cast<lapack_complex_double*>(x), 8);            
      #endif
      // real/imaginary components of correction to pressure and vel
      Complex* Cp = &(solver->corr_sol_hat[i-1][0]); 
      Complex* Cu = &(solver->corr_sol_hat[i-1][n]);
      Complex* Cv = &(solver->corr_sol_hat[i-1][2 * n]); 
      Complex* Cw = &(solver->corr_sol_hat[i-1][3 * n]);
      
      // e^(-kz), e^(k(z-H)), ze^(k(z-H)), e^(-kH), ze^(-kz)
      double enkz, ekzmh, zekzmh, zenkz;
      // compute correction 
      #pragma omp simd
      for (unsigned int j = 0; j < Nz; ++j) 
      {
        enkz = exp(-k * z[j]); ekzmh = enkh / enkz;
        zekzmh = z[j] * ekzmh; zenkz = z[j] * enkz;
        Cp[j] = x[0] * enkz + x[1] * ekzmh;
        // x vel correction
        Cu[j] = -I * x[0] * kx * zenkz * fac / k +
                I * x[1] * kx * zekzmh * fac / k +
                x[2] * enkz + x[3] * ekzmh;
        // y vel correction
        Cv[j] = -I * x[0] * ky * zenkz * fac / k +
                I * x[1] * ky * zekzmh * fac / k + 
                x[4] * enkz + x[5] * ekzmh;
        // z vel correction
        Cw[j] = x[0] * zenkz * fac + x[1] * zekzmh * fac +
                x[6] * enkz + x[7] * ekzmh;
      }
      // reflection for cheb transform
      #pragma omp simd
      for (unsigned int j = Nz; j < n; ++j)
      {
        Cp[j] = Cp[n - j];
        Cu[j] = Cu[n - j];
        Cv[j] = Cv[n - j];
        Cw[j] = Cw[n - j];
      }
      // forward transform in z to get cheb coeffs
      fftw_execute_dft(solver->fplan, reinterpret_cast<fftw_complex*>(Cp), 
                       reinterpret_cast<fftw_complex*>(Cp));
      fftw_execute_dft(solver->fplan, reinterpret_cast<fftw_complex*>(Cu), 
                       reinterpret_cast<fftw_complex*>(Cu));
      fftw_execute_dft(solver->fplan, reinterpret_cast<fftw_complex*>(Cv), 
                       reinterpret_cast<fftw_complex*>(Cv));
      fftw_execute_dft(solver->fplan, reinterpret_cast<fftw_complex*>(Cw), 
                     reinterpret_cast<fftw_complex*>(Cw));
      // truncate and copy output
      unsigned int offset = i * Nz;
      double* Cp_out_r = &(solver->P_hat_r[offset]);
      double* Cu_out_r = &(solver->U_hat_r[dof * offset]);
      double* Cv_out_r = &(solver->U_hat_r[1 + dof * offset]);
      double* Cw_out_r = &(solver->U_hat_r[2 + dof * offset]);
      double* Cp_out_i = &(solver->P_hat_i[offset]);
      double* Cu_out_i = &(solver->U_hat_i[dof * offset]);
      double* Cv_out_i = &(solver->U_hat_i[1 + dof * offset]);
      double* Cw_out_i = &(solver->U_hat_i[2 + dof * offset]);
      Cp_out_r[0] += creal(Cp[0]) / n;
      Cu_out_r[0] += creal(Cu[0]) / n;
      Cv_out_r[0] += creal(Cv[0]) / n;
      Cw_out_r[0] += creal(Cw[0]) / n;
      Cp_out_i[0] += cimag(Cp[0]) / n;
      Cu_out_i[0] += cimag(Cu[0]) / n;
      Cv_out_i[0] += cimag(Cv[0]) / n;
      Cw_out_i[0] += cimag(Cw[0]) / n;
      #pragma omp simd
      for (unsigned int j = 1; j < Nz - 1; ++j)
      {
        Cp_out_r[j]       += creal(Cp[j] + Cp[n - j]) / n;
        Cu_out_r[dof * j] += creal(Cu[j] + Cu[n - j]) / n;
        Cv_out_r[dof * j] += creal(Cv[j] + Cv[n - j]) / n;
        Cw_out_r[dof * j] += creal(Cw[j] + Cw[n - j]) / n;
        Cp_out_i[j]       += cimag(Cp[j] + Cp[n - j]) / n;
        Cu_out_i[dof * j] += cimag(Cu[j] + Cu[n - j]) / n;
        Cv_out_i[dof * j] += cimag(Cv[j] + Cv[n - j]) / n;
        Cw_out_i[dof * j] += cimag(Cw[j] + Cw[n - j]) / n;
      }
      Cp_out_r[Nz-1]         += creal(Cp[Nz-1]) / n; 
      Cu_out_r[dof * (Nz-1)] += creal(Cu[Nz-1]) / n; 
      Cv_out_r[dof * (Nz-1)] += creal(Cv[Nz-1]) / n; 
      Cw_out_r[dof * (Nz-1)] += creal(Cw[Nz-1]) / n; 
      Cp_out_i[Nz-1]         += cimag(Cp[Nz-1]) / n; 
      Cu_out_i[dof * (Nz-1)] += cimag(Cu[Nz-1]) / n; 
      Cv_out_i[dof * (Nz-1)] += cimag(Cv[Nz-1]) / n; 
      Cw_out_i[dof * (Nz-1)] += cimag(Cw[Nz-1]) / n; 
    }
  }
}
