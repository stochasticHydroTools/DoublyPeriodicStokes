#ifndef _DPTOOLS_H
#define _DPTOOLS_H
#include<complex.h> 
#include<fftw3.h>
#include<LinearSolvers.h>

// set mkl complex type
#define MKL_Complex16 double _Complex
// alias for c complex
typedef double _Complex Complex;

/* struct for data and methods of the DP no wall solver */
struct DoublyPeriodicStokes
{
  /* inputs - Cf - coeffs of forces, must be set externally by SetRHS in wrapper
            - Dx, Dy - wave numbers multiplied with i, and handling unmatched mode
            - C, Ginv, Ainv_B - components used in lin solve at each k
            - FIMat, SIMat - first and second cheb integral mappings
            - storage_kbp - precompute for kbpenta lin solver
            - K - sqrt(kx^2+ky^2)
            - H - Lz / 2
            - eta - viscosity
            - Nxy, Nz, dof - combined pts in xy, pts in z, degrees of freedom (3)
            - n_threads is num_threads for omp
  */       
  const Complex *Cf, *Dx, *Dy, *C, *Ginv, *Ainv_B, *FIMat, *SIMat;
  const double *storage_kbp, *K, *Cf_r, *Cf_i; double H, eta;

  unsigned int Nyx, Nz, dof, n_threads;
  /* container for cache aligned outputs per wave number
     stokes_out[i] has all output at wave number i, ordered as
     Df, P_RHS, Cp, Dp, U_RHS, Cu 
    
     See offsets in ctor for details
  */
  Complex** stokes_out;
  // containers for cache aligned intermediates per wave number
  Complex **X, **Y, **Xpenta, **U_BC_RHS;
  // external outputs - velocity and pressure coeffs split into real/imag
  double *U_hat_r, *U_hat_i, *P_hat_r, *P_hat_i;
  // forward transform plan for wall corrections
  fftw_plan fplan;
  // storage for correction sol ffts
  Complex** corr_sol_hat;
  // storage for evaltheta()
  double *phi_out_r, *phi_out_i;

  // strides for output at each k
  unsigned int str_df, str_cu, str_urhs;
  unsigned int str_cp, str_dp, str_prhs;
  unsigned int str_x, str_y, str_xp, str_ubcrhs, Nk; 
  // offsets for output at each k
  unsigned int offset_df, offset_prhs, offset_cp;
  unsigned int offset_dp, offset_urhs, offset_cu;
  // switch to compute terms for k0 
  // - k0 = 0 will do nothing
  // - k0 = 1 will populate the P_RHS part of stokes_out[0] (i.e. for k=0)
  // Note, treatment of solutions for k=0 is handled externally with the same switch
  // The switch here just prevents unecessary computation.
  unsigned int k0;
  // storage for outputs at k0 if requested
  double *p_rhs_k0_r, *p_rhs_k0_i, *u_rhs_k0_r, *u_rhs_k0_i;

  // read inputs, allocate buffers, set strides and offsets
  DoublyPeriodicStokes(const Complex* _Dx, const Complex* _Dy,
                       const double* _storage_kbp, const Complex* _C,
                       const Complex* _Ginv, const Complex* _Ainv_B,
                       const Complex* _FIMat, const Complex* _SIMat,
                       const double* _K, const double _H, const double _eta,
                       const unsigned int _Nyx, const unsigned int _Nz,
                       const unsigned int _dof, const unsigned int _k0, 
                       const int _n_threads);

  // fill writeable arrays with 0
  void zero_init();
  // create fftw forward plan
  void plan_fft();

  // solve the dp no wall problem
  void doublyPeriodicSolve_no_wall();
  // differentiate 1D Chebyshev series at wave number i 
  void chebCoeffDiff(unsigned int i, unsigned int _dof); 
  // compute p_RHS = Dx * Cf + Dy * Cg + Dh for wave number i
  void assemble_p_rhs(unsigned int i);
  /********************************************************************************

    Given a block linear system of the form 
    
    |A B||a |   |f|
    |C D||c0| = |0| 
         |d0|   |0|,  
    
    with the Schur complement of A : G = (C A^{-1} B - D),
                  
    this function takes the output from precomputeLinOps)
    and computes
      - (c0,d0) = inv(G)*(C*(inv(A)*f))
      - a = inv(A)*f - inv(A)*B*(c0,d0)
      - Dp = I1 * (a, c0, d0)
      - Cp = I2 * (a, c0, d0),

    where I1, I2 are the Chebyshev first and second integral mappings, and Cp, Dp are
    the Fourier-Chebysehv coefficients of pressure and its derivative (DP no walls).
    f is the RHS for the pressure equation
  */
  void schurSolve_p(unsigned int i);
  /* compute u_RHS =  (Dx * Cp - Cf,
                       Dy * Cp - Cg,
                       Dp - Ch) / eta  
  */
  void assemble_u_rhs(unsigned int i);
  // compute last 6 rows for velocity bvp (2 for each component)
  void assemble_u_bc_rhs(unsigned int i);
  /************************************************************************************
    Given a block linear system of the form 
    
    |A B||a |   |f|
    |C D||c0| = |alpha| 
         |d0|   |beta|,  
    
    with the Schur complement of A : G = (C A^{-1} B - D),
                  
    this function takes the output from precomputeLinOps
    and computes
      - (c0,d0) = inv(G)*[(C*(inv(A)*f)) - (alpha,beta)] 
      - a = inv(A)*f - inv(A)*B*(c0,d0)
      - Cu = I2 * (a, c0, d0),

    where I2 is the Chebyshev second integral mapping, and Cu has
    the Fourier-Chebysehv coefficients of a velocity component(DP no walls).
    f is the RHS for the velocity equation, and (alpha,beta) are the bc_rows
  */
  void schurSolve_u(unsigned int i);

  /*
    This function is used to evaluate a Chebyshev series at a given
    value of theta (point on the cheb grid)
    
    Parameters: 

      theta - determines the slice in z
            - eg) theta = pi is z = 0, theta = 0 is z = H 
      
    Side effects: overwrites
      phi_out_r/i - the output array (size (Nyx * dof, 1))
                  - these are the Fourier-Chebyshev coeffs on the x-y
                   plane at a given z value 
  */
  void evalTheta(double theta, const int ind=0);

  // free memory
  void clean();
};


// C wrapper for DoublyPeriodicStokes struct
extern "C"
{
  // initialize the dp no wall solver
  DoublyPeriodicStokes* Solver(const Complex* _Dx, const Complex* _Dy,
                               const double* _storage_kbp, const Complex* _C,
                               const Complex* _Ginv, const Complex* _Ainv_B,
                               const Complex* _FIMat, const Complex* _SIMat,
                               const double* _K, const double _H, const double _eta,
                               const unsigned int _Nyx, const unsigned int _Nz,
                               const unsigned int _dof, const unsigned int _k0,
                               const int _n_threads);
  // execute the dp no wall solver
  void DoublyPeriodicSolve_no_wall(DoublyPeriodicStokes* solver);
  // call zero_init method of solver
  void ZeroInit(DoublyPeriodicStokes* solver);
  // set the RHS for the solver - this must be done before calling solve()
  void SetRHS(DoublyPeriodicStokes* solver, const Complex* _Cf);
  // free memory
  void Clean(DoublyPeriodicStokes* solver);
  // get real/imag part of pressure/velocity coeffs
  double* GetU_real(DoublyPeriodicStokes* solver);
  double* GetU_imag(DoublyPeriodicStokes* solver);
  double* GetP_real(DoublyPeriodicStokes* solver);
  double* GetP_imag(DoublyPeriodicStokes* solver);
  double* GetP_RHS_real(DoublyPeriodicStokes* solver, unsigned int k);
  double* GetP_RHS_imag(DoublyPeriodicStokes* solver, unsigned int k);
  // compute u_rhs at k only if assemble = true
  // eg. first call GetU_RHS_real(s, k, true) and then GetU_RHS_imag(s, k, false)
  double* GetU_RHS_real(DoublyPeriodicStokes* solver, unsigned int k, bool assemble);
  double* GetU_RHS_imag(DoublyPeriodicStokes* solver, unsigned int k, bool assemble);
  void SetP(DoublyPeriodicStokes* solver, const Complex* cp,  unsigned int k);    
  void SetdP(DoublyPeriodicStokes* solver, const Complex* dp, unsigned int k);    
  void SetU(DoublyPeriodicStokes* solver, const Complex* cu,  unsigned int k); 
  double* GetPhi_out_real(DoublyPeriodicStokes* solver);
  double* GetPhi_out_imag(DoublyPeriodicStokes* solver);


  /*
    This function is used to evaluate a Chebyshev series at a given
    value of theta (point on the cheb grid)
    
    Parameters:
      solver - ptr to solver struct 
      theta - determines the slice in z
            - eg) theta = pi is z = 0, theta = 0 is z = H 
    Side effecits: overwrites
      phi_out_r/i - the output array (size (Nyx * dof, 1))
                  - these are the Fourier-Chebyshev coeffs on the x-y
                    plane at a given z value 
      Access by calling GetPhi_out_real() for eg.
  */
  void evalTheta(DoublyPeriodicStokes* solver, double theta);

  
  /*
    Evaluate the analytical correction to the DP solve to enforce no-slip 
    BCs at the bottom wall
 
    Parameters: 
      C(p,u,v,w)corr       - Fourier-Cheb coeffs of correction sol for  
                             pressure and velocity
                           - these are 0 at time of call
                           - overwritten during exectution
                           - the k = 0 element remains 0
      fhat - Fourier-Cheb coeffs at bottom wall
      z - Chebyshev points in z
      Kx, Ky - tensor prod of wave numbers in x,y
      eta - viscosity
      Nyx, Nz - Nyx = Ny * Nx for Ny,Nx points in x,y and Nz points in z
      dof - degrees of freedom   
      nthreads - num threads for omp loops
  */
  void evalCorrectionSol_bottomWall(DoublyPeriodicStokes* solver, 
                                    const double* Kx, const double* Ky,
                                    const double* z, double eta, unsigned int Nyx,
                                    unsigned int Nz, unsigned int dof, int n_threads);


  
  /*  
    Evaluate the analytical correction to the DP solve to enforce no-slip 
    BCs at the bottom and top wall
 
    Parameters: 
      C(p,u,v,w)  - Fourier-Cheb coeffs of correction sol for  
                    pressure and velocity
                  - these are 0 at time of call
                  - and overwritten during exectution
                  - the k = 0 element remains 0
      fbha  - Fourier-Cheb coeffs at bottom wall
      fthat - Fourier-Cheb coeffs at top wall
      z - Chebyshev points in z
      Kx, Ky - tensor prod of wave numbers in x,y
      Lz - extent of z grid
      eta - viscosity
      Nyx, Nz - Nyx = Ny * Nx for Ny,Nx points in x,y and Nz points in z
      dof - degrees of freedom   
      nthreads - num threads for omp loops
  */
  void evalCorrectionSol_slitChannel(DoublyPeriodicStokes* solver,
                                     const double* Kx, const double* Ky, const double* z, 
                                     double H, double eta, unsigned int Nyx, unsigned int Nz, 
                                     unsigned int dof, int n_threads);

  // flattened index into 4D array
  inline unsigned int at4D(unsigned int l, unsigned int i, unsigned int j,
                                 unsigned int k, const unsigned int Nl,
                                 const unsigned int Ni, const unsigned int Nj)
  {
    return l + Nl * (i + Ni * (j + Nj * k));
  }
  inline unsigned int at3D(unsigned int l, unsigned int ij, unsigned int k, 
                                 const unsigned int Nl, const unsigned int Nij)
  {
    return l + Nl * (ij + Nij * k);
  }
  
  inline unsigned int at2D(unsigned int i, unsigned int j, 
                                 const unsigned int Ni)
  {
    return i + Ni * j;
  }

  inline void mmul_transX(const Complex* A, const Complex* X, Complex* B, 
                   const unsigned int M, const unsigned int N, const unsigned int K)
  {
    Complex *b, x; const Complex *a;
    for (unsigned int k = 0; k < K; ++k)
    {
      b = &B[k * M]; 
      for (unsigned int n = 0; n < N; ++n)
      {
        a = &A[n * M]; x = X[k + n * K];
        #pragma omp simd
        for (unsigned int m = 0; m < M; ++m)
        {
          b[m] += a[m] * x; 
        }
      } 
    }
  }
  
  inline void mmul_no_transX(const Complex* A, const Complex* X, Complex* B, 
                             const unsigned int M, const unsigned int N, const unsigned int K)
  {
    Complex *b, x; const Complex *a;
    for (unsigned int k = 0; k < K; ++k)
    {
      b = &B[k * M]; 
      for (unsigned int n = 0; n < N; ++n)
      {
        a = &A[n * M]; x = X[n + k * N];
        #pragma omp simd
        for (unsigned int m = 0; m < M; ++m)
        {
          b[m] += a[m] * x; 
        }
      } 
    }
  }

  // 2 x 2 mat vec with some tweaks for pressure/velocity DP solve
  inline void Ginv_y(const Complex* ginv, Complex* y, const Complex* bc_rhs, unsigned int str)
  {
      Complex y0 = ginv[0] * (y[0] - bc_rhs[0]) + ginv[2] * (y[str] - bc_rhs[1]);
      Complex y1 = ginv[1] * (y[0] - bc_rhs[0]) + ginv[3] * (y[str] - bc_rhs[1]);
      y[0] = y0; y[str] = y1;
  }

}
#endif
