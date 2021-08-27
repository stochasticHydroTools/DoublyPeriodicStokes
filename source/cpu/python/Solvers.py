import numpy as np
from Chebyshev import * 
from config import num_threads
import ctypes
import scipy.io
from scipy.interpolate import interp1d, PchipInterpolator

# get precompute func
libLinSolve = ctypes.CDLL('liblinSolve.so')
# get doubly periodic tools, including main solver struct
libDPTools = ctypes.CDLL('libdpTools.so')
# see end of file for lib function signatures

class Stokes(object):
  """
  Unifying TP and DP solver routines.
  
  See main routines below, and DPTools.cpp.
  
  Attributes:
    Kx, Ky, K - meshgrid of wavenumbers and K = sqrt(Kx^2 + Ky^2)
    Dx, Dy - Fourier derivative operator iKx (with 0 middle mode if N is even)
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    FIMat - first Chebyshev integral matrix (Nz x Nz + 2)
    pints - precomputed integrals of Chebyshev polynomials for pressure correction
          - this is only used in DP no_wall if the k0 flag is 1
    uvints - precomputed integrals of Chebyshev polynomials for velocity correction
           - this is only used in DP no_wall if the k0 flag is 1
    A_kpb, Ainv_B, C, PIV, Ginv - Components of matrices in BVP for each k != 0 
                       - See precomputeLinOps and its side effects for details
    C_k0, Ginv_k0 - analog of C, Ginv above for k = 0
    BCR1, BCL1, BCR2, BCL2 - See DoublyPeriodicNoWallBCs for details
                           - this is only used in the bottom_wall and 
                           - slit_channel routines
    Nx, Ny, Nz - number of points in x,y,z
    Lx, Ly - extent of x and y grids
    H - half extent of z grid (Lz / 2)
    solverType - 'TP' is triply periodic
               - 'DP' is doubly periodic no wall
               - 'DPBW' is doubly periodic no slip bottom wall
               - 'DPSC' is doubly periodic no slip slit channel
    solver - pointer to c++ dp solver struct
    P_hat_r/i - real and imaginary part coeffs of pressure
    U_hat_r/i - real and imaginary part coeffs of velocity
    zpts - cheb pts in z, only needed for DPBW or DPSC
         - set with SetZ
  """
  def __init__(self, _Nx, _Ny, _Nz, _Lx, _Ly, _Lz, _dof, _eta, _k0, _solverType):
    """ 
    The constructor for the Solver class.
    
    Parameters:
      Nx, Ny, Nz - number of points in x,y,z
      Lx, Ly - extent of x and y grids
      H - half extent of z grid (Lz / 2)
      dof - deg of freedom (should always be 3)
      eta - viscosity
      k0  - if using a DP solver, choose whether to 0 the k=0 mode of the RHS for the no wall routine
            k0 = 0 - the k=0 mode of the RHS for pressure and velocity will be 0
                   - the k=0 mode of pressure and velocity will remain 0
            k0 = 1 - the k=0 mode of the RHS for pressure and velocity will not be 0
                   - there will be a bvp solve for the k=0 mode of pressure and velocity,
                     and there will be a correction after each solve
            For TP solverType, k0 is always 0.
    Side Effects: Populates precomputable attributes and initializes solver
                  Sets Cf, U/P_hat_r/i to None
    """
    self.solverType = _solverType; self.k0 = _k0; self.eta = _eta
    self.Nx = _Nx; self.Ny = _Ny; self.Nz = _Nz; 
    self.Lx = _Lx; self.Ly = _Ly; self.Lz = _Lz; self.H = _Lz / 2; self.dof = _dof
    self.Ntot = self.Nx * self.Ny * self.Nz * self.dof
    self.U_hat_i = None; self.U_hat_r = None; self.P_hat_i = None; self.P_hat_r = None
    self.Ut_hat_i = None; self.Ut_hat_r = None; self.solver = None 
    if _solverType == 'TP':
      self.k0 = 0
      self.Ksq, self.Dx, self.Dy, self.Dz\
       = TriplyPeriodicStokes_init(self.Nx, self.Ny, self.Nz, self.Lx, self.Ly, self.Lz)
    elif _solverType == 'DP' or _solverType == 'DPBW' or _solverType == 'DPSC': 
      # precompute what we can
      self.Kx, self.Ky, self.K, self.Dx, self.Dy,\
       self.FIMat, self.SIMat, self.pints, self.uvints, \
        self.A_kpb, self.Ainv_B, self.C, self.C_k0, self.Ginv,\
         self.Ginv_k0, self.BCR1, self.BCL1, self.BCR2, self.BCL2, \
          = DoublyPeriodicStokes_init(self.Nx, self.Ny, self.Nz, self.Lx, self.Ly, self.H)
      # initialize the dp no wall solver
      self.solver = libDPTools.Solver(self.Dx, self.Dy, self.A_kpb.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                      self.C, self.Ginv, self.Ainv_B, self.FIMat, self.SIMat,\
                                      self.K.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                      self.H, self.eta, self.Ny * self.Nx, self.Nz, self.dof,\
                                      self.k0, num_threads)
    else:
      raise ValueError("Solver type not recognized")
  
  def SetZ(self, _zpts):
    self.zpts = _zpts  

  def Solve(self, fG_hat_r, fG_hat_i, compute_curl=False):
    """
    Solve Stokes eq. in Fourier domain for triply periodic
    or Fourier-Chebyshev domain for doubly periodic

    If compute_curl=True, the fluid vorticity 
    is evaluated in Fourier space
    by taking the curl of the linear velocity
    """
    if self.solverType == 'TP':
      self.U_hat_r, self.U_hat_i\
        = TriplyPeriodicStokes(fG_hat_r, fG_hat_i, self.eta, self.Ksq,\
                               self.Dx, self.Dy, self.Dz, self.Ntot)
    elif self.solverType == 'DP':
      self.U_hat_r, self.U_hat_i, self.P_hat_r, self.P_hat_i\
        = DoublyPeriodicStokes_no_wall(self.solver, fG_hat_r, fG_hat_i,\
                                       self.FIMat, self.SIMat, self.pints, \
                                       self.uvints, self.C_k0, self.Ginv_k0,\
                                       self.k0, self.eta, self.Nx, self.Ny, self.Nz)
    elif self.solverType == 'DPBW':
      self.U_hat_r, self.U_hat_i, self.P_hat_r, self.P_hat_i\
        = DoublyPeriodicStokes_bottom_wall(self.solver, fG_hat_r, fG_hat_i,\
                                           self.zpts, self.Kx, self.Ky, self.FIMat,\
                                           self.SIMat, self.pints, self.uvints, self.C_k0,\
                                           self.Ginv_k0, self.BCR1, self.BCL2, self.k0, self.eta,\
                                           self.Nx, self.Ny, self.Nz, self.H)
    elif self.solverType == 'DPSC':
      self.U_hat_r, self.U_hat_i, self.P_hat_r, self.P_hat_i\
        = DoublyPeriodicStokes_slit_channel(self.solver, fG_hat_r, fG_hat_i,\
                                            self.zpts, self.Kx, self.Ky, self.FIMat,\
                                            self.SIMat, self.pints, self.uvints, self.C_k0,\
                                            self.Ginv_k0, self.BCR2, self.BCL2, self.k0, self.eta,\
                                            self.Nx, self.Ny, self.Nz, self.H)
    else:
      raise ValueError("Solver type not recognized")
    # we're done if no torque
    if not compute_curl:
      return self.U_hat_r, self.U_hat_i, self.P_hat_r, self.P_hat_i
    # evaluate angular velocity otherwise
    else:

      # use Fourier derivative operator in z for TP, and Chebyshev for DP
      if self.solverType == 'TP':
        Cu = self.U_hat_r[0::3] + 1j * self.U_hat_i[0::3]
        Cv = self.U_hat_r[1::3] + 1j * self.U_hat_i[1::3]
        Cw = self.U_hat_r[2::3] + 1j * self.U_hat_i[2::3]
        Cchi = self.Dz * Cu
        Cpsi = self.Dz * Cv 
        # vorticity 
        Cut = 0.5 * (self.Dy * Cw - Cpsi) 
        Cvt = 0.5 * (-self.Dx * Cw + Cchi) 
        Cwt = 0.5 * (self.Dx * Cv - self.Dy * Cu)
        self.Ut_hat_r = np.zeros((self.Ny * self.Nx * self.Nz * self.dof,))
        self.Ut_hat_i = np.zeros((self.Ny * self.Nx * self.Nz * self.dof,))
        self.Ut_hat_r[0::3] = np.real(Cut).copy()
        self.Ut_hat_r[1::3] = np.real(Cvt).copy()
        self.Ut_hat_r[2::3] = np.real(Cwt).copy()
        self.Ut_hat_i[0::3] = np.imag(Cut).copy()
        self.Ut_hat_i[1::3] = np.imag(Cvt).copy()
        self.Ut_hat_i[2::3] = np.imag(Cwt).copy()
      else:
        if self.solverType == 'DP':
          self.U_hat_r.shape = (self.Ny * self.Nx, self.Nz, self.dof)
          self.U_hat_i.shape = (self.Ny * self.Nx, self.Nz, self.dof)
        Cu = self.U_hat_r[:,:,0] + 1j * self.U_hat_i[:,:,0]
        Cv = self.U_hat_r[:,:,1] + 1j * self.U_hat_i[:,:,1]
        Cw = self.U_hat_r[:,:,2] + 1j * self.U_hat_i[:,:,2]
        Cchi = chebCoeffDiff_perm(Cu, self.Nx, self.Ny, self.Nz, 1, self.H).reshape((self.Ny * self.Nx, self.Nz))
        Cpsi = chebCoeffDiff_perm(Cv, self.Nx, self.Ny, self.Nz, 1, self.H).reshape((self.Ny * self.Nx, self.Nz))
        self.Dy.shape = self.Dx.shape = (self.Ny * self.Nx, 1) 
        # vorticity 
        Cut = 0.5 * (self.Dy * Cw - Cpsi) 
        Cvt = 0.5 * (-self.Dx * Cw + Cchi) 
        Cwt = 0.5 * (self.Dx * Cv - self.Dy * Cu)
        self.Ut_hat_r = np.zeros((self.Ny * self.Nx, self.Nz, self.dof,))
        self.Ut_hat_i = np.zeros((self.Ny * self.Nx, self.Nz, self.dof,))
        self.Ut_hat_r[:,:,0] = np.real(Cut).copy()
        self.Ut_hat_r[:,:,1] = np.real(Cvt).copy()
        self.Ut_hat_r[:,:,2] = np.real(Cwt).copy()
        self.Ut_hat_i[:,:,0] = np.imag(Cut).copy()
        self.Ut_hat_i[:,:,1] = np.imag(Cvt).copy()
        self.Ut_hat_i[:,:,2] = np.imag(Cwt).copy()
      return self.U_hat_r, self.U_hat_i, self.Ut_hat_r, self.Ut_hat_i, self.P_hat_r, self.P_hat_i
       

  def Clean(self):
    if self.solverType == 'DP' or self.solverType == 'DPBW' or self.solverType == 'DPSC': 
      if self.solver is not None:
        libDPTools.Clean(self.solver) 
    elif self.solverType == 'TP':
      return
    else:
      raise ValueError("Solver type not recognized")


###################################################################################
########################## Main solver routines ###################################
###################################################################################

def TriplyPeriodicStokes_init(Nx, Ny, Nz, Lx, Ly, Lz):
  """
  precompute wave number grids for triply periodic stokes solver
  """  
  # wave numbers
  kvec_x = 2*np.pi*np.concatenate((np.arange(0,np.floor(Nx/2)+1),\
                                   np.arange(-1*np.ceil(Nx/2)+1,0)), axis=None) / Lx
  
  kvec_y = 2*np.pi*np.concatenate((np.arange(0,np.floor(Ny/2)+1),\
                                   np.arange(-1*np.ceil(Ny/2)+1,0)), axis=None) / Ly
  
  kvec_z = 2*np.pi*np.concatenate((np.arange(0,np.floor(Nz/2)+1),\
                                   np.arange(-1*np.ceil(Nz/2)+1,0)), axis=None) / Lz

  #Kz, Ky, Kx = [a.flatten() for a in np.meshgrid(kvec_z, kvec_y, kvec_x, indexing='ij')]
  Ky, Kx, Kz = [a.flatten() for a in np.meshgrid(kvec_y, kvec_x, kvec_z, indexing='ij')]
  Ksq = Kx**2 + Ky**2 + Kz**2
  Dx = 1j * Kx; Dy = 1j * Ky; Dz = 1j * Kz;
  # zero unpaired mode
  if Nx % 2 == 0: 
    Dx.reshape((Ny,Nx,Nz))[:,int(Nx/2),:] = 0
  if Ny % 2 == 0:
    Dy.reshape((Ny,Nx,Nz))[int(Ny/2),:,:] = 0
  if Nz % 2 == 0:
    Dz.reshape((Ny,Nx,Nz))[:,:,int(Nz/2)] = 0
  return Ksq, Dx, Dy, Dz

def TriplyPeriodicStokes(fG_hat_r, fG_hat_i, eta, Ksq, Dx, Dy, Dz, Ntotal):
  """
  Solve triply periodic Stokes eq in Fourier domain given the Fourier
  coefficients of the forcing.
  
  Parameters:
    fG_hat_r, fG_hat_i - real and complex part of Fourier coefficients
                         of spread forces on the grid. These are both
                         arrays of doubles (not complex).
  
    eta - viscocity
    Ntotal - total number of points in x,y,z
    Kx, Ky, Kz (Dx..) - wave numbers and Fourier derivative operators 
  Returns:
    U_hat_r, U_hat_i - real and complex part of Fourier coefficients of
                       fluid velocity on the grid. 
  
  Note: We assume the net force on the unit cell is 0 by *ignoring* 
        the k = 0 mode. That is, the k=0 mode of the output solution
        will be 0.
  """
  # separate x,y,z components
  f_hat = (fG_hat_r[0::3] + 1j * fG_hat_i[0::3])
  g_hat = (fG_hat_r[1::3] + 1j * fG_hat_i[1::3])
  h_hat = (fG_hat_r[2::3] + 1j * fG_hat_i[2::3])

  # precompute parts of RHS
  p_hat = np.divide(-1.0 * (Dx * f_hat + Dy * g_hat + Dz * h_hat), Ksq, \
                  out = np.zeros_like(f_hat), where = Ksq != 0, dtype = np.complex) 
  I2 = np.divide(1, eta * Ksq, out = np.zeros_like(Ksq), where = Ksq != 0, dtype = np.double)
  # solve for Fourier coeffs of velocity
  u_hat = I2 * (f_hat - (Dx * p_hat))
  v_hat = I2 * (g_hat - (Dy * p_hat))
  w_hat = I2 * (h_hat - (Dz * p_hat))
  # ignore k = 0
  u_hat[0] = 0
  v_hat[0] = 0
  w_hat[0] = 0
  # interleave solution components and split
  # real/imaginary parts for passing back to c
  U_hat_r = np.zeros((Ntotal,), dtype = np.double)
  U_hat_r[0::3] = np.real(u_hat)
  U_hat_r[1::3] = np.real(v_hat)
  U_hat_r[2::3] = np.real(w_hat)
  U_hat_i = np.zeros((Ntotal,), dtype = np.double)
  U_hat_i[0::3] = np.imag(u_hat)
  U_hat_i[1::3] = np.imag(v_hat)
  U_hat_i[2::3] = np.imag(w_hat)
  return U_hat_r, U_hat_i

# Precomputations for all DP solvers
def DoublyPeriodicStokes_init(Nx, Ny, Nz, Lx, Ly, H):
  """
  Precompute the linear operators and boundary conditions
  for the doubly periodic no wall problem. The return
  values of this function are used in the no_wall,
  bottom_wall and slit_channel solvers

  This function is useful when several solves are done
  on the same grid, as it only needs to be called
  before the first solve eg) at t=0 in a time stepping context.

  Parameters:
    Nx, Ny, Nz - number of points in x,y,z
    Lx, Ly - extent of x and y grids
    H - half extent of z grid (Lz / 2)
  
  Side Effects: None
  Returns:
    Kx, Ky, K - meshgrid of wavenumbers and K = sqrt(Kx^2 + Ky^2)
    Dx, Dy - Fourier derivative operator iKx (with 0 middle mode if N is even)
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    FIMat - first Chebyshev integral matrix (Nz x Nz + 2)
    pints - precomputed integrals of Chebyshev polynomials for pressure correction
          - this is only used in DP no_wall if the k0 flag is 1
    uvints - precomputed integrals of Chebyshev polynomials for velocity correction
           - this is only used in DP no_wall if the k0 flag is 1
    LU, Ainv_B, C, PIV, Ginv - Components of matrices in BVP for each k != 0 
                       - See precomputeBandedLinOps and its side effects for details
    C_k0, Ginv_k0 - analog of C, Ginv above for k = 0
    BCR1, BCL1, BCR2, BCL2 - See DoublyPeriodicNoWallBCs for details
                           - this is only used in the bottom_wall and 
                           - slit_channel routines
  """
  # wave numbers
  kvec_x = 2*np.pi*np.concatenate((np.arange(0,np.floor(Nx/2)+1),\
                                   np.arange(-1*np.ceil(Nx/2)+1,0)), axis=None) / Lx
  
  kvec_y = 2*np.pi*np.concatenate((np.arange(0,np.floor(Ny/2)+1),\
                                   np.arange(-1*np.ceil(Ny/2)+1,0)), axis=None) / Ly
  Kx, Ky = [K.flatten() for K in np.meshgrid(kvec_x, kvec_y, indexing='xy')]
  Ksq = Kx**2 + Ky**2; K = np.sqrt(Ksq)
  Dx = 1j * Kx; Dy = 1j * Ky;
  # zero unpaired mode
  if Nx % 2 == 0: 
    Dx.reshape((Ny,Nx))[:,int(Nx/2)] = 0
  if Ny % 2 == 0:
    Dy.reshape((Ny,Nx))[int(Ny/2),:] = 0
  # get Chebyshev integration matrices
  FIMat = firstIntegralMatrix(Nz, H)
  SIMat = secondIntegralMatrix(Nz, H)
  diagsEye = getDiagsP2M2(np.eye(Nz))
  diagsSIMat = getDiagsP2M2(SIMat[:,:Nz])
  # precompute integrals of cheb polys for pressure/vel correction
  pints, uvints = precomputeInts(Nz, H)
  # get boundary conditions rows
  BCR1, BCR2, BCL1, BCL2 = DoublyPeriodic_no_wall_BCs(Nz)
  # assemble appropriate BCs for k = 0 and k != 0
  BCs_k0 = np.stack((BCR2, -BCL2), axis = 0)
  # (Nyx, 2, Nz + 2)
  # (Nyx, 2, Nz + 2)
  BCs_k = np.einsum('ij, k->kij', np.stack((BCR2, BCL2), axis = 0), H**2 * K) 
  BCs_k +=  BCs_k0 + np.stack((H * BCR1 - BCR2, H * BCL1 + BCL2), axis = 0) 
  # precompute system matrices for each k and permute to fortran layout
  A = diagsEye - np.einsum('i, k->ki', diagsSIMat, Ksq)
  A_kbp = np.zeros((Ny * Nx * (3 * Nz - 3),), dtype = np.float64) 
  # (Nyx, 2, Nz)
  B = np.einsum('ij, k->kij', SIMat[:,Nz::], -Ksq)
  B = np.transpose(B, (0,2,1)).copy().astype(np.float64)
  # (Nyx, Nz, 2)
  C = BCs_k[:,:,0:Nz]
  C = np.transpose(C, (0,2,1)).copy().astype(np.float64); C_k0 = BCs_k0[:,0:Nz]; 
  # (Nyx, 2, 2)
  D = BCs_k[:,:,Nz::]
  D = np.transpose(D, (0,2,1)).copy().astype(np.float64); D_k0 = -BCs_k0[:,Nz::];
  # precompute LU decompositions and precomputable linear solves for each k
  Ginv = np.zeros((4 * Ny * Nx,), dtype = np.float64)
  G = np.zeros((4 * Ny * Nx,), dtype = np.float64)
  # precomputations for inverting A at each k and computing Ainv * B
  # A_kbp is overwritten with pentadiagonal solver info for A,
  # B is overwritten with Ainv * B
  # Ginv is the inverse of the 2x2 schur complement
  precomputeLinOps(A, A_kbp, B, C, D, G, Ginv, Ny * Nx, Nz)
  # rename overwritten vars
  Ainv_B = B
  # get Ginv for k=0
  Ginv_k0 = 1.0 / (D_k0[0,0] * D_k0[1,1] - D_k0[0,1] * D_k0[1,0]) * \
                   np.array([[D_k0[1,1], -D_k0[0,1]], [-D_k0[1,0],D_k0[0,0]]])
  return Kx, Ky, K, Dx, Dy, FIMat.astype(np.complex128), SIMat.astype(np.complex128), pints, uvints,\
         A_kbp, Ainv_B.astype(np.complex128), C.astype(np.complex128), C_k0, Ginv.astype(np.complex128),\
         Ginv_k0, BCR1, BCL1, BCR2, BCL2

# DP no wall solver
def DoublyPeriodicStokes_no_wall(solver, fG_hat_r, fG_hat_i, FIMat, SIMat, pints, \
                                 uvints, C_k0, Ginv_k0, k0, eta, Nx, Ny, Nz):
  """
  Solve doubly periodic Stokes eq in Fourier-Chebyshev domain given the Fourier-Chebyshev
  coefficients of the forcing.
  
  Parameters:
    solver - Stokes solver object
    fG_hat_r, fG_hat_i - real and complex part of Fourier-Chebyshev coefficients
                         of spread forces on the grid. These are both
                         arrays of doubles (not complex).
  
    eta - viscocity
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    FIMat - first Chebyshev integral matrix (Nz x Nz + 2)
    pints - precomputed integrals of Chebyshev polynomials for pressure correction
          - this is only used in DP no_wall if the k0 flag is 1
    uvints - precomputed integrals of Chebyshev polynomials for velocity correction
           - this is only used in DP no_wall if the k0 flag is 1
    C_k0, Ginv_k0 - lin ops for k = 0
    k0 - switch determining how to handle k = 0 mode of pressure and velocity
       - if k0 = 0, k = 0 mode of RHS is assumed to be 0, and sol will be 0
       - if k0 = 1, the k = 0 mode of RHS is assumed to be non-zero. There
         will be a correction to the non-zero solution.
  
  Returns:
    U_hat_r, U_hat_i, P_hat_r, P_hat_i - real and complex part of 
                                       - Fourier-Chebyshev coefficients of
                                         fluid velocity and pressure on the grid. 
  """
  dof = 3; Nyx = Ny * Nx
  # separate x,y,z components
  # set rhs for solver
  libDPTools.SetRHS_split(solver, fG_hat_r, fG_hat_i)
  # solve for k > 0
  libDPTools.DoublyPeriodicSolve_no_wall(solver)
  # get split sol arrays
  U_hat_r = np.ctypeslib.as_array(libDPTools.GetU_real(solver), shape = (Nyx * Nz * dof,))
  U_hat_i = np.ctypeslib.as_array(libDPTools.GetU_imag(solver), shape = (Nyx * Nz * dof,))
  P_hat_r = np.ctypeslib.as_array(libDPTools.GetP_real(solver), shape = (Nyx * Nz,))
  P_hat_i = np.ctypeslib.as_array(libDPTools.GetP_imag(solver), shape = (Nyx * Nz,))

  # handle k = 0 for pressure and velocity
  if k0 == 1:
    # get real/imag part of p_rhs at k=0
    p_rhs_k0_r = np.ctypeslib.as_array(libDPTools.GetP_RHS_real(solver, 0), shape = (Nz,))
    p_rhs_k0_i = np.ctypeslib.as_array(libDPTools.GetP_RHS_imag(solver, 0), shape = (Nz,))
    p_rhs_k0 = p_rhs_k0_r + 1j * p_rhs_k0_i
    # get k=0 of forces
    ch_k0 = fG_hat_r[2:dof * Nz:dof] + 1j * fG_hat_i[2:dof * Nz:dof]
    # solve k0=1 type corrected bvp for pressure
    cp_k0, dp_k0 = DoublyPeriodicStokes_no_wall_solvePressureBVP_k0(\
                      p_rhs_k0, C_k0, Ginv_k0, SIMat, FIMat, ch_k0, pints, k0)
    # save to main split arrays
    P_hat_r[0:Nz] = np.real(cp_k0).copy()
    P_hat_i[0:Nz] = np.imag(cp_k0).copy()
    # save k=0 of pressure and deriv into solver (needed for urhs)
    libDPTools.SetP(solver, cp_k0, 0)
    libDPTools.SetdP(solver, dp_k0, 0)

    # compute u_rhs for k=0
    U_rhs_k0_r = np.ctypeslib.as_array(libDPTools.GetU_RHS_real(solver, 0, 1), shape = (Nz * dof,))
    U_rhs_k0_i = np.ctypeslib.as_array(libDPTools.GetU_RHS_imag(solver, 0, 0), shape = (Nz * dof,))
    U_rhs_k0 = U_rhs_k0_r + 1j * U_rhs_k0_i
    u_rhs_k0 = U_rhs_k0[0:dof * Nz:dof]
    v_rhs_k0 = U_rhs_k0[1:dof * Nz:dof]
    w_rhs_k0 = U_rhs_k0[2:dof * Nz:dof]
    # get k=0 of forces
    cf_k0 = fG_hat_r[0:dof * Nz:dof] + 1j * fG_hat_i[0:dof * Nz:dof]
    cg_k0 = fG_hat_r[1:dof * Nz:dof] + 1j * fG_hat_i[1:dof * Nz:dof]
    # solve k0=1 type corrected bvp for velocity
    cu_k0, cv_k0, cw_k0 = DoublyPeriodic_no_wall_solveVelocityBVP_k0(\
                                  u_rhs_k0, v_rhs_k0, w_rhs_k0, C_k0,\
                                  Ginv_k0, SIMat, cf_k0, cg_k0, uvints, eta, k0)
    cU_k0 = np.zeros((3 * Nz,), dtype = np.complex128)
    cU_k0[0::3] = cu_k0; cU_k0[1::3] = cv_k0; cU_k0[2::3] = cw_k0;
    # save to main split arrays
    U_hat_r[0:dof * Nz] = np.real(cU_k0).copy()
    U_hat_i[0:dof * Nz] = np.imag(cU_k0).copy()
    # save k=0 of vel into solver
    libDPTools.SetU(solver, cU_k0, 0)

  return U_hat_r, U_hat_i, P_hat_r, P_hat_i


# DP bottom wall solver
def DoublyPeriodicStokes_bottom_wall(solver, fG_hat_r, fG_hat_i, zpts, Kx, Ky, FIMat, SIMat, pints, 
                                     uvints, C_k0, Ginv_k0, BCR1, BCL2, k0, eta, Nx, Ny, Nz, H):
  """
  Solve Stokes eq in doubly periodic bottom wall in the Fourier-Chebyshev domain 
  given the Fourier-Chebyshev coefficients of the forcing. We first solve a DP
  subproblem (no walls), and then compute a correction to the DP solution to enforce
  the no-slip BCs.
  
  Parameters:
    solver - Stokes solver object
    fG_hat_r, fG_hat_i - real and complex part of Fourier-Chebyshev coefficients
                         of spread forces on the grid. These are both
                         arrays of doubles (not complex).
    zpts - Chebyshev grid pts 
    eta - viscocity
    Nx, Ny, Nz - number of points in x, y and z
    H - Lz / 2
    Kx, Ky -  meshgrid of wavenumbers and K = sqrt(Kx^2 + Ky^2)
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    FIMat - first Chebyshev integral matrix (Nz x Nz + 2)
    pints - precomputed integrals of Chebyshev polynomials for pressure correction
          - this is only used in DP no_wall if the k0 flag is 1
    uvints - precomputed integrals of Chebyshev polynomials for velocity correction
           - this is only used in DP no_wall if the k0 flag is 1
    C_k0, Ginv_k0 - lin ops for k = 0
    BCR1, BCL2 - See DoublyPeriodicNoWallBCs for details
  
  Returns:
    U_hat_r, U_hat_i, P_hat_r, P_hat_i - real and complex part of 
                                       - Fourier-Chebyshev coefficients of
                                         fluid velocity and pressure on the grid. 
  """
  dof = 3; Nyx = Ny * Nx;
  # set rhs for solver
  libDPTools.SetRHS_split(solver, fG_hat_r, fG_hat_i)
  # solve for k > 0
  libDPTools.DoublyPeriodicSolve_bottom_wall(solver, Kx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                          Ky.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                          zpts.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                          eta, Nyx, Nz, dof, num_threads)
  # get split sol arrays
  U_hat_r = np.ctypeslib.as_array(libDPTools.GetU_real(solver), shape = (Nyx, Nz, dof))
  U_hat_i = np.ctypeslib.as_array(libDPTools.GetU_imag(solver), shape = (Nyx,  Nz, dof))
  P_hat_r = np.ctypeslib.as_array(libDPTools.GetP_real(solver), shape = (Nyx,  Nz))
  P_hat_i = np.ctypeslib.as_array(libDPTools.GetP_imag(solver), shape = (Nyx,  Nz))

  # Cheb coeffs of forcing for k = 0
  Cf_k0 = fG_hat_r[0:dof * Nz:dof] + 1j * fG_hat_i[0:dof * Nz:dof]
  Cg_k0 = fG_hat_r[1:dof * Nz:dof] + 1j * fG_hat_i[1:dof * Nz:dof]
  Ch_k0 = fG_hat_r[2:dof * Nz:dof] + 1j * fG_hat_i[2:dof * Nz:dof]
  Dh_k0 = chebCoeffDiff(Ch_k0, 1, 1, Nz, 1, H).reshape((Nz,))
  # First cheb coeff of pressure for k = 0
  Cp_k0 = P_hat_r[0,0] + 1j * P_hat_i[0,0]
  Dp_k0 = chebCoeffDiff(P_hat_r[0,:] + 1j * P_hat_i[0,:], 1, 1, Nz, 1, H).reshape((Nz,)) 
  # RHS for k = 0 correction solve for pressure and vel
  p_RHS_k0 = Dh_k0
  u_RHS_k0 = -Cf_k0 / eta
  v_RHS_k0 = -Cg_k0 / eta
  w_RHS_k0 = (Dp_k0 - Ch_k0) / eta
  
  # correct pressure for k = 0
  Cpcorr_k0 = DoublyPeriodicStokes_wall_solvePressureBVP_k0(p_RHS_k0, C_k0, Ginv_k0, SIMat, Ch_k0, pints, Cp_k0)
  # correct velocities for k = 0
  BCs_k0_bw = np.asfortranarray(np.stack((BCR1, BCL2), axis = 0))
  C_k0_bw = BCs_k0_bw[:,0:Nz]
  D_k0_bw = -BCs_k0_bw[:,Nz::]
 
  Ginv_k0_bw = 1.0 / (D_k0_bw[0,0] * D_k0_bw[1,1] - D_k0_bw[0,1] * D_k0_bw[1,0]) * \
                   np.array([[D_k0_bw[1,1], -D_k0_bw[0,1]], [-D_k0_bw[1,0],D_k0_bw[0,0]]])
 
  Cucorr_k0, Cvcorr_k0 \
    = DoublyPeriodic_wall_solveVelocityBVP_k0(u_RHS_k0, v_RHS_k0, w_RHS_k0, C_k0_bw,\
                                            Ginv_k0_bw, SIMat, Cf_k0, Cg_k0)

  # add k=0 correction to sol
  P_hat_r[0,:] += np.real(Cpcorr_k0)
  P_hat_i[0,:] += np.imag(Cpcorr_k0)
  U_hat_r[0,:,0] += np.real(Cucorr_k0)
  U_hat_i[0,:,0] += np.imag(Cucorr_k0)
  U_hat_r[0,:,1] += np.real(Cvcorr_k0)
  U_hat_i[0,:,1] += np.imag(Cvcorr_k0)
  return U_hat_r, U_hat_i, P_hat_r, P_hat_i  

# DP slit channel solver
def DoublyPeriodicStokes_slit_channel(solver, fG_hat_r, fG_hat_i, zpts, Kx, Ky, FIMat, SIMat, pints, 
                                      uvints, C_k0, Ginv_k0, BCR2, BCL2, k0, eta, Nx, Ny, Nz, H):
  """
  Solve Stokes eq in doubly periodic slit channel in the Fourier-Chebyshev domain 
  given the Fourier-Chebyshev coefficients of the forcing. We first solve a DP
  subproblem (no walls), and then compute a correction to the DP solution to enforce
  the no-slip BCs.
  
  Parameters:
    solver - Stokes solver object
    fG_hat_r, fG_hat_i - real and complex part of Fourier-Chebyshev coefficients
                         of spread forces on the grid. These are both
                         arrays of doubles (not complex).
    zpts - Chebyshev grid pts 
    eta - viscocity
    Kx, Ky - meshgrid of wavenumbers 
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    FIMat - first Chebyshev integral matrix (Nz x Nz + 2)
    pints - precomputed integrals of Chebyshev polynomials for pressure correction
          - this is only used in DP no_wall if the k0 flag is 1
    uvints - precomputed integrals of Chebyshev polynomials for velocity correction
           - this is only used in DP no_wall if the k0 flag is 1
    C_k0, Ginv_k0 - lin ops for k = 0
    BCR2, BCL2 - See DoublyPeriodicNoWallBCs for details
  
  Returns:
    U_hat_r, U_hat_i, P_hat_r, P_hat_i - real and complex part of 
                                       - Fourier-Chebyshev coefficients of
                                         fluid velocity and pressure on the grid. 
  """
  dof = 3; Nyx = Ny * Nx; Lz = 2 * H
  # set rhs for solver
  libDPTools.SetRHS_split(solver, fG_hat_r, fG_hat_i)
  # solve for k > 0
  libDPTools.DoublyPeriodicSolve_slit_channel(solver, Kx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                           Ky.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                           zpts.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                           Lz, eta, Nyx, Nz, dof, num_threads)
  # get split sol arrays
  U_hat_r = np.ctypeslib.as_array(libDPTools.GetU_real(solver), shape = (Nyx, Nz, dof))
  U_hat_i = np.ctypeslib.as_array(libDPTools.GetU_imag(solver), shape = (Nyx,  Nz, dof))
  P_hat_r = np.ctypeslib.as_array(libDPTools.GetP_real(solver), shape = (Nyx,  Nz))
  P_hat_i = np.ctypeslib.as_array(libDPTools.GetP_imag(solver), shape = (Nyx,  Nz))
  # Cheb coeffs of forcing for k = 0
  Cf_k0 = fG_hat_r[0:dof * Nz:dof] + 1j * fG_hat_r[0:dof * Nz:dof]
  Cg_k0 = fG_hat_r[1:dof * Nz:dof] + 1j * fG_hat_r[1:dof * Nz:dof]
  Ch_k0 = fG_hat_r[2:dof * Nz:dof] + 1j * fG_hat_r[2:dof * Nz:dof]
  Dh_k0 = chebCoeffDiff(Ch_k0, 1, 1, Nz, 1, H).reshape((Nz,))
  # First cheb coeff of pressure for k = 0
  Cp_k0 = P_hat_r[0,0] + 1j * P_hat_i[0,0]
  Dp_k0 = chebCoeffDiff(P_hat_r[0,:] + 1j * P_hat_i[0,:], 1, 1, Nz, 1, H).reshape((Nz,)) 

  # RHS for k = 0 correction solve for pressure and vel
  p_RHS_k0 = Dh_k0
  u_RHS_k0 = -Cf_k0 / eta
  v_RHS_k0 = -Cg_k0 / eta
  w_RHS_k0 = (Dp_k0 - Ch_k0) / eta

  # correct pressure for k = 0
  Cpcorr_k0 = DoublyPeriodicStokes_wall_solvePressureBVP_k0(p_RHS_k0, C_k0, Ginv_k0, SIMat, Ch_k0, pints, Cp_k0)

  # correct velocities for k = 0
  BCs_k0_tw = np.asfortranarray(np.stack((BCR2, BCL2), axis = 0))
  C_k0_tw = BCs_k0_tw[:,0:Nz]
  D_k0_tw = -BCs_k0_tw[:,Nz::]
  
  Ginv_k0_tw = 1.0 / (D_k0_tw[0,0] * D_k0_tw[1,1] - D_k0_tw[0,1] * D_k0_tw[1,0]) * \
                   np.array([[D_k0_tw[1,1], -D_k0_tw[0,1]], [-D_k0_tw[1,0],D_k0_tw[0,0]]])
  
  Cucorr_k0, Cvcorr_k0 = DoublyPeriodic_wall_solveVelocityBVP_k0(\
                               u_RHS_k0, v_RHS_k0, w_RHS_k0, C_k0_tw,\
                               Ginv_k0_tw, SIMat, Cf_k0, Cg_k0)

  # add the correction solution for k=0 
  P_hat_r[0,:] += np.real(Cpcorr_k0)
  P_hat_i[0,:] += np.imag(Cpcorr_k0)
  U_hat_r[0,:,0] += np.real(Cucorr_k0)
  U_hat_i[0,:,0] += np.imag(Cucorr_k0)
  U_hat_r[0,:,1] += np.real(Cvcorr_k0)
  U_hat_i[0,:,1] += np.imag(Cvcorr_k0)
  return U_hat_r, U_hat_i, P_hat_r, P_hat_i  

# BCs for DP problem (not usually called externally)
def DoublyPeriodic_no_wall_BCs(N):
  """
  This function gives the following BCs for the 
  no-wall Doubly periodic Stokes BVPs in z at each wave number:
  BCR1 = first derivative (first integral) evaluated at x=1
  BCR2 = function (second integral) evaluated at x=1
  BCL1 = first derivative (first integral) evaluated at x=-1
  BCL2 = NEGATIVE OF function (second integral) evaluated at x=-1

  Parameters:
    N - number of Chebyshev points
  Returns:
    BCR1, BCR2, BCL1, BCL2 - BCs for DP no-wall BVP solve in z (to be assembled with k) 
  """
  BCR1 = np.zeros((N+2,), dtype = np.double)
  BCR2 = np.zeros((N+2,), dtype = np.double)
  BCL1 = np.zeros((N+2,), dtype = np.double)
  BCL2 = np.zeros((N+2,), dtype = np.double)
  # Special cases - right
  BCR1[N+1] = 1; BCR2[N] = 1; BCR1[0] = 1; BCR1[2] = -1/2;
  BCR2[N+1] += 1; BCR2[1] = -1/8; BCR2[3] = 1/8;
  BCR1[1] += 1/4; BCR1[3] -= 1/4; BCR2[0] += 1/4; BCR2[2] -= (1/8 + 1/24); 
  BCR2[4] += 1/24;
  # Special cases -left 
  BCL1[N+1] = 1; BCL2[N] = -1; BCL1[0] = -1; BCL1[2] = 1/2;  
  BCL2[N+1] += 1; BCL2[1] = -1/8; BCL2[3] = 1/8;  
  BCL1[1] += 1/4; BCL1[3] -= 1/4; BCL2[0] -= 1/4; BCL2[2] += 1/8 + 1/24;  
  BCL2[4] -= 1/24;  
  # Easy cases    
  jj = np.arange(3,N)
  BCR1[2:N-1] += 1 / (2 * jj)  
  BCL1[2:N-1] += np.power(-1, jj) / (2 * jj)  
  BCR1[4:N+1] -= 1 / (2 * jj) * (jj < N - 1)
  BCL1[4:N+1] -= np.power(-1, jj) / (2 * jj) * (jj < N - 1)    
  BCR2[1:N-2] += 1 / (2 * jj) * 1 / (2 * jj - 2)
  BCL2[1:N-2] -= 1 / (2 * jj) * 1 / (2 * jj - 2) * np.power(-1, jj)  
  BCR2[5:N+2] += 1 / (2 * jj) * 1 / (2 * jj + 2) * (jj < N - 2)
  BCL2[5:N+2] -= 1 / (2 * jj) * 1 / (2 * jj + 2) * np.power(-1, jj) * (jj < N - 2) 
  BCR2[3:N] -= (1 / (2 * jj) * 1 / (2 * jj - 2) + 1 / (2 * jj) * 1 / (2 * jj + 2) * (jj < N - 1))
  BCL2[3:N] += (1 / (2 * jj) * 1 / (2 * jj - 2) + 1 / (2 * jj) * 1 / (2 * jj + 2) * (jj < N - 1)) * np.power(-1, jj)
  return BCR1, BCR2, BCL1, BCL2 

###################################################################################
####################### Pressure solve subroutines for k = 0 ######################
###################################################################################

def DoublyPeriodicStokes_no_wall_solvePressureBVP_k0(p_RHS, C, Ginv, SIMat, FIMat, Ch, pints, k0):
  """
  This function solves the pressure poisson equation for k = 0, based on the switch k0.

  Parameters:
    p_RHS - Nz x 1 Fourier-Chebyshev coeffs of div(f) at k = 0
    C - 2 x Nz third block of solve matrix for k = 0 (see DoublyPeriodicStokes_no_wall() for details)
    Ginv - 2 x 2 inverse of schur complement of first block of solve matrix at k = 0
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    FIMat - first Chebyshev integral matrix (Nz x Nz + 2)
    Ch - Fourier-Chebyshev coeffs of z-component of f at k = 0
    pints - precomputed integrals of Chebyshev polynomials for pressure correction
    k0 - switch to handle k = 0 mode 
       - if k0 = 0, k = 0 mode of RHS is assumed to be 0, and sol will be 0
       - if k0 = 1, the k = 0 mode of RHS is assumed to be non-zero. There
         will be a correction to the non-zero solution.

  Returns:
    Cp, Dp - Fourier-Chebyshev coeffs of pressure and its derivative at k = 0
  """
  Nz = p_RHS.shape[0]
  # sol will be 0 if k=0 of RHS is 0
  if k0 == 0:
    Cp = np.zeros((Nz,))
    Dp = np.zeros((Nz,))
  # if k=0 of RHS is not 0, we solve the pressure poblem with homogenous BCs
  # and must correct the 2nd Chebyshev coefficient of that mode (see eq. 18 in Ondrej's report)
  elif k0 == 1:
    secD = np.concatenate((p_RHS, Ginv @ (C @ p_RHS)))
    Cp = SIMat @ secD
    Dp = FIMat @ secD
    Cp[1] += 0.5 * pints @ Ch
  return Cp, Dp

def DoublyPeriodicStokes_wall_solvePressureBVP_k0(p_RHS, C, Ginv, SIMat, Ch, pints, cp0):
  """
  This function gives the correction for k=0 of pressure in the DP+correction method.
  It works for either bottom wall or slit channel.

  Parameters:
    p_RHS - Nz x 1 Fourier-Chebyshev coeffs of div(f) at k = 0
    C - 2 x Nz third block of solve matrix for k = 0 (see DoublyPeriodicStokes_no_wall() for details)
    Ginv - 2 x 2 inverse of schur complement of first block of solve matrix at k = 0
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    Ch - Fourier-Chebyshev coeffs of z-component of f at k = 0
    pints - precomputed integrals of Chebyshev polynomials for pressure correction
    cp0 - coeff of first cheb poly for k = 0 component of the DP subproblem
  Returns:
    Cp, Dp - Fourier-Chebyshev coeffs of pressure and its derivative at k = 0
  """
  secD = np.concatenate((p_RHS, Ginv @ (C @ p_RHS)))
  Cp = SIMat @ secD
  Cp[1] += 0.5 * pints @ Ch
  Cp[0] = cp0 + 0.5 * pints @ Ch
  return Cp

###################################################################################
####################### Velocity solve subroutines for k = 0 ######################
###################################################################################

def DoublyPeriodic_no_wall_solveVelocityBVP_k0(u_RHS, v_RHS, w_RHS, C, Ginv, SIMat,\
                                               Cf, Cg, uvints, eta, k0):
  """
  This function solves the velocity BVPs for k = 0, based on the switch k0.

  Parameters:
    u_RHS, v_RHS, w_RHS - Nz x 1 Fourier-Chebyshev coeffs of div(f) at k = 0
    C - 2 x Nz third block of solve matrix for k = 0 (see DoublyPeriodicStokes_no_wall() for details)
    Ginv - 2 x 2 inverse of schur complement of first block of solve matrix at k = 0
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    FIMat - first Chebyshev integral matrix (Nz x Nz + 2)
    Cf, Cg - Fourier-Chebyshev coeffs of x and y components of f at k = 0
    uvints - precomputed integrals of Chebyshev polynomials for velocity correction
    k0 - switch to handle k = 0 mode 
       - if k0 = 0, k = 0 mode of RHS is assumed to be 0, and sol will be 0
       - if k0 = 1, the k = 0 mode of RHS is assumed to be non-zero. There
         will be a correction to the non-zero solution.

  Returns:
    Cu, Cv, Cw - Fourier-Chebyshev coeffs of velocity components at k = 0
  """
  Nz = u_RHS.shape[0]
  # sol will be 0 if k=0 of RHS is 0
  if k0 == 0:
    Cu = np.zeros((Nz,))
    Cv = np.zeros((Nz,))
    Cw = np.zeros((Nz,))
  # if k=0 of RHS is not 0, we solve with homogenous BCs and
  # add the linear in z term for the null space   
  elif k0 == 1:
    au = np.concatenate((u_RHS, Ginv @ (C @ u_RHS)))
    av = np.concatenate((v_RHS, Ginv @ (C @ v_RHS)))
    aw = np.concatenate((w_RHS, Ginv @ (C @ w_RHS)))
    Cu = SIMat @ au
    Cv = SIMat @ av
    Cw = SIMat @ aw
    Cu[1] += 0.5 / eta * uvints @ Cf
    Cv[1] += 0.5 / eta * uvints @ Cg
  return Cu, Cv, Cw

def DoublyPeriodic_wall_solveVelocityBVP_k0(u_RHS, v_RHS, w_RHS, C, Ginv, SIMat,Cf, Cg):
  """
  This function gives the correction for k=0 of the velocity in the DP+correction method.
  It works for either bottom wall or slit channel.

  Parameters:
    u_RHS, v_RHS, w_RHS - Nz x 1 Fourier-Chebyshev coeffs of div(f) at k = 0
    C - 2 x Nz third block of solve matrix for k = 0 (see DoublyPeriodicStokes_no_wall() for details)
    Ginv - 2 x 2 inverse of schur complement of first block of solve matrix at k = 0
    SIMat - second Chebyshev integral matrix (Nz x Nz + 2)
    FIMat - first Chebyshev integral matrix (Nz x Nz + 2)
    Cf, Cg - Fourier-Chebyshev coeffs of x and y components of f at k = 0

  Returns:
    Cu, Cv - Fourier-Chebyshev coeffs of velocity components at k = 0
  """
  Nz = u_RHS.shape[0]
  au = np.concatenate((u_RHS, Ginv @ (C @ u_RHS)))
  av = np.concatenate((v_RHS, Ginv @ (C @ v_RHS)))
  Cu = SIMat @ au
  Cv = SIMat @ av
  return Cu, Cv


###################################################################################
####################### DP wall correction subroutines ############################
###################################################################################

# define wrappers for dptools
#def evalTheta(phi_in, theta, Nyx, Nz, dof):
def evalTheta(solver, theta, Nyx, dof):
  """
    This function is used to evaluate a Chebyshev series at a given
    value of theta (point on the cheb grid)
    
    Parameters: 
      in - the input array (size (Nz * Nyx * dof,1))
         - these are the Fourier-Chebyshev coeffs on the grid
      theta - determines the slice in z
            - eg) theta = pi is z = 0, theta = 0 is z - Lz 
      Nyx - total number of points in x,y
      Nz  - number of points in z
      dof - degrees of freedom
    
    Side Effects: None
    
    Returns : Phi_out - the output array (size (Nyx * dof, 1))
                      - these are the Fourier-Chebyshev coeffs on the x-y
                      - plane at a given z value 
  """
  libDPTools.evalTheta(solver, theta);
  phi_out_r = np.ctypeslib.as_array(libDPTools.GetPhi_out_real(solver), shape = (Nyx * dof,))
  phi_out_i = np.ctypeslib.as_array(libDPTools.GetPhi_out_imag(solver), shape = (Nyx * dof,))
  return phi_out_r + 1j * phi_out_i

# correction sol for bottom wall or slit channel 
def evalCorrectionSol(solver, zpts, Kx, Ky, eta, Nx, Ny, Nz, dof):
  """ 
  Evaluate the analytical correction to the DP solve to enforce no-slip 
  BCs in either bottom wall or slit channel (calls c lib)
 
  Parameters :
    solver - Stokes solver object 
    zpts - Chebyshev points in z
    Kx, Ky - meshgrid of wave numbers in x,y
    eta - viscosity
    Nx, Ny, Nz - num points in x,y,z
    dof - degrees of freedom   
  
  Side Effects : None
  Returns : U[p]_hat_corr_r[i] - split real/complex Fourier-Chebyshev 
                               - correction for pressure and velocity
                               - the k = 0 element is 0
  """
  Nyx = Ny * Nx;
  U_hat_r = np.ctypeslib.as_array(libDPTools.GetU_real(solver), shape = (Nyx * Nz * dof,)).copy()
  U_hat_i = np.ctypeslib.as_array(libDPTools.GetU_imag(solver), shape = (Nyx * Nz * dof,)).copy()
  P_hat_r = np.ctypeslib.as_array(libDPTools.GetP_real(solver), shape = (Nyx * Nz,)).copy()
  P_hat_i = np.ctypeslib.as_array(libDPTools.GetP_imag(solver), shape = (Nyx * Nz,)).copy()
  if solver.solverType == 'DPBW':
    libDPTools.evalCorrectionSol_bottomWall(solver, Kx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                            Ky.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                            zpts.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                            eta, Nyx, Nz, dof, num_threads)
  elif solver.solverType == 'DPSC':
    libDPTools.evalCorrectionSol_slitChannel(solver, Kx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                             Ky.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                             zpts.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                                             Lz, eta, Nyx, Nz, dof, num_threads)
  else:
    raise ValueError('Only DPBW and DPSC have a correction solution') 
  U_hat_corr_r =  np.ctypeslib.as_array(libDPTools.GetU_real(solver), shape = (Nyx * Nz * dof,)) - U_hat_r
  U_hat_corr_i =  np.ctypeslib.as_array(libDPTools.GetU_imag(solver), shape = (Nyx * Nz * dof,)) - U_hat_i
  P_hat_corr_r =  np.ctypeslib.as_array(libDPTools.GetP_real(solver), shape = (Nyx * Nz,)) - U_hat_r
  P_hat_corr_i =  np.ctypeslib.as_array(libDPTools.GetP_imag(solver), shape = (Nyx * Nz,)) - U_hat_i
  return U_hat_corr_r, U_hat_corr_i, P_hat_corr_r, P_hat_corr_i

###################################################################################
####################### Wrappers and helpers ######################################
###################################################################################

# define wrappers for linear solvers
def precomputeLinOps(A, A_kbp, B, C, D, G, Ginv, Nyx, Nz):
  """
  Given a block linear system of the form 
  
  |A B||a |   |f|
  |C D||c0| = |alpha| 
       |d0|   |beta|,  
              
  
  where A is Nz x Nz and banded, B is Nz x 2, C is 2 x Nz and D is 2x2 (all real), 
  this function computes the inverse of the schur complement of A,
    i.e. Ginv = (C A^{-1} B - D)^{-1},
  as well as A^{-1}B, and info for later applications of A^{-1}. It uses the
  KBPENTA solver for pentadiagonal systems, with only 3 non-zero diags
  This does it for every k \in [0, Ny * Nx]

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
    Ginv is overwritten with G^{-1} explicitly computed using 2x2 inv formula

  """
  libLinSolve.precomputeLinOps(A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                               A_kbp.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                               B.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                               C.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                               D.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                               G.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                               Ginv.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),\
                               Nyx, Nz, num_threads)


def getDiagsP2M2(A):
  n = A.shape[1]
  d = np.diag(A)
  dp2 = np.diag(A, k = 2)
  dm2 = np.diag(A, k = -2)
  return np.concatenate((d, dp2, dm2), axis = None)

###################################################################################
##################### Library function declarations ###############################
###################################################################################

# declare lin solver lib funcs
libLinSolve.precomputeLinOps.argtypes = [ctypes.POINTER(ctypes.c_double),\
                                         ctypes.POINTER(ctypes.c_double),\
                                         ctypes.POINTER(ctypes.c_double),\
                                         ctypes.POINTER(ctypes.c_double),\
                                         ctypes.POINTER(ctypes.c_double),\
                                         ctypes.POINTER(ctypes.c_double),\
                                         ctypes.POINTER(ctypes.c_double),\
                                         ctypes.c_int, ctypes.c_int, ctypes.c_int]
libLinSolve.precomputeLinOps.restype = None

# declare dptools lib funcs
libDPTools.Solver.argtypes = [np.ctypeslib.ndpointer(np.complex128),\
                              np.ctypeslib.ndpointer(np.complex128),\
                              ctypes.POINTER(ctypes.c_double),\
                              np.ctypeslib.ndpointer(np.complex128),\
                              np.ctypeslib.ndpointer(np.complex128),\
                              np.ctypeslib.ndpointer(np.complex128),\
                              np.ctypeslib.ndpointer(np.complex128),\
                              np.ctypeslib.ndpointer(np.complex128),\
                              ctypes.POINTER(ctypes.c_double),\
                              ctypes.c_double, ctypes.c_double, 
                              ctypes.c_uint, ctypes.c_uint,\
                              ctypes.c_uint, ctypes.c_uint, ctypes.c_int]

libDPTools.Solver.restype = ctypes.c_void_p

libDPTools.SetRHS.argtypes = [ctypes.c_void_p, np.ctypeslib.ndpointer(np.complex128)]
libDPTools.SetRHS.restype = None
libDPTools.SetRHS_split.argtypes = [ctypes.c_void_p, np.ctypeslib.ndpointer(np.double),\
                                    np.ctypeslib.ndpointer(np.double)]
libDPTools.SetRHS_split.restype = None
  
libDPTools.ZeroInit.argtypes = [ctypes.c_void_p]
libDPTools.ZeroInit.restype = None

libDPTools.DoublyPeriodicSolve_no_wall.argtypes = [ctypes.c_void_p]
libDPTools.DoublyPeriodicSolve_no_wall.restype = None
  
libDPTools.GetU_real.argtypes = [ctypes.c_void_p]
libDPTools.GetU_real.restype = ctypes.POINTER(ctypes.c_double)
libDPTools.GetU_imag.argtypes = [ctypes.c_void_p]
libDPTools.GetU_imag.restype = ctypes.POINTER(ctypes.c_double)
libDPTools.GetP_real.argtypes = [ctypes.c_void_p]
libDPTools.GetP_real.restype = ctypes.POINTER(ctypes.c_double)
libDPTools.GetP_imag.argtypes = [ctypes.c_void_p]
libDPTools.GetP_imag.restype = ctypes.POINTER(ctypes.c_double)
libDPTools.GetP_RHS_real.argtypes = [ctypes.c_void_p, ctypes.c_uint]
libDPTools.GetP_RHS_real.restype = ctypes.POINTER(ctypes.c_double)
libDPTools.GetP_RHS_imag.argtypes = [ctypes.c_void_p, ctypes.c_uint]
libDPTools.GetP_RHS_imag.restype = ctypes.POINTER(ctypes.c_double)
libDPTools.GetU_RHS_real.argtypes = [ctypes.c_void_p, ctypes.c_uint, ctypes.c_bool]
libDPTools.GetU_RHS_real.restype = ctypes.POINTER(ctypes.c_double)
libDPTools.GetU_RHS_imag.argtypes = [ctypes.c_void_p, ctypes.c_uint, ctypes.c_bool]
libDPTools.GetU_RHS_imag.restype = ctypes.POINTER(ctypes.c_double)

libDPTools.GetPhi_out_real.argtypes = [ctypes.c_void_p]
libDPTools.GetPhi_out_real.restype = ctypes.POINTER(ctypes.c_double)
libDPTools.GetPhi_out_imag.argtypes = [ctypes.c_void_p]
libDPTools.GetPhi_out_imag.restype = ctypes.POINTER(ctypes.c_double)


libDPTools.SetP.argtypes = [ctypes.c_voidp, np.ctypeslib.ndpointer(np.complex128),\
                            ctypes.c_uint]
libDPTools.SetP.restype = None
libDPTools.SetdP.argtypes = [ctypes.c_void_p, np.ctypeslib.ndpointer(np.complex128),\
                            ctypes.c_uint]
libDPTools.SetdP.restype = None
libDPTools.SetU.argtypes = [ctypes.c_void_p, np.ctypeslib.ndpointer(np.complex128),\
                            ctypes.c_uint]
libDPTools.SetU.restype = None

  
libDPTools.Clean.argtypes = [ctypes.c_void_p]
libDPTools.Clean.restype = None

libDPTools.evalTheta.argtypes = [ctypes.c_void_p,\
                                 ctypes.c_double]
libDPTools.evalTheta.restype = None

libDPTools.evalCorrectionSol_bottomWall.argtypes = [ctypes.c_void_p,\
                                                    ctypes.POINTER(ctypes.c_double),\
                                                    ctypes.POINTER(ctypes.c_double),\
                                                    ctypes.POINTER(ctypes.c_double),\
                                                    ctypes.c_double, ctypes.c_uint,\
                                                    ctypes.c_uint, ctypes.c_uint, ctypes.c_int]
libDPTools.evalCorrectionSol_bottomWall.restype = None

libDPTools.evalCorrectionSol_slitChannel.argtypes = [ctypes.c_void_p,\
                                                     ctypes.POINTER(ctypes.c_double),\
                                                     ctypes.POINTER(ctypes.c_double),\
                                                     ctypes.POINTER(ctypes.c_double),\
                                                     ctypes.c_double, ctypes.c_double,
                                                     ctypes.c_uint, ctypes.c_uint,\
                                                     ctypes.c_uint, ctypes.c_int]
libDPTools.evalCorrectionSol_slitChannel.restype = None

libDPTools.DoublyPeriodicSolve_bottom_wall.argtypes = [ctypes.c_void_p,\
                                                    ctypes.POINTER(ctypes.c_double),\
                                                    ctypes.POINTER(ctypes.c_double),\
                                                    ctypes.POINTER(ctypes.c_double),\
                                                    ctypes.c_double, ctypes.c_uint,\
                                                    ctypes.c_uint, ctypes.c_uint, ctypes.c_int]
libDPTools.DoublyPeriodicSolve_bottom_wall.restype = None
 
libDPTools.DoublyPeriodicSolve_slit_channel.argtypes = [ctypes.c_void_p,\
                                                     ctypes.POINTER(ctypes.c_double),\
                                                     ctypes.POINTER(ctypes.c_double),\
                                                     ctypes.POINTER(ctypes.c_double),\
                                                     ctypes.c_double, ctypes.c_double,
                                                     ctypes.c_uint, ctypes.c_uint,\
                                                     ctypes.c_uint, ctypes.c_int]
libDPTools.DoublyPeriodicSolve_slit_channel.restype = None
