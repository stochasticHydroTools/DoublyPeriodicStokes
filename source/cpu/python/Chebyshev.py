import numpy as np
from scipy.sparse import diags

def clencurt(Np1, a, b):
  """
  Clenshaw-curtis nodes cpts and weights cwts,
  i.e. the Chebyshev points of the second kind, and
  associated weights. These are such that
  f(cpts) \cdot cwts = \int_a^b f(x) dx (they are rescaled to [a,b])
  The implementation follows that given in ATAP by Trefethen
  
  Parameters:
    Np1 - number of Chebyshev points and weights
    a, b - left and right endpoints for interval [a,b]

  Returns:
    cpts - Chebyshev points of the second kind on [a,b]
    cwts - clenshaw-curtis quadrature weights 
  """
  N = Np1 - 1; H = (b - a) / 2; bpad2 = (b + a) / 2;
  v = np.zeros((N-1,), dtype = np.double)
  cpts = np.zeros((Np1,), dtype = np.double)
  cwts = np.zeros((Np1,), dtype = np.double)
  if N % 2:
    w1 = H / (N * N); end = int((N - 1) / 2 + 1)
  else:
    w1 = H / (N * N - 1); end = int(N / 2)
  # compute cheb points
  cpts = H * np.cos(np.pi * np.arange(0,Np1) / N) + bpad2 
  # now do weights
  cwts[0] = w1; cwts[N] = w1;
  v[0:N] = 1.0; iis = np.arange(1, N)
  for k in range(1, end):
    v -= 2 * np.cos(2.0 * k * np.pi * iis / N) / (4.0 * k * k - 1.0)
  if N % 2 == 0:
    v -= np.cos(N * np.pi * iis / N) / (N * N -1)
  cwts[1:N] = 2 * H * v / N 
  return cpts, cwts 

def chebCoeffDiff(fhat, Nx, Ny, Nz, dof, H):
  """
  Compute the Fourier-Chebyshev coefficients of the 
  z-derivative of fhat.
  
  Parameters:
    fhat (complex) - Fourier-Chebyshev coefficients of grid
                     function f, given sampling on a Chebyshev 
                     grid in z, and uniform in x,y. This is the
                     combined real and complex part, i.e.
                     fhat = fhat_r + 1j * fhat_i
    Nx,Ny,Nz - number of points in x,y,z
    dof - degrees of freedom of fhat (num components of vec field)
    H - (b-a)/2 for z \in [a,b]
  
  Returns:
    Df (complex) - Fourier-Chebyshev coefficients of dfhat/dz. 
                   This the combined real and complex part
         
  """
  Df = np.zeros((Nz, Ny, Nx, dof), dtype = np.complex)  
  fhat_rs = np.reshape(fhat, (Nz, Ny, Nx, dof))
  Df[-2,:,:,:] = 2 / H * (Nz - 1) * fhat_rs[-1,:,:,:] 
  for j in range(2,Nz):
    Df[Nz-j-1,:,:,:] = Df[Nz-j+1,:,:,:] + 2 / H * (Nz - j) * fhat_rs[Nz-j,:,:,:]
  Df[0,:,:,:] /= 2
  return Df

def chebCoeffDiff_perm(fhat, Nx, Ny, Nz, dof, H):
  """
  Compute the Fourier-Chebyshev coefficients of the 
  z-derivative of fhat for permuted layout.
  
  Parameters:
    fhat (complex) - Fourier-Chebyshev coefficients of grid
                     function f, given sampling on a Chebyshev 
                     grid in z, and uniform in x,y. This is the
                     combined real and complex part, i.e.
                     fhat = fhat_r + 1j * fhat_i
    Nx,Ny,Nz - number of points in x,y,z
    dof - degrees of freedom of fhat (num components of vec field)
    H - (b-a)/2 for z \in [a,b]
  
  Returns:
    Df (complex) - Fourier-Chebyshev coefficients of dfhat/dz. 
                   This the combined real and complex part
         
  """
  Df = np.zeros((Ny, Nx, Nz, dof), dtype = np.complex)  
  fhat_rs = np.reshape(fhat, (Ny, Nx, Nz, dof))
  Df[:,:,-2,:] = 2 / H * (Nz - 1) * fhat_rs[:,:,-1,:] 
  for j in range(2,Nz):
    Df[:,:,Nz-j-1,:] = Df[:,:,Nz-j+1,:] + 2 / H * (Nz - j) * fhat_rs[:,:,Nz-j,:]
  Df[:,:,0,:] /= 2
  return Df
  
def precomputeInts(Nz, H):
  """
  Precompute the integrals for the Nz Chebyshev modes that come up in
  the null-space A_p, A_x, A_y values for the pressure and velocity
  DP BVP solve
  
  Parameters:
    Nz - number of points on the z grid
    H - (b-a)/2 for z \in [a,b]
  
  Returns:
    pints - precomputed integrals of Chebyshev polynomials  
            for pressure correction
    uvints - precomputed integrals of Chebyshev polynomials
            for velocity correctionj
  """
  theta = np.pi * np.arange(0,1000).reshape(-1, 1) / 999
  _, zwts = clencurt(1000, 0, 2 * H)
  zwts = np.reshape(zwts, (1, -1))
  pints = np.dot(zwts, (np.cos(np.arange(0, Nz) * theta)))
  uvints = H * np.dot(zwts, (np.cos(np.arange(0, Nz) * theta) * np.cos(theta)))
  return pints, uvints

def firstIntegralMatrix(N, H):
  """
  Compute the first Chebyshev integration matrix. This maps
  Chebyshev coefficients of a function to the coefficients
  of the integral of the function
  
  Parameters:
    N - number of Chebyshev points
    H - (b-a)/2 for a Chebyshev grid [a,b]

  Returns:
    FIMat - the first Chebyshev integral matrix, scaled by H
  """
  jj = np.arange(2,N)
  colm1 = np.concatenate(([1], 1 / (2 * jj)))
  colp1 = np.concatenate(([0, -1 / 2], -1 / (2 * jj) * (jj < N - 1)))
  FIMat = diags([colm1, colp1], [-1, 1], shape = (N, N+2)).toarray(order='F')
  FIMat[0,N+1] = 1
  return H * FIMat

def secondIntegralMatrix(N, H):
  """
  Compute the second Chebyshev integration matrix. This maps
  Chebyshev coefficients of a function to the coefficients
  of the second integral of the function
  
  Parameters:
    N - number of Chebyshev points
    H - (b-a)/2 for a Chebyshev grid [a,b]

  Returns:
    SIMat - the second Chebyshev integral matrix, scaled by H^2
  """

  jj = np.arange(3, N)
  colm2 = np.concatenate(([0.25], 1 / (2 * jj * (2 * jj - 2))))
  colp2 = np.concatenate(([0, 0.125, 1./24.], 1 / (2 * jj * (2 * jj + 2)) * (jj < N - 2)))
  col0  = np.concatenate(([0, -0.125, -1. / 8. - 1. / 24.], -1 / (2 * jj * (2 * jj - 2)) \
                          - 1 / (2 * jj * (2 * jj + 2)) * (jj < N - 1))) 
  SIMat = diags([colm2, col0, colp2], [-2, 0, 2], shape = (N,N+2)).toarray(order='F')
  SIMat[0,N] = 1; SIMat[1,N+1] = 1;
  return H**2 * SIMat
