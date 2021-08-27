import sys
import random
# import Python modules wrapping C libraries (and also numpy)
from Grid import *
from Particles import *
from SpreadInterp import *
from Transform import *
from Chebyshev import *
from Solvers import Stokes
from GridAndKernelConfig import configure_grid_and_kernels_xy, configure_grid_and_kernels_z
from timeit import default_timer as timer
"""
Script for benchmarking
"""
def mobility(xP, fP, tP, eta, Lx, Ly, minPz, maxPz, radP, kernTypes, solverType, has_torque, write=False, ref=False):
  """
  Compute the mobility of particles with positions xP, forces fP and torques tP. If has torque=True,
  both the linear and angular velocity will be returned. This function does not use kernel derivatives.

  Parameters:
    xP,fP,tP - positions, forces and torques(can be none) on particles
    eta - viscosity
    Lx, Ly - approximate length in x,y
    minPz,maxPz - min/max particle height, or approximate location of bottom/top boundary
    radP - particle radii
    kernTypes - array of types of kernel for each particle
             - 0 for ES 6 pt (both monopole and dipole if has_torque = true) 
             - 1 for ES 5 pt (both monopole and dipole if has_torque = true)
             - 2 for ES 4 pt (only monopole, dipole not supported) 
             - 3 for ES 6 pt monopole, ES 5 pt dipole
             - 4 for ES 6 pt monopole, ES 4 pt dipole
             - 5 for ES 5 pt monopole, ES 4 pt dipole 
    solverType :
     'TP' - triply periodic
     'DP' - doubly periodic
     'DPBW' - doubly periodic with bottom wall
     'DPSC - doubly periodic slit channel
    has_torque - boolean to specify whether problem involves torque (default is False)
    userInp - selects optimal kernel by default and toggles user input 
    write - write the particle positions, forces and velocities (for dipole as well if has_torque = True)
    ref - if true, doubles the support underneath the kernel (by doubling the grid) for a more resolved solution


  Returns: vP, omegaP - the linear and angular velocity of the particles (latter is 0 if has_torque=False)
  """

  # select x-y grid for the particles (if optInd=-1, optimal adjusted grid automatically chosen)
  Lx, Ly, hx, hy, Nx, Ny, wfP,\
  wdP, cbetam, cbetad, betafP, betatP\
    = configure_grid_and_kernels_xy(Lx, Ly, radP, kernTypes, 0, has_torque)
  # set the z grid 
  Lz, hz, Nz, z0\
    = configure_grid_and_kernels_z(minPz, maxPz, hx, wfP, wdP, solverType, has_torque)
  # chebyshev grid and weights for z
  zpts, zwts = clencurt(Nz, 0, Lz)
  nP = kernTypes.size; dof = 3

  vP = np.zeros((dof * nP,), dtype = np.double); omegaP = 0
  if has_torque:
    curlT = np.zeros((Nx * Ny * Nz * dof,), dtype = np.complex)
    omegaP = np.zeros((nP * dof))
  else:
    omegaP = 0
  # grid periodicity
  periodic_x = periodic_y = periodic_z = True;
  # boundary conditions specified for ends of each axis 
  # 0 - mirror wall, 1 - inverse mirror wall , 2 - none 
  # start out with none (works for TP, DP)
  BCs = 2 * np.ones(dof * 6, dtype = np.uintc)
  # correct BCs based on domain
  if solverType != 'TP':
    periodic_z = False
  if solverType == 'DPBW':
    # apply mirror-inv on bottom wall only for each solution component
    BCs[5 * dof] = BCs[5 * dof + 1] = BCs[5 * dof + 2] = 1
  elif solverType == 'DPSC':
    # apply mirror-inv on both walls for each solution component
    BCs[4 * dof] = BCs[4 * dof + 1] = BCs[4 * dof + 2] = 1
    BCs[5 * dof] = BCs[5 * dof + 1] = BCs[5 * dof + 2] = 1
  elif solverType != 'TP' and solverType != 'DP':
    raise ValueError("solverType not supported")
  # choose what to do with k=0 in the no wall solver
  # k0 = 0 - the k=0 mode of the pressure and velocity will be 0
  # k0 = 1 - the k=0 mode of the pressure and velocity will not be 0
  #        - there will be a bvp solve and correction to the k=0 mode after each solve
  k0 = 0
  solver = Stokes(Nx, Ny, Nz, Lx, Ly, Lz, dof, eta, k0, solverType)
  if solverType != 'TP':
    solver.SetZ(zpts)
    # plan fftws
    transformer = Transformer(Nx, Ny, Nz, dof, 1)
    if has_torque:
      transformerT = Transformer(Nx, Ny, Nz, dof, 1)
  else:
    transformer = Transformer(Nx, Ny, Nz, dof, 0)
    if has_torque:
      transformerT = Transformer(Nx, Ny, Nz, dof, 0)
    zpts = zwts = None
 
  # make the grid and set up the particles on it
  grid = GridGen(Lx, Ly, Lz, 0.0, 0.0, z0, hx, hy, hz, Nx, Ny, Nz, dof, periodic_x, periodic_y, periodic_z, BCs, zpts, zwts)
  grid.Make()
  particles = ParticlesGen(nP, dof, xP, fP, radP, wfP, cbetam, betafP)
  particles.UseCbeta()
  particles.Make()
  particles.Setup(grid)
  if write:
    particles.WriteParticles('monopole_position_and_force.txt')

  # spread forces and transform to Fourier-Cheb
  F = Spread(particles, grid); 
  transformer.Ftransform(F)
  F_hat_r = transformer.out_real; tG_hat_r = 0;
  F_hat_i = transformer.out_imag; tG_hat_i = 0;
 
  if has_torque:
    # make torque grid and set up torque particles on it
    tgrid = GridGen(Lx, Ly, Lz, 0.0, 0.0, z0, hx, hy, hz, Nx, Ny, Nz, dof, periodic_x, periodic_y, periodic_z, BCs, zpts, zwts)
    tgrid.Make()
    tparticles = ParticlesGen(nP, dof, xP, tP, radP, wdP, cbetad, betatP)
    tparticles.UseCbeta()
    tparticles.Make()
    tparticles.Setup(tgrid)
    if write:
      tparticles.WriteParticles('dipole_position_and_torque.txt')
    # spread torques, transform to Fourier-cheb and compute curl
    tG = Spread(tparticles, tgrid)
    transformerT.Ftransform(tG)
    if solverType == 'TP':
      tG_hat = transformerT.out_real + 1j * transformerT.out_imag
      dyTz = solver.Dy * tG_hat[2::3]
      dxTz = solver.Dx * tG_hat[2::3]
      dxTy = solver.Dx * tG_hat[1::3]
      dyTx = solver.Dy * tG_hat[0::3]
      dzTy = solver.Dz * tG_hat[1::3]
      dzTx = solver.Dz * tG_hat[0::3]
    else: 
      tG_hat = (transformerT.out_real + 1j * transformerT.out_imag).reshape((Ny * Nx, Nz, dof))
      solver.Dy.shape = solver.Dx.shape = (Ny * Nx, 1) 
      dyTz = (solver.Dy * tG_hat[:,:,2]).reshape((Nx * Ny * Nz,))
      dxTz = (solver.Dx * tG_hat[:,:,2]).reshape((Nx * Ny * Nz,))
      dxTy = (solver.Dx * tG_hat[:,:,1]).reshape((Nx * Ny * Nz,))
      dyTx = (solver.Dy * tG_hat[:,:,0]).reshape((Nx * Ny * Nz,))
      dzTy = chebCoeffDiff_perm(tG_hat[:,:,1], Nx, Ny, Nz, 1, Lz / 2).reshape((Nx * Ny * Nz,)) 
      dzTx = chebCoeffDiff_perm(tG_hat[:,:,0], Nx, Ny, Nz, 1, Lz / 2).reshape((Nx * Ny * Nz,)) 
    curlT[0::3] = 0.5 * (dyTz - dzTy)
    curlT[1::3] = 0.5 * (dzTx - dxTz)
    curlT[2::3] = 0.5 * (dxTy - dyTx)
    tG_hat_r = np.real(curlT); tG_hat_i = np.imag(curlT) 
  # combine force/torque RHS
  F_hat_r += tG_hat_r
  F_hat_i += tG_hat_i
  # solve pde and backtransform to get linear velocities on grid
  if not has_torque:
    U_hat_r, U_hat_i, _, _ = solver.Solve(F_hat_r, F_hat_i, has_torque) 
    transformer.Btransform(U_hat_r, U_hat_i)
    uG_r = transformer.out_real
    uG_i = transformer.out_imag
    # interpolate linear velocities on particles
    grid.SetData(uG_r)
    Interpolate(particles, grid, vP)
    if write:
      particles.WriteParticles('monopole_position_and_velocity.txt')
    grid.Clean()
    particles.Clean()
    transformer.Clean()
  else:
    U_hat_r, U_hat_i, Ut_hat_r, Ut_hat_i, _, _ = solver.Solve(F_hat_r, F_hat_i, has_torque) 
    transformer.Btransform(U_hat_r, U_hat_i)
    transformerT.Btransform(Ut_hat_r, Ut_hat_i)
    uG_r = transformer.out_real
    uG_i = transformer.out_imag
    utG_r = transformerT.out_real
    utG_i = transformerT.out_imag
    # interpolate linear velocities on particles
    grid.SetData(uG_r)
    Interpolate(particles, grid, vP)
    # interpolate angular velocities on particles
    tgrid.SetData(utG_r)
    Interpolate(tparticles, tgrid, omegaP)
    if write:
      particles.WriteParticles('monopole_position_and_velocity.txt')
      tparticles.WriteParticles('dipole_position_and_velocity.txt')
    grid.Clean()
    particles.Clean()
    tgrid.Clean()
    tparticles.Clean()
    transformer.Clean()
    transformerT.Clean()
  solver.Clean()
  return vP, omegaP


def mobility_derivkern(xP, fP, tP, eta, Lx, Ly, minPz, maxPz, radP, kernTypes, solverType, has_torque, write=False, ref=False, bench=False, warmup=10,trials=20):
  """
  Compute the mobility of particles with positions xP, forces fP and torques tP. If has torque=True,
  both the linear and angular velocity will be returned. The latter will use the derivatives of the kernel for spreading and interpolation.

  Parameters:
    xP,fP,tP - positions, forces and torques(can be none) on particles
    eta - viscosity
    Lx, Ly - approximate length in x,y
    minPz,maxPz - min/max particle height, or approximate location of bottom/top boundary
    radP - particle radii
    kernTypes - array of types of kernel for each particle
             - 0 for ES 6 pt (both monopole and dipole if has_torque = true) 
             - 1 for ES 5 pt (both monopole and dipole if has_torque = true)
             - 2 for ES 4 pt (only monopole, dipole not supported) 
             - 3 for ES 6 pt monopole, ES 5 pt dipole
             - 4 for ES 6 pt monopole, ES 4 pt dipole
             - 5 for ES 5 pt monopole, ES 4 pt dipole 
    solverType :
     'TP' - triply periodic
     'DP' - doubly periodic
     'DPBW' - doubly periodic with bottom wall
     'DPSC - doubly periodic slit channel
    has_torque - boolean to specify whether problem involves torque (default is False)

  Returns: vP, omegaP - the linear and angular velocity of the particles (latter is 0 if has_torque=False)
  """
  # select x-y grid for the particles (if optInd=-1, optimal adjusted grid automatically chosen)
  Lx, Ly, hx, hy, Nx, Ny, wfP,\
  wdP, cbetam, cbetad, betafP, betatP\
    = configure_grid_and_kernels_xy(Lx, Ly, radP, kernTypes, 0, has_torque)
  # set the z grid 
  Lz, hz, Nz, z0\
    = configure_grid_and_kernels_z(minPz, maxPz, hx, wfP, wdP, solverType, has_torque)
  # chebyshev grid and weights for z
  zpts, zwts = clencurt(Nz, 0, Lz)
  ####### print final settings #######
  dispstr = '\nFinal grid settings for '
  if solverType == 'TP':
    dispstr = dispstr + 'Triply periodic domain:\n'
  elif solverType == 'DP':
    dispstr = dispstr + 'No wall domain:\n'
  elif solverType == 'DPBW':
    dispstr = dispstr + 'Bottom wall domain:\n'
  else:
    dispstr = dispstr + 'Two wall domain:\n'
  print(dispstr)
  print('\t (Nx, Ny, Nz) = (%d, %d, %d)' % (Nx, Ny, Nz)) 
  print('\t hxy = %.16f' % (hx))
  print('\t Lx = %.16f' % (Lx)) 
  print('\t Ly = %.16f' % (Ly))
  print('\t z_a = %.16f \n \t z_b = %.16f' % (0, Lz))
  nP = kernTypes.size; dof = 3
  # grid periodicity
  periodic_x = periodic_y = periodic_z = True;
  # boundary conditions specified for ends of each axis 
  # 0 - mirror wall, 1 - inverse mirror wall , 2 - none 
  # start out with none (works for TP, DP)
  BCs = 2 * np.ones(dof * 6, dtype = np.uintc)
  # correct BCs based on domain
  if solverType != 'TP':
    periodic_z = False
  if solverType == 'DPBW':
    # apply mirror-inv on bottom wall only for each solution component
    BCs[5 * dof] = BCs[5 * dof + 1] = BCs[5 * dof + 2] = 1
  elif solverType == 'DPSC':
    # apply mirror-inv on both walls for each solution component
    BCs[4 * dof] = BCs[4 * dof + 1] = BCs[4 * dof + 2] = 1
    BCs[5 * dof] = BCs[5 * dof + 1] = BCs[5 * dof + 2] = 1
  elif solverType != 'TP' and solverType != 'DP':
    raise ValueError("solverType not supported")
  # choose what to do with k=0 in the no wall solver
  # k0 = 0 - the k=0 mode of the pressure and velocity will be 0
  # k0 = 1 - the k=0 mode of the pressure and velocity will not be 0
  #        - there will be a bvp solve and correction to the k=0 mode after each solve
  k0 = 0
  solver = Stokes(Nx, Ny, Nz, Lx, Ly, Lz, dof, eta, k0, solverType)
  if solverType != 'TP':
    solver.SetZ(zpts)
    # plan fftws
    transformer = Transformer(Nx, Ny, Nz, dof, 1)
  else:
    transformer = Transformer(Nx, Ny, Nz, dof, 0)
    zpts = zwts = None 
  # make the grid and set up the particles on it
  grid = GridGen(Lx, Ly, Lz, 0.0, 0.0, z0, hx, hy, hz, Nx, Ny, Nz, dof, periodic_x, periodic_y, periodic_z, BCs, zpts, zwts)
  grid.Make()
  particles = ParticlesGen(nP, dof, xP, fP, radP, wfP, cbetam, betafP)
  particles.UseCbeta()
  particles.Make()
  particles.Setup(grid)
  vP = np.zeros((dof * nP,), dtype = np.double); omegaP = 0
  if has_torque:
    # make torque grid and set up torque particles on it
    tgrid = GridGen(Lx, Ly, Lz, 0.0, 0.0, z0, hx, hy, hz, Nx, Ny, Nz, dof, periodic_x, periodic_y, periodic_z, BCs, zpts, zwts)
    tgrid.Make()
    tparticles = ParticlesGen(nP, dof, xP, tP, radP, wdP, cbetad, betatP)
    tparticles.UseCbeta()
    tparticles.IsDipole() # need to specify these are dipoles
    tparticles.Make()
    tparticles.Setup(tgrid)
    # define  mirror BC for spreading z deriv  
    BCst_z = 2 * np.ones(dof * 6, dtype = np.uintc)
    if solverType != 'TP' and solverType != 'DP':
      BCst_z[5 * dof] = BCst_z[5 * dof + 1] = BCst_z[5 * dof + 2] = 0
      if solverType == 'DPSC':
        BCst_z[4 * dof] = BCst_z[4 * dof + 1] = BCst_z[4 * dof + 2] = 0
    omegaP = np.zeros((nP * dof), dtype = np.double)
    dxvP = np.zeros((nP * dof), dtype = np.double)
    dyvP = np.zeros((nP * dof), dtype = np.double)
    dzvP = np.zeros((nP * dof), dtype = np.double)
    

  if not bench:
    trials = 1 

  ttot = tSI = tSI_t = tT = tS = 0

  for j in range(0,trials):

    t0 = timer()

    # spread forces
    F = Spread(particles, grid); 
    
    t1 = timer()
    
    if j >= warmup:
      ttot += t1 - t0
      tSI += t1 - t0

    if has_torque:
     
      t0 = timer()

      # spread derivative of torques, compute curl and add to F 
      dxT = SpreadDx(tparticles, tgrid)
      libSpreadInterp.addDx(grid.grid, tgrid.grid)
      dyT = SpreadDy(tparticles, tgrid)
      libSpreadInterp.addDy(grid.grid, tgrid.grid)
      tgrid.SetBCs(BCst_z)
      dzT = SpreadDz(tparticles, tgrid)
      libSpreadInterp.addDz(grid.grid, tgrid.grid)
      
      t1 = timer()
     
      if j >= warmup:
        ttot += t1 - t0
        tSI_t += t1 - t0

 
    t0 = timer() 
    
    # forward transform
    transformer.Ftransform(F)
    
    t1 = timer()
    
    if j >= warmup:
      ttot += t1 - t0
      tT += t1 - t0
    
    F_hat_r = transformer.out_real; 
    F_hat_i = transformer.out_imag; 



    t0 = timer()

    # solve pde 
    U_hat_r, U_hat_i, _, _ = solver.Solve(F_hat_r, F_hat_i)

    t1 = timer()

    if j >= warmup:
      ttot += t1-t0
      tS += t1-t0
     
    t0 = timer()
  
    # back transform
    transformer.Btransform(U_hat_r, U_hat_i)

    t1 = timer()

    if j >= warmup:
      ttot += t1 - t0
      tT += t1 - t0

    uG_r = transformer.out_real

    t0 = timer()

    # interpolate linear velocities on particles
    grid.SetData(uG_r)
    Interpolate(particles, grid, vP)

    t1 = timer()

    if j >= warmup:
      ttot += t1 - t0
      tSI += t1 - t0

    # reset forces
    particles.SetData(fP)

    if has_torque:

      t0 = timer()
      
      # interpolate derivative of linear velocities on particles
      tgrid.SetData(uG_r)
      # set BC for x-y deriv
      tgrid.SetBCs(BCs)
      InterpolateDx(tparticles, tgrid, dxvP)
      InterpolateDy(tparticles, tgrid, dyvP)
      # set BC for z deriv
      tgrid.SetBCs(BCst_z)
      InterpolateDz(tparticles, tgrid, dzvP)
      # reset boundary condition, since we changed it for z
      tgrid.SetBCs(BCs)
      # compute angular velocity of particles
      omegaP[0::3] = -1/2 * (dyvP[2::3] - dzvP[1::3])
      omegaP[1::3] = -1/2 * (dzvP[0::3] - dxvP[2::3]) 
      omegaP[2::3] = -1/2 * (dxvP[1::3] - dyvP[0::3]) 
      t1 = timer()
     
      if j >= warmup:
        ttot += t1 - t0
        tSI_t += t1 - t0
        print("Average time: ", ttot / (j+1-warmup))
      
      # reset torques
      tparticles.SetData(tP)

  grid.Clean()
  particles.Clean()
  transformer.Clean()
  solver.Clean()
  if has_torque:
      tgrid.Clean()
      tparticles.Clean()

  if bench:
    print("\tSpread Interp Force     : ", tSI / (trials-warmup))
    print("\tSpread Interp Torque    : ", tSI_t / (trials-warmup))
    print("\tFCT+iFCT                : ", tT / (trials-warmup))
    print("\tSolve                   : ", tS / (trials-warmup))
    print("\tAverage time              : ", ttot / (trials-warmup))
  return vP, omegaP


def mobility_matrix(xP, eta, Lx, Ly, minPz, maxPz, radP, kernTypes, solverType, ref=False, has_torque=True):
  """
  Compute the mobility matrix for each configuration in xP(nConfig, 3 * nP)

  Parameters:
    xP - particle positions (3*nP)
    eta - viscosity
    Lx,Ly - length of unit cell in x,y
    minPz, maxPz - min and max possible particle height
    radP - partice radii
    kernTypes - array of types of kernel for each particle
             - 0 for ES 6 pt (both monopole and dipole if has_torque = true) 
             - 1 for ES 5 pt (both monopole and dipole if has_torque = true)
             - 2 for ES 4 pt (only monopole, dipole not supported) 
             - 3 for ES 6 pt monopole, ES 5 pt dipole
             - 4 for ES 6 pt monopole, ES 4 pt dipole
             - 5 for ES 5 pt monopole, ES 4 pt dipole 
    solverType :
     'TP' - triply periodic
     'DP' - doubly periodic
     'DPBW' - doubly periodic with bottom wall
     'DPSC - doubly periodic slit channel
    ref - (boolean) if true, computation will be a reference on a doubled grid 
   
  """

  nConfig, nP3 = xP.shape; nP = int(nP3 / 3)
  # construct mobility matrix at each height
  if has_torque:
    M = np.zeros((nConfig, 6 * nP, 6 * nP), dtype = np.double)
    F = np.eye(6 * nP, dtype = np.double)
    NormMat = np.zeros((6 * nP, 6 * nP))
    NormMat[0:3*nP, 0:3*nP] = 1 / (6 * np.pi * eta * np.repeat(radP, 3)).reshape((1, 3 * nP))
    NormMat[0:3*nP, 3*nP::] = NormMat[3*nP::, 0:3*nP] = 1 / (6 * np.pi * eta * np.repeat(np.power(radP, 2), 3)).reshape(1, 3 * nP)
    NormMat[3*nP::, 3*nP::] = 1 / (8 * np.pi * eta * np.repeat(np.power(radP, 3), 3)).reshape(1, 3 * nP)
  else:
    M = np.zeros((nConfig, 3 * nP, 3 * nP), dtype = np.double)
    F = np.eye(3 * nP, dtype = np.double)
    NormMat = np.zeros((3 * nP, 3 * nP))
    NormMat[0:3*nP, 0:3*nP] = 1 / (6 * np.pi * eta * np.repeat(radP, 3)).reshape((1, 3 * nP))
  for iC in range(0, nConfig):
    print("config #", iC)
    if has_torque:
      for j in range(0, 6 * nP):
        mob, mobt = mobility(xP[iC,:].copy(), F[0:3*nP,j].copy(), F[3*nP::,j].copy(), eta, Lx, Ly, minPz, maxPz, radP, kernTypes, solverType, has_torque, False, ref)
        M[iC,:,j] = np.concatenate((mob, mobt))
    else:
      for j in range(0, 3 * nP):
        mob, _ = mobility(xP[iC,:].copy(), F[0:3*nP,j].copy(), None, eta, Lx, Ly, minPz, maxPz, radP, kernTypes, solverType, has_torque, False, ref)
        M[iC,:,j] = mob
    M[iC,:,:] /= NormMat

  return M, NormMat
