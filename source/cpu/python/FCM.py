from Grid import *
from Particles import *
from SpreadInterp import *
from Transform import *
from Chebyshev import clencurt
from GridAndKernelConfig import configure_grid_and_kernels_xy, configure_grid_and_kernels_z
from Solvers import Stokes

class FCM(object):
  """
  Interface for Stokes mobility problems using the Force Coupling Method with the ES kernel 

  Attributes:
    domType - domain type for the problem
            - 'TP' for triply periodic
            - 'DP' for doubly periodic
            - 'DPBW' for doubly periodic with bottom wall
            - 'DPSC' for doubly periodic slit channel
    radP - array of particle radii
    kernTypes - array of types of kernel for each particle
             - 0 for ES 6 pt (both monopole and dipole if has_torque = true) 
             - 1 for ES 5 pt (both monopole and dipole if has_torque = true)
             - 2 for ES 4 pt (only monopole, dipole not supported) 
             - 3 for ES 6 pt monopole, ES 5 pt dipole
             - 4 for ES 6 pt monopole, ES 4 pt dipole
             - 5 for ES 5 pt monopole, ES 4 pt dipole 
    has_torque - if True, problem involves torque (default is False)
    min[max]X[Y,Z] - Intervals specifying the unit cell [minX,maxX]x[minY,maxY]x[minZ,maxZ]
    Lx[y,z] - maxX[Y,Z]-minX[Y,Z] side lengths of the unit cell
    xP - particle positions (x1,y1,z1,x2,y2,z2,..)
    hx[y,z] - grid spacing in x,y (z only for TP, otherwise hz=0)
    Nx[y,z] - number of points in x, y, z
    wm[d]   - number of points underneath monopole/dipole kernels of each particle
            - if has_torque = false, wd = 0
            - if has_torque = true, wm and wd are filled
    cbetam[d] - dimensionless radius for each kernelType as function of beta
    betam[d]P - ES kernel beta for monopole and dipole (betadP = 0 if has_torque=False)
    zpts - Chebyshev grid [0,Lz] (in reverse) (= 0 if domType == 'TP')
    zwts - Clenshaw-curtis wiehgts (= 0 if domType == 'TP')
    z0    - actual floor of z grid [a,b] which we shift to [0,b-a]
    periodic_x[yz] - periodicity in each direction
    BCs - boundary conditions
    BCst_z - boundary condition for z deriv ES kernel
    [t]grid - grid object (Python obj managing c pointer - see Grid.py) for force/torque
    [t]particles - particles object (Python obj managing c pointer - see Particles.py) for force/torque
    vP - linear particle velocity
    omegaP - rotational particle velocity
    d[xyz]vP - derivative of linear velocity (used to compute omegaP)
    F - spread forces/torque curl on grid
    solver - solver object (see Solvers.py)
    transformer - object for fourier/cheb transforms (see Transforms.py)
    F_hat_r[i] - real/imaginary part of spread force
    U_hat_r[i],P_hat_r[i] - coefs of transform of velocity and pressure on the grid
    uG_r - velocity on the grid

  Example Usage:
    
    nP = 2048
    domType = 'DPBW'; eta = 0.957e-3; has_torque = True
    minX = 0.0; maxX = 128.7923; minY = 0.0; maxZ = 128.7923
    minZ = 0.0; maxZ = 20.0
    radP = 1.0155 * np.ones(nP, dtype = np.double)
    kernTypes = np.zeros(nP, dtype = np.int)
    
    problem = FCM(radP, kernTypes, domType, has_torque)
    problem.SetUnitCell([minX,maxX], [minY,maxY], [minZ,maxZ])
    problem.Initialize(eta, 0)
    problem.SetPositions(xP)
    
    # compute lin/rot vel
    V = problem.Mdot(F)
    # update pos 
    xP += 0.1 * np.abs(V[0:3 * nP])
    problem.SetPositions(xP)
    # compute lin/rot vel 
    V = problem.Mdot(F)
    # deallocate
    problem.Clean()
  """  
  def __init__(self, radP, kernTypes, domType, has_torque = False):
    """
    Constructor for FCM class 
  
    Inputs:
      radP - array of radii 
      kernTypes - array of kernel types
      domType - domain type
      has_torque - bool specifying if problem has torque

    """
     
    if domType != 'TP' and domType != 'DP' and\
       domType != 'DPSC' and domType != 'DPBW':
      raise ValueError('Domain type is not supported')
    
    self.domType = domType
    self.nP, = kernTypes.shape
    
    if self.nP != radP.shape[0]:
      raise ValueError('Radii list and kernel types list must have the same length')  
    
    self.has_torque = has_torque  
 
    for j in range(0, self.nP): 
      if (kernTypes[j] > 5 or kernTypes[j] < 0) or (has_torque and kernTypes[j] == 2):
        print('Invalid kernel type for particle #', j)
        print('Valid types are:\n\
              - 0 for ES 6 pt (both monopole and dipole if has_torque = true)\n\
              - 1 for ES 5 pt (both monopole and dipole if has_torque = true)\n\
              - 2 for ES 4 pt (only monopole, dipole not supported)\n\
              - 3 for ES 6 pt monopole, ES 5 pt dipole\n\
              - 4 for ES 6 pt monopole, ES 4 pt dipole\n\
              - 5 for ES 5 pt monopole, ES 4 pt dipole\n')
        raise ValueError('Kernel types list is invalid.')
      if radP[j] < 0:
        print('Negative radius provided for particle #', j)
        raise ValueError('Radii list is invalid.')

    self.kernTypes = kernTypes
    self.radP = radP
    self.dof = 3
    self.posinit = False 
    self.particles = None
    self.tparticles = None
    self.grid = None
    self.tgrid = None
    self.transformer = None
    self.solver = None 
 
  def SetUnitCell(self, X, Y, Z):
    """
    Set the unit cell
    
    Inputs:
      X,Y,Z - intervals [minX,maxX] etc.
    """
    print('\nSetting the unit cell...')
    if X[0] > X[1] or Y[0] > Y[1] or Z[0] > Z[1]:
      raise ValueError('One of the frame intervals is invalid')
    self.minX = X[0]; self.maxX = X[1]
    self.minY = Y[0]; self.maxY = Y[1]
    self.minZ = Z[0]; self.maxZ = Z[1]
    self.Lx = X[1] - X[0] 
    self.Ly = Y[1] - Y[0] 
    self.Lz = Z[1] - Z[0]
    print('unit cell := [%.4f,%.4f]x[%.4f,%.4f]x[%.4f,%.4f]\n' % (X[0], X[1], Y[0], Y[1], Z[0], Z[1])) 


  def Initialize(self, viscosity, optInd, fac=1.5):
    """
    Initialize the grids, precompute the solver, and plan FFTs
    
    Inputs:
      viscosity - for the Stokes problem
      optInd - which of the candidate grids to select (see stdout)
      fac - buffer factor for DP problems (default/recommended is 1.5)
    Stdout: 
          - The optimal *adjusted* grid is displayed (optInd = -1)
          - Several candidate fft-friendly grids are displayed (optInd = 0,1,..) 
          - The final kernel settings for each particle size are displayed
    
    Side Effects: The following members are defined:
      hx[y,z] - grid spacing in x,y (z only for TP, otherwise hz=0)
      Nx[y,z] - number of points in x, y, z
      wm[d]   - number of points underneath monopole/dipole kernels of each particle
              - if has_torque = false, wd = 0
              - if has_torque = true, wm and wd are filled
      cbetam[d] - dimensionless radius for each kernelType as function of beta
      betam[d]P - ES kernel beta for monopole and dipole (betadP = 0 if has_torque=False)
      zpts - Chebyshev grid [0,Lz] (in reverse) (= 0 if domType == 'TP')
      zwts - Clenshaw-curtis wiehgts (= 0 if domType == 'TP')
      z0    - actual floor of z grid [a,b] which we shift to [0,b-a]
      periodic_x[yz] - periodicity in each direction
      BCs - boundary conditions
      BCst_z - boundary condition for z deriv ES kernel
      [t]grid - grid object (Python obj managing c pointer - see Grid.py) for force/torque
      [t]particles - particles object (Python obj managing c pointer - see Particles.py) for force/torque
    
    Side Effects:
      grid members are initialized
      solver member is initialized and memory allocated
      transformer members are initialized, memory allocated, and fftw plans created/loaded
    """
    print('Determining the grid and configuring kernels...\n')
    # select x-y grid for the particles (if optInd=-1, optimal adjusted grid automatically chosen)
    self.Lx, self.Ly, self.hx, self.hy, self.Nx, self.Ny, self.wm,\
    self.wd, self.cbetam, self.cbetad, self.betamP, self.betadP\
      = configure_grid_and_kernels_xy(self.Lx, self.Ly, self.radP, self.kernTypes, optInd, self.has_torque)
    # set the z grid 
    self.Lz, self.hz, self.Nz, self.z0\
      = configure_grid_and_kernels_z(self.minZ, self.maxZ, self.hx, self.wm, self.wd, self.domType, self.has_torque, fac)
    # chebyshev grid and weights for z
    self.zpts, self.zwts = clencurt(self.Nz, 0, self.Lz)
    ####### print final settings #######
    dispstr = '\nFinal grid settings for '
    if self.domType == 'TP':
      dispstr = dispstr + 'Triply periodic domain:\n'
    elif self.domType == 'DP':
      dispstr = dispstr + 'No wall domain:\n'
    elif self.domType == 'DPBW':
      dispstr = dispstr + 'Bottom wall domain:\n'
    else:
      dispstr = dispstr + 'Two wall domain:\n'
    print(dispstr)
    print('\t (Nx, Ny, Nz) = (%d, %d, %d)' % (self.Nx, self.Ny, self.Nz)) 
    print('\t hxy = %.16f' % (self.hx))
    print('\t Lx = %.16f' % (self.Lx)) 
    print('\t Ly = %.16f' % (self.Ly))
    print('\t z_a = %.16f \n \t z_b = %.16f\n' % (0, self.Lz))
    # grid periodicity
    self.periodic_x = self.periodic_y = self.periodic_z = True;
    # boundary conditions specified for ends of each axis 
    # 0 - mirror wall, 1 - inverse mirror wall , 2 - none 
    # start out with none (works for TP, DP)
    self.BCs = 2 * np.ones(self.dof * 6, dtype = np.uintc)
    # correct self.BCs based on domain
    if self.domType != 'TP':
      self.periodic_z = False
    if self.domType == 'DPBW':
      # apply mirror-inv on bottom wall only for each solution component
      self.BCs[5 * self.dof] = self.BCs[5 * self.dof + 1] = self.BCs[5 * self.dof + 2] = 1
    elif self.domType == 'DPSC':
      # apply mirror-inv on both walls for each solution component
      self.BCs[4 * self.dof] = self.BCs[4 * self.dof + 1] = self.BCs[4 * self.dof + 2] = 1
      self.BCs[5 * self.dof] = self.BCs[5 * self.dof + 1] = self.BCs[5 * self.dof + 2] = 1
    elif self.domType != 'TP' and self.domType != 'DP':
      raise ValueError("domType not supported")
    print('Initializing force grid...\n')
    # make the grid and set up the particles on it
    self.grid = GridGen(self.Lx, self.Ly, self.Lz, self.minX, self.minY, self.z0, self.hx, self.hy, self.hz,\
                        self.Nx, self.Ny, self.Nz, self.dof, self.periodic_x,\
                        self.periodic_y, self.periodic_z, self.BCs, self.zpts, self.zwts)
    self.grid.Make()
    if self.has_torque:
      # make torque grid and set up torque particles on it
      print('Initializing torque grid...\n')
      self.tgrid = GridGen(self.Lx, self.Ly, self.Lz, self.minX, self.minY, self.z0, self.hx, self.hy, self.hz,\
                           self.Nx, self.Ny, self.Nz, self.dof, self.periodic_x,\
                           self.periodic_y, self.periodic_z, self.BCs, self.zpts, self.zwts)
      self.tgrid.Make()
      # define  mirror BC for spreading z deriv  
      self.BCst_z = 2 * np.ones(self.dof * 6, dtype = np.uintc)
      if self.domType != 'TP' and self.domType != 'DP':
        self.BCst_z[5 * self.dof] = self.BCst_z[5 * self.dof + 1] = self.BCst_z[5 * self.dof + 2] = 0
        if self.domType == 'DPSC':
          self.BCst_z[4 * self.dof] = self.BCst_z[4 * self.dof + 1] = self.BCst_z[4 * self.dof + 2] = 0
    print('Precomputing solver and planning FFTs...\n')
    self.viscosity = viscosity
    # choose what to do with k=0 in the no wall solver
    # k0 = 0 - the k=0 mode of the pressure and velocity will be 0
    # k0 = 1 - the k=0 mode of the pressure and velocity will not be 0
    #        - there will be a bvp solve and correction to the k=0 mode after each solve
    k0 = 0
    self.solver = Stokes(self.Nx, self.Ny, self.Nz, self.Lx, self.Ly, self.Lz, self.dof, self.viscosity, k0, self.domType)
    if self.domType != 'TP':
      self.solver.SetZ(self.zpts)
      # plan fftws
      self.transformer = Transformer(self.Nx, self.Ny, self.Nz, self.dof, 1)
    else:
      self.transformer = Transformer(self.Nx, self.Ny, self.Nz, self.dof, 0)
    print('Solver is ready\n')

  def SetPositions(self, xP):
    """
    Set the positions of the particles 
    On the first call, this builds grid-particle search structures
    and allocates memory for the grids
    On subsequent calls, the search structures are re-used to 
    map the new positions to the grid.
   
    If a position is outside of the z-extent (for DP problems),
    the code will exit with an error.
 
    Inputs:
      xP - positions (x1,y1,z1,x2,y2,z2,..)
    """
    self.nP = int(xP.shape[0] / 3)
    # build search structures on first call and check validity of configuration
    if not self.posinit:
      print('Setting initial positions...')
      for j in range(0, 3 * self.nP, 3):
        if xP[j] < self.minX or xP[j] > self.maxX or\
           xP[j+1] < self.minY or xP[j+1] > self.maxY or\
           xP[j+2] < self.minZ or xP[j+2] > self.maxZ:
          print('Position is outside of the unit cell for particle # ', j / 3)
          raise ValueError('Initial positions must be within unit cell')
      self.xP = xP
      print('Building monopole particle-grid search structures...')
      self.particles = ParticlesGen(self.nP, self.dof, self.xP, np.zeros((3 * self.nP,)),\
                                    self.radP, self.wm, self.cbetam, self.betamP)
      self.particles.UseCbeta()
      self.particles.Make()
      self.particles.Setup(self.grid)
      self.vP = np.zeros((self.dof * self.nP,), dtype = np.double); self.omegaP = 0
      if self.has_torque:
        print('Building dipole particle-grid search structures...')
        self.tparticles = ParticlesGen(self.nP, self.dof, self.xP, np.zeros((3 * self.nP,)),\
                                       self.radP, self.wd, self.cbetad, self.betadP)
        self.tparticles.UseCbeta()
        self.tparticles.IsDipole() # need to specify these are dipoles
        self.tparticles.Make()
        self.tparticles.Setup(self.tgrid)
        self.omegaP = np.zeros((self.nP * self.dof), dtype = np.double)
        self.dxvP = np.zeros((self.nP * self.dof), dtype = np.double)
        self.dyvP = np.zeros((self.nP * self.dof), dtype = np.double)
        self.dzvP = np.zeros((self.nP * self.dof), dtype = np.double)
      self.posinit = True 
    else: 
      if xP.shape[0] != 3 * self.nP:
        raise ValueError('Positions array has incorrect length')
    self.xP = xP
    # use existing search structures to update positions
    self.particles.Update(self.xP, self.grid)
    if self.has_torque:
      self.tparticles.Update(self.xP, self.tgrid)
  
  def Mdot(self, forces, torques=None):
    """
    Compute V = M*F

    Inputs:
      forces/torques (fx1,fy1,fz1,fx2,fy2,fz2,...), (tx1,ty1,tz1,tx2,ty2,tz2,..)

    Output:
      vP, omegaP    - lin/rot vel on the particles
                    - (vx1,vy1,vz1,...vxn,vyn,vzn) if has_torque=False
                    - (vx1,vy1,vz1,...vxn,vyn,vzn,wx1,wy1,wz1,...,wxn,wyn,wzn) if has_torque=True
    """
    nf, = forces.shape; nt, = torques.shape
    if (self.has_torque and nf + nt != 6 * self.nP) or\
       (not self.has_torque and nf != 3 * self.nP):
      raise ValueError('Dimension mismatch in particle force array')


    self.particles.SetData(forces)    

    # spread forces
    self.F = Spread(self.particles, self.grid); 
    if self.has_torque and torques is not None:
     
      self.tparticles.SetData(torques)    

      # spread derivative of torques, compute curl and add to F 
      dxT = SpreadDx(self.tparticles, self.tgrid)
      libSpreadInterp.addDx(self.grid.grid, self.tgrid.grid)
      dyT = SpreadDy(self.tparticles, self.tgrid)
      libSpreadInterp.addDy(self.grid.grid, self.tgrid.grid)
      # set BC for z deriv
      self.tgrid.SetBCs(self.BCst_z)
      dzT = SpreadDz(self.tparticles, self.tgrid)
      libSpreadInterp.addDz(self.grid.grid, self.tgrid.grid)
 
    # forward transform
    self.transformer.Ftransform(self.F)
    self.F_hat_r = self.transformer.out_real; 
    self.F_hat_i = self.transformer.out_imag; 

    # solve pde 
    self.U_hat_r, self.U_hat_i, self.P_hat_r, self.P_hat_i = self.solver.Solve(self.F_hat_r, self.F_hat_i)
  
    # back transform
    self.transformer.Btransform(self.U_hat_r, self.U_hat_i)
    self.uG_r = self.transformer.out_real

    # interpolate linear velocities on particles
    self.grid.SetData(self.uG_r)
    Interpolate(self.particles, self.grid, self.vP)

    if self.has_torque and torques is not None:

      # interpolate derivative of linear velocities on particles
      self.tgrid.SetData(self.uG_r)
      # set BC for x-y deriv
      self.tgrid.SetBCs(self.BCs)
      InterpolateDx(self.tparticles, self.tgrid, self.dxvP)
      InterpolateDy(self.tparticles, self.tgrid, self.dyvP)
      # set BC for z deriv
      self.tgrid.SetBCs(self.BCst_z)
      InterpolateDz(self.tparticles, self.tgrid, self.dzvP)
      # reset boundary condition, since we changed it for z
      self.tgrid.SetBCs(self.BCs)
      # compute angular velocity of particles
      self.omegaP[0::3] = -1/2 * (self.dyvP[2::3] - self.dzvP[1::3])
      self.omegaP[1::3] = -1/2 * (self.dzvP[0::3] - self.dxvP[2::3]) 
      self.omegaP[2::3] = -1/2 * (self.dxvP[1::3] - self.dyvP[0::3]) 

    if self.has_torque:
      return self.vP, self.omegaP
    else:
      return self.vP

  def Clean(self):
    """
    Deallocate memory - must be called when finished with FCM module
    """
    print('Deallocating memory')
    self.grid.Clean()
    self.particles.Clean()
    self.transformer.Clean()
    self.solver.Clean()
    if self.has_torque:
        self.tgrid.Clean()
        self.tparticles.Clean()
    self.posinit = False 
    self.particles = None
    self.tparticles = None
    self.grid = None
    self.tgrid = None
    self.transformer = None
    self.solver = None 


