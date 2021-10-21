#Sachin Natesh and Raul P. Pelaez 2021.
#The class FCMJoint in this source code joins both the GPU and CPU python interfaces into a common one. An usage example is available in dpstokes_common.py.
#See dpstokesCPU.py or dpstokesGPU.py for usage examples for each individual interface.
import numpy as np
import os
import sys
if "FCM_CPUCONFIG_LAUNCHED" not in os.environ:
    print("ERROR: Please source cpuconfig.sh before using this module")
    sys.exit()
try:
  import uammd # GPU solver
except:
  print('Could not load uammd GPU Python interface: Have you run make in the root directory?')

try:
  from FCM import *
except:
  print('Could not load CPU Python interface: Have you run make in the root directory?')
  from GridAndKernelConfig import *

class FCMJoint:
    """A common interface for the CPU and GPU triply and doubly periodic FCM solvers."""
    def __init__(self, device = 'cpu'):
        self.device = device
        self.gpuCreated = False
        self.cpuCreated = False
        if device == 'cpu':
          self.precision = np.float64
        else:
          self.precision = np.float32 if uammd.getPrecision() == 'single' else np.float64

    def Initialize(self, numberParticles, hydrodynamicRadius, viscosity,
                   domType, has_torque,
                   xmax, ymax, zmin, zmax,
                   xmin=0, ymin=0, kernType = 0, optInd=0, ref=False):
      """
      Initialize the DPStokes module, can be called on an already initialize module to change the parameters.
        Inputs:
          numberParticles
          hydrodynamicRadius
          xmin, xmax, ymin, ymax, zmin, zmax - The simulation box limits
          viscosity - for the Stokes problem
          optInd - which of the candidate grids to select (see stdout). Default is 0.
          radP - particle radius 
          kernTypes - kernel type (default=0):
                    - 0 for ES 6 pt (both monopole and dipole if has_torque = true) 
                    - 1 for ES 5 pt (both monopole and dipole if has_torque = true)
                    - 2 for ES 4 pt (only monopole, dipole not supported) 
                    - 3 for ES 6 pt monopole, ES 5 pt dipole
                    - 4 for ES 6 pt monopole, ES 4 pt dipole
                    - 5 for ES 5 pt monopole, ES 4 pt dipole 
          domType - domain type:
                  - 'TP' for triply periodic
                  - 'DP' for doubly periodic
                  - 'DPBW' for doubly periodic with bottom wall
                  - 'DPSC' for doubly periodic slit channel
          has_torque - bool specifying if problem has torque
          ref - doubly resolved grid (2x pts in each direction, for reference computation. default=False)

        Stdout: 
          - The optimal *adjusted* grid is displayed (optInd = -1)
          - Several candidate fft-friendly grids are displayed (optInd = 0,1,..)
          - The final kernel settings for each particle size are displayed
      """
      self.Clean()
      self.has_torque = has_torque
      self.viscosity = viscosity
      self.radP = hydrodynamicRadius*np.ones(numberParticles)
      self.kernTypes = kernType*np.ones(numberParticles, dtype=int)
      if self.device == 'cpu':
          self.cpusolver = FCM(self.radP, self.kernTypes, domType, has_torque)
          self.cpusolver.SetUnitCell([xmin,xmax], [ymin,ymax], [zmin,zmax])
          self.cpusolver.Initialize(viscosity, optInd=optInd, ref=ref)
          self.cpuCreated = True
      elif self.device == 'gpu':
          Lx, Ly, hx, hy, nx, ny, w, w_d, cbeta, cbeta_d, beta, beta_d\
              = configure_grid_and_kernels_xy(xmax-xmin, ymax-ymin, self.radP, self.kernTypes, optInd, has_torque, ref)
          zmax, hz, nz, zmin = configure_grid_and_kernels_z(zmin, zmax, hx, w, w_d, domType, fac=1.5, ref=ref)
          if domType == 'TP':
              mode = 'periodic'
          elif domType == 'DP':
              mode = 'nowall'
          elif domType == 'DPBW':
              mode = 'bottom'
          elif domType == 'DPSC':
              mode = 'slit'
          self.par = uammd.StokesParameters(viscosity=viscosity,
                                       Lx=Lx, Ly=Ly,
                                       zmin=zmin, zmax=zmax,
                                       w=w[0], w_d=w_d[0],
                                       beta=beta[0]*w[0], beta_d=beta_d[0]*w_d[0],
                                       nx=nx, ny=ny, nz=nz, mode=mode)
          print(self.par)
          if not self.gpuCreated:
              self.gpusolver = uammd.DPStokes()
              self.gpuCreated = True
          self.gpusolver.initialize(self.par, numberParticles)

    def Clean(self):
        """Release all memory allocated by the module"""
        if self.device == 'cpu' and self.cpuCreated:
            self.cpusolver.Clean()
            self.cpuCreated = False
        elif self.device == 'gpu' and self.gpuCreated:
            self.gpuCreated = False
            self.gpusolver.clear()

    def SetPositions(self, positions):
        """Set the positions to compute the mobility matrix"""
        if self.device == 'cpu':
            self.cpusolver.SetPositions(positions)
        elif self.device == 'gpu':
            self.gpusolver.setPositions(self.precision(positions))

    def Mdot(self, forces, torques=np.array(0)):
        """Computes the product of the Mobility tensor with the provided forces and torques. 
           If torques are not present, they are assumed to be zero and angular displacements will not be computed
        """
        if self.device == 'cpu':
            MF, MT  = self.cpusolver.Mdot(forces, torques)
        elif self.device == 'gpu':
            MF = np.zeros(forces.size, self.precision)
            MT = np.zeros(torques.size, self.precision)
            self.gpusolver.Mdot(forces=self.precision(forces), torques=self.precision(torques),
                                velocities=MF, angularVelocities=MT)
        
        return MF, MT
