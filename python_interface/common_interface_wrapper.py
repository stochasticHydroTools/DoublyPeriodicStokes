import numpy as np
import os
import sys
if "FCM_CPUCONFIG_LAUNCHED" not in os.environ:
    print("ERROR: Please source cpuconfig.sh before using this module")
    sys.exit()
try:
  import uammd # GPU solver
except:
  print('Could not load uammd GPU Python interface')

try:
  from FCM import *
except:
  print('Could not load CPU Python interface')
  from GridAndKernelConfig import *

class FCMJoint:

    def __init__(self, device = 'cpu'):
        self.__device = device
        self.__gpuCreated = False
        self.__cpuCreated = False
        if device == 'cpu':
          # CPU solver
          self.precision = np.float64
        else:
          self.precision = np.float32

    def Initialize(self, numberParticles, hydrodynamicRadius, viscosity,
                   kernType, domType,
                   has_torque,
                   xmax, ymax, zmin, zmax,
                   xmin=0, ymin=0, optInd=0):        
        self.Clean()
        self.__has_torque = has_torque
        radP = hydrodynamicRadius*np.ones(numberParticles)
        kernTypes = kernType*np.ones(numberParticles, dtype=int)
        if self.__device == 'cpu':
            self.cpusolver = FCM(radP, kernTypes, domType, has_torque)
            self.cpusolver.SetUnitCell([xmin,xmax], [ymin,ymax], [zmin,zmax])
            self.cpusolver.Initialize(viscosity, optInd=optInd)
            self.__cpuCreated = True
        elif self.__device == 'gpu':
            Lx, Ly, hx, hy, nx, ny, w, w_d, cbeta, cbeta_d, beta, beta_d\
                = configure_grid_and_kernels_xy(xmax-xmin, ymax-ymin, radP, kernTypes, optInd, has_torque)
            zmax, hz, nz, zmin = configure_grid_and_kernels_z(zmin, zmax, hx, w, w_d, domType, has_torque)
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
            if not self.__gpuCreated:
                self.gpusolver = uammd.DPStokes()
                self.__gpuCreated = True
            self.gpusolver.initialize(self.par, numberParticles)

    def Clean(self):
        if self.__device == 'cpu' and self.__cpuCreated:
            self.cpusolver.Clean()
            self.__cpuCreated = False
        elif self.__device == 'gpu' and self.__gpuCreated:
            self.__gpuCreated = False
            self.gpusolver.clear()

    def SetPositions(self, positions):
        if self.__device == 'cpu':
            self.cpusolver.SetPositions(positions)
        elif self.__device == 'gpu':
            self.gpusolver.setPositions(positions)

    def Mdot(self, forces, torques=np.array(0)):
        if self.__device == 'cpu':
            if self.__has_torque:
                F = np.concatenate((forces, torques))
            else:
                F = forces
            V = self.cpusolver.Mdot(F)
            MF = V[0:forces.size]
            MT = V[forces.size::]
            return MF, MT
        elif self.__device == 'gpu':
            __MT = np.zeros(torques.size, self.precision)
            __MF = np.zeros(forces.size, self.precision)
            self.gpusolver.Mdot(forces=forces, torques=torques,
                                velocities=__MF, angularVelocities=__MT)
            return __MF, __MT


