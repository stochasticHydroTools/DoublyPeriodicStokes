import numpy as np
import uammd # GPU solver
import os
import sys
if "FCM_CPUCONFIG_LAUNCHED" not in os.environ:
    print("ERROR: Please source cpuconfig.sh before using this module")
    sys.exit()
from FCM import *  # CPU solver

# Donev: Explain to me how changing the domain size only (for open domains) is supposed to be done here -- calling Initialize again?
#A, Raul: Correct, I would need to zoom discuss with you both to separate it. I am not sure about what exactly Sachin's setUnitCell does (a.i. does it adapt the input size?)

# I don't get where the CPU method clean is called etc. but I am sure I am just not getting exactly what some lines do...
# A, Raul: It is called in Clean(self):. If you are asking "why" is it called in the constructor, it is because AFAIK, Sachin's impl requires Clean before reinitialization.
class FCMJoint:

    def __init__(self, device = 'cpu'):
        self.__device = device
        self.__gpuCreated = False
        self.__cpuCreated = False
        if device == 'cpu':
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
            Lx, Ly, _, _, nx, ny, w, w_d, cbeta, cbeta_d, beta, beta_d\
                = configure_grid_and_kernels_xy(xmax-xmin, ymax-ymin, radP, kernTypes, optInd, has_torque)
            zmax, _, nz, zmin = configure_grid_and_kernels_z(zmin, zmax, hx, w, w_d, domType, has_torque)
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
        elif self.__device == 'gpu' and self.__gpuCreated:
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


