import sys
import random
from FCM import *
np.random.seed(1234)

# example usage of the FCM CPU module


# first set the domain type, extents, viscosity and number of paticles 

# domType can be 'DPBW', 'DPSC', 'DP', 'TP'
# for bottom wall, slit channel, no wall, and triply periodic
domType = 'DPBW'
nP = 2048; 
has_torque = True
viscosity = 0.957e-3
xmin = 0.0; xmax = 128.7923
ymin = 0.0; ymax = 128.7923
zmin = 0.0; zmax = 9.1740639106166668
# set the radii and kernels for each particle
radP = 1.0155 * np.ones(nP, dtype = np.double)
kernTypes = np.zeros(nP, dtype = np.int)

# create the stokes mobility problem
problem = FCM(radP, kernTypes, domType, has_torque)
problem.SetUnitCell([xmin,xmax], [ymin,ymax], [zmin,zmax])
problem.Initialize(viscosity, optInd=0)


# now define some random positions
xP = np.zeros(3 * nP, dtype = np.double) # Arbitrary xP inside the Z domain
xP[0::3] = np.random.uniform(xmin, xmax, nP) #X
xP[1::3] = np.random.uniform(ymin, ymax, nP) #Y
xP[2::3] = np.random.uniform(zmin, zmax, nP) #Z
# pass the positions to the problem
problem.SetPositions(xP)

# set some random torques and forces
torques = np.random.normal(0, 1, 3 * nP) 
forces = np.random.normal(0, 1, 3 * nP)

if has_torque:
  F = np.concatenate((forces,torques))
else:
  F = forces

# solve the mobility problem
V = problem.Mdot(F)

# clear memory
problem.Clean()

#Print something
print("Positions of the first three particles: ")
xP = np.reshape(xP, (nP, 3));
print(xP[0:3])

MF = np.reshape(V[0:3*nP], (nP, 3));
print("Linear velocity of the first three particles: ")
print(MF[0:3])

MT = np.reshape(V[3*nP::], (nP, 3));
print("Angular velocity of the first 3 particles: ")
print(MT[0:3])
