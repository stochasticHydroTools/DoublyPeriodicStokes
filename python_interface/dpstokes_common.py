#Raul P. Pelaez 2021. Common interface for the FCM (Triply periodic + DPStokes) CPU (Sachin) and GPU (Raul) solvers.
#Examples on the interfaces for the CPU and GPU versions are available in dpstokesCPU.py and dpstokesGPU.py respectively.
#import common_interface_wrapper and run help(FCMJoint) for detailed usage information.
import numpy as np
from common_interface_wrapper import FCMJoint
np.random.seed(1234)

# Device to run the code, can be either 'cpu' or 'gpu'
device = 'cpu'

# first set the domain type, extents, viscosity and number of paticles 
# domType can be 'DPBW', 'DPSC', 'DP', 'TP'
# for bottom wall, slit channel, no wall, and triply periodic
domType = 'DPBW'
nP = 2048 #Number of particles
has_torque = True #Set to True if angular displacements are needed 
viscosity = 0.957e-3
#Simulation domain limits
xmin = 0.0; xmax = 128.7923
ymin = 0.0; ymax = 128.7923
zmin = 0.0; zmax = 9.1740639106166668
# Initialize the solver with all the parameters
solver = FCMJoint(device)
#Initialize can be called several times in order to change the parameters
solver.Initialize(numberParticles=nP, hydrodynamicRadius=1.0155, kernType=0,
                  domType=domType, has_torque=has_torque,
                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  viscosity=viscosity, optInd=0)

# now define some random positions
xP = np.zeros(3 * nP, dtype = np.double) # Arbitrary xP inside the Z domain
xP[0::3] = np.random.uniform(xmin, xmax, nP) #X
xP[1::3] = np.random.uniform(ymin, ymax, nP) #Y
xP[2::3] = np.random.uniform(zmin, zmax, nP) #Z
# pass the positions to the problem
solver.SetPositions(xP)

# set some random torques and forces
torques = np.random.normal(0, 1, 3 * nP) 
forces = np.random.normal(0, 1, 3 * nP)

# solve the mobility problem
V, T = solver.Mdot(forces, torques)

# clear memory
solver.Clean()

#Print something
print("Positions of the first three particles: ")
xP = np.reshape(xP, (nP, 3));
print(xP[0:3])

MF = np.reshape(V, (nP, 3));
print("Linear velocity of the first three particles: ")
print(MF[0:3])

MT = np.reshape(T, (nP, 3));
print("Angular velocity of the first 3 particles: ")
print(MT[0:3])

