import sys
import random
from FCM import *
np.random.seed(0)

# example usage of the FCM module

pos  = np.loadtxt('./Test_Data_For_Rollers/Const_Torque_t_15.clones', skiprows=1, usecols=[0,1,2])
data = np.loadtxt('./Test_Data_For_Rollers/One_Blob/N_Images_64.txt')

nP = 2048; 
domType = 'DPBW'
eta = 0.957e-3
has_torque = False
minX = 0.0; maxX = 128.7923
minY = 0.0; maxY = 128.7923
minZ = 0.0; maxZ = 20.0
xP = np.reshape(pos, (3 * nP,)).copy()
if has_torque:
  F = np.concatenate((data[:,0],data[:,1]))
else:
  F = data[:,0].copy()
radP = 1.0155 * np.ones(nP, dtype = np.double)
kernTypes = np.zeros(nP, dtype = np.int)

problem = FCM(radP, kernTypes, domType, has_torque)
problem.SetUnitCell([minX,maxX], [minY,maxY], [minZ,maxZ])
problem.Initialize(eta, 0)
problem.SetPositions(xP)

V = problem.Mdot(F)
print(V[0:3*nP])
# dummy update
xP += 0.1 * np.abs(V[0:3 * nP])
problem.SetPositions(xP)

V = problem.Mdot(F)

problem.Clean()
