from mobility_bench import *
import scipy.io
import matplotlib.pyplot as plt
np.random.seed(0)
# compute the mobility of a 2048 particle configuration 
# with random forces and torques in DPBW

########################## SOLVER INFO ####################################3
# solver/domain Type ('TP', 'DP', 'DPBW', 'DPSC')
#solverType = 'DPSC'
solverType = 'DPBW'
# viscosity
eta = 0.957e-3# solver degrees of freedom (should be 3)
# solver degrees of freedom (should be 3)
dof = 3; 
# switch for torque
has_torque = True

########################## GRID INFO ####################################3
# grid extent in x,y
Lx = 128.7923; Ly = 128.7923

########################## PARTICLE INFO ####################################3
# number of particles
nBlobs = 1
#nBlobs = 12#nBlobs = 42
#nBlobs = 42
nP = 2048 * nBlobs
# particle configuration
xP = np.loadtxt('./Test_Data_For_Rollers/Const_Torque_t_15.clones', skiprows=1, usecols=[0,1,2])
#xP = np.loadtxt('./Test_Data_For_Rollers/MultiBlob/Cfg_12_blobs_per.txt', skiprows=1, usecols=[0,1,2])
#xP = np.loadtxt('./Test_Data_For_Rollers/MultiBlob/Cfg_42_blobs_per.txt', skiprows=1, usecols=[0,1,2])
if nBlobs == 1:
  data = np.loadtxt('./Test_Data_For_Rollers/One_Blob/N_Images_64.txt')
# min/max particle height for determining z extent based on solverType
minPz = 0; maxPz = np.max(xP[:,2])
#minPz = 0; maxPz = 9.1740639106166668#np.max(xP[:,2])
# forces
F = data[:,0].copy()
#F = np.random.normal(0.0, 1.0, 3*nP)
# torques
if has_torque:
  T = data[:,1].copy()
  #T = np.random.normal(0.0, 1.0, 3*nP)
else:
  T = None
xP = np.reshape(xP, (nP * 3,)).copy()
# radii   
radP = 1.0155 * np.ones(nP, dtype = np.double)
#radP = 0.42287520345118434 * np.ones(nP, dtype = np.double)
#radP = 0.247328128441116 * np.ones(nP, dtype = np.double)# kernel types (see config_grid_and_kernels())
# kernel types (see config_grid_and_kernels())
if has_torque:
  kernTypes = np.zeros(nP, dtype = np.int) # use w=6 for monopole (and dipole if has_torque=True)
else:
  kernTypes = 2 * np.ones(nP, dtype = np.int) # use w=4 for monopole 

############################ SOLVE #####################################
mob, mobt = mobility_derivkern(xP, F, T, eta, Lx, Ly, minPz, maxPz, radP, kernTypes, solverType, has_torque, write=False, ref=False, bench=True)


print(mob)
######################### SCATTER PLOT #######################################

#Mtt_x_f = data[:,2]
#Mtr_x_T = data[:,3]
#Mrt_x_F = data[:,4]
#Mrr_x_T = data[:,5]
#
#if has_torque:
#  mob_multiblob = Mtt_x_f + Mtr_x_T
#  mobt_multiblob = Mrt_x_F + Mrr_x_T
#  fig, ax = plt.subplots(1,2)
#  ax[0].scatter(mob_multiblob, mob,c='r',label='Linear Velocity')
#  ax[0].set_xlabel('RigidMultiBlob', fontsize=15)
#  ax[0].set_ylabel('FCM', fontsize=15)
#  ax[0].legend(fontsize=15)
#  ax[1].scatter(mobt_multiblob, mobt,c='b', label='Rotational Velocity')
#  ax[1].set_xlabel('RigidMultiBlob', fontsize=15)
#  ax[1].set_ylabel('FCM',fontsize=15)
#  ax[1].legend(fontsize=15)
#  plt.show()
#else:
#  mob_multiblob = Mtt_x_f
#  fig, ax = plt.subplots(1,1)
#  ax.scatter(mob_multiblob, mob,c='r',label='Linear Velocity')
#  ax.set_xlabel('RigidMultiBlob',fontsize=15)
#  ax.set_ylabel('FCM',fontsize=15)
#  ax.legend(fontsize=15)
#  plt.show()
