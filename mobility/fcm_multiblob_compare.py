from common_interface_wrapper import FCMJoint
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
# compute the mobility of a 2048 particle configuration 
# with random forces and torques in DPBW

device = 'cpu'
domType = 'DPBW'
has_torque = True
eta = 0.957e-3
nBlobs = 12 #nBlobs = 1 #nBlobs = 42
nP = 2048 * nBlobs
if nBlobs == 1:
  radP = 1.0155 
  xP = np.loadtxt('./Test_Data_For_Rollers/Const_Torque_t_15.clones', skiprows=1, usecols=[0,1,2])
  data = np.loadtxt('./Test_Data_For_Rollers/One_Blob/N_Images_64.txt')
  F = data[:,0].copy()
  if has_torque:
    T = data[:,1].copy()
elif nBlobs == 12:
  radP = 0.42287520345118434 
  xP = np.loadtxt('./Test_Data_For_Rollers/MultiBlob/Cfg_12_blobs_per.txt', skiprows=1, usecols=[0,1,2])
  F = np.random.normal(0.0, 1.0, 3*nP)
  if has_torque:
    T = np.random.normal(0.0, 1.0, 3*nP)
elif nBlobs == 42:
  radP = 0.247328128441116 
  xP = np.loadtxt('./Test_Data_For_Rollers/MultiBlob/Cfg_42_blobs_per.txt', skiprows=1, usecols=[0,1,2])
  F = np.random.normal(0.0, 1.0, 3*nP)
  if has_torque:
    T = np.random.normal(0.0, 1.0, 3*nP)
if not has_torque:
  T = None

xmin = 0.0; xmax = 128.7923
ymin = 0.0; ymax = 128.7923
zmin = 0.0; zmax = np.max(xP[:,2])
xP = np.reshape(xP, (nP * 3,)).copy()

kernType = 0 * has_torque + 2 * (not has_torque)

solver = FCMJoint(device)
solver.Initialize(numberParticles=nP, hydrodynamicRadius=radP, kernType=kernType,
                  domType=domType, has_torque=has_torque,
                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  viscosity=eta, optInd=0)

solver.SetPositions(xP)
mob, mobt = solver.Mdot(F, T)

print(mob)

if nBlobs == 1:
  Mtt_x_f = data[:,2]
  Mtr_x_T = data[:,3]
  Mrt_x_F = data[:,4]
  Mrr_x_T = data[:,5]
  
  if has_torque:
    mob_multiblob = Mtt_x_f + Mtr_x_T
    mobt_multiblob = Mrt_x_F + Mrr_x_T
    fig, ax = plt.subplots(1,2)
    ax[0].scatter(mob_multiblob, mob,c='r',label='Linear Velocity')
    ax[0].set_xlabel('RigidMultiBlob', fontsize=15)
    ax[0].set_ylabel('FCM', fontsize=15)
    ax[0].legend(fontsize=15)
    ax[1].scatter(mobt_multiblob, mobt,c='b', label='Rotational Velocity')
    ax[1].set_xlabel('RigidMultiBlob', fontsize=15)
    ax[1].set_ylabel('FCM',fontsize=15)
    ax[1].legend(fontsize=15)
    plt.show()
  else:
    mob_multiblob = Mtt_x_f
    fig, ax = plt.subplots(1,1)
    ax.scatter(mob_multiblob, mob,c='r',label='Linear Velocity')
    ax.set_xlabel('RigidMultiBlob',fontsize=15)
    ax.set_ylabel('FCM',fontsize=15)
    ax.legend(fontsize=15)
    plt.show()
