import numpy as np
from mobility_matrix import *
import scipy.io

device = 'cpu'
domType = 'DPBW'
kernType = 0

eta = 1/4/np.sqrt(np.pi)
xmin = 0.0; xmax = 76.8
ymin = 0.0; ymax = 76.8
zmin = 0.0; zmax = 19.2
radP = 1.0
nP = 2
nTrials = 100

if kernType == 0:
  has_torque = True
elif kernType == 2:
  has_torque = False

solver = FCMJoint(device); nP = 2
solver.Initialize(numberParticles=nP, hydrodynamicRadius=radP, kernType=0,
                  domType=domType, has_torque=has_torque,
                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  viscosity=eta, optInd=0, ref=False, useRegKernel=False)

xPs = np.zeros((nTrials, 3 * nP), dtype = np.double)

for iTrial in range(1, nTrials):
  xPs[iTrial, 0] = np.random.random() * xmax
  xPs[iTrial, 1] = np.random.random() * ymax
  xPs[iTrial, 2] = np.random.random() * zmax
  xPs[iTrial, 3] = np.random.random() * xmax
  xPs[iTrial, 4] = np.random.random() * ymax
  xPs[iTrial, 5] = np.random.random() * zmax

Ms, NormMat = mobility_matrix(xPs, solver)
solver.Clean()

Masym = np.zeros(nTrials)
Posdef = np.zeros(nTrials)
for iTrial in range(1, nTrials):
  M = Ms[iTrial,:,:]
  Masym[iTrial] = np.linalg.norm(M - M.T) / np.linalg.norm(M)
  w, v = np.linalg.eig(M); w = np.sort(w)
  Posdef[iTrial] = w[0]

if domType == 'DPBW':
  if kernType == 0:
    outname = 'mobility_asymm_posdef_bw_w6.mat'
  elif kernType == 2:
    outname = 'mobility_asymm_posdef_bw_w4.mat'
elif domType == 'DPSC':
  if kernType == 0:
    outname = 'mobility_asymm_posdef_sc_w6.mat'
  elif kernType == 2:
    outname = 'mobility_asymm_posdef_sc_w4.mat'
scipy.io.savemat(outname, dict(Masym=Masym, Posdef=Posdef))

