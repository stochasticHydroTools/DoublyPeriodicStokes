import numpy as np
from mobility_matrix import *
import scipy.io

device = 'cpu'
domType = 'DPSC'
nP = 2
has_torque = True
eta = 1/4/np.sqrt(np.pi)
xmin = 0.0; xmax = 76.8
ymin = 0.0; ymax = 76.8
zmin = 0.0; zmax = 19.2
radP = 1.0
ref = True

solver = FCMJoint(device)
solver.Initialize(numberParticles=nP, hydrodynamicRadius=radP, kernType=0,
                  domType=domType, has_torque=has_torque,
                  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                  viscosity=eta, optInd=0, ref=ref)

nHeights = 74
heights1 = np.linspace(0,1,20)
heights2 = np.linspace(1,1.3,19)
heights3 = np.linspace(1.3,zmax/4,19)
heights4 = np.linspace(zmax/4,zmax/2,19)
heights = np.union1d(heights1,np.union1d(heights2,np.union1d(heights3,heights4))).reshape((nHeights,1))

seps = np.array([3 * radP, 4 * radP, 8 * radP])
nSep = seps.shape[0]
xPs = np.zeros((nHeights, 3 * nP), dtype = np.double)

M = np.zeros((nSep, nHeights, 6 * nP, 6 * nP), dtype = np.double)
#M = np.zeros((nSep, nHeights, 3 * nP, 3 * nP), dtype = np.double)
xPs = np.zeros((nHeights, 3 * nP), dtype = np.double)
for iSep in range(0, nSep):
  for iH in range(1, nHeights):
    xPs[iH, 0] = xmax / 2 + seps[iSep] / 2
    xPs[iH, 1] = xmax / 2
    xPs[iH, 3] = xmax / 2 - seps[iSep] / 2
    xPs[iH, 4] = xmax / 2
    xPs[iH, 2] = xPs[iH, 5] = heights[iH]

  M[iSep,:,:,:], NormMat = mobility_matrix(xPs, solver)

solver.Clean()

if domType == 'DPBW':
  if ref:
    outname = 'pair_mobility_bw_ref.mat'
  else:
    outname = 'pair_mobility_bw.mat'

elif domType == 'DPSC':
  if ref:
    outname = 'pair_mobility_sc_ref.mat'
  else:
    outname = 'pair_mobility_sc.mat'
scipy.io.savemat(outname, dict(M=M, heights=heights/radP))

