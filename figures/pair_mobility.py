import numpy as np
from mobility_matrix import *
import scipy.io

device = 'cpu'
domType = 'DPBW'
kernType = 0
ref2x = False
ref_rpy = False
useimg = False


if kernType == 0:
  has_torque = True
elif kernType == 2:
  has_torque = False


eta = 1/4/np.sqrt(np.pi)
xmin = 0.0; xmax = 76.8
ymin = 0.0; ymax = 76.8
zmin = 0.0; zmax = 19.2
radP = 1.0

if useimg:
  num_imgs = 200
else:
  num_imgs = 0

nHeights = 26
heights1 = np.linspace(0,1,10)
heights2 = np.linspace(1,zmax/4,9)
heights3 = np.linspace(zmax/4,zmax/2,9)
heights = np.union1d(heights1,np.union1d(heights2,heights3)).reshape((nHeights,1)) 


solver = FCMJoint(device); nP = 2
if (domType == 'DPBW' and not ref_rpy) or domType == 'DPSC':
  solver.Initialize(numberParticles=nP, hydrodynamicRadius=radP, kernType=0,
                    domType=domType, has_torque=has_torque,
                    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                    viscosity=eta, optInd=0, ref=ref2x, useRegKernel=False)

  seps = np.array([3 * radP, 4 * radP, 8 * radP])
  nSep = seps.shape[0]
  xPs = np.zeros((nHeights, 3 * nP), dtype = np.double)
  if kernType == 0:
    M = np.zeros((nSep, nHeights, 6 * nP, 6 * nP), dtype = np.double)
  elif kernType == 2:
    M = np.zeros((nSep, nHeights, 3 * nP, 3 * nP), dtype = np.double)
     
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
elif domType == 'DPBW' and ref_rpy:
  L = np.array([xmax, ymax, 0])
  heights = np.linspace(1,zmax/2,nHeights)
  seps = np.array([3 * radP, 4 * radP, 8 * radP])
  nSep = seps.shape[0]
  xPs = np.zeros((nHeights, 3 * nP), dtype = np.double)
  M = np.zeros((nSep, nHeights, 6 * nP, 6 * nP), dtype = np.double)
  for iSep in range(0, nSep):
    for iH in range(0, nHeights):
      xPs[iH,0:3] = [0, 0, heights[iH]]
      xPs[iH,3::] = [seps[iSep], 0, heights[iH]]
    M[iSep,:,:,:], NormMat = mobility_matrix_numba(xPs, eta, L, radP, num_imgs, self=False)


if domType == 'DPBW':
  if ref_rpy:
    if useimg:
      outname = 'pair_mobility_bw_ref.mat'
    else:
      outname = 'pair_mobility_bw_ref_noimg.mat'
  elif ref2x:
    outname = 'pair_mobility_bw_ref2x.mat'
  elif kernType == 0:
    outname = 'pair_mobility_bw.mat'
  elif kernType == 2:
    outname = 'pair_mobility_bw_w4.mat'
elif domType == 'DPSC':
  if ref2x:
    outname = 'pair_mobility_sc_ref.mat'
  elif kernType == 0:
    outname = 'pair_mobility_sc.mat'
  elif kernType == 2:
    outname = 'pair_mobility_sc_w4.mat'
scipy.io.savemat(outname, dict(M=M, heights=heights/radP))

