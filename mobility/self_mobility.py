import numpy as np
from mobility_matrix import *
import scipy.io

device = 'cpu'
domType = 'DPBW'
kernType = 0
ref2x = True
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

nHeights = 60
heights = np.linspace(0, zmax/4, 50)
heights = np.concatenate((heights, np.linspace(zmax/4+heights[2]-heights[1],zmax/2,10)))

solver = FCMJoint(device)
if (domType == 'DPBW' and not ref_rpy) or domType == 'DPSC':
  xPs = np.zeros((nHeights, 3), dtype = np.double)
  for iH in range(0, nHeights):
    xPs[iH, 0] = xPs[iH, 1] = xmax / 2
    xPs[iH, 2] = heights[iH]
  
  solver.Initialize(numberParticles=1, hydrodynamicRadius=radP, kernType=kernType,
                    domType=domType, has_torque=has_torque,
                    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax,
                    viscosity=eta, optInd=0, ref=ref2x, useRegKernel=False)
  M, NormMat = mobility_matrix(xPs, solver)
  solver.Clean()

elif domType == 'DPBW' and ref_rpy:
  L = np.array([xmax, ymax, 0])
  nHeights = 70
  xPs = np.zeros((nHeights, 6), dtype = np.double)
  heights = np.linspace(1,zmax/4,60);
  heights = np.concatenate((heights, np.linspace(zmax/4+heights[2]-heights[1],zmax/2,10)))

  for iH in range(0, nHeights):
    xPs[iH,0:3] = [0, 0, heights[iH]]
    xPs[iH,3::] = [1e-20, 0, heights[iH]]
  
  M, NormMat = mobility_matrix_numba(xPs, eta, L, radP, num_imgs=num_imgs, self=True)


if domType == 'DPBW':
  if ref_rpy:
    if useimg:
      outname = 'self_mobility_bw_ref.mat'
    else:
      outname = 'self_mobility_bw_ref_noimg.mat'
  elif ref2x:
    outname = 'self_mobility_bw_ref2x.mat'
  elif kernType == 0:
    outname = 'self_mobility_bw.mat'
  elif kernType == 2:
    outname = 'self_mobility_bw_w4.mat'
elif domType == 'DPSC':
  if ref2x:
    outname = 'self_mobility_sc_ref.mat'
  elif kernType == 0:
    outname = 'self_mobility_sc.mat'
  elif kernType == 2:
    outname = 'self_mobility_sc_w4.mat'

scipy.io.savemat(outname, dict(M=M, heights=heights/radP))
