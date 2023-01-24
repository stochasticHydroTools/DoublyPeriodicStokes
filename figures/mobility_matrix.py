import numpy as np
from common_interface_wrapper import *
from mobility_numba import *

def mobility_matrix(xPs, solver):
  """
  Compute the mobility matrix for each configuration in xPs(nConfig, 3 * nP)

  Parameters:
    xPs - particle positions (nConfig, 3 * nP)
    solver - FCMJoint solver object which has been initialized
  """

  nConfig, nP3 = xPs.shape; nP = int(nP3 / 3)
  eta = solver.viscosity
  has_torque = solver.has_torque
  radP = solver.radP
  # construct mobility matrix at each height
  if has_torque:
    M = np.zeros((nConfig, 6 * nP, 6 * nP), dtype = np.double)
    F = np.eye(6 * nP, dtype = np.double)
    NormMat = np.zeros((6 * nP, 6 * nP))
    NormMat[0:3*nP, 0:3*nP] = 1 / (6 * np.pi * eta * np.repeat(radP, 3)).reshape((1, 3 * nP))
    NormMat[0:3*nP, 3*nP::] = NormMat[3*nP::, 0:3*nP] = 1 / (6 * np.pi * eta * np.repeat(np.power(radP, 2), 3)).reshape(1, 3 * nP)
    NormMat[3*nP::, 3*nP::] = 1 / (8 * np.pi * eta * np.repeat(np.power(radP, 3), 3)).reshape(1, 3 * nP)
  else:
    M = np.zeros((nConfig, 3 * nP, 3 * nP), dtype = np.double)
    F = np.eye(3 * nP, dtype = np.double)
    NormMat = np.zeros((3 * nP, 3 * nP))
    NormMat[0:3*nP, 0:3*nP] = 1 / (6 * np.pi * eta * np.repeat(radP, 3)).reshape((1, 3 * nP))
  
  for iC in range(0, nConfig):
    print("config #", iC)
    xP = xPs[iC,:].copy()
    solver.SetPositions(xP)
    if has_torque:
      for j in range(0, 6 * nP):
        forces = F[0:3*nP,j].copy(); torques = F[3*nP::,j].copy();
        mob, mobt = solver.Mdot(forces, torques)
        M[iC,:,j] = np.concatenate((mob, mobt))
    else:
      for j in range(0, 3 * nP):
        forces = F[0:3*nP,j].copy()
        mob, _ = solver.Mdot(forces)
        M[iC,:,j] = mob
    M[iC,:,:] /= NormMat

  return M, NormMat

# periodized or aperiodic RPY approximation for bottom wall only
def mobility_matrix_numba(xPs, eta, L, a, num_imgs, self):
  nConfig, nP3 = xPs.shape; nP = int(nP3 / 3); 
  if self:
    nP = 1
  radP = a * np.ones(nP)
  # construct mobility matrix at each height
  M = np.zeros((nConfig, 6 * nP, 6 * nP), dtype = np.double)
  F = np.eye(6 * nP, dtype = np.double)
  NormMat = np.zeros((6 * nP, 6 * nP))
  NormMat[0:3*nP, 0:3*nP] = 1 / (6 * np.pi * eta * np.repeat(radP, 3)).reshape((1, 3 * nP))
  NormMat[0:3*nP, 3*nP::] = NormMat[3*nP::, 0:3*nP] = 1 / (6 * np.pi * eta * np.repeat(np.power(radP, 2), 3)).reshape(1, 3 * nP)
  NormMat[3*nP::, 3*nP::] = 1 / (8 * np.pi * eta * np.repeat(np.power(radP, 3), 3)).reshape(1, 3 * nP)
  
  for iC in range(0, nConfig):
    print("config #", iC)
    xP = xPs[iC,:].copy()
    for j in range(0, 6 * nP):
      forces = F[0:3*nP,j].copy(); torques = F[3*nP::,j].copy();
      if self:
        forces = np.concatenate((forces,[0,0,0]))
        torques = np.concatenate((torques,[0,0,0]))
      mob = single_wall_mobility_trans_times_force_numba(xP, forces, eta, a, L, num_imgs) +\
            single_wall_mobility_trans_times_torque_numba(xP, torques, eta, a, L, num_imgs)
      mobt = single_wall_mobility_rot_times_force_numba(xP, forces, eta, a, L, num_imgs) +\
             single_wall_mobility_rot_times_torque_numba(xP, torques, eta, a, L, num_imgs)
      if self:
        mob = mob[3::]; mobt = mobt[3::]
      M[iC,:,j] = np.concatenate((mob, mobt))
    M[iC,:,:] /= NormMat
  return M, NormMat
