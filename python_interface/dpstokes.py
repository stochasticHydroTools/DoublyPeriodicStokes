#Raul P. Pelaez 2021. Example of the python interface to DPStokes
#This code can use two uammd hydrodynamics modules:
# -DPStokes for doubly periodic hydrodynamics with either no walls (nowall), a bottom wall (bottom) or two walls (slit)
# -FCM for triply periodic hydrodynamics
#The parameters related with the spread/interp kernel are a little precarious now, waiting for a better unification
#For FCM, I recomend setting Lxy, Lz, hydrodynamic radius and tolerance and let the module compute everything else accordingly.

#About periodic wrapping
#Periodic wrapping for the particles is handled internally, so the origin of the particles in the plane is irrelevant.
#For the z direction we have to make a distinction between TP and DP:
# -The triply periodic case will simply interpret H=zmax-zmin, so the actual origin of the particles is also irrelevant.
# -For the doubly periodic version, zmax and zmin should correspond to the actual maximum and minimum height of the particles in the domain. A particle having z>zmax or z<zmin results in undefined behavior.
#   If a bottom wall geometry is used, zmin will correspond to the height of the bottom wall. Similarly, an additional wall will be placed at zmax in the slit channel mode.
import numpy as np
import uammd
mode='slit'
#mode='bottom'   #Bottom wall geometry
#mode='slit'     #Slit channel
#mode='nowall'   #Doubly Periodic without walls
#mode='periodic' #Triply Periodic using FCM


viscosity = 1/(6*np.pi);
#About the parameters below:
#  An user should only worry about the parameters here (Lxy, zmin, zmax, hydroRadius, w, and w_d).
#  The rest of the parameters should be considered implementation details.
#  For testing purposes and convenience the user can also tweak all the other parameters
#  However, note that some of them are incompatible.
#    In particular, if present, the hydroRadius will control every other parameter (nxyz, betasand alphas). Manually setting these should be considered an advanced usage of the interface.
# The default values for the parameters that can be autocomputed is -1, if a combination of incompatible parameters is provided a warning, or error, will be issued (for instance setting hydroRadius and beta)
#NOTE: As of today, the code is not actually able to autocompute the parameters from the hydroRadius in DP mode.
#hydrodynamicRadius = 1.0; 
Lx = Ly = 32.0
zmin = -16.0
zmax = 16.0
#Width of ES kernel
w = 6.0;
w_d = 6.0;
#NOTE: The TP interface as of today ignores all parameters related with the kernel, in that case you can choose and hydrodynamic radius and (optionally) a number of cells (which will be defaulted to enforce a certain default tolerance). Note that in this case the hydrodynamic radius might not be posible to enforce exactly, in which case a warning is issued with the closest possible chosen.

#In DP (and i the future also TP), instead of choosing the hydroRadius and letting the code compute the rest internally, 
# the user can leave hydroRadius at the default and set all the parameters below instead.
# Now, this is here to ease testing while we figure out how to autocompute them.
#Number of grid points in the plane
nx = ny = 48;
#Number of grid points in the chebyshev grid (just nz, not (2*nz-2))
nz = 92;
#Beta parameter of ES kernel (just beta, not beta*w)
beta = 1.8;
beta_d = 1.8;
#Alpha parameter of ES kernel.
# If alpha<w*h*0.5 an error will be thrown, since it could result in the kernel trying to evaluate sqrt of a negative number.  
h=Lx/nx
alpha = w*h*0.5;
alpha_d = w_d*h*0.5;


par = uammd.StokesParameters(viscosity=viscosity,
                             Lx=Lx,Ly=Ly,
                             zmin=zmin, zmax=zmax,
                             w=w, w_d=w_d,
                             alpha=alpha, alpha_d=alpha_d,
                             beta=beta, beta_d=beta_d,
#                             hydrodynamicRadius=hydrodynamicRadius,
                             nx=nx, ny=ny, nz=nz, mode=mode)
print(par)
numberParticles = 2
dpstokes = uammd.DPStokes(par, numberParticles)
precision = np.float64;
positions = np.array([0, 0, -3, 0, 0, 3], precision);
forces = np.array([1, 1, 1, -1, -1, -1], precision);
torques = np.array([1, 1, 1, 1, 1, 1], precision);
#Linear displacements
MF=np.zeros(3*numberParticles, precision);
#Angular displacements
MT=np.zeros(3*numberParticles, precision);
dpstokes.Mdot(positions,forces, torques, MF, MT)

MF = np.reshape(MF,(numberParticles,3));
print("MF: ")
print(MF)

MT = np.reshape(MT,(numberParticles,3));
print("MT: ")
print(MT)
