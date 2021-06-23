#Raul P. Pelaez 2021. Example of the python interface to DPStokes
#This code can use two uammd hydrodynamics modules:
# -DPStokes for doubly periodic hydrodynamics with either no walls (nowall), a bottom wall (bottom) or two walls (slit)
# -FCM for triply periodic hydrodynamics
#The parameters related with the spread/interp kernel are a little precarious now, waiting for a better unification
#For DPStokes, I recomment setting nxy, nz and fix support=6. gw is ignored
#For FCM, I recomend setting nxy, nz and tolerance and let the module compute the support accordingly.
import numpy as np
import uammd
Lx = Ly = 32.0
H=32.0
mode='slit'
#mode='bottom'   #Bottom wall geometry
#mode='slit'     #Slit channel
#mode='nowall'   #Doubly Periodic without walls
#mode='periodic' #Triply Periodic using FCM

hydrodynamicRadius = 1.0;
viscosity = 1/(6*np.pi);

#Number of grid points in the plane
nx = ny = 48;
#Number of grid points in the chebyshev grid (just nz, not (2*nz-2))
nz = 92;
#Width of ES kernel
w = 6.0;
w_d = 6.0;
#Beta parameter of ES kernel (just beta, not beta*w)
beta =2*np.pi;
beta_d = 2*np.pi;
#Alpha parameter of ES kernel
h=Lx/nx
alpha = w*h*0.5;
alpha_d = w_d*h*0.5;
par = uammd.StokesParameters(viscosity=viscosity,
                             Lx=Lx,Ly=Ly,
                             zmin=-H*0.5, zmax=H*0.5,
                             w=w, w_d=w_d,
                             alpha=alpha, alpha_d=alpha_d,
                             beta=beta, beta_d=beta_d,
                             hydrodynamicRadius=hydrodynamicRadius,
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
