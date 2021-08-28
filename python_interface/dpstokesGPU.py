#Raul P. Pelaez 2021. Example of the python interface to some UAMMD hydrodynamic modules
#This interface allows to compute the hydrodynamic displacements of a group of particles due to
# forces and/or torques acting on them:
# [Displacements; AngularDisplacements] = Mobility·[Forces; Torques]
#Hydrodynamic displacements can be computed in several domain geometries.
#In particular, this interface allows to access two uammd hydrodynamics modules:
# -DPStokes for doubly periodic (DP) hydrodynamics with either no walls (nowall), a bottom wall (bottom) or two walls (slit)
# -FCM for triply periodic (TP) hydrodynamics

#Both algorithms are grid based. Communication of the particles' forces and torques is carried out via the ES kernel.
#In TP, the grid is regular in all directions, while in DP the grid is defined in the Chebishev basis in the third direction.

#USAGE:
#0- Compile the library. It is required to pass arrays to the interface with the same precision as the library was compiled in. By default, the Makefile will compile UAMMD in single precision (corresponding to np.float32).
##1- Import the uammd python wrapper:
# import uammd
##You will need numpy for passing data to uammd:
# import numpy as np
##2- Create the DPStokes object:
# dpstokes = uammd.DPStokes()
##3- Encapsulate the parameters in the uammd.StokesParameters object:
## A list of parameters and their meaning can be found below
# par = uammd.StokesParameters(viscosity=viscosity,
#                             Lx=Lx,Ly=Ly,
#                             zmin=zmin, zmax=zmax,
#                             w=w, w_d=w_d,
#                             alpha=alpha, alpha_d=alpha_d,
#                             beta=beta, beta_d=beta_d,
#                             nx=nx, ny=ny, nz=nz, mode=mode)
##4- Initialize DPStokes with the parameters and the number of particles:
# dpstokes.initialize(par, numberParticles)
## If some parameters need to change, simply call initialize again with the new ones.
## Mind you, initialization is in general a really slow operation.
##5- Set the positions to construct the mobility matrix with:
## All arrays have the format [x0 y0 z0 x1 y1 z1 ...] (interleaved)
## For instance:
# positions = np.array([0, 0, -3, 0, 0, 3], np.float32)
# Corresponds to two particles, the first one located at (0,0,-3) and the second at (0,0,3)
# dpstokes.setPositions(positions)
##6- Compute hydrodynamic displacements for a group of forces and/or torques acting on the previously set positions
# dpstokes.Mdot(forces=forces, torques=torques, velocities=MF, angularVelocities=MT)
## The torques and angularVelocities can be omited, in which case they are assumed to be zero:
# dpstokes.Mdot(forces=forces, velocities=MF)
## Both MF and MT are output arguemnts and must be numpy arrays with the precision for which the library was compiled for.
## Otherwise the output will be ignored by python.
## The contents of both arrays will be overwritten.
##7- Clean up any memory allocated by the module, which will remain in an unusable state until initialization is called again
# dpstokes.clear()



#LIST OF PARAMETERS:

#viscosity: The viscosity of the solvent

#mode:      The domain geometry, can be any of the following:
#       -"periodic": TP FCM
#       -"nowall":   DP, open boundaries in the third direction
#       -"bottom":   DP, there is a wall at the bottom of the domain
#       -"slit":     DP, there are walls at both the bottom and top of the domain

#Lx,Ly:     The domain size in the plane (which is always periodic). Periodic wrapping is carried out by UAMMD, so the positions' origin in the plane is irrelevant.

#zmin, zmax: The domain size in the Z direction. The meaning of this parameter changes with each mode:
#             -In TP the domain is periodic in z with a size Lz=zmax-zmin (the origin of the particles' positions is irrelevant
#             -In DP zmin and zmax denote the allowed range for the heights of the particles. In the DP modes particles must always be contained between zmin and zmax.
#                *For "nowall" this range is not physical and can be considered a requirement of the implementation (and similarly for zmax in the "bottom" wall case).
#                *For walled modes ("bottom" and "slit") zmin/zmax correspond to the locations of the bottom/top walls respectively.

#w, alpha, beta: Parameters of the ES kernel to spread forces
#w_d, alpha_d, beta_d: Parameters of the ES kernel to spread torques.
#     In both TP and DP torques are spreaded by performing the curl in fourier space, allowing to spread the torques using the same kernel as for the forces (instead of its derivative).
#     In all cases alpha defaults to w*h*0.5, where h is the grid size (Lx/nx).

#nx, ny, nz: Number of grid cells in each direction.

import numpy as np
import uammd
np.random.seed(0)
#Precision for the data arrays, ensure it is the same as the compiled library (single by default)
precision = np.float32

#Some arbitrary parameters

mode='bottom'    #Doubly Periodic with a bottom wall
#mode='slit'     #Doubly Periodic slit channel
#mode='nowall'   #Doubly Periodic without walls
#mode='periodic' #Triply Periodic

numberParticles = 2048
viscosity = 0.957e-3
Lx = Ly = 128.7923000000000116
nx = ny = 216
zmax =  9.1740639106166668
zmin = 0
nz = 26
#The parameters of the ES kernel will be related to the hydrodynamic radius
w = w_d = 6
beta = 8.264036224425126
beta_d = 13.779780233768633

#Create the handle to the module. Several instances can be created if necessary.
#They will share memory, so two identical DPStokes instances will take approx. the same space as one.
dpstokes = uammd.DPStokes()

#Pack all the parameters
par = uammd.StokesParameters(viscosity=viscosity,
                             Lx=Lx,Ly=Ly,
                             zmin=zmin, zmax=zmax,
                             w=w, w_d=w_d,
#                             alpha=alpha, alpha_d=alpha_d,
                             beta=beta, beta_d=beta_d,
                             nx=nx, ny=ny, nz=nz, mode=mode)
print(par)

#Initialize the module
dpstokes.initialize(par, numberParticles)

#Set positions to construct the mobility matrix with
positions = np.zeros(numberParticles*3, precision) # Arbitrary positions inside the Z domain
positions[::3] = np.random.uniform(0, Lx, numberParticles) #X
positions[1::3] = np.random.uniform(0, Lx, numberParticles) #Y
positions[2::3] = np.random.uniform(zmin, zmax, numberParticles) #Z
dpstokes.setPositions(positions)

#Compute displacements and store them in MF, MT.
torques = np.random.normal(0,1,3*numberParticles) #Random torques
forces = np.random.normal(0,1,3*numberParticles) #Random forces

MT = np.zeros(3*numberParticles, precision)
MF = np.zeros(3*numberParticles, precision)
dpstokes.Mdot(forces=forces, torques=torques, velocities=MF, angularVelocities=MT)
#dpstokes.Mdot(forces=forces, velocities=MF) #Alternative without torques

#Clear all memory
dpstokes.clear()



#Print something
print("Positions of the first three particles: ")
positions = np.reshape(positions,(numberParticles,3));
print(positions[0:3])

MF = np.reshape(MF,(numberParticles,3));
print("Displacements of the first three particles: ")
print(MF[0:3])

MT = np.reshape(MT,(numberParticles,3));
print("Angular displacements of the first 3 particles: ")
print(MT[0:3])
