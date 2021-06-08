#Raul P. Pelaez 2021. Example of the python interface to DPStokes
#This code can use two uammd hydrodynamics modules:
# -DPStokes for doubly periodic hydrodynamics with either no walls (nowall), a bottom wall (bottom) or two walls (slit)
# -FCM for triply periodic hydrodynamics
#The parameters related with the spread/interp kernel are a little precarious now, waiting for a better unification
#For DPStokes, I recomment setting nxy, nz and fix support=6. gw is ignored
#For FCM, I recomend setting gw and tolerance and let the module compute nxy, nz and support accordingly.
import numpy as np
import uammd
Lxy = 32
H=32
mode='slit'
#mode='bottom'   #Bottom wall geometry
#mode='slit'     #Slit channel
#mode='nowall'   #Doubly Periodic without walls
#mode='periodic' #Triply Periodic using FCM
par = uammd.StokesParameters(Lxy=Lxy, H=H, gw = 1.0, Nxy = 48,
                               nz=92, support=6, viscosity=1/(6*np.pi), mode=mode)
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
