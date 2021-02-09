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
                               nz=92, support=13, viscosity=1/(6*np.pi), mode=mode)
print(par)
numberParticles = 2
dpstokes = uammd.DPStokes(par, numberParticles)
precision = np.float64;
positions = np.array([0, 0, -3, 0, 0, 3], precision);
forces = np.array([1, 1, 1, -1, -1, -1], precision);
MF=np.zeros(3*numberParticles, precision);
dpstokes.Mdot(positions,forces, MF)
MF = np.reshape(MF,(numberParticles,3));

print("MF: ")
print(MF)
