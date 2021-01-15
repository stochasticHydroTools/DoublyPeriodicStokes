import numpy as np
import uammd
Lxy = 32
H=32
par = uammd.DPStokesParameters(Lxy=Lxy, H=H, gw = 1.0, Nxy=32, nz=32, support=13, viscosity=1/(6*np.pi))
print(par)
numberParticles = 2

dpstokes = uammd.DPStokes(par, numberParticles)

precision = np.float32;

positions = np.array([0, 0, -1, 0, 0, 1], precision);
forces = np.array([1, 1, 1, -1, -1, -1], precision);
MF=np.zeros(3*numberParticles, precision);

dpstokes.Mdot(positions,forces, MF)
MF = np.reshape(MF,(numberParticles,3));

print("MF: ")
print(MF)
                
