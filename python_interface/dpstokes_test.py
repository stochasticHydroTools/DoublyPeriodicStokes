#Raul P. Pelaez 2021. Test for the DP algorithm, compares with Sachin's data.
# Donev: Fix comments to say "compares GPU and CPU versions" and not "Sachin's data" if you intend this to be in the public release
# For instructions on usage see dpstokes.py
# This script will not work without Sachin's comparison data files (called fcm_multiblob_compare.tgz, not included in the repo)

import numpy as np
import uammd
import sys
from timeit import default_timer as timer
precision = np.float32

mode='bottom'   #Bottom wall geometry
#mode='slit'     #Slit channel
#mode='nowall'   #Doubly Periodic without walls
#mode='periodic' #Triply Periodic using FCM

#These two sets of parameters are for comparing accuracy
numberParticles = 2048
folder="fcm_multiblob_compare/trans_rot_w6/"
inputFile = folder+"monopole_position_and_force.txt"
dotorques = True
viscosity = 0.957e-3
Lx = Ly = 128.7923000000000116
nx = ny = 216
zmax =  9.1740639106166668
nz = 26 
w = w_d = 6
beta = 8.264036224425126
beta_d = 13.779780233768633
randomFT=False


#numberParticles = 100
#folder="compare_data/bottom_wall/trans_rot_w6/"
#inputFile = folder+"monopole_position_and_force.txt"
#dotorques = True
#viscosity = 1/(4.0*np.sqrt(np.pi))
#Lx = Ly = 82.0572088991620916
#nx = ny = 142
#zmax =  24.6963437569977771
#nz = 68
#w = w_d = 6
#beta = 1.3267*6
#beta_d = 2.216*6
#randomFT=False

#These two are for measuring performance
#numberParticles = 24576
#dotorques = True
#viscosity = 0.957e-3
#Lx = Ly = 128.7923000000000116
#nx = ny = 528
#zmax =  9.1740639106166668 #67.5920828561666696 
#nz = 55 #181 
#w = w_d = 6
#beta = 8.024707442720228
#beta_d = 13.419437243071116
#randomFT=True

#numberParticles = 86016
#dotorques = True
#viscosity = 0.957e-3
#Lx = Ly = 128.7923000000000116
#nx = ny = 900
#zmax =  9.1740639106166668 #67.5920828561666696 
#nz = 91
#w = w_d = 6
#beta = 7.973081016984452
#beta_d = 13.341819112745952
#randomFT=True
#inputFile = "multiblob/monopole_position.42txt"


dpstokes = uammd.DPStokes()
par = uammd.StokesParameters(viscosity=viscosity,
                             Lx=Lx,Ly=Ly,
                             zmin=0, zmax=zmax,
                             w=w, w_d=w_d,
                             beta=beta, beta_d=beta_d,
                             nx=nx, ny=ny, nz=nz, mode=mode)
print(par)

positions = np.zeros(numberParticles*3, precision)
forces = np.zeros(numberParticles*3, precision)
torques = np.zeros(numberParticles*3, precision)

#Read positions and forces
with open(inputFile) as f:
    content = f.readlines()
    lineNumber = 0
    for line in content:
        lineNumber += 1
        i = lineNumber
        if i > numberParticles:
            break
        numbers = np.array([np.float64(x) for x in line.split()])
        numbers[0] -= par.Lx*0.5
        numbers[1] -= par.Ly*0.5
        positions[3*(i-1):3*i] = numbers[0:3]
        if randomFT:
            forces[3*(i-1):3*i] = np.random.normal(0,1,3)
        else:
            forces[3*(i-1):3*i] = numbers[3:6]

#Read torques if necessary
if dotorques:
    if randomFT:
        torques = np.random.normal(0,1,3*numberParticles)
    else:
        inputFile = folder+"dipole_position_and_torque.txt"
        with open(inputFile) as f:
            content = f.readlines()
            lineNumber = 0
            for line in content:
                lineNumber += 1
                i = lineNumber
                if i > numberParticles:
                    break
                numbers = np.array([np.float64(x) for x in line.split()])        
                torques[3*(i-1):3*i] = numbers[3:6]



dpstokes.initialize(par, numberParticles)


#Compute the same several times to measure performance
ttot = timer()
Ntest=10
for i in range(0,Ntest):
    MF = np.zeros(3*numberParticles, precision)
    dpstokes.setPositions(positions)
    if dotorques:
        MT = np.zeros(3*numberParticles, precision)
        t0 = timer()
        dpstokes.Mdot(forces=forces, torques=torques, velocities=MF, angularVelocities=MT)
        t1 = timer()
        print("Time : ", t1-t0)
    else:
        dpstokes.Mdot(forces=forces, velocities=MF)
tfin = timer()
print("Time total : ", (tfin- ttot)/Ntest*1e3, "ms")
dpstokes.clear()


#If measuring accuracy compare results to reference
if not randomFT:
    inputFile = folder + "monopole_position_and_velocity.txt"
    with open(inputFile) as f:
        content = f.readlines()
        vel_reference = np.zeros(numberParticles*3, precision)
        lineNumber = 0
        for line in content:
            lineNumber += 1
            i = lineNumber
            if i > numberParticles:
                break
            numbers = np.array([np.float64(x) for x in line.split()])
            vel_reference[3*(i-1):3*i] = numbers[3:6]    
        vel_reference = np.reshape(vel_reference, (numberParticles, 3))
        MF = np.reshape(MF, (numberParticles, 3))
        print("MF: ")
        print(MF[0:3])
            
        print("reference: ")
        print(vel_reference[0:3])
            
        err = (MF - vel_reference)/(np.max(np.abs(vel_reference)))
        print("err: ")
        print(err[0:3])
        print("fac : ")
        print(MF[0:3]/vel_reference[0:3])
        print("max err: ")
        print(np.max(np.abs(err[0::3])), np.sqrt(np.var(err)))
        print(np.max(np.abs(err[1::3])), np.sqrt(np.var(err)))
        print(np.max(np.abs(err[2::3])), np.sqrt(np.var(err)))	
        if dotorques:
            inputFile = folder + "dipole_position_and_velocity.txt"
            angvel_reference = np.zeros(numberParticles*3, precision)
            with open(inputFile) as f:
                content = f.readlines()
                lineNumber = 0
                for line in content:
                    lineNumber += 1
                    i = lineNumber
                    if i > numberParticles:
                        break
                    numbers = np.array([np.float64(x) for x in line.split()])
                    angvel_reference[3*(i-1):3*i] = numbers[3:6]    
                angvel_reference = np.reshape(angvel_reference, (numberParticles, 3))
                MT = np.reshape(MT, (numberParticles, 3))
                print("MT: ")
                print(MT[0:3])
                print("reference: ")
                print(angvel_reference[0:3])
                
                err = (MT - angvel_reference)/np.max(np.abs(angvel_reference))
                print("err: ")
                print(err[0:3])
                print("max err: ")
                print(np.max(np.abs(err)), np.sqrt(np.var(err)))
