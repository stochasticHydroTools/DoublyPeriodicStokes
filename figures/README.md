This directory contains scripts for running mobility tests with the python interface to 
the `FCMJoint` module.

## Scripts 
self_mobility.py - compute the self mobility of a single particle at various heights
                   above the bottom wall in DPBW or DPSC domains
                 - reference computations on a doubled grid
                 - reference computation in DPBW using periodized RPY 
                 
pair_mobility.py - compute the pair mobility of two particles at various heights
                   above the bottom wall, and separations from each other in DPBW or DPSC.
                   
mobility_matrix.py - compute the mobility matrix for several particle configurations using
                     the `FCMJoint` interface
                     
mobility_numba.py - periodized (and regular) RPY implementation with bottomw wall corrections

pair_mobility_matrix_asym.py - check the asymmetry and positive definiteness of the mobility
                               matrix in DPBW or DPSC over several trials



Also included are scripts and raw data to regenerate all of the figures.