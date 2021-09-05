#!/bin/bash
#Define this env that informs the module of the execution of this script
export FCM_CPUCONFIG_LAUNCHED=1
# NOTE: users should only modify num_threads, OMP_STACKSIZE and thread pinning settings
#       do not change OMP_NUM_THREADS or the nested parallelism settings

# number of threads for OpenMP
num_threads=6
# let the shell use the maximum amount of stack memory
ulimit -s unlimited
# set mem for thread stack
export OMP_STACKSIZE=256m
# set num threads env var to 1 and disable nested parallelism
export OMP_NUM_THREADS=1
export OMP_NESTED=false
export OMP_MAX_ACTIVE_LEVELS=1 # for newer omp versions
# thread pinning settings
#export OMP_PLACES="{0}:6:1"
#export OMP_PROC_BIND=true
#export OMP_DISPLAY_ENV=true

# uncomment below if using Intel compiler

#source /opt/intel/mkl/bin/mklvars.sh intel64  
#export MKL_THREADING_LAYER=sequential
#export KMP_AFFINITY="verbose,proclist=[0,1,2,3,4,5,6,7,8,9],explicit"

sed -i "/num_threads/c num_threads=${num_threads}" /home/raul/Dropbox/Trabajo/uammd_python/DPStokesTests/source/cpu/python/config.py
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/raul/Dropbox/Trabajo/uammd_python/DPStokesTests/source/cpu/lib
export PYTHONPATH=${PYTHONPATH}:/home/raul/Dropbox/Trabajo/uammd_python/DPStokesTests/source/cpu/python
