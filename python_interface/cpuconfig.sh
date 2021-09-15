#!/bin/bash

DPSTOKES_ROOT=/home/sachin/Research/NYU/ParticleSuspension/DoublyPeriodicStokes
DPSTOKES_INSTALL=/home/sachin/Research/NYU/ParticleSuspension/DoublyPeriodicStokes/python_interface/lib
CPU=Intel

#################### BEGIN USER EDIT ####################################

# number of threads for OpenMP
num_threads=6
# let the shell use the maximum amount of stack memory
ulimit -s unlimited
# set mem for thread stack
export OMP_STACKSIZE=256m
pin_threads=false

# thread pinning settings
if [ "$CPU" == "Intel" -a "$pin_threads" == "true" ]; then
  export MKL_THREADING_LAYER=sequential
  export KMP_AFFINITY="verbose,proclist=[0,1,2,3,4,5,6],explicit"
elif [ "$pin_threads" == "true" ]; then
  export OMP_PLACES="{0}:6:1"
  export OMP_PROC_BIND=true
  export OMP_DISPLAY_ENV=true
fi

##################### END USER EDIT #######################################
sed -i "/num_threads/c num_threads=${num_threads}" ${DPSTOKES_ROOT}/source/cpu/python/config.py
export LD_LIBRARY_PATH=${DPSTOKES_INSTALL}:${LD_LIBRARY_PATH}
export PYTHONPATH=${DPSTOKES_ROOT}/source/cpu/python:${DPSTOKES_ROOT}/python_interface:${DPSTOKES_INSTALL}:${PYTHONPATH}
export FCM_CPUCONFIG_LAUNCHED=1
