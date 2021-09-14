#!/bin/bash

DPSTOKES_ROOT=
DPSTOKES_INSTALL=
CPU=GNU

#################### BEGIN USER EDIT ####################################

# number of threads for OpenMP
num_threads=10
# let the shell use the maximum amount of stack memory
ulimit -s unlimited
# set mem for thread stack
export OMP_STACKSIZE=256m

# thread pinning settings
CPU=GNU
if [ "$CPU" == "Intel" ]; then
  export MKL_THREADING_LAYER=sequential
  export KMP_AFFINITY="verbose,proclist=[0,1,2,3,4,5,6,7,8,9],explicit"
else
  export OMP_PLACES="{0}:10:1"
  export OMP_PROC_BIND=true
  export OMP_DISPLAY_ENV=true
fi

##################### END USER EDIT #######################################
sed -i "/num_threads/c num_threads=${num_threads}" ${DPSTOKES_ROOT}/source/cpu/python/config.py
export LD_LIBRARY_PATH=${DPSTOKES_INSTALL}:${LD_LIBRARY_PATH}
export PYTHONPATH=${DPSTOKES_ROOT}/source/cpu/python:${DPSTOKES_ROOT}/python_interface:${DPSTOKES_INSTALL}:${PYTHONPATH}
export FCM_CPUCONFIG_LAUNCHED=1
