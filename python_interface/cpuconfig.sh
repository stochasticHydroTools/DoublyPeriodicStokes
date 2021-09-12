#!/bin/bash

SRC_DIR=/home/srn324/DoublyPeriodicStokes
INSTALL_DIR=/home/srn324/DoublyPeriodicStokes/python_interface/lib
CPU=Intel

#################### BEGIN USER EDIT ####################################

# number of threads for OpenMP
num_threads=10
# let the shell use the maximum amount of stack memory
ulimit -s unlimited
# set mem for thread stack
export OMP_STACKSIZE=256m

# thread pinning settings
if [ "$CPU" == "GNU" ]; then
  export OMP_PLACES="{0}:10:1"
  export OMP_PROC_BIND=true
  export OMP_DISPLAY_ENV=true
elif [ "$CPU" == "Intel" ]; then
  export MKL_THREADING_LAYER=sequential
  export KMP_AFFINITY="verbose,proclist=[0,1,2,3,4,5,6,7,8,9],explicit"
fi

##################### END USER EDIT #######################################
sed -i "/num_threads/c num_threads=${num_threads}" ${SRC_DIR}/source/cpu/python/config.py
export LD_LIBRARY_PATH=${INSTALL_DIR}:${LD_LIBRARY_PATH}
export PYTHONPATH=${SRC_DIR}/source/cpu/python:${SRC_DIR}/python_interface:${INSTALL_DIR}:${PYTHONPATH}
export FCM_CPUCONFIG_LAUNCHED=1
