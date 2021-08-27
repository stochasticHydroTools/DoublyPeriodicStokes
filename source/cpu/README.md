# FCM for Stokes mobility problems with the ES kernel
A spectral immersed boundary method for particle suspensions in Stokes flow on doubly and triply periodic domains.

### Build Dependencies ###
You will need to `apt install` at least the following dependencies:

* build-essential
* cmake3
* libomp-dev
* gcc 7.5.0 or later (eg. module load gcc-9.2 on cims machines)
* FFTW
* LAPACK and the C interface LAPACKE 
    (https://www.assistedcoding.eu/2017/11/04/how-to-install-lapacke-ubuntu/)
* Python 3+ with NumPy/SciPy 

### Build Instructions - GNU Compiler ###
Now, execute the following from the top of the source tree: 
```
$ mkdir build && cd build
$ cmake3 -Dfftw_wisdom=on -Duse_stack=on ..
$ make -j6 (or however many threads you'd like to use)
$ make install
```
Executing the commands above will build all libraries and executables. The libraries are
installed in `$INSTALL_PATH/lib`. Executables are installed in `$INSTALL_PATH/bin`. 
By default, `$INSTALL_PATH` is the top of the source tree. Setting the `fftw_wisdom`
flag to off will build the libraries using fftw without the wisdom utility, as
described below. Build with the `use_stack` flag set to off if you run into
stack overflow errors during calls to spread/interp.

### Build Instructions - Intel Compiler ###
If using the Intel compilers, there is an additional step of downloading and
installing FFTW locally.

First, load the Intel compilers with something like
```
$ module load intel-2019
```
or
```
$ source /opt/intel/bin/compilervars.sh intel64 
$ source /opt/intel/mkl/bin/mklvars.sh intel64
```
Then, download and install FFTW with

```
$ wget ftp://ftp.fftw.org/pub/fftw/fftw-3.3.9.tar.gz
$ tar xvzf fftw-3.3.9.tar.gz
$ cd fftw-3.3.9
$ sed -i 's/fopenmp/qopenmp/g' configure
$ CC=icc F77=ifort ./configure --prefix=$SRC_DIR/fftw_install --enable-shared --enable-openmp --enable-sse2 --enable-avx --enable-avx2
$ make
$ make install
```
where `$SRC_DIR` is the top level of the project source tree (eg. `/path/to/DoublyPeriodicStokes/source/cpu`)

Lastly, build the libraries with 
```
$ mkdir build && cd build
$ cmake3 -Dfftw_wisdom=on -Duse_stack=on -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc ..
$ make -j6 (or however many threads you'd like to use)
$ make install
```

### Usage ###

#### Environment Settings ####
Before using any function (python or c++), set the following environment
variables for OpenMP (eg. in bash):
```
$ ulimit -s unlimited
$ export OMP_STACKSIZE=256m
$ export OMP_NESTED=false (or OMP_MAX_ACTIVE_LEVELS=1)
$ export OMP_NUM_THREADS=1
$ export OMP_MAX_ACTIVE_LEVELS=1
```

The `PYTHONPATH` and `LD_LIBRARY_PATH` have to be expanded with the location
of the Python files (in `python` folder), and the shared libraries (in `lib` folder):
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/path/to/DoublyPeriodicStokes/source/cpu/lib
export PYTHONPATH=${PYTHONPATH}:/path/to/DoublyPeriodicStokes/source/cpu/python

```
Then, set the `num_threads` variable in the file `python/config.py` to
the number of threads you want to use, or accomplish this with `sed`.
```
$ num_threads=10
$ sed -i "/num_threads/c num_threads=${num_threads}" /path/to/DoublyPeriodicStokes/source/cpu/python/config.py
```
All OpenMP parallel loops throughout the library will use this number of 
threads, and nested parallelism is disabled. 

A convenience script `config.sh` is available in the `examples` folder after installation. 
Executing the following will accomplish the above snippets, and the file can be 
edited as needed. It must be sourced every time a new shell is used for a run.
```
$ source config.sh
```
For thread pinning, one must use different settings depending on how the library was built. These
can be added to config.sh. For example, for a 10 core single socket machine with GNU compilation
```
$ export OMP_PLACES="{0}:10:1"
$ export OMP_PROC_BIND=true

```
or with Intel compilation
```
source /opt/intel/mkl/bin/mklvars.sh intel64  
export MKL_THREADING_LAYER=sequential
export KMP_AFFINITY="verbose,proclist=[0,1,2,3,4,5,6,7,8,9],explicit"
```

#### FFTW Settings ####
Each solver works on Fourier or Fourier-Chebyshev coefficients of
3D vector fields computed with FFTW. If the library is built with
`fftw_wisdom=on` passed to `CMake`, the first time a transform of a 
given size is requested, the FFT forward and backward plans are 
saved to disk in the folder `./fftw_wisdom`.

The first call incurs some startup cost, since we use FFTW_MEASURE - one 
of the slower flags for fftw which produces a more optimal algorithm.

The next time a transform of the same size is requested, the existing
wisdom is read from disk, signficantly reducing planning time, but giving
the same optimal algorithm.

Optimality can be further improved by switching FFTW_MEASURE to FFTW_PATIENT
by adding `-Dfftw_wisdom_patient=on` to the `cmake` command

For example, one way to generate wisdom beforehand whithout actually running 
a solver is to execute a python script with something like the following:

```
from Transform import *

# size of grid
Nx = Ny = 128; Nz = 65; dof = 3
# type of transform (0 = Fourier, 1 = Fourier-Cheb)
tType = 1
# Initialize plans, generate wisdom if it doesn't already exist.
# Assuming num_threads = 6 in config.py, the call will search for, 
# or create files in fftw_wisdom/ called:
# fftw_forward_wisdom_Nx128_Ny128_Nz126_dof3_nthr6_measure
# fftw_backward_wisdom_Nx128_Ny128_Nz126_dof3_nthr6_measure
transformer = Transformer(Nx, Ny, Nz, dof, tType)
# clean memory
transformer.Clean()
```

The next time a transform of the same size and type is requested, the wisdom
will be loaded from disk.

NOTE: Generated wisdom is dependent on
```
1) The number of threads set in config.py
2) The OpenMP environment settings
    - Eg) OMP_NESTED, OMP_PROC_BIND, OMP_PLACES
```
Trying to use wisdom generated with a different number of threads, or 
even the correct number of threads but a different OpenMP env config
will cause the planner to fail/stall! It is also recommended to 
regenerate wisdom each time the library is recompiled.

You can disable the use of wisdom by setting `fftw_wisdom=off` in the 
command line input to `CMake`. In this case, FFTW_ESTIMATE is the
planner mode used for all transforms. There is minimal cost
in computing such plans, though the resulting transform will be
musch slower for larger transform sizes.


#### Using the Python Interface ####
Each function in the Python interface (in the `python` folder) is equipped with a doc string, 
and calling `help(module_name)` in Python will print the documentation for all functions in the module. 
In general, each file in the `python` folder exposes a standalone C library which implements a
step in the FCM solution process:
 
- defining particles and grids (`Particles.py`,`Grid.py`, `Chebyshev.py`)
- configuring particle kernels and grids (`GridAndKernelConfig.py`)
- creating FFTW plans (`Transform.py`)
- spreading and interpolation (`SpreadInterp.py`)
- spectral Stokes solver (`Solvers.py`)
- FCM interface (`FCM.py`)

Below is an example script for using the `FCM` module, also provided in `examples/fcm_example.py`.

```
import sys
from FCM import *

# example usage of the FCM module

pos  = np.loadtxt('./Test_Data_For_Rollers/Const_Torque_t_15.clones', skiprows=1, usecols=[0,1,2])
data = np.loadtxt('./Test_Data_For_Rollers/One_Blob/N_Images_64.txt')

nP = 2048; 
domType = 'DPBW'
eta = 0.957e-3
has_torque = True
minX = 0.0; maxX = 128.7923
minY = 0.0; maxY = 128.7923
minZ = 0.0; maxZ = 20.0
xP = np.reshape(pos, (3 * nP,)).copy()
F = np.concatenate((data[:,0],data[:,1]))
radP = 1.0155 * np.ones(nP, dtype = np.double)
kernTypes = np.zeros(nP, dtype = np.int)

problem = FCM(radP, kernTypes, domType, has_torque)
problem.SetUnitCell([minX,maxX], [minY,maxY], [minZ,maxZ])
problem.Initialize(eta, 0)
problem.SetPositions(xP)

V = problem.Mdot(F)

# dummy update
xP += 0.1 * np.abs(V[0:3 * nP])
problem.SetPositions(xP)

V = problem.Mdot(F)

problem.Clean()

```
An example workflow for using this script in Bash from a newly created `run` 
directory is provided below:
```
$ cd && mkdir run && cd run
$ cp /path/to/DoublyPeriodicStokes/source/cpu/examples/config.sh . 
$ cp /path/to/DoublyPeriodicStokes/source/cpu/examples/fcm_example.py .
$ cp /path/to/DoublyPeriodicStokes/source/cpu/examples/Test_Data_For_Rollers.tgz .
$ tar xvzf Test_Data_For_Rollers.tgz
$ source config.sh
$ python fcm_example.py
```

### Organization ###
The C++ library files are in `src` and `include`. The C wrapper is in `wrapper`.

The Python wrappers to the exposed C library, as well as currently available solvers are in `python`.

Examples using the python wrappers and solvers can be found in `examples`.
