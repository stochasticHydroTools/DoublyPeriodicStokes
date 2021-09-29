This folder contains a script for benchmarking the cpu and gpu solvers, along with data needed for the script.
The data archive must be unarchived before executing the script, like 
```shell
tar xvzf Test_Data_For_Rollers.tgz
python3 fcm_multiblob_compare cpu (or gpu)
```
The `fftw_wisdom` folder contains precomputed FFTW plans (using the PATIENT planner) for the cpu solver, 
which are only valid on Courant's `blob` machine at `blob.cims.nyu.edu`. The plans use 10 threads pinned
to each physical core, and each thread is allowed to migrate within the two threads per physical core.
