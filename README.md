# Introduction

Grid-point Multiple-Conservation dynamical cORE

Check barotropic test results [here](https://github.com/gmcore-project/gmcore/wiki/Test-Archive).

# Status

- [ ] Parallelization using MPI:
  - [X] 1D latitudional decomposition (done)
  - [X] 2D decomposition (partially done)
  - [ ] Optimize for X86 (~2021.04)
- [ ] Nesting at middle and low latitudes (~2021.11).
- [ ] Acceleration using GPU (~?).
- [ ] Baroclinic version (~2021.02).
  - [X] Hydrostatic baroclinic version (done)
    - [X] Rossby-Haurwitz wave test
    - [X] Mountain induced wave test
    - [X] Steady state test
    - [X] Baroclinic wave test
    - [X] Held-Suarez test 
  - [X] Nonhydrostatic baroclinic version (~2021.02)
    - [X] X-Z version (done)
    - [X] Quasi-2D mountain wave on reduced sphere (done)
    - [X] Circular mountain wave on reduced sphere (done)
    - [X] Internal gravity wave (done)
- [ ] Advection module (~2021.04).
- [ ] Incorporation with physics parameterisation (2021.05-2021.12).
- [ ] Data assimilation (~?).

# Usage

First make sure you have installed netCDF library, and set `NETCDF_ROOT` environment variable to it. Then clone the repository:
```
$ git clone https://github.com/LASG-GAMIL/GMCORE gmcore
```
There is a Python script `run_tests.py`, which will clone the testbed repository, and run several tests, but it assumes MPI to be installed or SLURM job manager is available:
```
$ ./run_tests.py -w <work_directory> --slurm -q <job_queue> -n <process_number> --ntasks-per-node <n>
```
It will take some time to run the tests. When the tests are finished, cd to `<work_directory>`, and use some visualization tools, such as Panoply, to view the results.

# Authors

- Li Dong <dongli@lasg.iap.ac.cn>
- Jianghao Li

You are welcome to join our team to develop a robust global model!
