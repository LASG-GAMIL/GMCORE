#!/bin/bash

np=10

# Fetch submodules if needed.
git submodule update --init

# Build the model.
cd build
FC=gfortran cmake ..
make
cd ..

# Run all the tests.
if (( $# == 1 )); then
  work_dir=$1
  if [[ ! -d $work_dir ]]; then
    mkdir -p $work_dir
  fi
  cp src/tests/swm/namelist.swm.rh4.360x180 $work_dir
  cp src/tests/swm/namelist.swm.mz.360x180  $work_dir
  cp src/tests/swm/namelist.swm.sg.360x180  $work_dir
  cp src/tests/swm/namelist.swm.jz.360x180  $work_dir
  cp src/tests/swm/namelist.swm.cp.360x180  $work_dir
  cp src/tests/swm/namelist.swm.sw.360x180  $work_dir
  cp build/gmcore_swm_driver.exe $work_dir
  for namelist in $(ls $work_dir/namelist.swm.*); do
    sed -i "s/num_proc_lat = [0-9]*/num_proc_lat = $np/" $namelist
  done
else
  echo '[Error]: You should set a work directory to run the model!'
  exit 1
fi

cd $work_dir

mpirun="mpiexec -np $np"

echo '=========================================================================='
echo 'Rossby-Haurwitz wave test'
time $mpirun ./gmcore_swm_driver.exe namelist.swm.rh4.360x180

echo '=========================================================================='
echo 'Zonal mountain wave'
time $mpirun ./gmcore_swm_driver.exe namelist.swm.mz.360x180

echo '=========================================================================='
echo 'Jet zonal flow test'
time $mpirun ./gmcore_swm_driver.exe namelist.swm.jz.360x180

echo '=========================================================================='
echo 'Cross pole flow test'
time $mpirun ./gmcore_swm_driver.exe namelist.swm.cp.360x180

echo '=========================================================================='
echo 'Steady geostrophic flow test'
time $mpirun ./gmcore_swm_driver.exe namelist.swm.sg.360x180

echo '=========================================================================='
echo 'Shallow water waves test'
time $mpirun ./gmcore_swm_driver.exe namelist.swm.sw.360x180
