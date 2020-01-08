#!/bin/bash

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
  cp build/gmcore_swm_driver.exe $work_dir
else
  echo '[Error]: You should set a work directory to run the model!'
  exit 1
fi

cd $work_dir

echo '=========================================================================='
echo 'Rossby-Haurwitz wave test'
time ./gmcore_swm_driver.exe namelist.swm.rh4.360x180

echo '=========================================================================='
echo 'Zonal mountain wave'
time ./gmcore_swm_driver.exe namelist.swm.mz.360x180

echo '=========================================================================='
echo 'Jet zonal flow test'
time ./gmcore_swm_driver.exe namelist.swm.jz.360x180

echo '=========================================================================='
echo 'Cross pole flow test'
time ./gmcore_swm_driver.exe namelist.swm.cp.360x180

echo '=========================================================================='
echo 'Steady geostrophic flow test'
time ./gmcore_swm_driver.exe namelist.swm.sg.360x180
