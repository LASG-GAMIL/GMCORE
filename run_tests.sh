#!/bin/bash

np=2
if (( $# == 1 )); then
  if [[ ! -d $1 ]]; then
    mkdir -p $1
  fi
  work_dir=$(cd $1 && pwd)
else
  echo '[Error]: You should set a work directory to run the model!'
  exit 1
fi

gmcore_root=$(cd $(dirname $BASH_SOURCE) && pwd)

# Fetch submodules if needed.
cd $gmcore_root
git submodule update --init

# Build the model.
if [[ ! -d $work_dir/build ]]; then
  mkdir $work_dir/build
fi
cd $work_dir/build
FC=gfortran cmake $gmcore_root
make
cd ..

# Run all the tests.
cp $gmcore_root/src/tests/swm/namelist.swm.rh4.360x180 $work_dir
cp $gmcore_root/src/tests/swm/namelist.swm.mz.360x180  $work_dir
cp $gmcore_root/src/tests/swm/namelist.swm.sg.360x180  $work_dir
cp $gmcore_root/src/tests/swm/namelist.swm.jz.360x180  $work_dir
cp $gmcore_root/src/tests/swm/namelist.swm.cp.360x180  $work_dir
cp $gmcore_root/src/tests/swm/namelist.swm.sw.360x180  $work_dir
for namelist in $(ls $work_dir/namelist.swm.*); do
  sed -i '.bak' "s/num_proc_lat = [0-9]*/num_proc_lat = $np/" $namelist
done
rm *.bak

cd $work_dir

mpirun="mpiexec -np $np"

echo '=========================================================================='
echo 'Rossby-Haurwitz wave test'
time $mpirun build/gmcore_swm_driver.exe namelist.swm.rh4.360x180

echo '=========================================================================='
echo 'Zonal mountain wave'
time $mpirun build/gmcore_swm_driver.exe namelist.swm.mz.360x180

echo '=========================================================================='
echo 'Jet zonal flow test'
time $mpirun build/gmcore_swm_driver.exe namelist.swm.jz.360x180

echo '=========================================================================='
echo 'Cross pole flow test'
time $mpirun build/gmcore_swm_driver.exe namelist.swm.cp.360x180

echo '=========================================================================='
echo 'Steady geostrophic flow test'
time $mpirun build/gmcore_swm_driver.exe namelist.swm.sg.360x180

echo '=========================================================================='
echo 'Shallow water waves test'
time $mpirun build/gmcore_swm_driver.exe namelist.swm.sw.360x180
