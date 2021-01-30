cmake_minimum_required(VERSION 3.0)

project(gmcore LANGUAGES Fortran)

option(V_POLE "Stagger V on the two Poles" OFF)
option(R4 "Use single precision" OFF)
option(R16 "Use quadruple precision" OFF)
option(ENSURE_ORDER "Ensure computation order in parallel" OFF)

if (V_POLE)
  add_definitions(-DV_POLE)
	set(CPPFLAGS ${CPPFLAGS} -DV_POLE)
	message(STATUS "Place v at Poles.")
endif ()
if (R4)
  add_definitions(-DREAL_KIND=4)
	set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=4)
  message(STATUS "Use single precision.")
elseif (R16)
  add_definitions(-DREAL_KIND=16)
	set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=16)
  message(STATUS "Use quadruple precision.")
endif ()
if (ENSURE_ORDER)
  add_definitions(-DENSURE_ORDER)
  set(CPPFLAGS ${CPPFLAGS} -DENSURE_ORDER)
  message(STATUS "Ensure computation order in parallel.")
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffpe-trap=invalid,zero,overflow")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -fcheck=all -g -fbacktrace")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast")
  endif ()
  # FIXME: OpenMP causes sigmentation faults.
  # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O0 -check bounds -g -traceback -mkl")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O3 -align array64byte -mkl")
  endif ()
endif ()

find_package(MPI)

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

if (DEFINED ENV{MKLROOT})
  include_directories("$ENV{MKLROOT}/include")
  link_directories("$ENV{MKLROOT}/lib")
elseif (DEFINED ENV{MKL_INC} AND DEFINED ENV{MKL_LIB})
  include_directories("$ENV{MKL_INC}")
  link_directories("$ENV{MKL_LIB}")
else ()
  message (FATAL_ERROR "Unable to find MKL library!")
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(MKL_LIBRARIES "lapack64")
endif ()

add_subdirectory(lib/container)
add_subdirectory(lib/datetime)
add_subdirectory(lib/string)
add_subdirectory(lib/fiona)
add_subdirectory(lib/flogger)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/reduce_mod.F90
  COMMAND cpp ${CPPFLAGS} -P -C ${PROJECT_SOURCE_DIR}/src/core/reduce_mod.F90 > ${CMAKE_CURRENT_BINARY_DIR}/reduce_mod.F90
  MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/src/core/reduce_mod.F90
  VERBATIM
)

set(gmcore_files
  src/utils/const_mod.F90
  src/utils/formula_mod.F90
  src/utils/namelist_mod.F90
  src/utils/time_mod.F90
  src/utils/math_mod.F90
  src/utils/history_mod.F90
  src/utils/initial_mod.F90
  src/utils/restart_mod.F90
  src/utils/time_mod.F90
  src/utils/sphere_geometry_mod.F90
  src/utils/block_mod.F90
  src/utils/process_mod.F90
  src/utils/halo_mod.F90
  src/utils/allocator_mod.F90
  src/utils/parallel_types_mod.F90
  src/utils/parallel_mod.F90
  src/utils/tridiag_mod.F90
  src/utils/debug_mod.F90
  src/core/time_schemes_mod.F90
  src/core/pgf/pgf_mod.F90
  src/core/pgf/pgf_swm_mod.F90
  src/core/pgf/pgf_sb81_mod.F90
  src/core/pgf/pgf_lin97_mod.F90
  src/core/pgf/pgf_dflx_mod.F90
  src/core/mesh_mod.F90
  src/core/vert_coord/sigma_coord_mod.F90
  src/core/vert_coord/hybrid_coord_test_mod.F90
  src/core/vert_coord/hybrid_coord_ecmwf_mod.F90
  src/core/vert_coord/hybrid_coord_mod.F90
  src/core/vert_coord/vert_coord_mod.F90
  src/core/refer_state/refer_state_types_mod.F90
  src/core/refer_state/refer_state_mod.F90
  src/core/refer_state/refer_state_wrf_mod.F90
  src/core/state_mod.F90
  src/core/static_mod.F90
  src/core/tend_mod.F90
  src/core/pv_mod.F90
  src/core/ke_mod.F90
  src/core/nh_mod.F90
  src/core/interp_mod.F90
  src/core/operators_mod.F90
  src/core/reduced_types_mod.F90
  ${CMAKE_CURRENT_BINARY_DIR}/reduce_mod.F90
  src/core/nest_mod.F90
  src/core/damp/zonal_damp_mod.F90
  src/core/damp/div_damp_mod.F90
  src/core/damp/vor_damp_mod.F90
  src/core/damp/rayleigh_damp_mod.F90
  src/core/damp/smag_damp_mod.F90
  src/core/damp/damp_mod.F90
  src/core/gmcore_mod.F90
  src/physics/test_forcing_mod.F90
  src/tests/swm/rossby_haurwitz_wave_test_mod.F90
  src/tests/swm/mountain_zonal_flow_test_mod.F90
  src/tests/swm/jet_zonal_flow_test_mod.F90
  src/tests/swm/steady_geostrophic_flow_test_mod.F90
  src/tests/swm/cross_pole_flow_test_mod.F90
  src/tests/swm/shallow_water_waves_test_mod.F90
  src/tests/steady_state_test_mod.F90
  src/tests/rossby_haurwitz_wave_3d_test_mod.F90
  src/tests/mountain_wave_test_mod.F90
  src/tests/baroclinic_wave_test_mod.F90
  src/tests/held_suarez_test_mod.F90
  src/tests/steady_state_pgf_test_mod.F90
  src/tests/ksp15_test_mod.F90
  lib/quadpack.f90
)

add_library(gmcore ${gmcore_files})
target_link_libraries(gmcore netcdff fortran_container fortran_datetime fiona flogger ${MKL_LIBRARIES})

set(gmcore_prepare_files
  src/prepare/topo_mod.F90
  src/prepare/latlon_interp_mod.F90
  src/prepare/vert_interp_mod.F90
  src/prepare/era5_reader_mod.F90
  src/prepare/mpas_reader_mod.F90
  src/prepare/bkg_mod.F90
  src/prepare/gmcore_prepare.F90
)

add_executable(gmcore_prepare.exe ${gmcore_prepare_files})
target_link_libraries(gmcore_prepare.exe gmcore)

add_executable(gmcore_swm_driver.exe src/drivers/gmcore_swm_driver.F90)
target_link_libraries(gmcore_swm_driver.exe gmcore)

add_executable(gmcore_driver.exe src/drivers/gmcore_driver.F90)
target_link_libraries(gmcore_driver.exe gmcore)
