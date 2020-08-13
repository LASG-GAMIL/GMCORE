cmake_minimum_required(VERSION 3.0)

project(gmcore LANGUAGES Fortran)

option(V_POLE "Stagger V on the two Poles" OFF)
option(SIGMA "Use sigma vertical coordinate" OFF)
option(R4 "Use single precision" OFF)
option(R16 "Use quadruple precision" OFF)

if (V_POLE)
  add_definitions(-DV_POLE)
endif ()
if (SIGMA)
  add_definitions(-DSIGMA)
  message(STATUS "Use sigma vertical coordinate.")
endif ()
if (R4)
  add_definitions(-DREAL_KIND=4)
  message(STATUS "Use single precision.")
elseif (R16)
  add_definitions(-DREAL_KIND=16)
  message(STATUS "Use quadruple precision.")
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
endif ()

find_package(MPI)
include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})

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

add_subdirectory(lib/container)
add_subdirectory(lib/datetime)
add_subdirectory(lib/string)
add_subdirectory(lib/fiona)
add_subdirectory(lib/flogger)

add_custom_command(
  OUTPUT ${PROJECT_SOURCE_DIR}/build/reduce_mod.F90
  COMMAND cpp -P -C ${PROJECT_SOURCE_DIR}/src/core/reduce_mod.F90 > ${PROJECT_SOURCE_DIR}/build/reduce_mod.F90
  MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/src/core/reduce_mod.F90
  VERBATIM
)

set(source_files
  src/utils/const_mod.F90
  src/utils/formula_mod.F90
  src/utils/namelist_mod.F90
  src/utils/time_mod.F90
  src/utils/math_mod.F90
  src/utils/history_mod.F90
  src/utils/restart_mod.F90
  src/utils/time_mod.F90
  src/utils/sphere_geometry_mod.F90
  src/utils/block_mod.F90
  src/utils/process_mod.F90
  src/utils/halo_mod.F90
  src/utils/allocator_mod.F90
  src/utils/parallel_types_mod.F90
  src/utils/parallel_mod.F90
  src/utils/debug_mod.F90
  src/core/mesh_mod.F90
  src/core/sigma_coord_mod.F90
  src/core/hybrid_coord_test_mod.F90
  src/core/hybrid_coord_ecmwf_mod.F90
  src/core/hybrid_coord_mod.F90
  src/core/vert_coord_mod.F90
  src/core/state_mod.F90
  src/core/static_mod.F90
  src/core/tend_mod.F90
  src/core/pv_mod.F90
  src/core/ke_mod.F90
  src/core/interp_mod.F90
  src/core/operators_mod.F90
  src/core/reduced_types_mod.F90
  build/reduce_mod.F90
  src/core/nest_mod.F90
  src/core/damp_mod.F90
  src/core/gmcore_mod.F90
  src/tests/swm/rossby_haurwitz_wave_test_mod.F90
  src/tests/swm/mountain_zonal_flow_test_mod.F90
  src/tests/swm/jet_zonal_flow_test_mod.F90
  src/tests/swm/steady_geostrophic_flow_test_mod.F90
  src/tests/swm/cross_pole_flow_test_mod.F90
  src/tests/swm/shallow_water_waves_test_mod.F90
  src/tests/steady_state_test_mod.F90
  src/tests/rossby_haurwitz_wave_3d_test_mod.F90
  src/tests/mountain_wave_test_mod.F90
  src/tests/held_suarez_test_mod.F90
  lib/quadpack.f90
)

add_library(gmcore ${source_files})
target_link_libraries(gmcore netcdff fortran_container fortran_datetime fiona flogger ${MPI_Fortran_LIBRARIES})

add_executable(gmcore_swm_driver.exe src/drivers/gmcore_swm_driver.F90)
target_link_libraries(gmcore_swm_driver.exe gmcore)

add_executable(gmcore_driver.exe src/drivers/gmcore_driver.F90)
target_link_libraries(gmcore_driver.exe gmcore)
