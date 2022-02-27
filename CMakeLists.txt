cmake_minimum_required(VERSION 3.0)

project(gmcore LANGUAGES Fortran)

option(R4 "Use single precision" OFF)
option(R16 "Use quadruple precision" OFF)
option(ENSURE_ORDER "Ensure computation order in parallel" OFF)

if (R4)
  add_definitions(-DREAL_KIND=4)
  set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=4)
  message(STATUS "Use single precision.")
elseif (R16)
  add_definitions(-DREAL_KIND=16)
  set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=16)
  message(STATUS "Use quadruple precision.")
else ()
  add_definitions(-DREAL_KIND=8)
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
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O0 -check bounds -g -traceback -mkl=sequential -fp-model precise")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O3 -align array64byte -mkl=sequential")
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
set(EXTERNAL_LIBS netcdff netcdf)

if (DEFINED ENV{MKLROOT})
  add_definitions(-DHAS_MKL)
  include_directories("$ENV{MKLROOT}/include")
  link_directories("$ENV{MKLROOT}/lib")
elseif (DEFINED ENV{MKL_INC} AND DEFINED ENV{MKL_LIB})
  add_definitions(-DHAS_MKL)
  include_directories("$ENV{MKL_INC}")
  link_directories("$ENV{MKL_LIB}")
else ()
  message (WARNING "Unable to find MKL library!")
endif ()

if (HAS_MKL)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} lapack64)
  endif ()
endif ()

add_subdirectory(lib/container)
add_subdirectory(lib/datetime)
add_subdirectory(lib/string)
add_subdirectory(lib/fiona)
add_subdirectory(lib/flogger)

set(gmcore_files
  src/utils/const_mod.F90
  src/utils/formula_mod.F90
  src/utils/namelist_mod.F90
  src/utils/time_mod.F90
  src/utils/math_mod.F90
  src/utils/diag_state_mod.F90
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
  src/utils/tridiag_mkl_mod.F90
  src/utils/tridiag_spk_mod.F90
  src/utils/tridiag_mod.F90
  src/utils/debug_mod.F90
  src/core/time_schemes_mod.F90
  src/core/pgf/pgf_mod.F90
  src/core/pgf/pgf_swm_mod.F90
  src/core/pgf/pgf_lin97_mod.F90
  src/core/mesh_mod.F90
  src/core/vert_coord/sigma_coord_mod.F90
  src/core/vert_coord/hybrid_coord_test_mod.F90
  src/core/vert_coord/hybrid_coord_ecmwf_mod.F90
  src/core/vert_coord/hybrid_coord_mars_mod.F90
  src/core/vert_coord/hybrid_coord_mod.F90
  src/core/vert_coord/vert_coord_mod.F90
  src/core/state_mod.F90
  src/core/static_mod.F90
  src/core/tend_mod.F90
  src/core/nh_mod.F90
  src/core/interp_mod.F90
  src/core/upwind_mod.F90
  src/core/weno_mod.F90
  src/core/operators_mod.F90
  src/core/filter_types_mod.F90
  src/core/filter_mod.F90
  src/core/damp/div_damp_mod.F90
  src/core/damp/vor_damp_mod.F90
  src/core/damp/laplace_damp_mod.F90
  src/core/damp/smag_damp_mod.F90
  src/core/damp/damp_mod.F90
  src/core/transport/tracer_mod.F90
  src/core/transport/ffsl_mod.F90
  src/core/gmcore_mod.F90
  src/prepare/prepare_mod.F90
  src/prepare/topo_mod.F90
  src/prepare/latlon_interp_mod.F90
  src/prepare/vert_interp_mod.F90
  src/prepare/era5_reader_mod.F90
  src/prepare/mpas_reader_mod.F90
  src/prepare/waccm_reader_mod.F90
  src/prepare/openmars_reader_mod.F90
  src/prepare/bkg_mod.F90
  src/physics/test_forcing_mod.F90
  src/tests/swm/rossby_haurwitz_wave_test_mod.F90
  src/tests/swm/mountain_zonal_flow_test_mod.F90
  src/tests/swm/jet_zonal_flow_test_mod.F90
  src/tests/swm/steady_geostrophic_flow_test_mod.F90
  src/tests/swm/cross_pole_flow_test_mod.F90
  src/tests/swm/shallow_water_waves_test_mod.F90
  src/tests/swm/vortex_erosion_test_mod.F90
  src/tests/steady_state_test_mod.F90
  src/tests/rossby_haurwitz_wave_3d_test_mod.F90
  src/tests/mountain_wave_test_mod.F90
  src/tests/baroclinic_wave_test_mod.F90
  src/tests/held_suarez_test_mod.F90
  src/tests/steady_state_pgf_test_mod.F90
  src/tests/ksp15_test_mod.F90
  src/tests/dcmip31_test_mod.F90
  src/tests/mars_cold_run_mod.F90
  lib/quadpack.f90
)

add_library(gmcore ${gmcore_files})
target_link_libraries(gmcore fortran_container fortran_datetime fiona flogger ${EXTERNAL_LIBS})

add_executable(gmcore_prepare.exe src/drivers/gmcore_prepare.F90)
target_link_libraries(gmcore_prepare.exe gmcore)

add_executable(gmcore_swm_driver.exe src/drivers/gmcore_swm_driver.F90)
target_link_libraries(gmcore_swm_driver.exe gmcore)

add_executable(gmcore_driver.exe src/drivers/gmcore_driver.F90)
target_link_libraries(gmcore_driver.exe gmcore)
