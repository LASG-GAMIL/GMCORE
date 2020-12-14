module namelist_mod

  use const_mod
  use time_mod, start_time => start_time_array, end_time => end_time_array

  implicit none

  character(256) :: case_desc = 'N/A'
  character(256) :: case_name = 'N/A'
  character(30) :: test_case = 'N/A'
  character(30) :: history_interval(1) = 'N/A'
  character(30) :: restart_interval    = 'N/A'
  character(30) :: print_interval = '1 hours'
  character(256) :: initial_file = 'N/A'
  character(256) :: restart_file = 'N/A'

  logical :: restart = .false.

  integer num_lon
  integer num_lat
  integer :: num_lev = 1

  logical :: baroclinic = .false.
  logical :: hydrostatic = .true.

  integer :: num_proc_lon(20) = 0
  integer :: num_proc_lat(20) = 0

  character(30) :: tangent_wgt_scheme = 'classic'

  character(30) :: vert_coord_scheme = 'hybrid'
  character(30) :: refer_state_scheme = 'wrf'

  integer :: ke_scheme = 1
  real(r8) :: ke_cell_wgt = 3.0_r8 / 8.0_r8

  integer :: pv_scheme = 1
  logical :: pv_pole_stokes = .true.
  character(8) :: pgf_scheme = 'lin97'
  integer :: coriolis_scheme = 1
  integer :: upwind_order = -1 ! -1, 1, 3

  integer :: fast_cycles = 1
  character(30) :: split_scheme = ''
  character(30) :: time_scheme = 'pc2'

  real(r8) :: coarse_polar_lats_exp = 0
  real(r8) :: coarse_polar_lats_mul = 2

  ! Reduce settings
  integer :: reduce_factors(100) = 0
  integer :: color_proc_lat_reduce(100) = 0
  integer :: num_proc_lon_reduce(100) = 0
  logical :: reduce_pv_directly = .true.

  ! Damping settings
  logical :: use_polar_damp = .false.
  integer :: polar_damp_order = 4
  integer :: polar_damp_cycles = 1
  logical :: use_div_damp = .false.
  integer :: div_damp_order = 2
  integer :: div_damp_j0 = 0
  real(r8) :: div_damp_exp = 0.01_r8
  real(r8) :: div_damp_mul = 0.5_r8
  real(r8) :: div_damp_coef2 = 1.0_r8 / 128.0_r8
  real(r8) :: div_damp_coef4 = 0.01_r8
  logical :: use_vor_damp = .false.
  integer :: vor_damp_order = 2
  real(r8) :: vor_damp_lat0 = 50.0_r8
  real(r8) :: vor_damp_decay = 0.2_r8
  real(r8) :: vor_damp_coef2 = 0.001_r8
  real(r8) :: vor_damp_coef4 = 0.01_r8

  ! Nest settings
  character(30) :: nest_time_scheme   = 'pc2'
  integer       :: nest_max_dom       = 1
  integer       :: nest_parent_id(20) = 1
  real(r8)      :: nest_lon_beg(20) = inf
  real(r8)      :: nest_lon_end(20) = inf
  real(r8)      :: nest_lat_beg(20) = inf
  real(r8)      :: nest_lat_end(20) = inf

  namelist /gmcore_control/     &
    case_name                 , &
    test_case                 , &
    case_desc                 , &
    num_lon                   , &
    num_lat                   , &
    num_lev                   , &
    hydrostatic               , &
    num_proc_lon              , &
    num_proc_lat              , &
    start_time                , &
    end_time                  , &
    dt_in_seconds             , &
    run_hours                 , &
    run_days                  , &
    history_interval          , &
    restart_interval          , &
    print_interval            , &
    initial_file              , &
    restart_file              , &
    restart                   , &
    tangent_wgt_scheme        , &
    vert_coord_scheme         , &
    refer_state_scheme        , &
    ke_scheme                 , &
    ke_cell_wgt               , &
    pv_scheme                 , &
    pv_pole_stokes            , &
    pgf_scheme                , &
    coriolis_scheme           , &
    upwind_order              , &
    split_scheme              , &
    time_scheme               , &
    fast_cycles               , &
    coarse_polar_lats_exp     , &
    coarse_polar_lats_mul     , &
    reduce_factors            , &
    color_proc_lat_reduce     , &
    num_proc_lon_reduce       , &
    reduce_pv_directly        , &
    use_polar_damp            , &
    polar_damp_order          , &
    polar_damp_cycles         , &
    use_div_damp              , &
    div_damp_order            , &
    div_damp_j0               , &
    div_damp_exp              , &
    div_damp_mul              , &
    div_damp_coef2            , &
    div_damp_coef4            , &
    use_vor_damp              , &
    vor_damp_order            , &
    vor_damp_lat0             , &
    vor_damp_decay            , &
    vor_damp_coef2            , &
    nest_time_scheme          , &
    nest_max_dom              , &
    nest_parent_id            , &
    nest_lon_beg              , &
    nest_lon_end              , &
    nest_lat_beg              , &
    nest_lat_end

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gmcore_control)
    close(10)

    if (num_lev > 1) then
      baroclinic = .true.
    end if

  end subroutine parse_namelist

  subroutine print_namelist()

    write(*, *) '====================== Parameters ======================'
    write(*, *) 'num_lon =', num_lon
    write(*, *) 'num_lat =', num_lat
    write(*, *) 'dt_in_seconds =', dt_in_seconds
    write(*, *) 'pgf_scheme =', pgf_scheme
    write(*, *) 'ke_scheme =', ke_scheme
    write(*, *) 'pv_scheme =', pv_scheme
    write(*, *) 'pv_pole_stokes =', pv_pole_stokes
    write(*, *) 'time_scheme =', time_scheme
    write(*, *) 'upwind_order =', upwind_order
    write(*, *) 'use_div_damp =', use_div_damp
    write(*, *) 'div_damp_coef2 =', div_damp_coef2
    write(*, *) 'use_vor_damp =', use_vor_damp
    write(*, *) 'vor_damp_lat0 =', vor_damp_lat0
    write(*, *) 'vor_damp_decay =', vor_damp_decay
    write(*, *) 'vor_damp_coef2 =', vor_damp_coef2
    write(*, *) 'use_polar_damp =', use_polar_damp
    write(*, *) '========================================================'

  end subroutine print_namelist

end module namelist_mod
