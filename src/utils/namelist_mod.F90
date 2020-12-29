module namelist_mod

  use string
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
  logical :: nonhydrostatic = .false.

  integer :: num_proc_lon(20) = 0
  integer :: num_proc_lat(20) = 0

  character(30) :: tangent_wgt_scheme = 'classic'

  character(30) :: vert_coord_scheme = 'hybrid'
  character(30) :: vert_coord_template = 'N/A'
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

  real(r8) :: coarse_polar_lat0 = 0
  real(r8) :: coarse_polar_decay = 0.2

  ! Reduce settings
  integer :: reduce_factors(100) = 0
  logical :: reduce_pv_directly = .true.
  logical :: do_reduce_ke = .true.

  ! Damping settings
  logical :: use_polar_damp = .false.
  integer :: polar_damp_order = 4
  integer :: polar_damp_cycles = 1
  logical :: use_div_damp = .false.
  integer :: div_damp_order = 2
  integer :: div_damp_j0 = 0
  real(r8) :: div_damp_upper = 2.0_r8
  real(r8) :: div_damp_polar = 0.5_r8
  real(r8) :: div_damp_exp = 0.01_r8
  real(r8) :: div_damp_coef2 = 1.0_r8 / 128.0_r8
  real(r8) :: div_damp_coef4 = 0.01_r8
  logical :: use_vor_damp = .false.
  integer :: vor_damp_order = 2
  real(r8) :: vor_damp_lat0 = 50.0_r8
  real(r8) :: vor_damp_decay = 0.2_r8
  real(r8) :: vor_damp_coef2 = 0.001_r8
  real(r8) :: vor_damp_coef4 = 0.01_r8
  logical :: use_rayleigh_damp = .false.
  logical :: use_smag_damp = .false.
  real(r8) :: smag_damp_coef = 0.2

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
    vert_coord_template       , &
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
    coarse_polar_lat0         , &
    coarse_polar_decay        , &
    reduce_factors            , &
    reduce_pv_directly        , &
    do_reduce_ke              , &
    use_polar_damp            , &
    polar_damp_order          , &
    polar_damp_cycles         , &
    use_div_damp              , &
    div_damp_order            , &
    div_damp_polar            , &
    div_damp_upper            , &
    div_damp_j0               , &
    div_damp_exp              , &
    div_damp_coef2            , &
    div_damp_coef4            , &
    use_vor_damp              , &
    vor_damp_order            , &
    vor_damp_lat0             , &
    vor_damp_decay            , &
    vor_damp_coef2            , &
    use_rayleigh_damp         , &
    use_smag_damp             , &
    smag_damp_coef            , &
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

    nonhydrostatic = .not. hydrostatic

  end subroutine parse_namelist

  subroutine print_namelist()

    write(*, *) '=================== GMCORE Parameters ==================='
    write(*, *) 'num_lon             = ', to_string(num_lon)
    write(*, *) 'num_lat             = ', to_string(num_lat)
    write(*, *) 'num_lev             = ', to_string(num_lev)
    write(*, *) 'vert_coord_scheme   = ', trim(vert_coord_scheme)
    write(*, *) 'vert_coord_template = ', trim(vert_coord_template)
    write(*, *) 'dt_in_seconds       = ', to_string(int(dt_in_seconds))
    write(*, *) 'pgf_scheme          = ', trim(pgf_scheme)
    write(*, *) 'ke_scheme           = ', to_string(ke_scheme)
    write(*, *) 'pv_scheme           = ', to_string(pv_scheme)
    write(*, *) 'pv_pole_stokes      = ', to_string(pv_pole_stokes)
    write(*, *) 'time_scheme         = ', trim(time_scheme)
    write(*, *) 'upwind_order        = ', to_string(upwind_order)
    write(*, *) 'use_div_damp        = ', to_string(use_div_damp)
    write(*, *) 'div_damp_coef2      = ', to_string(div_damp_coef2, 3)
    write(*, *) 'use_vor_damp        = ', to_string(use_vor_damp)
    write(*, *) 'vor_damp_lat0       = ', to_string(vor_damp_lat0, 1)
    write(*, *) 'vor_damp_decay      = ', to_string(vor_damp_decay, 1)
    write(*, *) 'vor_damp_coef2      = ', to_string(vor_damp_coef2, 3)
    write(*, *) 'use_polar_damp      = ', to_string(use_polar_damp)
    write(*, *) 'use_rayleigh_damp   = ', to_string(use_rayleigh_damp)
    write(*, *) 'use_smag_damp       = ', to_string(use_smag_damp)
    write(*, *) 'smag_damp_coef      = ', to_string(smag_damp_coef, 1)
    write(*, *) '========================================================='

  end subroutine print_namelist

end module namelist_mod
