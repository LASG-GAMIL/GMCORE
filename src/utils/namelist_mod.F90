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
  character(256) :: restart_file = 'N/A'

  logical :: restart = .false.

  integer num_lon
  integer num_lat

  integer :: num_proc_lon(20) = 0
  integer :: num_proc_lat(20) = 0

  character(30) :: tangent_wgt_scheme = 'classic'

  integer :: pv_scheme = 2
  logical :: pv_pole_stokes = .true.

  integer :: coriolis_scheme = 1

  integer :: fast_cycles = 1
  character(30) :: split_scheme = ''
  character(30) :: time_scheme = 'pc2'

  real(r8) :: coarse_polar_lats_exp = 0
  real(r8) :: coarse_polar_lats_mul = 2

  integer :: reduce_factors(100) = 0
  integer :: damp_orders(100) = 0
  logical :: adaptive_damp = .false.
  
  ! Nest settings
  character(30) :: nest_time_scheme   = 'pc2'
  integer       :: nest_max_dom       = 1
  integer       :: nest_parent_id(20) = 1
  real(r8)      :: nest_lon_beg(20) = inf
  real(r8)      :: nest_lon_end(20) = inf
  real(r8)      :: nest_lat_beg(20) = inf
  real(r8)      :: nest_lat_end(20) = inf
  
  namelist /gmcore_swm_control/ &
    case_name                 , &
    test_case                 , &
    case_desc                 , &
    num_lon                   , &
    num_lat                   , &
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
    restart_file              , &
    restart                   , &
    tangent_wgt_scheme        , &
    pv_scheme                 , &
    pv_pole_stokes            , &
    coriolis_scheme           , &
    split_scheme              , &
    time_scheme               , &
    fast_cycles               , &
    coarse_polar_lats_exp     , &
    coarse_polar_lats_mul     , &
    reduce_factors            , &
    damp_orders               , &
    adaptive_damp             , &
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
    read(10, nml=gmcore_swm_control)
    close(10)

  end subroutine parse_namelist

end module namelist_mod
