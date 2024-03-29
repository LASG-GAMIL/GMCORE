module namelist_mod

  use string
  use flogger
  use const_mod, only: r8, const_init
  use time_mod, start_time => start_time_array, end_time => end_time_array

  implicit none

  character(30)   :: planet               = 'earth'

  real(8)         :: dt_dyn               = 0
  real(8)         :: dt_phys              = 0
  real(8)         :: dt_adv               = 0

  character(256)  :: case_desc            = 'N/A'
  character(256)  :: case_name            = 'N/A'
  character(30 )  :: test_case            = 'N/A'
  character(30)   :: initial_time         = '2000-01-01T00:00:00'
  character(30 )  :: history_interval(1)  = 'N/A'
  character(30 )  :: restart_interval     = 'N/A'
  character(30 )  :: print_interval       = '1 hours'
  character(256)  :: initial_file         = 'N/A'
  character(256)  :: restart_file         = 'N/A'
  character(256)  :: topo_file            = 'N/A'
  character(30 )  :: topo_type            = 'etopo1' ! etopo1, mola32
  character(256)  :: bkg_file             = 'N/A'
  character(30 )  :: bkg_type             = 'era5'

  logical         :: restart              = .false.

  integer         :: num_lon
  integer         :: num_lat
  integer         :: num_lev              = 1

  logical         :: baroclinic           = .false.
  logical         :: hydrostatic          = .true.
  logical         :: nonhydrostatic       = .false.
  logical         :: advection            = .false.

  integer         :: num_proc_lon(20)     = 0
  integer         :: num_proc_lat(20)     = 0
  integer         :: lon_halo_width       = 2


  character(30)   :: tangent_wgt_scheme   = 'classic'

  real(r8)        :: implicit_w_wgt       = 0.5_r8

  character(30)   :: vert_coord_scheme    = 'hybrid'
  character(30)   :: vert_coord_template  = 'N/A'
  character(30)   :: refer_state_scheme   = 'wrf'
  real(r8)        :: ptop                 = 2.194e2_r8
  real(r8)        :: hybrid_coord_p0      = 1.0e5_r8

  integer         :: ke_scheme            = 2
  real(r8)        :: ke_cell_wgt          = 0.5_r8

  character(30)   :: pv_scheme            = 'upwind' ! midpoint, upwind, ffsl
  logical         :: pv_pole_stokes       = .true.
  integer         :: upwind_order_pv      = 3
  real(r8)        :: upwind_wgt_pv        = 1

  character(8)    :: pgf_scheme           = 'lin97'
  integer         :: coriolis_scheme      = 1

  character(8)    :: adv_scheme           = 'ffsl'
  character(8)    :: limiter_type         = 'mono'
  character(8)    :: ffsl_flux_type       = 'ppm'
  character(8)    :: tvd_limiter_type     = 'van_leer'

  character(8)    :: zonal_tridiag_solver = 'spk' ! mkl, spk

  integer         :: weno_order           = -1 ! -1, 3
  integer         :: upwind_order         = 3  ! -1, 1, 3
  real(r8)        :: upwind_wgt           = 1.0_r8
  real(r8)        :: upwind_wgt_pt        = 0.25_r8

  integer         :: vert_weno_order      = -1 ! -1, 3
  integer         :: vert_upwind_order    = 3  ! -1, 1, 3
  real(r8)        :: vert_upwind_wgt      = 1.0_r8

  character(30)   :: time_scheme          = 'wrfrk3'

  real(r8)        :: coarse_pole_mul      = 0
  real(r8)        :: coarse_pole_decay    = 100.0

  ! Filter settings
  real(r8)        :: max_wave_speed       = 300
  real(r8)        :: max_cfl              = 0.5
  real(r8)        :: filter_coef_a        = 1.0
  real(r8)        :: filter_coef_b        = 0.4
  real(r8)        :: filter_coef_c        = 0.2
  real(r8)        :: filter_coef_d        = 0.375
  real(r8)        :: filter_coef_e        = 0.16

  ! Damping settings
  logical         :: use_topo_smooth      = .false.
  integer         :: topo_smooth_cycles   = 1
  logical         :: use_div_damp         = .false.
  integer         :: div_damp_order       = 2
  integer         :: div_damp_k0          = 3
  real(r8)        :: div_damp_imp_lat0    = 90
  real(r8)        :: div_damp_top         = 3.0_r8
  real(r8)        :: div_damp_pole        = 0.0_r8
  real(r8)        :: div_damp_lat0        = 90
  real(r8)        :: div_damp_coef2       = 1.0_r8 / 128.0_r8
  real(r8)        :: div_damp_coef4       = 0.01_r8
  real(r8)        :: rayleigh_damp_w_coef = 0.2
  real(r8)        :: rayleigh_damp_top    = 10.0d3 ! m
  logical         :: use_smag_damp        = .false.
  real(r8)        :: smag_damp_coef       = 0.1

  ! Output settings
  character(8)    :: output_i0_dtype      = 'r8'
  logical         :: output_h0            = .true.
  character(8)    :: output_h0_dtype      = 'r4'
  logical         :: output_h1            = .false.
  logical         :: split_h0             = .false.
  character(30)   :: output_h0_new_file   = ''
  character(8)    :: output_h0_vars(100)  = ''
  integer         :: output_group_size    = 0

  namelist /gmcore_control/     &
    planet                    , &
    case_name                 , &
    test_case                 , &
    case_desc                 , &
    num_lon                   , &
    num_lat                   , &
    num_lev                   , &
    nonhydrostatic            , &
    advection                 , &
    num_proc_lon              , &
    num_proc_lat              , &
    lon_halo_width            , &
    initial_time              , &
    start_time                , &
    end_time                  , &
    dt_dyn                    , &
    dt_phys                   , &
    dt_adv                    , &
    run_years                 , &
    run_hours                 , &
    run_days                  , &
    history_interval          , &
    restart_interval          , &
    print_interval            , &
    initial_file              , &
    output_i0_dtype           , &
    restart_file              , &
    restart                   , &
    topo_file                 , &
    topo_type                 , &
    bkg_file                  , &
    bkg_type                  , &
    tangent_wgt_scheme        , &
    implicit_w_wgt            , &
    vert_coord_scheme         , &
    vert_coord_template       , &
    refer_state_scheme        , &
    ptop                      , &
    hybrid_coord_p0           , &
    ke_scheme                 , &
    ke_cell_wgt               , &
    pv_scheme                 , &
    pv_pole_stokes            , &
    upwind_order_pv           , &
    upwind_wgt_pv             , &
    pgf_scheme                , &
    coriolis_scheme           , &
    adv_scheme                , &
    limiter_type              , &
    ffsl_flux_type            , &
    tvd_limiter_type          , &
    zonal_tridiag_solver      , &
    weno_order                , &
    upwind_order              , &
    upwind_wgt                , &
    upwind_wgt_pt             , &
    vert_weno_order           , &
    vert_upwind_order         , &
    vert_upwind_wgt           , &
    time_scheme               , &
    max_wave_speed            , &
    max_cfl                   , &
    filter_coef_a             , &
    filter_coef_b             , &
    filter_coef_c             , &
    filter_coef_d             , &
    filter_coef_e             , &
    coarse_pole_mul           , &
    coarse_pole_decay         , &
    use_topo_smooth           , &
    topo_smooth_cycles        , &
    use_div_damp              , &
    div_damp_order            , &
    div_damp_imp_lat0         , &
    div_damp_k0               , &
    div_damp_top              , &
    div_damp_pole             , &
    div_damp_lat0             , &
    div_damp_coef2            , &
    div_damp_coef4            , &
    rayleigh_damp_w_coef      , &
    rayleigh_damp_top         , &
    use_smag_damp             , &
    smag_damp_coef            , &
    output_h0                 , &
    output_h0_dtype           , &
    output_h1                 , &
    split_h0                  , &
    output_h0_new_file        , &
    output_h0_vars            , &
    output_group_size

contains

  subroutine parse_namelist(file_path)

    character(*), intent(in) :: file_path

    open(10, file=file_path, status='old')
    read(10, nml=gmcore_control)
    close(10)

    ! Here we set baroclinic according to levels.
    baroclinic = num_lev > 1
    if (.not. baroclinic) then
      hydrostatic    = .false.
      nonhydrostatic = .false.
      ke_scheme      = 1
    else
      hydrostatic = .not. nonhydrostatic
    end if

    if (advection) then
      hydrostatic    = .false.
      baroclinic     = .false.
      nonhydrostatic = .false.
    end if

    call const_init(planet)

    if (dt_dyn == 0) dt_dyn = dt_adv
    if (dt_adv == 0) dt_adv = dt_dyn

  end subroutine parse_namelist

  subroutine print_namelist()

      write(*, *) '=================== GMCORE Parameters ==================='
      write(*, *) 'case_name           = ', trim(case_name)
      write(*, *) 'num_lon             = ', to_str(num_lon)
      write(*, *) 'num_lat             = ', to_str(num_lat)
      write(*, *) 'num_lev             = ', to_str(num_lev)
    if (coarse_pole_mul /= 0) then
      write(*, *) 'coarse_pole_mul     = ', to_str(coarse_pole_mul, 3)
      write(*, *) 'coarse_pole_decay   = ', to_str(coarse_pole_decay, 3)
    end if
      write(*, *) 'hydrostatic         = ', to_str(hydrostatic)
      write(*, *) 'nonhydrostatic      = ', to_str(nonhydrostatic)
      write(*, *) 'vert_coord_scheme   = ', trim(vert_coord_scheme)
      write(*, *) 'vert_coord_template = ', trim(vert_coord_template)
      write(*, *) 'ptop                = ', to_str(ptop, 4)
      write(*, *) 'hybrid_coord_p0     = ', hybrid_coord_p0
      write(*, *) 'dt_dyn              = ', to_str(dt_dyn, 2)
      write(*, *) 'dt_adv              = ', to_str(dt_adv, 2)
      write(*, *) 'max_wave_speed      = ', max_wave_speed
      write(*, *) 'max_cfl             = ', max_cfl
      write(*, *) 'filter_coef_a       = ', filter_coef_a
      write(*, *) 'filter_coef_b       = ', filter_coef_b
      write(*, *) 'filter_coef_c       = ', filter_coef_c
      write(*, *) 'filter_coef_d       = ', filter_coef_d
      write(*, *) 'filter_coef_e       = ', filter_coef_e
      write(*, *) 'pgf_scheme          = ', trim(pgf_scheme)
      write(*, *) 'adv_scheme          = ', trim(adv_scheme)
      write(*, *) 'limiter_type        = ', trim(limiter_type)
    if (adv_scheme == 'ffsl') then
      write(*, *) 'ffsl_flux_type      = ', trim(ffsl_flux_type)
    end if
      write(*, *) 'ke_scheme           = ', to_str(ke_scheme)
    if (ke_scheme == 2) then
      write(*, *) 'ke_cell_wgt         = ', to_str(ke_cell_wgt, 2)
    end if
      write(*, *) 'pv_scheme           = ', trim(pv_scheme)
      write(*, *) 'pv_pole_stokes      = ', to_str(pv_pole_stokes)
    if (pv_scheme == 'upwind') then
      write(*, *) 'upwind_order_pv     = ', to_str(upwind_order_pv)
      write(*, *) 'upwind_wgt_pv       = ', to_str(upwind_wgt_pv, 2)
    end if
      write(*, *) 'time_scheme         = ', trim(time_scheme)
      write(*, *) 'upwind_order        = ', to_str(upwind_order)
      write(*, *) 'use_topo_smooth     = ', to_str(use_topo_smooth)
    if (use_topo_smooth) then
      write(*, *) 'topo_smooth_cycles  = ', to_str(topo_smooth_cycles)
    end if
      write(*, *) 'use_div_damp        = ', to_str(use_div_damp)
    if (use_div_damp) then
      write(*, *) 'div_damp_coef2      = ', div_damp_coef2
      write(*, *) 'div_damp_top        = ', to_str(div_damp_top, 3)
      write(*, *) 'div_damp_pole       = ', div_damp_pole
      write(*, *) 'div_damp_lat0       = ', div_damp_lat0
    end if
    if (nonhydrostatic) then
      write(*, *) 'implicit_w_wgt      = ', to_str(implicit_w_wgt, 3)
      write(*, *) 'rayleigh_damp_w_coef= ', to_str(rayleigh_damp_w_coef, 2)
      write(*, *) 'rayleigh_damp_top   = ', to_str(rayleigh_damp_top   , 2)
    end if
      write(*, *) 'use_smag_damp       = ', to_str(use_smag_damp)
    if (use_smag_damp) then
      write(*, *) 'smag_damp_coef      = ', to_str(smag_damp_coef, 1)
    end if
      write(*, *) '========================================================='

  end subroutine print_namelist

end module namelist_mod
