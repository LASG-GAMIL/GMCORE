module history_mod

  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use allocator_mod
  use time_mod, dt => dt_in_seconds
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use block_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write_state
  public history_write_debug

  ! A-grid velocity
  real(r8), allocatable, dimension(:,:) :: u
  real(r8), allocatable, dimension(:,:) :: v
  real(r8), allocatable, dimension(:,:) :: h
  real(r8), allocatable, dimension(:,:) :: hs

contains

  subroutine history_init()

    character(10) time_value, time_units
    real(r8) seconds

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(history_interval(1), ' ', 1)
    time_units = split_string(history_interval(1), ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid history interval ' // trim(history_interval(1)) // '!')
    end select

    call fiona_init(time_units, start_time_str)

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name))
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%num_full_lon, add_var=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%num_full_lat, add_var=.true.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%num_half_lon, add_var=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%num_half_lat, add_var=.true.)
    call fiona_add_var('h0', 'u'    , long_name='u wind component'         , units='m s-1' , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h0', 'v'    , long_name='v wind component'         , units='m s-1' , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h0', 'h'    , long_name='height'                   , units='m'     , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h0', 'hs'   , long_name='surface height'           , units='m'     , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h0', 'pv'   , long_name='potential vorticity'      , units='s-1'   , dim_names=['ilon', 'ilat', 'time'])
    call fiona_add_var('h0', 'tm'   , long_name='total mass'               , units='m'     , dim_names=['time'])
    call fiona_add_var('h0', 'te'   , long_name='total energy'             , units='m4 s-4', dim_names=['time'], data_type='real(8)')
    call fiona_add_var('h0', 'tpe'  , long_name='total potential enstrophy', units='m2 s-5', dim_names=['time'], data_type='real(8)')

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name))
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time'     , add_var=.true.)
    call fiona_add_dim('h1', 'lon'      , size=global_mesh%num_full_lon, add_var=.true.)
    call fiona_add_dim('h1', 'lat'      , size=global_mesh%num_full_lat, add_var=.true.)
    call fiona_add_dim('h1', 'ilon'     , size=global_mesh%num_half_lon, add_var=.true.)
    call fiona_add_dim('h1', 'ilat'     , size=global_mesh%num_half_lat, add_var=.true.)
    call fiona_add_var('h1', 'qhv'      , long_name='nonliear zonal Coriolis force'                , units='m s-2'  , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'qhu'      , long_name='nonliear meridional Coriolis force'           , units='m s-2'  , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'dpedlon'  , long_name='zonal geopotential energy gradient force'     , units=''       , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'dkedlon'  , long_name='zonal kinetic energy gradient force'          , units=''       , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'dpedlat'  , long_name='meridional geopotential energy gradient force', units=''       , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'dkedlat'  , long_name='meridional kinetic energy gradient force'     , units=''       , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'dmfdlon'  , long_name='zonal mass flux divergence'                   , units=''       , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h1', 'dmfdlat'  , long_name='meridional mass flux divergence'              , units=''       , dim_names=['lon ', 'lat ', 'time'])
    call fiona_add_var('h1', 'us'       , long_name='staggered u wind component'                   , units='m s-1'  , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'vs'       , long_name='staggered v wind component'                   , units='m s-1'  , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'hv'       , long_name='mass on vertices'                             , units='m'      , dim_names=['ilon', 'ilat', 'time'])
    call fiona_add_var('h1', 'mf_lon_n' , long_name='normal mass flux on U grid'                   , units='m2 s-1' , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'mf_lon_t' , long_name='tangent mass flux on U grid'                  , units='m2 s-1' , dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'mf_lat_n' , long_name='normal mass flux on V grid'                   , units='m2 s-1' , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'mf_lat_t' , long_name='tangent mass flux on V grid'                  , units='m2 s-1' , dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'pv_lon'   , long_name='pv on U grid'                                 , units='m-1 s-1', dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'pv_lat'   , long_name='pv on V grid'                                 , units='m-1 s-1', dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'dpv_lon_t', long_name='meridional dpv on U grid'                     , units='m-1 s-1', dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'dpv_lon_n', long_name='zonal dpv on U grid'                          , units='m-1 s-1', dim_names=['ilon', 'lat ', 'time'])
    call fiona_add_var('h1', 'dpv_lat_t', long_name='zonal dpv on V grid'                          , units='m-1 s-1', dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'dpv_lat_n', long_name='meridional dpv on V grid'                     , units='m-1 s-1', dim_names=['lon ', 'ilat', 'time'])
    call fiona_add_var('h1', 'ke'       , long_name='kinetic energy on cell grid'                  , units='',        dim_names=['lon ', 'lat ', 'time'])

    if (.not. allocated(u )) call allocate_array(global_mesh, u )
    if (.not. allocated(v )) call allocate_array(global_mesh, v )
    if (.not. allocated(h )) call allocate_array(global_mesh, h )
    if (.not. allocated(hs)) call allocate_array(global_mesh, hs)

    call time_add_alert('history_write', seconds=seconds)

  end subroutine history_init

  subroutine history_final()

    if (allocated(u )) deallocate(u )
    if (allocated(v )) deallocate(v )
    if (allocated(h )) deallocate(h )
    if (allocated(hs)) deallocate(hs)

  end subroutine history_final

  subroutine history_write_state(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    integer i, j

    mesh => blocks(1)%mesh
    state => blocks(1)%state(itime)
    static => blocks(1)%static

    ! Convert wind from C grid to A grid.
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        v(i,j) = 0.5_r8 * (state%v(i,j) + state%v(i,j+1))
#else
        v(i,j) = 0.5_r8 * (state%v(i,j) + state%v(i,j-1))
#endif
        u(i,j) = 0.5_r8 * (state%u(i,j) + state%u(i-1,j))
        h(i,j) = (static%ghs(i,j) + state%gd(i,j)) / g
        hs(i,j) = static%ghs(i,j) / g
      end do
    end do

    call fiona_start_output('h0', elapsed_seconds, new_file=time_step == 0)
    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    call fiona_output('h0', 'u'   , u                (1:global_mesh%num_full_lon,1:global_mesh%num_full_lat))
    call fiona_output('h0', 'v'   , v                (1:global_mesh%num_full_lon,1:global_mesh%num_full_lat))
    call fiona_output('h0', 'h'   , h                (1:global_mesh%num_full_lon,1:global_mesh%num_full_lat))
    call fiona_output('h0', 'hs'  , hs               (1:global_mesh%num_full_lon,1:global_mesh%num_full_lat))
    call fiona_output('h0', 'pv'  , state%pv         (1:global_mesh%num_half_lon,1:global_mesh%num_half_lat))
    call fiona_output('h0', 'tm'  , state%total_m)
    call fiona_output('h0', 'te'  , state%total_e)
    call fiona_output('h0', 'tpe' , state%total_pe)
    call fiona_end_output('h0')

  end subroutine history_write_state

  subroutine history_write_debug(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    type(state_type), pointer :: state
    type(tend_type), pointer :: tend

    state => blocks(1)%state(itime)
    tend => blocks(1)%tend(itime)

    call fiona_start_output('h1', elapsed_seconds, new_file=time_step == 0)
    call fiona_output('h1', 'lon'      , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h1', 'lat'      , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h1', 'ilon'     , global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h1', 'ilat'     , global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    call fiona_output('h1', 'qhv'      , tend%qhv         (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'qhu'      , tend%qhu         (1:global_mesh%num_full_lon,1:global_mesh%num_half_lat))
    call fiona_output('h1', 'dpedlon'  , tend%dpedlon     (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'dkedlon'  , tend%dkedlon     (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'dpedlat'  , tend%dpedlat     (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'dkedlat'  , tend%dkedlat     (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'dmfdlon'  , tend%dmfdlon     (1:global_mesh%num_full_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'dmfdlat'  , tend%dmfdlat     (1:global_mesh%num_full_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'us'       , state%u          (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'vs'       , state%v          (1:global_mesh%num_full_lon,1:global_mesh%num_half_lat))
    call fiona_output('h1', 'hv'       , state%m_vtx      (1:global_mesh%num_half_lon,1:global_mesh%num_half_lat))
    call fiona_output('h1', 'mf_lon_n' , state%mf_lon_n   (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'mf_lon_t' , state%mf_lon_t   (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'mf_lat_n' , state%mf_lat_n   (1:global_mesh%num_full_lon,1:global_mesh%num_half_lat))
    call fiona_output('h1', 'mf_lat_t' , state%mf_lat_t   (1:global_mesh%num_full_lon,1:global_mesh%num_half_lat))
    call fiona_output('h1', 'pv_lon'   , state%pv_lon     (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'pv_lat'   , state%pv_lat     (1:global_mesh%num_full_lon,1:global_mesh%num_half_lat))
    call fiona_output('h1', 'dpv_lon_t', state%dpv_lon_t  (1:global_mesh%num_full_lon,1:global_mesh%num_half_lat))
    call fiona_output('h1', 'dpv_lon_n', state%dpv_lon_n  (1:global_mesh%num_full_lon,1:global_mesh%num_half_lat))
    call fiona_output('h1', 'dpv_lat_t', state%dpv_lat_t  (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'dpv_lat_n', state%dpv_lat_n  (1:global_mesh%num_half_lon,1:global_mesh%num_full_lat))
    call fiona_output('h1', 'ke'       , state%ke         (1:global_mesh%num_full_lon,1:global_mesh%num_full_lat))
    call fiona_end_output('h1')

  end subroutine history_write_debug

end module history_mod
