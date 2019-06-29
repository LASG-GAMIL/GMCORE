module history_mod

  use fiona
  use flogger
  use const_mod
  use namelist_mod
  use string_mod
  use parallel_mod
  use allocator_mod
  use time_mod, dt => dt_in_seconds
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write

  ! A-grid velocity
  real(real_kind), allocatable :: u(:,:)
  real(real_kind), allocatable :: v(:,:)
  real(real_kind), allocatable :: h(:,:)

  interface history_write
    module procedure history_write_state
  end interface history_write

contains

  subroutine history_init()

    character(10) time_value, time_units
    real(real_kind) seconds

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = string_split(history_interval(1), 1)
    time_units = string_split(history_interval(1), 2)
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
      call log_error('Invalid history invalid ' // trim(history_interval(1)) // '!')
    end select

    call io_init(time_units, start_time_str)

    call io_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name))
    call io_add_att('h0', 'time_step_size', dt)
    call io_add_dim('h0', 'time',  add_var=.true.)
    call io_add_dim('h0', 'lon',   size=mesh%num_full_lon, add_var=.true.)
    call io_add_dim('h0', 'lat',   size=mesh%num_full_lat, add_var=.true.)
    call io_add_dim('h0', 'ilon',  size=mesh%num_half_lon, add_var=.true.)
    call io_add_dim('h0', 'ilat',  size=mesh%num_half_lat, add_var=.true.)
    call io_add_var('h0', 'u',     long_name='u wind component',    units='m s-1', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('h0', 'v',     long_name='v wind component',    units='m s-1', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('h0', 'h',     long_name='height',              units='m'    , dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('h0', 'hs',    long_name='surface height',      units='m'    , dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('h0', 'pv',    long_name='potential vorticity', units='s-1'  , dim_names=['ilon', 'ilat', 'time'])


    call io_add_var('h0', 'us'      , long_name='staggered u wind component' , units='m s-1'  , dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('h0', 'vs'      , long_name='staggered v wind component' , units='m s-1'  , dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('h0', 'hv'      , long_name='mass on vertices'           , units='m'      , dim_names=['ilon', 'ilat', 'time'])
    call io_add_var('h0', 'mf_lon_n', long_name='normal zonal mass flux'     , units='m2 s-1' , dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('h0', 'mf_lat_t', long_name='tangent merdional mass flux', units='m2 s-1' , dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('h0', 'mf_lat_n', long_name='normal merdional mass flux' , units='m2 s-1' , dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('h0', 'mf_lon_t', long_name='tangent zonal mass flux'    , units='m2 s-1' , dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('h0', 'pv_lon'  , long_name='pv on U grid'               , units='m-1 s-1', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('h0', 'pv_lat'  , long_name='pv on V grid'               , units='m-1 s-1', dim_names=['lon ', 'ilat', 'time'])

    if (.not. allocated(u)) call allocate_array(mesh, u)
    if (.not. allocated(v)) call allocate_array(mesh, v)
    if (.not. allocated(h)) call allocate_array(mesh, h)

    call time_add_alert('history_write', seconds=seconds)

  end subroutine history_init

  subroutine history_final()

    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)
    if (allocated(h)) deallocate(h)

  end subroutine history_final

  subroutine history_write_state(static, state)

    type(static_type), intent(in) :: static
    type(state_type ), intent(in) :: state

    integer i, j

    ! Convert wind from C grid to A grid.
    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        v(i,j) = 0.5d0 * (state%v(i,j) + state%v(i,j+1))
#else
        v(i,j) = 0.5d0 * (state%v(i,j) + state%v(i,j-1))
#endif
        u(i,j) = 0.5d0 * (state%u(i,j) + state%u(i-1,j))
        h(i,j) = (state%gd(i,j) + static%ghs(i,j)) / g
      end do
    end do

    call io_start_output('h0', elapsed_seconds, tag=curr_time_str, new_file=.false.)
    call io_output('h0', 'lon' , mesh%full_lon_deg(1:mesh%num_full_lon))
    call io_output('h0', 'lat' , mesh%full_lat_deg(1:mesh%num_full_lat))
    call io_output('h0', 'ilon', mesh%half_lon_deg(1:mesh%num_half_lon))
    call io_output('h0', 'ilat', mesh%half_lat_deg(1:mesh%num_half_lat))
    call io_output('h0', 'u'   , u         (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('h0', 'v'   , v         (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('h0', 'h'   , h         (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('h0', 'hs'  , static%ghs(1:mesh%num_full_lon,1:mesh%num_full_lat) / g)
    call io_output('h0', 'pv'  , state%pv  (1:mesh%num_half_lon,1:mesh%num_half_lat))

    call io_output('h0', 'us'      , state%u              (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call io_output('h0', 'vs'      , state%v              (1:mesh%num_full_lon,1:mesh%num_half_lat))
    call io_output('h0', 'hv'      , state%mass_vertex    (1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('h0', 'mf_lon_n', state%mass_flux_lon_n(1:mesh%num_half_lon,1:mesh%num_full_lat))
    call io_output('h0', 'mf_lat_t', state%mass_flux_lat_t(1:mesh%num_half_lon,1:mesh%num_full_lat))
    call io_output('h0', 'mf_lat_n', state%mass_flux_lat_n(1:mesh%num_full_lon,1:mesh%num_half_lat))
    call io_output('h0', 'mf_lon_t', state%mass_flux_lon_t(1:mesh%num_full_lon,1:mesh%num_half_lat))
    call io_output('h0', 'pv_lon'  , state%pv_lon         (1:mesh%num_half_lon,1:mesh%num_full_lat))
    call io_output('h0', 'pv_lat'  , state%pv_lat         (1:mesh%num_full_lon,1:mesh%num_half_lat))

    call io_end_output('h0')

  end subroutine history_write_state

end module history_mod
