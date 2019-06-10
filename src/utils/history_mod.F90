module history_mod

  use io_mod
  use log_mod
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

    call io_create_dataset('h0', desc=case_desc, file_prefix=case_name // '.h0')
    call io_add_att('h0', 'time_step_size', dt)
    call io_add_dim('h0', 'time', add_var=.true.)
    call io_add_dim('h0', 'lon',  size=mesh%num_full_lon, add_var=.true.)
    call io_add_dim('h0', 'lat',  size=mesh%num_full_lat, add_var=.true.)
    call io_add_dim('h0', 'ilon', size=mesh%num_half_lon, add_var=.true.)
    call io_add_dim('h0', 'ilat', size=mesh%num_half_lat, add_var=.true.)
    call io_add_var('h0', 'u',    long_name='u wind component',     units='m s-1',  dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('h0', 'v',    long_name='v wind component',     units='m s-1',  dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('h0', 'h',    long_name='height',               units='m',      dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('h0', 'hs',   long_name='surface height',       units='m',      dim_names=['lon ', 'lat ', 'time'])

    if (.not. allocated(u)) call allocate_array(mesh, u)
    if (.not. allocated(v)) call allocate_array(mesh, v)
    if (.not. allocated(h)) call allocate_array(mesh, h)

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
    case ('steps')
      seconds = seconds * dt
    case default
      call log_error('Invalid history invalid ' // trim(history_interval(1)) // '!')
    end select


    call time_add_alert('history_write', seconds=seconds)

  end subroutine history_init

  subroutine history_final()

    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)
    if (allocated(h)) deallocate(h)

  end subroutine history_final

  subroutine history_write_state(state, static)

    type(state_type), intent(in) :: state
    type(static_type), intent(in) :: static

    integer i, j

    ! Convert wind from C grid to A grid.
    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        u(i,j) = 0.5 * (state%u(i,j) + state%u(i-1,j))
        v(i,j) = 0.5 * (state%v(i,j) + state%v(i,j-1))
        h(i,j) = (state%gd(i,j) + static%ghs(i,j)) / g
      end do
    end do

    call io_start_output('h0')
    call io_output('h0', 'lon',  mesh%full_lon_deg(:))
    call io_output('h0', 'lat',  mesh%full_lat_deg(:))
    call io_output('h0', 'ilon', mesh%half_lon_deg(:))
    call io_output('h0', 'ilat', mesh%half_lat_deg(:))
    call io_output('h0', 'u',    u(1:mesh%num_half_lon,1:mesh%num_full_lat))
    call io_output('h0', 'v',    v(1:mesh%num_full_lon,1:mesh%num_half_lat))
    call io_output('h0', 'h',    h(1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('h0', 'hs',   static%ghs(1:mesh%num_full_lon,1:mesh%num_full_lat) / g)
    call io_end_output('h0')

  end subroutine history_write_state

end module history_mod
