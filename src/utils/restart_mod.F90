module restart_mod

  use fiona
  use string
  use flogger
  use datetime
  use const_mod
  use namelist_mod
  use time_mod
  use block_mod
  use parallel_mod

  implicit none

  private

  public restart_init
  public restart_write
  public restart_read
  public restart_final

contains

  subroutine restart_init()

    character(10) time_value, time_units
    real(8) seconds

    if (restart_interval == 'N/A') then
      if (is_root_proc()) call log_warning('Parameter restart_interval is not set, so no restart file outputted.')
      return
    end if
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(restart_interval, ' ', 1)
    time_units = split_string(restart_interval, ' ', 2)
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
      call log_error('Invalid restart interval ' // trim(restart_interval) // '!')
    end select

    call time_add_alert('restart_write', seconds=seconds)

  end subroutine restart_init

  subroutine restart_write(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)
    character(4) lon_dims(4), lat_dims(4), cell_dims(4)
    character(4) cell_dims_2d(3)

    if (baroclinic) then
       lon_dims(1) = 'ilon';  lon_dims(2) =  'lat';  lon_dims(3) =  'lev';  lon_dims(4) = 'time'
       lat_dims(1) =  'lon';  lat_dims(2) = 'ilat';  lat_dims(3) =  'lev';  lat_dims(4) = 'time'
      cell_dims(1) =  'lon'; cell_dims(2) =  'lat'; cell_dims(3) =  'lev'; cell_dims(4) = 'time'
    else

    end if
    cell_dims_2d(1) = 'lon'; cell_dims_2d(2) = 'lat'; cell_dims_2d(3) = 'time'

    call fiona_create_dataset('r0', desc=case_desc, file_prefix=trim(case_name) // '.' // trim(curr_time_str), mpi_comm=proc%comm)
    call fiona_add_att('r0', 'time_step_size', dt_in_seconds)
    call fiona_add_att('r0', 'restart_interval', restart_interval)
    call fiona_add_dim('r0', 'time', add_var=.true.)
    call fiona_add_dim('r0', 'lon' , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'lat' , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'ilon', size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('r0', 'ilat', size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    if (baroclinic) then
      call fiona_add_dim('r0', 'lev' , size=global_mesh%num_full_lev, add_var=.true.)
      call fiona_add_dim('r0', 'ilev', size=global_mesh%num_half_lev, add_var=.true.)
    end if
    call fiona_add_var('r0', 'u'   , long_name='u wind component'           , units='m s-1' , dim_names=lon_dims    , data_type='real(8)')
    call fiona_add_var('r0', 'v'   , long_name='v wind component'           , units='m s-1' , dim_names=lat_dims    , data_type='real(8)')
    call fiona_add_var('r0', 'gzs' , long_name='surface geopotential height', units='m2 s-2', dim_names=cell_dims_2d, data_type='real(8)')
    if (baroclinic) then
      call fiona_add_var('r0', 'phs', long_name='hydrostatic surface pressure', units='Pa'    , dim_names=cell_dims_2d, data_type='real(8)')
      call fiona_add_var('r0', 'pt' , long_name='potential temperature'       , units='K'     , dim_names=cell_dims   , data_type='real(8)')
    else
      call fiona_add_var('r0', 'gz' , long_name='geopotential height'         , units='m2 s-2', dim_names=cell_dims   , data_type='real(8)')
    end if

    call fiona_start_output('r0', elapsed_seconds, new_file=.true.)
    call fiona_output('r0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('r0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('r0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('r0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    do iblk = 1, size(blocks)
      mesh   => blocks(iblk)%mesh
      state  => blocks(iblk)%state(itime)
      static => blocks(iblk)%static

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_output('r0', 'gzs', static%gzs(is:ie,js:je), start=start, count=count)
      if (baroclinic) then
        call fiona_output('r0', 'phs', state%phs(is:ie,js:je      ), start=start, count=count)
        call fiona_output('r0', 'pt' , state%pt (is:ie,js:je,ks:ke), start=start, count=count)
      else
        call fiona_output('r0', 'gz' , state%gz (is:ie,js:je,ks:ke), start=start, count=count)
      end if

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_output('r0', 'u'  , state %u(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

      call fiona_output('r0', 'v'  , state %v(is:ie,js:je,ks:ke), start=start, count=count)
    end do
    call fiona_end_output('r0')

  end subroutine restart_write

  subroutine restart_read(blocks, itime)

    type(block_type), intent(inout), target :: blocks(:)
    integer, intent(in) :: itime

    type(block_type), pointer :: block
    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    type(datetime_type) time
    integer iblk, time_step, is, ie, js, je, ks, ke
    integer start(3), count(3)
    real(8) time_value
    character(30) time_units

    if (restart_file == 'N/A') then
      call log_error('Parameter restart_file is needed to restart!')
    end if

    call fiona_open_dataset('r0', file_path=restart_file, mpi_comm=proc%comm)
    call fiona_start_input('r0')

    time_step = 1

    call fiona_input('r0', 'time', time_value, time_step=time_step)
    call fiona_get_att('r0', 'time', 'units', time_units)
    do iblk = 1, size(blocks)
      block  => blocks(iblk)
      mesh   => block%mesh
      state  => block%state(itime)
      static => block%static

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_input('r0', 'gzs', static%gzs(is:ie,js:je), start=start, count=count, time_step=time_step)
      call fill_halo(block, static%gzs, full_lon=.true., full_lat=.true.)
      if (baroclinic) then
        call fiona_input('r0', 'phs', state%phs(is:ie,js:je      ), start=start, count=count, time_step=time_step)
        call fill_halo(block, state%phs, full_lon=.true., full_lat=.true.)
        call fiona_input('r0', 'pt' , state%pt (is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
        call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      else
        call fiona_input('r0', 'gz' , state%gz (is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
        call fill_halo(block, state%gz, full_lon=.true., full_lat=.true., full_lev=.true.)
      end if

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_input('r0', 'u'  , state%u(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
      call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

      call fiona_input('r0', 'v'  , state%v(is:ie,js:je,ks:ke), start=start, count=count, time_step=time_step)
      call fill_halo(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end do
    call fiona_end_input('r0')

    call time_fast_forward(time_value, time_units)
    if (is_root_proc()) call log_notice('Restart to ' // trim(curr_time_str))

  end subroutine restart_read

  subroutine restart_final()

  end subroutine restart_final

end module restart_mod
