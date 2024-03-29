module initial_mod

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

  public initial_write
  public initial_read

contains

  subroutine initial_write(initial_file, initial_time)

    character(*), intent(in) :: initial_file
    character(*), intent(in) :: initial_time

    character(4) cell_dims(4), cell_dims_2d(3)
    character(4) lon_dims(4)
    character(4) lat_dims(4)
    character(4) lev_dims(4)
    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)

    cell_dims   (1) =  'lon';  cell_dims   (2) =  'lat';  cell_dims   (3) =  'lev';     cell_dims(4) = 'time'
     lon_dims   (1) = 'ilon';   lon_dims   (2) =  'lat';   lon_dims   (3) =  'lev';      lon_dims(4) = 'time'
     lat_dims   (1) =  'lon';   lat_dims   (2) = 'ilat';   lat_dims   (3) =  'lev';      lat_dims(4) = 'time'
     lev_dims   (1) =  'lon';   lev_dims   (2) =  'lat';   lev_dims   (3) = 'ilev';      lev_dims(4) = 'time'
    cell_dims_2d(1) =  'lon';  cell_dims_2d(2) =  'lat';  cell_dims_2d(3) = 'time'

    if (is_root_proc()) call log_notice('Write ' // trim(initial_file) // '.')

    call fiona_create_dataset('i0', file_path=initial_file, start_time=initial_time, time_units='hours', mpi_comm=proc%comm, group_size=output_group_size)
    call fiona_add_dim('i0', 'time', add_var=.true.)
    call fiona_add_dim('i0',  'lon', size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0',  'lat', size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0', 'ilon', size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('i0', 'ilat', size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    if (baroclinic) then
      call fiona_add_dim('i0',  'lev', size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
      call fiona_add_dim('i0', 'ilev', size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
      call fiona_add_var('i0', 'pt'  , long_name='potential temperature'       , units='K'    , dtype=output_i0_dtype, dim_names=cell_dims)
      call fiona_add_var('i0', 'phs' , long_name='surface hydrostatic pressure', units='Pa'   , dtype=output_i0_dtype, dim_names=cell_dims_2d)
      call fiona_add_var('i0', 'u'   , long_name='u wind component'            , units='m s-1', dtype=output_i0_dtype, dim_names=lon_dims)
      call fiona_add_var('i0', 'v'   , long_name='v wind component'            , units='m s-1', dtype=output_i0_dtype, dim_names=lat_dims)
      call fiona_add_var('i0', 'p'   , long_name='hydrostatic pressure'        , units='Pa'   , dtype=output_i0_dtype, dim_names=lev_dims)
    end if
    if (nonhydrostatic) then
      call fiona_add_var('i0', 'z', long_name='height', units='m', dtype=output_i0_dtype, dim_names=lev_dims)
    end if
    call fiona_add_var('i0', 'zs'    , long_name='surface height', units='m' , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    call fiona_add_var('i0', 'zs_std', long_name='surface height subgrid standard deviation', units='m2' , dtype=output_i0_dtype, dim_names=cell_dims_2d)
    call fiona_add_var('i0', 'landmask', long_name='land mask', units='-', dtype=output_i0_dtype, dim_names=cell_dims_2d)

    call fiona_start_output('i0', 0.0d0)
    call fiona_output('i0',  'lon', global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('i0',  'lat', global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('i0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('i0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    if (baroclinic) then
      call fiona_output('i0',  'lev', global_mesh%full_lev)
      call fiona_output('i0', 'ilev', global_mesh%half_lev)
    end if

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh,     &
                 state  => blocks(iblk)%state(1), &
                 static => blocks(iblk)%static)
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

      if (baroclinic) then
        call fiona_output('i0', 'pt' , state%pt (is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output('i0', 'phs', state%phs(is:ie,js:je      ), start=start, count=count)
      end if
      call fiona_output('i0', 'zs'      , static%gzs     (is:ie,js:je) / g, start=start, count=count)
      call fiona_output('i0', 'zs_std'  , static%zs_std  (is:ie,js:je)    , start=start, count=count)
      call fiona_output('i0', 'landmask', static%landmask(is:ie,js:je)    , start=start, count=count)

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_output('i0', 'u', state%u_lon(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

      call fiona_output('i0', 'v', state%v_lat(is:ie,js:je,ks:ke), start=start, count=count)

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_half_lev]

      if (baroclinic) then
        call fiona_output('i0', 'p', state%ph_lev(is:ie,js:je,ks:ke)    , start=start, count=count)
      end if
      if (nonhydrostatic) then
        call fiona_output('i0', 'z', state%gz_lev(is:ie,js:je,ks:ke) / g, start=start, count=count)
      end if
      end associate
    end do
    call fiona_end_output('i0')

    if (is_root_proc()) call log_notice('Done write.')

  end subroutine initial_write

  subroutine initial_read(initial_file_)

    character(*), intent(in), optional :: initial_file_

    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)

    if (present(initial_file_)) then
      call fiona_open_dataset('i0', file_path=initial_file_)
      if (is_root_proc()) call log_notice('Read initial data from ' // trim(initial_file_) // '.')
    else
      call fiona_open_dataset('i0', file_path=initial_file)
      if (is_root_proc()) call log_notice('Read initial data from ' // trim(initial_file) // '.')
    end if
    call fiona_start_input('i0')

    do iblk = 1, size(blocks)
      associate (block  => blocks(iblk)                    , &
                 mesh   => blocks(iblk)%mesh               , &
                 state  => blocks(iblk)%state(old_time_idx), &
                 static => blocks(iblk)%static)
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_input('i0', 'zs', static%gzs(is:ie,js:je), start=start, count=count)
      static%gzs = static%gzs * g
      call fill_halo(block, static%gzs, full_lon=.true., full_lat=.true.)
      call fiona_input('i0', 'zs_std', static%zs_std(is:ie,js:je), start=start, count=count)
      call fill_halo(block, static%zs_std, full_lon=.true., full_lat=.true.)
      call fiona_input('i0', 'landmask', static%landmask(is:ie,js:je), start=start, count=count)
      call fill_halo(block, static%landmask, full_lon=.true., full_lat=.true.)
      if (baroclinic) then
        call fiona_input('i0', 'phs', state%phs(is:ie,js:je      ), start=start, count=count)
        call fill_halo(block, state%phs, full_lon=.true., full_lat=.true.)
        call fiona_input('i0', 'pt' , state%pt (is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      else
        call fiona_input('i0', 'z' , state%gz (is:ie,js:je,ks:ke), start=start, count=count)
        state%gz = state%gz * g
        call fill_halo(block, state%gz, full_lon=.true., full_lat=.true., full_lev=.true.)
      end if

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_input('i0', 'u'  , state%u_lon(is:ie,js:je,ks:ke), start=start, count=count)
      call fill_halo(block, state%u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

      call fiona_input('i0', 'v'  , state%v_lat(is:ie,js:je,ks:ke), start=start, count=count)
      call fill_halo(block, state%v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_half_lev]

      if (baroclinic) then
        call fiona_input('i0', 'p', state%ph_lev(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block, state%ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      end if
      if (nonhydrostatic) then
        call fiona_input('i0', 'z', state%gz_lev(is:ie,js:je,ks:ke), start=start, count=count)
        state%gz_lev = state%gz_lev * g
        call fill_halo(block, state%gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      end if
      end associate
    end do
    call fiona_end_input('i0')

    if (is_root_proc()) call log_notice('Done read initial file.')

  end subroutine initial_read

end module initial_mod
