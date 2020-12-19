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

  public initial_read

contains

  subroutine initial_read(initial_file_)

    character(*), intent(in), optional :: initial_file_

    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)

    if (present(initial_file_)) then
      call fiona_open_dataset('i0', file_path=initial_file_, mpi_comm=proc%comm)
      if (is_root_proc()) call log_notice('Read initial data from ' // trim(initial_file_) // '.')
    else
      call fiona_open_dataset('i0', file_path=initial_file, mpi_comm=proc%comm)
      if (is_root_proc()) call log_notice('Read initial data from ' // trim(initial_file) // '.')
    end if
    call fiona_start_input('i0')

    do iblk = 1, size(proc%blocks)
      associate (block  => proc%blocks(iblk)                    , &
                 mesh   => proc%blocks(iblk)%mesh               , &
                 state  => proc%blocks(iblk)%state(old_time_idx), &
                 static => proc%blocks(iblk)%static)
        is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
        js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

        call fiona_input('i0', 'gzs', static%gzs(is:ie,js:je), start=start, count=count)
        call fill_halo(block, static%gzs, full_lon=.true., full_lat=.true.)
        if (baroclinic) then
          call fiona_input('i0', 'phs', state%phs(is:ie,js:je      ), start=start, count=count)
          call fill_halo(block, state%phs, full_lon=.true., full_lat=.true.)
          call fiona_input('i0', 'pt' , state%pt (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
        else
          call fiona_input('i0', 'gz' , state%gz (is:ie,js:je,ks:ke), start=start, count=count)
          call fill_halo(block, state%gz, full_lon=.true., full_lat=.true., full_lev=.true.)
        end if

        is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
        js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

        call fiona_input('i0', 'u'  , state%u(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

        is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
        js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

        call fiona_input('i0', 'v'  , state%v(is:ie,js:je,ks:ke), start=start, count=count)
        call fill_halo(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
      end associate
    end do
    call fiona_end_input('i0')

  end subroutine initial_read

end module initial_mod
