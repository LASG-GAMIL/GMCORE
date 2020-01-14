module process_mod

  use mpi
  use namelist_mod
  use mesh_mod
  use block_mod

  implicit none

  private

  public process_init
  public process_stop
  public process_final
  public proc

  type process_type
    integer comm
    integer dims(2)
    integer id
    integer :: ngb(4) = MPI_PROC_NULL
    type(block_type), allocatable :: blocks(:)
  end type process_type

  type(process_type) proc

contains

  subroutine process_init()

    integer ierr, i, j
    integer nproc, proc_coords(2)
    integer num_lon, num_lat, res_num, half_num
    integer lon_ibeg, lon_iend, lat_ibeg, lat_iend

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

    ! Decompose domain as a whole.
    if (num_proc_lon(1) /= 0 .and. num_proc_lat(1) /= 0) then
      proc%dims(1) = num_proc_lon(1)
      proc%dims(2) = num_proc_lat(1)
    else
      proc%dims = 0
      call MPI_DIMS_CREATE(nproc, 2, proc%dims, ierr)
    end if
    call MPI_CART_CREATE(MPI_COMM_WORLD, 2, proc%dims, [.true.,.false.], .true., proc%comm, ierr)
    call MPI_COMM_RANK(proc%comm, proc%id, ierr)
    call MPI_CART_COORDS(proc%comm, proc%id, 2, proc_coords, ierr)
    call MPI_CART_SHIFT(proc%comm, 0, 1, proc%ngb(1), proc%ngb(2), ierr)
    call MPI_CART_SHIFT(proc%comm, 1, 1, proc%ngb(3), proc%ngb(4), ierr)

    res_num = mod(global_mesh%num_full_lon, proc%dims(1))
    lon_ibeg = 1
    do i = 0, proc_coords(1) - 1
      if (res_num /= 0 .and. i < res_num - 1) then
        num_lon = global_mesh%num_full_lon / proc%dims(1) + 1
      else
        num_lon = global_mesh%num_full_lon / proc%dims(1)
      end if
      lon_ibeg = lon_ibeg + num_lon
    end do
    if (res_num /= 0 .and. proc_coords(1) < res_num) then
      num_lon = global_mesh%num_full_lon / proc%dims(1) + 1
    else
      num_lon = global_mesh%num_full_lon / proc%dims(1)
    end if
    lon_iend = lon_ibeg + num_lon - 1
#ifdef V_POLE
#else
    res_num = mod(global_mesh%num_full_lat, proc%dims(2))
    half_num = proc%dims(2) - proc%dims(2) / 2
    lat_ibeg = 1
    do j = 0, proc_coords(2) - 1
      if (res_num /= 0 .and. j >= half_num .and. j < half_num + res_num) then
        num_lat = global_mesh%num_full_lat / proc%dims(2) + 1
      else
        num_lat = global_mesh%num_full_lat / proc%dims(2)
      end if
      lat_ibeg = lat_ibeg + num_lat
    end do
    if (res_num /= 0 .and. proc_coords(2) >= half_num .and. proc_coords(2) < half_num + res_num) then
      num_lat = global_mesh%num_full_lat / proc%dims(2) + 1
    else
      num_lat = global_mesh%num_full_lat / proc%dims(2)
    end if
    lat_iend = lat_ibeg + num_lat - 1
#endif

    if (.not. allocated(proc%blocks)) allocate(proc%blocks(1))

    call proc%blocks(1)%init(proc%id, max(global_mesh%lon_halo_width, maxval(reduce_factors)), &
      global_mesh%lat_halo_width, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    allocate(proc%blocks(1)%halo(4))
    call proc%blocks(1)%halo(1)%init(proc%blocks(1)%mesh, proc_id=proc%ngb(1), west_lat_ibeg=lat_ibeg, west_lat_iend=lat_iend)
    call proc%blocks(1)%halo(2)%init(proc%blocks(1)%mesh, proc_id=proc%ngb(2), east_lat_ibeg=lat_ibeg, east_lat_iend=lat_iend)
    call proc%blocks(1)%halo(3)%init(proc%blocks(1)%mesh, proc_id=proc%ngb(3), south_lon_ibeg=lon_ibeg, south_lon_iend=lon_iend)
    call proc%blocks(1)%halo(4)%init(proc%blocks(1)%mesh, proc_id=proc%ngb(4), north_lon_ibeg=lon_ibeg, north_lon_iend=lon_iend)

  end subroutine process_init

  subroutine process_stop(code)

    integer, intent(in) :: code

    integer ierr

    call MPI_BARRIER(proc%comm, ierr)
    call MPI_ABORT(proc%comm, code, ierr)

  end subroutine process_stop

  subroutine process_final()

    integer ierr

    if (allocated(proc%blocks)) deallocate(proc%blocks)

    call MPI_FINALIZE(ierr)

  end subroutine process_final

end module process_mod
