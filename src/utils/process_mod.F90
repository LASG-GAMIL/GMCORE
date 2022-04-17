module process_mod

  use mpi
  use flogger
  use string
  use const_mod
  use namelist_mod
  use mesh_mod
  use block_mod
  use parallel_types_mod

  implicit none

  private

  public process_init
  public process_create_blocks
  public process_stop
  public process_final
  public proc
  public is_root_proc
  public zonal_circle_type

  integer, public, parameter :: decomp_1d_lat        = 1
  integer, public, parameter :: decomp_2d_simple     = 2

  integer, public, parameter :: decomp_normal_region = 5

contains

  subroutine process_init(comm)

    integer, intent(in), optional :: comm

    integer ierr

    if (present(comm)) then
      proc%comm = comm
    else
      call MPI_INIT(ierr)
      proc%comm = MPI_COMM_WORLD
    end if
    call MPI_COMM_GROUP(proc%comm, proc%group, ierr)
    call MPI_COMM_SIZE(proc%comm, proc%np, ierr)
    call MPI_COMM_RANK(proc%comm, proc%id, ierr)

    call setup_mpi_simple()
    call decompose_domains()
    call setup_zonal_comm()

  end subroutine process_init

  subroutine process_stop(code)

    integer, intent(in) :: code

    integer ierr

    call MPI_ABORT(proc%comm, code, ierr)
    call proc%zonal_circle%clear()

  end subroutine process_stop

  subroutine process_final()

    integer i, ierr

    do i = 1, size(blocks)
      call blocks(i)%clear()
    end do

    if (allocated(proc%ngb)) deallocate(proc%ngb)
    if (allocated(blocks  )) deallocate(blocks  )
    if (proc%group      /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%group     , ierr)
    if (proc%cart_group /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%cart_group, ierr)

    call MPI_FINALIZE(ierr)

  end subroutine process_final

  pure logical function is_root_proc()

    is_root_proc = proc%id == 0

  end function is_root_proc

  subroutine setup_mpi_simple()

    integer ierr, np, tmp_comm, i
    logical periods(2)

    proc%decomp_type = decomp_1d_lat
    proc%decomp_loc  = decomp_normal_region

    if (num_proc_lon(1) * num_proc_lat(1) == proc%np) then
      ! Check if process topology in namelist is compatible with MPI runtime.
      np = 0
      do i = 1, 1
        np = np + num_proc_lon(i) * num_proc_lat(i)
      end do
      if (proc%np /= np .and. is_root_proc()) then
        call log_notice('Namelist num_proc_lon and num_proc_lat are not compatible with MPI runtime. Reset to MPI runtime.')
        num_proc_lat(1) = proc%np
      end if
      ! Set the process topology into proc object.
      np = 0
      do i = 1, 1
        np = np + num_proc_lon(i) * num_proc_lat(i)
        if (proc%id + 1 <= np) then
          proc%cart_dims(1) = num_proc_lon(i)
          proc%cart_dims(2) = num_proc_lat(i)
          proc%idom = i
          exit
        end if
      end do
    else
      proc%cart_dims = [1, proc%np]
      proc%idom = 1
    end if
    periods = [.true.,.false.]
    ! Set MPI process topology.
    call MPI_COMM_SPLIT(proc%comm, proc%idom, proc%id, tmp_comm, ierr)
    call MPI_CART_CREATE(tmp_comm, 2, proc%cart_dims, periods, .true., proc%cart_comm, ierr)
    call MPI_COMM_GROUP(proc%cart_comm, proc%cart_group, ierr)
    call MPI_COMM_FREE(tmp_comm, ierr)
    call MPI_COMM_RANK(proc%cart_comm, proc%cart_id, ierr)
    call MPI_CART_COORDS(proc%cart_comm, proc%cart_id, 2, proc%cart_coords, ierr)

  end subroutine setup_mpi_simple

  subroutine decompose_domains()

    integer ierr, tmp_id(1), i, j

    ! Set neighborhood of the process.
    if (allocated(proc%ngb)) deallocate(proc%ngb)
    select case (proc%decomp_loc)
    case (decomp_normal_region)
      allocate(proc%ngb(4))
      call MPI_CART_SHIFT(proc%cart_comm, 0, 1, proc%ngb(west )%cart_id, proc%ngb(east )%cart_id, ierr)
      call MPI_CART_SHIFT(proc%cart_comm, 1, 1, proc%ngb(south)%cart_id, proc%ngb(north)%cart_id, ierr)
    end select

    ! Translate Cartesian ID of neighbors to global ID.
    do i = 1, size(proc%ngb)
      if (proc%ngb(i)%id == MPI_PROC_NULL) then
        call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [proc%ngb(i)%cart_id], proc%group, tmp_id, ierr)
        proc%ngb(i)%id = tmp_id(1)
      end if
    end do

    ! Handle processes at poles.
    if (proc%ngb(south)%id == MPI_PROC_NULL) then
      i = proc%id + proc%cart_dims(1) / 2 * proc%cart_dims(2)
      if (i >= proc%np) i = i - proc%np
      proc%ngb(south)%id = i
      proc%at_south_pole = .true.
    end if
    if (proc%ngb(north)%id == MPI_PROC_NULL) then
      i = proc%id + proc%cart_dims(1) / 2 * proc%cart_dims(2)
      if (i >= proc%np) i = i - proc%np
      proc%ngb(north)%id = i
      proc%at_north_pole = .true.
    end if

    ! Set initial values for num_lon, num_lat, lon_ibeg, lat_ibeg.
    proc%num_lon = global_mesh%num_full_lon
    select case (proc%decomp_loc)
    case (decomp_normal_region)
      proc%num_lat = global_mesh%num_full_lat
    end select

    call round_robin(proc%cart_dims(1), proc%cart_coords(1), proc%num_lon, proc%lon_ibeg, proc%lon_iend)
    call round_robin(proc%cart_dims(2), proc%cart_coords(2), proc%num_lat, proc%lat_ibeg, proc%lat_iend)

  end subroutine decompose_domains

  subroutine setup_zonal_comm()

    if (proc%idom == 1) then ! Only root domain has polar region.
      call proc%zonal_circle%init(proc)
    end if

  end subroutine setup_zonal_comm

  subroutine process_create_blocks()

    integer hw, i, j, dtype
    integer max_hw, lon_halo_width
    integer ierr, status(MPI_STATUS_SIZE)

    if (.not. allocated(blocks)) allocate(blocks(1))

    ! Use lon_halo_width in global_mesh.
    call blocks(1)%init_stage_1(proc%id, global_mesh%lon_halo_width, global_mesh%lat_halo_width, &
                                proc%lon_ibeg, proc%lon_iend, proc%lat_ibeg, proc%lat_iend)

    ! Each process calculate lon_halo_width from its filter%ngrid_lat(:).
    max_hw = 2
    do j = blocks(1)%mesh%half_lat_ibeg, blocks(1)%mesh%half_lat_iend
      max_hw = max(max_hw, (blocks(1)%filter%ngrid_lat(j) - 1) / 2)
    end do
    lon_halo_width = max(max_hw, global_mesh%lon_halo_width)

    ! Get lon_halo_width from southern and northern neighbors.
    if (proc%ngb(south)%id /= MPI_PROC_NULL) then
      call MPI_SENDRECV(lon_halo_width, 1, MPI_INT, proc%ngb(south)%id, 100, &
                        proc%ngb(south)%lon_halo_width, 1, MPI_INT, proc%ngb(south)%id, 100, &
                        proc%comm, status, ierr)
    else
      proc%ngb(south)%lon_halo_width = 0
    end if
    if (proc%ngb(north)%id /= MPI_PROC_NULL) then
      call MPI_SENDRECV(lon_halo_width, 1, MPI_INT, proc%ngb(north)%id, 100, &
                        proc%ngb(north)%lon_halo_width, 1, MPI_INT, proc%ngb(north)%id, 100, &
                        proc%comm, status, ierr)
    else
      proc%ngb(north)%lon_halo_width = 0
    end if

    call global_mesh%reinit(max(lon_halo_width, proc%ngb(south)%lon_halo_width, proc%ngb(north)%lon_halo_width))
    call blocks(1)%init_stage_2(global_mesh%lon_halo_width)

    call proc%ngb(west )%init(west , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
    call proc%ngb(east )%init(east , lat_ibeg=proc%lat_ibeg, lat_iend=proc%lat_iend)
    hw = max(lon_halo_width, proc%ngb(south)%lon_halo_width)
    call proc%ngb(south)%init(south, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw)
    hw = max(lon_halo_width, proc%ngb(north)%lon_halo_width)
    call proc%ngb(north)%init(north, lon_ibeg=proc%lon_ibeg-hw, lon_iend=proc%lon_iend+hw)

    select case (r8)
    case (4)
      dtype = MPI_REAL
    case (8)
      dtype = MPI_DOUBLE
    case (16)
      dtype = MPI_REAL16
    case default
      call log_error('Unsupported parameter r8!')
    end select

    ! Setup halos (only normal halos for the time being).
    allocate(blocks(1)%halo(size(proc%ngb)))
    do i = 1, size(proc%ngb)
      select case (proc%ngb(i)%orient)
      case (west, east)
        call blocks(1)%halo(i)%init(blocks(1)%mesh, proc%ngb(i)%orient, dtype,                    &
                                    host_id=proc%id, ngb_proc_id=proc%ngb(i)%id,                  &
                                    lat_ibeg=proc%ngb(i)%lat_ibeg, lat_iend=proc%ngb(i)%lat_iend)
      case (south, north)
        call blocks(1)%halo(i)%init(blocks(1)%mesh, proc%ngb(i)%orient, dtype,                    &
                                    host_id=proc%id, ngb_proc_id=proc%ngb(i)%id,                  &
                                    lon_ibeg=proc%ngb(i)%lon_ibeg, lon_iend=proc%ngb(i)%lon_iend, &
                                    at_south_pole=proc%at_south_pole, at_north_pole=proc%at_north_pole)
      end select
    end do

  end subroutine process_create_blocks

end module process_mod
