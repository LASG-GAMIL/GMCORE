module process_mod

  use mpi
  use flogger
  use string
  use namelist_mod
  use mesh_mod
  use block_mod

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

  type process_neighbor_type
    integer :: id       = MPI_PROC_NULL
    integer :: cart_id  = MPI_PROC_NULL
    integer :: orient   = 0
    integer :: lon_ibeg = inf_i4
    integer :: lon_iend = inf_i4
    integer :: lat_ibeg = inf_i4
    integer :: lat_iend = inf_i4
    integer :: lon_halo_width = 0
  contains
    procedure :: init => process_neighbor_init
  end type process_neighbor_type

  type zonal_circle_type
    integer :: group = MPI_GROUP_NULL
    integer :: comm  = MPI_COMM_NULL
    integer :: np    = 0
    integer :: id    = MPI_PROC_NULL
    integer :: west_ngb_id    = MPI_PROC_NULL
    integer :: east_ngb_id    = MPI_PROC_NULL
    integer, allocatable :: recv_type_r4(:,:) ! 0: one level, 1: full_lev, 2: half_lev
    integer, allocatable :: recv_type_r8(:,:) ! 0: one level, 1: full_lev, 2: half_lev
  contains
    procedure :: init => zonal_circle_init
    final :: zonal_circle_final
  end type zonal_circle_type

  type process_type
    integer :: comm           = MPI_COMM_NULL
    integer :: cart_comm      = MPI_COMM_NULL
    integer :: group          = MPI_GROUP_NULL
    integer :: cart_group     = MPI_GROUP_NULL
    integer :: cart_dims(2)   = 0
    integer :: cart_coords(2) = 0
    integer :: id             = MPI_PROC_NULL          ! MPI process ID
    integer :: cart_id        = MPI_PROC_NULL          ! MPI process ID in cart_comm
    integer idom                                       ! Nest domain index (root domain is 1)
    integer np
    integer num_lon
    integer num_lat
    integer lon_ibeg
    integer lon_iend
    integer lat_ibeg
    integer lat_iend
    logical :: at_south_pole = .false.
    logical :: at_north_pole = .false.
    type(zonal_circle_type) zonal_circle
    type(process_neighbor_type), allocatable :: ngb(:) ! Neighbor processes
    type(block_type), allocatable :: blocks(:)

    integer decomp_type
    integer decomp_loc
  end type process_type

  type(process_type) proc

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

  end subroutine process_stop

  subroutine process_final()

    integer ierr

    if (allocated(proc%ngb   )) deallocate(proc%ngb   )
    if (allocated(proc%blocks)) deallocate(proc%blocks)
    if (proc%group       /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%group      , ierr)
    if (proc%cart_group  /= MPI_GROUP_NULL) call MPI_GROUP_FREE(proc%cart_group , ierr)

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
      call proc%zonal_circle%init()
    end if

  end subroutine setup_zonal_comm

  subroutine process_create_blocks()

    integer hw, i, j, dtype
    integer max_hw, lon_halo_width
    integer ierr, status(MPI_STATUS_SIZE)

    if (.not. allocated(proc%blocks)) allocate(proc%blocks(1))

    ! Use lon_halo_width in global_mesh.
    call proc%blocks(1)%init_stage_1(proc%id, global_mesh%lon_halo_width, global_mesh%lat_halo_width, &
                                     proc%lon_ibeg, proc%lon_iend, proc%lat_ibeg, proc%lat_iend)

    ! Each process calculate lon_halo_width from its filter%ngrid_lat(:).
    max_hw = 2
    do j = proc%blocks(1)%mesh%half_lat_ibeg, proc%blocks(1)%mesh%half_lat_iend
      max_hw = max(max_hw, (proc%blocks(1)%filter%ngrid_lat(j) - 1) / 2)
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
    call proc%blocks(1)%init_stage_2(global_mesh%lon_halo_width)

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
    allocate(proc%blocks(1)%halo(size(proc%ngb)))
    do i = 1, size(proc%ngb)
      select case (proc%ngb(i)%orient)
      case (west, east)
        call proc%blocks(1)%halo(i)%init(proc%blocks(1)%mesh, proc%ngb(i)%orient, dtype,               &
                                         host_id=proc%id, ngb_proc_id=proc%ngb(i)%id,                  &
                                         lat_ibeg=proc%ngb(i)%lat_ibeg, lat_iend=proc%ngb(i)%lat_iend)
      case (south, north)
        call proc%blocks(1)%halo(i)%init(proc%blocks(1)%mesh, proc%ngb(i)%orient, dtype,               &
                                         host_id=proc%id, ngb_proc_id=proc%ngb(i)%id,                  &
                                         lon_ibeg=proc%ngb(i)%lon_ibeg, lon_iend=proc%ngb(i)%lon_iend)
      end select
    end do
    ! Initialize async objects.
    do i = 1, size(proc%blocks(1)%state)
      do j = 1, size(proc%blocks(1)%state(i)%async)
        call proc%blocks(1)%state(i)%async(j)%init(size(proc%ngb))
      end do
    end do

  end subroutine process_create_blocks

  subroutine process_neighbor_init(this, orient, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    class(process_neighbor_type), intent(inout) :: this
    integer, intent(in) :: orient
    integer, intent(in), optional :: lon_ibeg
    integer, intent(in), optional :: lon_iend
    integer, intent(in), optional :: lat_ibeg
    integer, intent(in), optional :: lat_iend

    this%orient = orient

    select case (orient)
    case (west, east)
      this%lat_ibeg = lat_ibeg
      this%lat_iend = lat_iend
    case (south, north)
      this%lon_ibeg = lon_ibeg
      this%lon_iend = lon_iend
    end select

  end subroutine process_neighbor_init

  subroutine zonal_circle_init(this)

    class(zonal_circle_type), intent(inout) :: this

    integer ierr, i, num_lon, ibeg, iend
    integer west_cart_id, east_cart_id, tmp_id(1)
    integer, allocatable :: zonal_proc_id(:)

    allocate(zonal_proc_id(proc%cart_dims(1)))
    do i = 1, proc%cart_dims(1)
      call MPI_CART_RANK(proc%cart_comm, [i-1,proc%cart_coords(2)], zonal_proc_id(i), ierr)
    end do
    call MPI_GROUP_INCL(proc%cart_group, size(zonal_proc_id), zonal_proc_id, this%group, ierr)
    call MPI_COMM_CREATE_GROUP(proc%cart_comm, this%group, sum(zonal_proc_id), this%comm, ierr)
    call MPI_COMM_SIZE(this%comm, this%np, ierr)
    call MPI_COMM_RANK(this%comm, this%id, ierr)
    deallocate(zonal_proc_id)

    ! Get IDs of the west and east neighbors in zonal circle comm.
    call MPI_CART_SHIFT(proc%cart_comm, 0, 1, west_cart_id, east_cart_id, ierr)
    call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [west_cart_id], this%group, tmp_id, ierr); this%west_ngb_id = tmp_id(1)
    call MPI_GROUP_TRANSLATE_RANKS(proc%cart_group, 1, [east_cart_id], this%group, tmp_id, ierr); this%east_ngb_id = tmp_id(1)

    if (this%id == 0) then
      ! Single precision
      allocate(this%recv_type_r4(this%np,0:2))
      do i = 1, this%np
        num_lon = global_mesh%num_full_lon
        call round_robin(this%np, i - 1, num_lon, ibeg, iend)
        call MPI_TYPE_CREATE_SUBARRAY(1, [global_mesh%num_full_lon], &
                                         [                 num_lon], &
                                         [ibeg-1], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,0), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,0), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%num_full_lon,global_mesh%num_full_lev], &
                                         [                 num_lon,global_mesh%num_full_lev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,1), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,1), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%num_full_lon,global_mesh%num_half_lev], &
                                         [                 num_lon,global_mesh%num_half_lev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_REAL, &
                                         this%recv_type_r4(i,2), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r4(i,2), ierr)
      end do
      ! Double precision
      allocate(this%recv_type_r8(this%np,0:2))
      do i = 1, this%np
        num_lon = global_mesh%num_full_lon
        call round_robin(this%np, i - 1, num_lon, ibeg, iend)
        call MPI_TYPE_CREATE_SUBARRAY(1, [global_mesh%num_full_lon], &
                                         [                 num_lon], &
                                         [ibeg-1], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,0), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,0), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%num_full_lon,global_mesh%num_full_lev], &
                                         [                 num_lon,global_mesh%num_full_lev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,1), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,1), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, [global_mesh%num_full_lon,global_mesh%num_half_lev], &
                                         [                 num_lon,global_mesh%num_half_lev], &
                                         [ibeg-1,0], MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                         this%recv_type_r8(i,2), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_r8(i,2), ierr)
      end do
    end if

  end subroutine zonal_circle_init

  subroutine zonal_circle_final(this)

    type(zonal_circle_type), intent(inout) :: this

    integer i, k, ierr

    if (allocated(this%recv_type_r4)) then
      do k = 0, 2
        do i = 1, this%np
          call MPI_TYPE_FREE(this%recv_type_r4(i,k), ierr)
        end do
        deallocate(this%recv_type_r4)
      end do
    end if

    if (allocated(this%recv_type_r8)) then
      do k = 0, 2
        do i = 1, this%np
          call MPI_TYPE_FREE(this%recv_type_r8(i,k), ierr)
        end do
        deallocate(this%recv_type_r8)
      end do
    end if

    if (this%group /= MPI_GROUP_NULL) call MPI_GROUP_FREE(this%group, ierr)

  end subroutine zonal_circle_final

  subroutine round_robin(dim, coord, num, ibeg, iend)

    integer, intent(in) :: dim
    integer, intent(in) :: coord
    integer, intent(inout) :: num
    integer, intent(out) :: ibeg ! Start from 1.
    integer, intent(out) :: iend ! Start from 1.

    integer res_num, tmp_num, i

    res_num = mod(num, dim)
    ibeg = 1
    do i = 0, coord - 1
      if (res_num /= 0 .and. i < res_num) then
        tmp_num = num / dim + 1
      else
        tmp_num = num / dim
      end if
      ibeg = ibeg + tmp_num
    end do
    if (res_num /= 0 .and. coord < res_num) then
      num = num / dim + 1
    else
      num = num / dim
    end if
    iend = ibeg + num - 1

  end subroutine round_robin

end module process_mod
