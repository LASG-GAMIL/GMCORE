module parallel_types_mod

  use mpi
  use const_mod
  use mesh_mod

  implicit none

  private

  public zonal_circle_type
  public process_type
  public round_robin
  public proc

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
    procedure :: init  => zonal_circle_init
    procedure :: clear => zonal_circle_clear
    final :: zonal_circle_final
  end type zonal_circle_type

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
    character(30) :: hostname = ''

    integer decomp_type
    integer decomp_loc
  end type process_type

  type(process_type) proc

contains

  subroutine zonal_circle_init(this, proc)

    class(zonal_circle_type), intent(inout) :: this
    type(process_type), intent(in) :: proc

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

  subroutine zonal_circle_clear(this)

    class(zonal_circle_type), intent(inout) :: this

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

  end subroutine zonal_circle_clear

  subroutine zonal_circle_final(this)

    type(zonal_circle_type), intent(inout) :: this

    call this%clear()

  end subroutine zonal_circle_final

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

end module parallel_types_mod
