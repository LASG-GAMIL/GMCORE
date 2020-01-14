module halo_mod

  use mpi
  use const_mod
  use mesh_mod

  implicit none

  private

  public halo_type

  type halo_type
    integer :: proc_id
    integer :: iblk
    integer :: west_lat_ibeg  = -999
    integer :: west_lat_iend  = -999
    integer :: east_lat_ibeg  = -999
    integer :: east_lat_iend  = -999
    integer :: south_lon_ibeg = -999
    integer :: south_lon_iend = -999
    integer :: north_lon_ibeg = -999
    integer :: north_lon_iend = -999
    integer :: full_send_type = MPI_DATATYPE_NULL
    integer :: half_send_type = MPI_DATATYPE_NULL
    integer :: full_recv_type = MPI_DATATYPE_NULL
    integer :: half_recv_type = MPI_DATATYPE_NULL
  contains
    procedure :: init => halo_init
    final :: halo_final
  end type halo_type

contains

  subroutine halo_init(this, mesh, proc_id, iblk, west_lat_ibeg, west_lat_iend, &
                       east_lat_ibeg, east_lat_iend, south_lon_ibeg, south_lon_iend, &
                       north_lon_ibeg, north_lon_iend)

    class(halo_type), intent(out) :: this
    type(mesh_type), intent(in) :: mesh
    integer, intent(in), optional :: proc_id
    integer, intent(in), optional :: iblk
    integer, intent(in), optional :: west_lat_ibeg
    integer, intent(in), optional :: west_lat_iend
    integer, intent(in), optional :: east_lat_ibeg
    integer, intent(in), optional :: east_lat_iend
    integer, intent(in), optional :: south_lon_ibeg
    integer, intent(in), optional :: south_lon_iend
    integer, intent(in), optional :: north_lon_ibeg
    integer, intent(in), optional :: north_lon_iend

    integer lon_offset, lat_offset, ierr
    integer full_array_size(2), half_array_size(2)
    integer subarray_size(2), subarray_start(2)

    lon_offset = mesh%lon_halo_width - 1 ! MPI start index from zero, so minus one.
    lat_offset = mesh%lat_halo_width - 1
    full_array_size = [mesh%num_full_lon+2*mesh%lon_halo_width,mesh%num_full_lat+2*mesh%lat_halo_width]
    half_array_size = [mesh%num_half_lon+2*mesh%lon_halo_width,mesh%num_half_lat+2*mesh%lat_halo_width]

    !                          wx                          nx                          wx      
    !                          |                           |                           |
    !                  |---------------|---------------------------------------|---------------|
    !                  |_______________|_______________________________________|_______________|__ 
    !                  |       |       |///////|///////|///////|///////|///////|       |       |  |
    !     1 + wy + ny -|       |       |///////|///////|///////|///////|///////|       |       |  |- wy
    !                  |       |       |///////|///////|///////|///////|///////|       |       |  |
    !                  |_______|_______|_______________________________________|_______________| _|
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    ! 1 + wy + ny - 1 -|///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |- ny
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !          1 + wy -|///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|__|
    !                  |       |       |///////|///////|///////|///////|///////|       |       |  |
    !               1 -|       |       |///////|///////|///////|///////|///////|       |       |  |- wy
    !                  |       |       |///////|///////|///////|///////|///////|       |       |  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|__|
    !                      |               |                       |       |       |
    !                      1             1 + wx                    |1 + wx + nx - 1|
    !                                                              |          1 + wx + nx
    !                                                         nx - wx + 1

    if (present(west_lat_ibeg) .and. present(west_lat_iend)) then
      this%west_lat_ibeg = west_lat_ibeg
      this%west_lat_iend = west_lat_iend
      subarray_size = [mesh%lon_halo_width,west_lat_iend-west_lat_ibeg+1]
      subarray_start = [1+lon_offset,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, full_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%full_send_type, ierr)
      call MPI_TYPE_COMMIT(this%full_send_type, ierr)
      subarray_start = [0,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, full_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%full_recv_type, ierr)
      call MPI_TYPE_COMMIT(this%full_recv_type, ierr)
      subarray_size = [mesh%lon_halo_width,mesh%num_half_lat]
      subarray_start = [1+lon_offset,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, half_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%half_send_type, ierr)
      call MPI_TYPE_COMMIT(this%half_send_type, ierr)
      subarray_start = [0,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, half_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%half_recv_type, ierr)
      call MPI_TYPE_COMMIT(this%half_recv_type, ierr)
    else if (present(east_lat_ibeg) .and. present(east_lat_iend)) then
      this%east_lat_ibeg = east_lat_ibeg
      this%east_lat_iend = east_lat_iend
      subarray_size = [mesh%lon_halo_width,east_lat_iend-east_lat_ibeg+1]
      subarray_start = [mesh%num_full_lon-mesh%lon_halo_width+1+lon_offset,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, full_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%full_send_type, ierr)
      call MPI_TYPE_COMMIT(this%full_send_type, ierr)
      subarray_start = [mesh%num_full_lon+1+lon_offset,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, full_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%full_recv_type, ierr)
      call MPI_TYPE_COMMIT(this%full_recv_type, ierr)
      subarray_size = [mesh%lon_halo_width,mesh%num_half_lat]
      subarray_start = [mesh%num_half_lon-mesh%lon_halo_width+1+lon_offset,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, half_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%half_send_type, ierr)
      call MPI_TYPE_COMMIT(this%half_send_type, ierr)
      subarray_start = [mesh%num_half_lon+1+lon_offset,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, half_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%half_recv_type, ierr)
      call MPI_TYPE_COMMIT(this%half_recv_type, ierr)
    else if (present(south_lon_ibeg) .and. present(south_lon_iend)) then
      this%south_lon_ibeg = south_lon_ibeg
      this%south_lon_iend = south_lon_iend
      subarray_size = [south_lon_iend-south_lon_ibeg+1,mesh%lat_halo_width]
      subarray_start = [1+lon_offset,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, full_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%full_send_type, ierr)
      call MPI_TYPE_COMMIT(this%full_send_type, ierr)
      subarray_start = [1+lon_offset,0]
      call MPI_TYPE_CREATE_SUBARRAY(2, full_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%full_recv_type, ierr)
      call MPI_TYPE_COMMIT(this%full_recv_type, ierr)
      subarray_size = [mesh%num_half_lon,mesh%lat_halo_width]
      subarray_start = [1+lon_offset,1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, half_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%half_send_type, ierr)
      call MPI_TYPE_COMMIT(this%half_send_type, ierr)
      subarray_start = [1+lon_offset,0]
      call MPI_TYPE_CREATE_SUBARRAY(2, half_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%half_recv_type, ierr)
      call MPI_TYPE_COMMIT(this%half_recv_type, ierr)
    else if (present(north_lon_ibeg) .and. present(north_lon_iend)) then
      this%north_lon_ibeg = north_lon_ibeg
      this%north_lon_iend = north_lon_iend
      subarray_size = [north_lon_iend-north_lon_ibeg+1,mesh%lat_halo_width]
      subarray_start = [1+lon_offset,mesh%num_full_lat-mesh%lat_halo_width+1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, full_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%full_send_type, ierr)
      call MPI_TYPE_COMMIT(this%full_send_type, ierr)
      subarray_start = [1+lon_offset,mesh%num_full_lat+1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, full_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%full_recv_type, ierr)
      call MPI_TYPE_COMMIT(this%full_recv_type, ierr)
      subarray_size = [mesh%num_half_lon,mesh%lat_halo_width]
      subarray_start = [1+lon_offset,mesh%num_half_lat-mesh%lat_halo_width+1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, half_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%half_send_type, ierr)
      call MPI_TYPE_COMMIT(this%half_send_type, ierr)
      subarray_start = [1+lon_offset,mesh%num_half_lat+1+lat_offset]
      call MPI_TYPE_CREATE_SUBARRAY(2, half_array_size, subarray_size, subarray_start, &
                                    MPI_ORDER_FORTRAN, MPI_DOUBLE, this%half_recv_type, ierr)
      call MPI_TYPE_COMMIT(this%half_recv_type, ierr)
    end if

  end subroutine halo_init

  subroutine halo_final(this)

    type(halo_type), intent(inout) :: this

    integer ierr

    if (this%full_send_type /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%full_send_type, ierr)
    if (this%full_recv_type /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%full_recv_type, ierr)
    if (this%half_send_type /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%half_send_type, ierr)
    if (this%half_recv_type /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%half_recv_type, ierr)

  end subroutine halo_final

end module halo_mod
