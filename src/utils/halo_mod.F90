module halo_mod

  use mpi
  use flogger
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
    ! (1,1): full_lon,full_lat (1,2): full_lon,half_lat
    ! (2,1): half_lon,full_lat (2,2): half_lon,half_lat
    integer :: array_size(2,2,2) = 0
    integer :: send_subarray_size(2,2,2)
    integer :: recv_subarray_size(2,2,2)
    integer :: send_subarray_start(2,2,2)
    integer :: recv_subarray_start(2,2,2)
    integer :: send_type(2,2) = MPI_DATATYPE_NULL
    integer :: recv_type(2,2) = MPI_DATATYPE_NULL
  contains
    procedure :: init_normal => halo_init_normal
    procedure :: init_nest => halo_init_nest
    final :: halo_final
  end type halo_type

contains

  subroutine halo_init_normal(this, mesh, ngb_proc_id, iblk, west_lat_ibeg, west_lat_iend, &
                              east_lat_ibeg, east_lat_iend, south_lon_ibeg, south_lon_iend, &
                              north_lon_ibeg, north_lon_iend)

    class(halo_type), intent(out) :: this
    type(mesh_type), intent(in) :: mesh
    integer, intent(in), optional :: ngb_proc_id
    integer, intent(in), optional :: iblk
    integer, intent(in), optional :: west_lat_ibeg
    integer, intent(in), optional :: west_lat_iend
    integer, intent(in), optional :: east_lat_ibeg
    integer, intent(in), optional :: east_lat_iend
    integer, intent(in), optional :: south_lon_ibeg
    integer, intent(in), optional :: south_lon_iend
    integer, intent(in), optional :: north_lon_ibeg
    integer, intent(in), optional :: north_lon_iend

    integer i, j, ierr

    if (present(ngb_proc_id)) then
      this%proc_id = ngb_proc_id
    else if (present(iblk)) then
      call log_error('Handle internal halo!', __FILE__, __LINE__)
    end if

    this%array_size(:,1,1) = [mesh%num_full_lon+2*mesh%lon_halo_width,mesh%num_full_lat+2*mesh%lat_halo_width]
    this%array_size(:,2,1) = [mesh%num_half_lon+2*mesh%lon_halo_width,mesh%num_full_lat+2*mesh%lat_halo_width]
    this%array_size(:,1,2) = [mesh%num_full_lon+2*mesh%lon_halo_width,mesh%num_half_lat+2*mesh%lat_halo_width]
    this%array_size(:,2,2) = [mesh%num_half_lon+2*mesh%lon_halo_width,mesh%num_half_lat+2*mesh%lat_halo_width]

    !                          wx                          nx                          wx      
    !                          |                           |                           |
    !                  |---------------|---------------------------------------|---------------|
    !                  |_______________|_______________________________________|_______________|__ 
    !                  |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|  |
    !         wy + ny -|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|  |- wy
    !                  |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|  |
    !                  |_______|_______|_______________________________________|_______________|__|
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !     wy + ny - 1 -|///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |- ny
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !              wy -|///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |///////|///////|       |       |       |       |       |///////|///////|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|__|
    !                  |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|  |
    !               0 -|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|  |- wy
    !                  |\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|\\\\\\\|  |
    !                  |_______|_______|_______|_______|_______|_______|_______|_______|_______|__|
    !                      |               |                       |       |       |
    !                      0               wx                      | wx + nx - 1   |
    !                                                              |            wx + nx
    !                                                              nx

    ! NOTE: MPI array index starts from zero.
    if (present(west_lat_ibeg) .and. present(west_lat_iend)) then
      this%west_lat_ibeg = west_lat_ibeg
      this%west_lat_iend = west_lat_iend
      ! full_lon + full_lat
      this%send_subarray_size (:,1,1) = [mesh%lon_halo_width,mesh%num_full_lat]
      this%recv_subarray_size (:,1,1) = this%send_subarray_size(:,1,1)
      this%send_subarray_start(:,1,1) = [mesh%lon_halo_width,mesh%lat_halo_width]
      this%recv_subarray_start(:,1,1) = [                  0,mesh%lat_halo_width]
      ! half_lon + full_lat
      this%send_subarray_size (:,2,1) = [mesh%lon_halo_width,mesh%num_full_lat]
      this%recv_subarray_size (:,2,1) = this%send_subarray_size(:,2,1)
      this%send_subarray_start(:,2,1) = [mesh%lon_halo_width,mesh%lat_halo_width]
      this%recv_subarray_start(:,2,1) = [                  0,mesh%lat_halo_width]
      ! full_lon + half_lat
      this%send_subarray_size (:,1,2) = [mesh%lon_halo_width,mesh%num_half_lat]
      this%recv_subarray_size (:,1,2) = this%send_subarray_size(:,1,2)
      this%send_subarray_start(:,1,2) = [mesh%lon_halo_width,mesh%lat_halo_width]
      this%recv_subarray_start(:,1,2) = [                  0,mesh%lat_halo_width]
      ! half_lon + half_lat
      this%send_subarray_size (:,2,2) = [mesh%lon_halo_width,mesh%num_half_lat]
      this%recv_subarray_size (:,2,2) = this%send_subarray_size(:,2,2)
      this%send_subarray_start(:,2,2) = [mesh%lon_halo_width,mesh%lat_halo_width]
      this%recv_subarray_start(:,2,2) = [                  0,mesh%lat_halo_width]
    else if (present(east_lat_ibeg) .and. present(east_lat_iend)) then
      this%east_lat_ibeg = east_lat_ibeg
      this%east_lat_iend = east_lat_iend
      ! full_lon + full_lat
      this%send_subarray_size (:,1,1) = [mesh%lon_halo_width,mesh%num_full_lat]
      this%recv_subarray_size (:,1,1) = this%send_subarray_size(:,1,1)
      this%send_subarray_start(:,1,1) = [mesh%num_full_lon                    ,mesh%lat_halo_width]
      this%recv_subarray_start(:,1,1) = [mesh%num_full_lon+mesh%lon_halo_width,mesh%lat_halo_width]
      ! half_lon + full_lat
      this%send_subarray_size (:,2,1) = [mesh%lon_halo_width,mesh%num_full_lat]
      this%recv_subarray_size (:,2,1) = this%send_subarray_size(:,2,1)
      this%send_subarray_start(:,2,1) = [mesh%num_half_lon                    ,mesh%lat_halo_width]
      this%recv_subarray_start(:,2,1) = [mesh%num_half_lon+mesh%lon_halo_width,mesh%lat_halo_width]
      ! full_lon + half_lat
      this%send_subarray_size (:,1,2) = [mesh%lon_halo_width,mesh%num_half_lat]
      this%recv_subarray_size (:,1,2) = this%send_subarray_size(:,1,2)
      this%send_subarray_start(:,1,2) = [mesh%num_full_lon                    ,mesh%lat_halo_width]
      this%recv_subarray_start(:,1,2) = [mesh%num_full_lon+mesh%lon_halo_width,mesh%lat_halo_width]
      ! half_lon + half_lat
      this%send_subarray_size (:,2,2) = [                  mesh%lon_halo_width,mesh%num_half_lat  ]
      this%recv_subarray_size (:,2,2) = this%send_subarray_size(:,2,2)
      this%send_subarray_start(:,2,2) = [mesh%num_half_lon                    ,mesh%lat_halo_width]
      this%recv_subarray_start(:,2,2) = [mesh%num_half_lon+mesh%lon_halo_width,mesh%lat_halo_width]
    else if (present(south_lon_ibeg) .and. present(south_lon_iend)) then
      this%south_lon_ibeg = south_lon_ibeg
      this%south_lon_iend = south_lon_iend
      ! full_lon + full_lat
      this%send_subarray_size (:,1,1) = [this%array_size(1,1,1),mesh%lat_halo_width]
      this%recv_subarray_size (:,1,1) = this%send_subarray_size(:,1,1)
      this%send_subarray_start(:,1,1) = [0,mesh%lat_halo_width]
      this%recv_subarray_start(:,1,1) = [0,0]
      ! half_lon + full_lat
      this%send_subarray_size (:,2,1) = [this%array_size(1,2,1),mesh%lat_halo_width]
      this%recv_subarray_size (:,2,1) = this%send_subarray_size(:,2,1)
      this%send_subarray_start(:,2,1) = [0,mesh%lat_halo_width]
      this%recv_subarray_start(:,2,1) = [0,0]
      ! full_lon + half_lat
      this%send_subarray_size (:,1,2) = [this%array_size(1,1,2),mesh%lat_halo_width]
      this%recv_subarray_size (:,1,2) = this%send_subarray_size(:,1,2)
      this%send_subarray_start(:,1,2) = [0,mesh%lat_halo_width]
      this%recv_subarray_start(:,1,2) = [0,0]
      ! half_lon + half_lat
      this%send_subarray_size (:,2,2) = [this%array_size(1,2,2),mesh%lat_halo_width]
      this%recv_subarray_size (:,2,2) = this%send_subarray_size(:,2,2)
      this%send_subarray_start(:,2,2) = [0,mesh%lat_halo_width]
      this%recv_subarray_start(:,2,2) = [0,0]
    else if (present(north_lon_ibeg) .and. present(north_lon_iend)) then
      this%north_lon_ibeg = north_lon_ibeg
      this%north_lon_iend = north_lon_iend
      ! full_lon + full_lat
      this%send_subarray_size (:,1,1) = [this%array_size(1,1,1),mesh%lat_halo_width]
      this%recv_subarray_size (:,1,1) = this%send_subarray_size(:,1,1)
      this%send_subarray_start(:,1,1) = [0,mesh%num_full_lat]
      this%recv_subarray_start(:,1,1) = [0,mesh%num_full_lat+mesh%lat_halo_width]
      ! half_lon + full_lat
      this%send_subarray_size (:,2,1) = [this%array_size(1,2,1),mesh%lat_halo_width]
      this%recv_subarray_size (:,2,1) = this%send_subarray_size(:,2,1)
      this%send_subarray_start(:,2,1) = [0,mesh%num_full_lat]
      this%recv_subarray_start(:,2,1) = [0,mesh%num_full_lat+mesh%lat_halo_width]
      ! full_lon + half_lat
      this%send_subarray_size (:,1,2) = [this%array_size(1,1,2),mesh%lat_halo_width]
      this%recv_subarray_size (:,1,2) = this%send_subarray_size(:,1,2)
      this%send_subarray_start(:,1,2) = [0,mesh%num_half_lat]
      this%recv_subarray_start(:,1,2) = [0,mesh%num_half_lat+mesh%lat_halo_width]
      ! half_lon + half_lat
      this%send_subarray_size (:,2,2) = [this%array_size(1,2,2),mesh%lat_halo_width]
      this%recv_subarray_size (:,2,2) = this%send_subarray_size(:,2,2)
      this%send_subarray_start(:,2,2) = [0,mesh%num_half_lat]
      this%recv_subarray_start(:,2,2) = [0,mesh%num_half_lat+mesh%lat_halo_width]
    end if

    do j = 1, 2
      do i = 1, 2
        call MPI_TYPE_CREATE_SUBARRAY(2, this%array_size(:,i,j), this%send_subarray_size(:,i,j), &
                                      this%send_subarray_start(:,i,j), MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                      this%send_type(i,j), ierr)
        call MPI_TYPE_COMMIT(this%send_type(i,j), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, this%array_size(:,i,j), this%recv_subarray_size(:,i,j), &
                                      this%recv_subarray_start(:,i,j), MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                      this%recv_type(i,j), ierr)
        call MPI_TYPE_COMMIT(this%recv_type(i,j), ierr)
      end do
    end do

  end subroutine halo_init_normal

  subroutine halo_init_nest(this, parent_mesh, parent_proc_id)

    class(halo_type), intent(inout) :: this
    type(mesh_type), intent(in) :: parent_mesh
    integer, intent(in) :: parent_proc_id

  end subroutine halo_init_nest

  subroutine halo_final(this)

    type(halo_type), intent(inout) :: this

    integer i, j
    integer ierr

    do j = 1, 2
      do i = 1, 2
        if (this%send_type(i,j) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%send_type(i,j), ierr)
        if (this%recv_type(i,j) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%recv_type(i,j), ierr)
      end do
    end do

  end subroutine halo_final

end module halo_mod
