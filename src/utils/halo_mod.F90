module halo_mod

  use mpi
  use flogger
  use const_mod
  use mesh_mod

  implicit none

  private

  public halo_type

  integer, parameter :: cross_proc_halo = 1
  integer, parameter :: cross_comm_halo = 2
  integer, parameter :: inner_halo = 3
  integer, parameter :: nest_halo = 4

  type halo_type
    integer :: comm = MPI_COMM_NULL
    integer :: proc_id = MPI_PROC_NULL
    integer :: iblk = 0
    integer :: orient = 0
    integer :: type = 0
    ! (1,1): full_lon,full_lat (1,2): full_lon,half_lat
    ! (2,1): half_lon,full_lat (2,2): half_lon,half_lat
    integer :: send_type_2d(2,2) = MPI_DATATYPE_NULL
    integer :: recv_type_2d(2,2) = MPI_DATATYPE_NULL
    ! (1,1,1): full_lon,full_lat,full_lev (1,2,1): full_lon,half_lat,full_lev
    ! (2,1,1): half_lon,full_lat,full_lev (2,2,1): half_lon,half_lat,full_lev
    ! (1,1,2): full_lon,full_lat,half_lev (1,2,2): full_lon,half_lat,half_lev
    ! (2,1,2): half_lon,full_lat,half_lev (2,2,2): half_lon,half_lat,half_lev
    integer :: send_type_3d(2,2,2) = MPI_DATATYPE_NULL
    integer :: recv_type_3d(2,2,2) = MPI_DATATYPE_NULL
  contains
    procedure :: init_normal => halo_init_normal
    procedure :: init_nest => halo_init_nest
    final :: halo_final
  end type halo_type

contains

  subroutine halo_init_normal(this, mesh, orient, ngb_proc_id, iblk)

    class(halo_type), intent(out) :: this
    type(mesh_type), intent(in) :: mesh
    integer, intent(in) :: orient
    integer, intent(in), optional :: ngb_proc_id
    integer, intent(in), optional :: iblk

    integer array_size(3,2,2)
    integer send_subarray_size(3,2,2)
    integer recv_subarray_size(3,2,2)
    integer send_subarray_start(3,2,2)
    integer recv_subarray_start(3,2,2)
    integer num_lev(2)
    integer i, j, k, ierr

    if (present(ngb_proc_id)) then
      this%proc_id = ngb_proc_id
    else if (present(iblk)) then
      call log_error('Handle internal halo!', __FILE__, __LINE__)
    end if

    ! NOTE: MPI array index starts from zero.

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

    this%orient = orient
    this%type = cross_proc_halo
    num_lev = [mesh%num_full_lev,mesh%num_half_lev]

    do k = 1, 2
      array_size(:,1,1) = [mesh%num_full_lon+2*mesh%lon_halo_width,mesh%num_full_lat+2*mesh%lat_halo_width,num_lev(k)]
      array_size(:,2,1) = [mesh%num_half_lon+2*mesh%lon_halo_width,mesh%num_full_lat+2*mesh%lat_halo_width,num_lev(k)]
      array_size(:,1,2) = [mesh%num_full_lon+2*mesh%lon_halo_width,mesh%num_half_lat+2*mesh%lat_halo_width,num_lev(k)]
      array_size(:,2,2) = [mesh%num_half_lon+2*mesh%lon_halo_width,mesh%num_half_lat+2*mesh%lat_halo_width,num_lev(k)]
      select case (orient)
      case (west)
        ! full_lon + full_lat
        send_subarray_size (:,1,1) = [mesh%lon_halo_width,mesh%num_full_lat,num_lev(k)]
        recv_subarray_size (:,1,1) = send_subarray_size(:,1,1)
        send_subarray_start(:,1,1) = [mesh%lon_halo_width,mesh%lat_halo_width,0]
        recv_subarray_start(:,1,1) = [                  0,mesh%lat_halo_width,0]
        ! half_lon + full_lat
        send_subarray_size (:,2,1) = [mesh%lon_halo_width,mesh%num_full_lat,num_lev(k)]
        recv_subarray_size (:,2,1) = send_subarray_size(:,2,1)
        send_subarray_start(:,2,1) = [mesh%lon_halo_width,mesh%lat_halo_width,0]
        recv_subarray_start(:,2,1) = [                  0,mesh%lat_halo_width,0]
        ! full_lon + half_lat
        send_subarray_size (:,1,2) = [mesh%lon_halo_width,mesh%num_half_lat,num_lev(k)]
        recv_subarray_size (:,1,2) = send_subarray_size(:,1,2)
        send_subarray_start(:,1,2) = [mesh%lon_halo_width,mesh%lat_halo_width,0]
        recv_subarray_start(:,1,2) = [                  0,mesh%lat_halo_width,0]
        ! half_lon + half_lat
        send_subarray_size (:,2,2) = [mesh%lon_halo_width,mesh%num_half_lat,num_lev(k)]
        recv_subarray_size (:,2,2) = send_subarray_size(:,2,2)
        send_subarray_start(:,2,2) = [mesh%lon_halo_width,mesh%lat_halo_width,0]
        recv_subarray_start(:,2,2) = [                  0,mesh%lat_halo_width,0]
      case (east)
        ! full_lon + full_lat
        send_subarray_size (:,1,1) = [mesh%lon_halo_width,mesh%num_full_lat,num_lev(k)]
        recv_subarray_size (:,1,1) = send_subarray_size(:,1,1)
        send_subarray_start(:,1,1) = [mesh%num_full_lon                    ,mesh%lat_halo_width,0]
        recv_subarray_start(:,1,1) = [mesh%num_full_lon+mesh%lon_halo_width,mesh%lat_halo_width,0]
        ! half_lon + full_lat
        send_subarray_size (:,2,1) = [mesh%lon_halo_width,mesh%num_full_lat,num_lev(k)]
        recv_subarray_size (:,2,1) = send_subarray_size(:,2,1)
        send_subarray_start(:,2,1) = [mesh%num_half_lon                    ,mesh%lat_halo_width,0]
        recv_subarray_start(:,2,1) = [mesh%num_half_lon+mesh%lon_halo_width,mesh%lat_halo_width,0]
        ! full_lon + half_lat
        send_subarray_size (:,1,2) = [mesh%lon_halo_width,mesh%num_half_lat,num_lev(k)]
        recv_subarray_size (:,1,2) = send_subarray_size(:,1,2)
        send_subarray_start(:,1,2) = [mesh%num_full_lon                    ,mesh%lat_halo_width,0]
        recv_subarray_start(:,1,2) = [mesh%num_full_lon+mesh%lon_halo_width,mesh%lat_halo_width,0]
        ! half_lon + half_lat
        send_subarray_size (:,2,2) = [mesh%lon_halo_width,mesh%num_half_lat,num_lev(k)]
        recv_subarray_size (:,2,2) = send_subarray_size(:,2,2)
        send_subarray_start(:,2,2) = [mesh%num_half_lon                    ,mesh%lat_halo_width,0]
        recv_subarray_start(:,2,2) = [mesh%num_half_lon+mesh%lon_halo_width,mesh%lat_halo_width,0]
      case (south)
        ! full_lon + full_lat
        send_subarray_size (:,1,1) = [array_size(1,1,1),mesh%lat_halo_width,num_lev(k)]
        recv_subarray_size (:,1,1) = send_subarray_size(:,1,1)
        send_subarray_start(:,1,1) = [0,mesh%lat_halo_width,0]
        recv_subarray_start(:,1,1) = [0,0,0]
        ! half_lon + full_lat
        send_subarray_size (:,2,1) = [array_size(1,2,1),mesh%lat_halo_width,num_lev(k)]
        recv_subarray_size (:,2,1) = send_subarray_size(:,2,1)
        send_subarray_start(:,2,1) = [0,mesh%lat_halo_width,0]
        recv_subarray_start(:,2,1) = [0,0,0]
        ! full_lon + half_lat
        send_subarray_size (:,1,2) = [array_size(1,1,2),mesh%lat_halo_width,num_lev(k)]
        recv_subarray_size (:,1,2) = send_subarray_size(:,1,2)
        send_subarray_start(:,1,2) = [0,mesh%lat_halo_width,0]
        recv_subarray_start(:,1,2) = [0,0,0]
        ! half_lon + half_lat
        send_subarray_size (:,2,2) = [array_size(1,2,2),mesh%lat_halo_width,num_lev(k)]
        recv_subarray_size (:,2,2) = send_subarray_size(:,2,2)
        send_subarray_start(:,2,2) = [0,mesh%lat_halo_width,0]
        recv_subarray_start(:,2,2) = [0,0,0]
      case (north)
        ! full_lon + full_lat
        send_subarray_size (:,1,1) = [array_size(1,1,1),mesh%lat_halo_width,num_lev(k)]
        recv_subarray_size (:,1,1) = send_subarray_size(:,1,1)
        send_subarray_start(:,1,1) = [0,mesh%num_full_lat,0]
        recv_subarray_start(:,1,1) = [0,mesh%num_full_lat+mesh%lat_halo_width,0]
        ! half_lon + full_lat
        send_subarray_size (:,2,1) = [array_size(1,2,1),mesh%lat_halo_width,num_lev(k)]
        recv_subarray_size (:,2,1) = send_subarray_size(:,2,1)
        send_subarray_start(:,2,1) = [0,mesh%num_full_lat,0]
        recv_subarray_start(:,2,1) = [0,mesh%num_full_lat+mesh%lat_halo_width,0]
        ! full_lon + half_lat
        send_subarray_size (:,1,2) = [array_size(1,1,2),mesh%lat_halo_width,num_lev(k)]
        recv_subarray_size (:,1,2) = send_subarray_size(:,1,2)
        send_subarray_start(:,1,2) = [0,mesh%num_half_lat,0]
        recv_subarray_start(:,1,2) = [0,mesh%num_half_lat+mesh%lat_halo_width,0]
        ! half_lon + half_lat
        send_subarray_size (:,2,2) = [array_size(1,2,2),mesh%lat_halo_width,num_lev(k)]
        recv_subarray_size (:,2,2) = send_subarray_size(:,2,2)
        send_subarray_start(:,2,2) = [0,mesh%num_half_lat,0]
        recv_subarray_start(:,2,2) = [0,mesh%num_half_lat+mesh%lat_halo_width,0]
      end select
      do j = 1, 2
        do i = 1, 2
          call MPI_TYPE_CREATE_SUBARRAY(3, array_size(:,i,j), send_subarray_size(:,i,j), &
                                        send_subarray_start(:,i,j), MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                        this%send_type_3d(i,j,k), ierr)
          call MPI_TYPE_COMMIT(this%send_type_3d(i,j,k), ierr)
          call MPI_TYPE_CREATE_SUBARRAY(3, array_size(:,i,j), recv_subarray_size(:,i,j), &
                                        recv_subarray_start(:,i,j), MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                        this%recv_type_3d(i,j,k), ierr)
          call MPI_TYPE_COMMIT(this%recv_type_3d(i,j,k), ierr)
        end do
      end do
    end do

    do j = 1, 2
      do i = 1, 2
        call MPI_TYPE_CREATE_SUBARRAY(2, array_size(1:2,i,j), send_subarray_size(1:2,i,j), &
                                      send_subarray_start(1:2,i,j), MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                      this%send_type_2d(i,j), ierr)
        call MPI_TYPE_COMMIT(this%send_type_2d(i,j), ierr)
        call MPI_TYPE_CREATE_SUBARRAY(2, array_size(1:2,i,j), recv_subarray_size(1:2,i,j), &
                                      recv_subarray_start(1:2,i,j), MPI_ORDER_FORTRAN, MPI_DOUBLE, &
                                      this%recv_type_2d(i,j), ierr)
        call MPI_TYPE_COMMIT(this%recv_type_2d(i,j), ierr)
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

    integer i, j, k
    integer ierr

    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          if (this%send_type_3d(i,j,k) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%send_type_3d(i,j,k), ierr)
          if (this%recv_type_3d(i,j,k) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%recv_type_3d(i,j,k), ierr)
        end do
      end do
    end do

    do j = 1, 2
      do i = 1, 2
        if (this%send_type_2d(i,j) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%send_type_2d(i,j), ierr)
        if (this%recv_type_2d(i,j) /= MPI_DATATYPE_NULL) call MPI_TYPE_FREE(this%recv_type_2d(i,j), ierr)
      end do
    end do

  end subroutine halo_final

end module halo_mod
