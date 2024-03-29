module parallel_mod

  use mpi
  use flogger
  use const_mod
  use mesh_mod
  use block_mod
  use process_mod

  implicit none

  private

  public proc
  public process_init
  public process_stop
  public process_final
  public is_root_proc
  public fill_zonal_halo
  public fill_halo
  public zero_halo
  public zonal_sum
  public zonal_max
  public global_sum
  public global_max
  public barrier
  public gather_zonal_array
  public scatter_zonal_array

  interface fill_zonal_halo
    module procedure fill_zonal_halo_1d_r4
    module procedure fill_zonal_halo_1d_r8
    module procedure fill_zonal_halo_2d_r4
    module procedure fill_zonal_halo_2d_r8
  end interface fill_zonal_halo

  interface fill_halo
    module procedure fill_halo_2d_r4
    module procedure fill_halo_2d_r8
    module procedure fill_halo_3d_r4
    module procedure fill_halo_3d_r8
  end interface fill_halo

  interface zero_halo
    module procedure zero_halo_1d_r4
    module procedure zero_halo_1d_r8
  end interface zero_halo

  interface zonal_sum
    module procedure zonal_sum_0d_i4
    module procedure zonal_sum_0d_r4
    module procedure zonal_sum_0d_r8
    module procedure zonal_sum_1d_r4
    module procedure zonal_sum_1d_r8
  end interface zonal_sum

  interface zonal_max
    module procedure zonal_max_r4
    module procedure zonal_max_r8
  end interface zonal_max

  interface global_sum
    module procedure global_sum_0d_r4
    module procedure global_sum_0d_r8
  end interface global_sum

  interface global_max
    module procedure global_max_0d_r4
    module procedure global_max_0d_r8
    module procedure global_max_0d_i4
  end interface global_max

  interface gather_zonal_array
    module procedure gather_zonal_array_1d_i4
    module procedure gather_zonal_array_1d_r4
    module procedure gather_zonal_array_1d_r8
    module procedure gather_zonal_array_2d_r4
    module procedure gather_zonal_array_2d_r8
  end interface gather_zonal_array

  interface scatter_zonal_array
    module procedure scatter_zonal_array_1d_r4
    module procedure scatter_zonal_array_1d_r8
    module procedure scatter_zonal_array_2d_r4
    module procedure scatter_zonal_array_2d_r8
  end interface scatter_zonal_array

contains

  subroutine fill_zonal_halo_1d_r4(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(4), intent(inout)  :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i1, i2, i3, i4, ierr

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                              n - 2 w + 1
      i1 = size(array) - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(i1:i2), halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 1, &
                        array(i3:i4), halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 1, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                             n - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = size(array) - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(i1:i2), halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 2, &
                        array(i3:i4), halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 2, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine fill_zonal_halo_1d_r4

  subroutine fill_zonal_halo_1d_r8(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i1, i2, i3, i4, ierr

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                              n - 2 w + 1
      i1 = size(array) - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(i1:i2), halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 1, &
                        array(i3:i4), halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 1, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                             n - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = size(array) - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(i1:i2), halo_width, block%halo(west)%dtype, block%halo(west)%proc_id, 2, &
                        array(i3:i4), halo_width, block%halo(east)%dtype, block%halo(east)%proc_id, 2, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine fill_zonal_halo_1d_r8

  subroutine fill_zonal_halo_2d_r4(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(4), intent(inout)  :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer nx, nz, i1, i2, i3, i4, ierr

    nz = size(array, 1)
    nx = size(array, 2)

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                             nx - 2 w + 1
      i1 = nx - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nz, block%halo(east)%dtype, block%halo(east)%proc_id, 1, &
                        array(:,i3:i4), halo_width * nz, block%halo(west)%dtype, block%halo(west)%proc_id, 1, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                            nx - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = nx - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nz, block%halo(west)%dtype, block%halo(west)%proc_id, 2, &
                        array(:,i3:i4), halo_width * nz, block%halo(east)%dtype, block%halo(east)%proc_id, 2, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine fill_zonal_halo_2d_r4

  subroutine fill_zonal_halo_2d_r8(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer nx, nz, i1, i2, i3, i4, ierr

    nz = size(array, 1)
    nx = size(array, 2)

    if (merge(west_halo, .true., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i3  | i4  |     |     |     |     | i1  | i2  |     |  n  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |                                   |
      !    1                             nx - 2 w + 1
      i1 = nx - 2 * halo_width + 1
      i2 = i1 + halo_width - 1
      i3 = 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nz, block%halo(east)%dtype, block%halo(east)%proc_id, 1, &
                        array(:,i3:i4), halo_width * nz, block%halo(west)%dtype, block%halo(west)%proc_id, 1, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if
    if (merge(east_halo, .true., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i1  | i2  |     |     |     |     | i3  | i4  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                            nx - w + 1
      i1 = 1 + halo_width
      i2 = i1 + halo_width - 1
      i3 = nx - halo_width + 1
      i4 = i3 + halo_width - 1
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nz, block%halo(west)%dtype, block%halo(west)%proc_id, 2, &
                        array(:,i3:i4), halo_width * nz, block%halo(east)%dtype, block%halo(east)%proc_id, 2, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine fill_zonal_halo_2d_r8

  subroutine fill_halo_2d_r4(block, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(4), intent(inout) :: array(:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer i, j, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_2d(i,j), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_2d(i,j), block%halo(west)%proc_id, 3, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_2d(i,j), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_2d(i,j), block%halo(east)%proc_id, 7, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 15, &
                        array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 15, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine fill_halo_2d_r4

  subroutine fill_halo_2d_r8(block, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer t1, t2, i, j, nx, ny, mx, hx, hy, ierr
    integer send_req, recv_req
    real(8) tmp(size(array,1),block%mesh%lat_halo_width)

    t1 = merge(1, 2, full_lon)
    t2 = merge(1, 2, full_lat)
    nx = size(array, 1)
    mx = size(array, 1) / 2
    ny = size(array, 2)
    hx = block%mesh%lon_halo_width
    hy = block%mesh%lat_halo_width

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_2d(t1,t2), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_2d(t1,t2), block%halo(west)%proc_id, 3, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_2d(t1,t2), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_2d(t1,t2), block%halo(east)%proc_id, 7, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(array, 1, block%halo(north)%send_type_2d(t1,t2), block%halo(north)%proc_id, 9, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(array, 1, block%halo(south)%recv_type_2d(t1,t2), block%halo(south)%proc_id, 9, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(array, 1, block%halo(south)%send_type_2d(t1,t2), block%halo(south)%proc_id, 10, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(array, 1, block%halo(north)%recv_type_2d(t1,t2), block%halo(north)%proc_id, 10, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo)) .and. proc%at_south_pole) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_2d(t1,t2), block%halo(south)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_2d(t1,t2), block%halo(south)%proc_id, 11, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:) = array(:,1:hy)
      if (block%halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,hy-j+1) = tmp(mx+1:nx-hx,j)
          array(mx+1:nx-hx,hy-j+1) = tmp(hx+1:mx   ,j)
        end do
      else
        do j = 1, hy
          array(:,hy-j+1) = tmp(:,j)
        end do
      end if
    end if

    if (merge(north_halo, .true., present(north_halo)) .and. proc%at_north_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_2d(t1,t2), block%halo(north)%proc_id, 12, &
                        array, 1, block%halo(north)%recv_type_2d(t1,t2), block%halo(north)%proc_id, 12, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:) = array(:,ny-hy+1:ny)
      if (block%halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,ny-hy+j) = tmp(mx+1:nx-hx,hy-j+1)
          array(mx+1:nx-hx,ny-hy+j) = tmp(hx+1:mx   ,hy-j+1)
        end do
      else
        do j = 1, hy
          array(:,ny-hy+j) = tmp(:,hy-j+1)
        end do
      end if
    end if

  end subroutine fill_halo_2d_r8

  subroutine fill_halo_3d_r4(block, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(4), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer i, j, k, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)
    k = merge(1, 2, merge(full_lev, .true., present(full_lev)))

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_3d(i,j,k), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_3d(i,j,k), block%halo(west)%proc_id, 3, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_3d(i,j,k), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_3d(i,j,k), block%halo(east)%proc_id, 7, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15, &
                        array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine fill_halo_3d_r4

  subroutine fill_halo_3d_r8(block, array, full_lon, full_lat, full_lev, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: full_lev
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer t1, t2, t3, i, j, nx, ny, mx, hx, hy, ierr
    integer send_req, recv_req
    real(8) tmp(size(array,1),block%mesh%lat_halo_width,size(array,3))

    t1 = merge(1, 2, full_lon)
    t2 = merge(1, 2, full_lat)
    t3 = merge(1, 2, merge(full_lev, .true., present(full_lev)))
    nx = size(array, 1)
    mx = size(array, 1) / 2
    ny = size(array, 2)
    hx = block%mesh%lon_halo_width
    hy = block%mesh%lat_halo_width

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_3d(t1,t2,t3), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_3d(t1,t2,t3), block%halo(west)%proc_id, 3, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_3d(t1,t2,t3), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_3d(t1,t2,t3), block%halo(east)%proc_id, 7, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_north_pole) then
        call MPI_ISEND(array, 1, block%halo(north)%send_type_3d(t1,t2,t3), block%halo(north)%proc_id, 9, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_south_pole) then
        call MPI_IRECV(array, 1, block%halo(south)%recv_type_3d(t1,t2,t3), block%halo(south)%proc_id, 9, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      send_req = MPI_REQUEST_NULL; recv_req = MPI_REQUEST_NULL
      if (.not. proc%at_south_pole) then
        call MPI_ISEND(array, 1, block%halo(south)%send_type_3d(t1,t2,t3), block%halo(south)%proc_id, 10, &
                       proc%comm, send_req, ierr)
      end if
      if (.not. proc%at_north_pole) then
        call MPI_IRECV(array, 1, block%halo(north)%recv_type_3d(t1,t2,t3), block%halo(north)%proc_id, 10, &
                       proc%comm, recv_req, ierr)
      end if
      call MPI_WAIT(send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(recv_req, MPI_STATUS_IGNORE, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo)) .and. proc%at_south_pole) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_3d(t1,t2,t3), block%halo(south)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_3d(t1,t2,t3), block%halo(south)%proc_id, 11, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      ! Reverse array order.
      tmp(:,:,:) = array(:,1:hy,:)
      if (block%halo(south)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,hy-j+1,:) = tmp(mx+1:nx-hx,j,:)
          array(mx+1:nx-hx,hy-j+1,:) = tmp(hx+1:mx   ,j,:)
        end do
      else
        do j = 1, hy
          array(:,hy-j+1,:) = tmp(:,j,:)
        end do
      end if
    end if

    if (merge(north_halo, .true., present(north_halo)) .and. proc%at_north_pole) then
      send_req = MPI_REQUEST_NULL; recv_req  = MPI_REQUEST_NULL
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_3d(t1,t2,t3), block%halo(north)%proc_id, 12, &
                        array, 1, block%halo(north)%recv_type_3d(t1,t2,t3), block%halo(north)%proc_id, 12, &
                        proc%comm, MPI_STATUS_IGNORE, ierr)
      tmp(:,:,:) = array(:,ny-hy+1:ny,:)
      if (block%halo(north)%proc_id == proc%id) then ! 1D decompostion, also reverse in lon
        do j = 1, hy
          array(hx+1:mx   ,ny-hy+j,:) = tmp(mx+1:nx-hx,hy-j+1,:)
          array(mx+1:nx-hx,ny-hy+j,:) = tmp(hx+1:mx   ,hy-j+1,:)
        end do
      else
        do j = 1, hy
          array(:,ny-hy+j,:) = tmp(:,hy-j+1,:)
        end do
      end if
    end if

  end subroutine fill_halo_3d_r8

  subroutine zero_halo_1d_r4(block, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    real(4), intent(inout) :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i1, i2

    if (merge(west_halo, .false., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i1  | i2  |     |     |     |     |     |     |     |     |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |     |
      !    1  1 + w - 1
      i1 = 1
      i2 = 1 + block%mesh%lon_halo_width - 1
      array(i1:i2) = 0.0d0
    end if

    if (merge(east_halo, .false., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     |     |     |     |     |     |     | i1  | i2  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                                                    |     |
      !                                                n - w + 1 n
      i1 = size(array) - block%mesh%lon_halo_width + 1
      i2 = size(array)
      array(i1:i2) = 0.0d0
    end if

  end subroutine zero_halo_1d_r4

  subroutine zero_halo_1d_r8(block, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i1, i2

    if (merge(west_halo, .false., present(west_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! | i1  | i2  |     |     |     |     |     |     |     |     |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !    |     |
      !    1  1 + w - 1
      i1 = 1
      i2 = 1 + block%mesh%lon_halo_width - 1
      array(i1:i2) = 0.0d0
    end if

    if (merge(east_halo, .false., present(east_halo))) then
      !   west halo |                                   | east_halo
      !  ___________|___________________________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     |     |     |     |     |     |     | i1  | i2  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                                                    |     |
      !                                                n - w + 1 n
      i1 = size(array) - block%mesh%lon_halo_width + 1
      i2 = size(array)
      array(i1:i2) = 0.0d0
    end if

  end subroutine zero_halo_1d_r8

  subroutine zonal_sum_0d_i4(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    integer, intent(in) :: work(:)
    integer, intent(out) :: value

#ifdef ENSURE_ORDER
    integer allvalue(global_mesh%num_full_lon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue)
      call MPI_BCAST(value, 1, MPI_INT, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work), value, 1, MPI_INT, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_i4

  subroutine zonal_sum_0d_r4(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: work(:)
    real(4), intent(out) :: value

#ifdef ENSURE_ORDER
    real(4) allvalue(global_mesh%num_full_lon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue)
      call MPI_BCAST(value, 1, MPI_REAL, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work), value, 1, MPI_REAL, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_r4

  subroutine zonal_sum_0d_r8(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: work(:)
    real(8), intent(out) :: value

#ifdef ENSURE_ORDER
    real(8) allvalue(global_mesh%num_full_lon)
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue)
      call MPI_BCAST(value, 1, MPI_DOUBLE, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work), value, 1, MPI_DOUBLE, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_0d_r8

  subroutine zonal_sum_1d_r4(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: work(:,:)
    real(4), intent(out) :: value(:)

#ifdef ENSURE_ORDER
    real(4) allvalue(global_mesh%num_full_lon,size(value))
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work, dim=1)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue, dim=1)
      call MPI_BCAST(value, size(value), MPI_REAL, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work, dim=1), value, size(value), MPI_REAL, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_1d_r4

  subroutine zonal_sum_1d_r8(zonal_circle, work, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: work(:,:)
    real(8), intent(out) :: value(:)

#ifdef ENSURE_ORDER
    real(8) allvalue(global_mesh%num_full_lon,size(value))
#endif
    integer ierr

#ifdef ENSURE_ORDER
    if (zonal_circle%np == 1) then
      value = sum(work, dim=1)
    else
      call gather_zonal_array(zonal_circle, work, allvalue)
      if (zonal_circle%id == 0) value = sum(allvalue, dim=1)
      call MPI_BCAST(value, size(value), MPI_DOUBLE, 0, zonal_circle%comm, ierr)
    end if
#else
    call MPI_ALLREDUCE(sum(work, dim=1), value, size(value), MPI_DOUBLE, MPI_SUM, zonal_circle%comm, ierr)
#endif

  end subroutine zonal_sum_1d_r8

  subroutine zonal_max_r4(zonal_circle, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(inout) :: value

    integer ierr
    real(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_MAX, zonal_circle%comm, ierr)
    value = res

  end subroutine zonal_max_r4

  subroutine zonal_max_r8(zonal_circle, value)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_MAX, zonal_circle%comm, ierr)
    value = res

  end subroutine zonal_max_r8

  subroutine global_sum_0d_r4(comm, value)

    integer, intent(in) :: comm
    real(4), intent(inout) :: value

    integer ierr
    real(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_0d_r4

  subroutine global_sum_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_0d_r8

  subroutine global_max_0d_r4(comm, value)

    integer, intent(in) :: comm
    real(4), intent(inout) :: value

    integer ierr
    real(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_REAL, MPI_MAX, comm, ierr)
    value = res

  end subroutine global_max_0d_r4

  subroutine global_max_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_MAX, comm, ierr)
    value = res

  end subroutine global_max_0d_r8

  subroutine global_max_0d_i4(comm, value)

    integer, intent(in) :: comm
    integer(4), intent(inout) :: value

    integer ierr
    integer(4) res

    call MPI_ALLREDUCE(value, res, 1, MPI_INT, MPI_MAX, comm, ierr)
    value = res

  end subroutine global_max_0d_i4

  subroutine barrier()

    integer ierr

    call MPI_BARRIER(proc%comm, ierr)

  end subroutine barrier

  subroutine gather_zonal_array_1d_i4(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    integer, intent(in) :: local_array(:)
    integer, intent(out) :: array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      array(1:size(local_array)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_i4(i,0), i - 1, 30, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_INT, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_i4

  subroutine gather_zonal_array_1d_r4(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: local_array(:)
    real(4), intent(out) :: array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      array(1:size(local_array)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r4(i,0), i - 1, 30, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_REAL, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_r4

  subroutine gather_zonal_array_1d_r8(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: local_array(:)
    real(8), intent(out) :: array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      array(1:size(local_array)) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r8(i,0), i - 1, 30, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_DOUBLE, 0, 30, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_1d_r8

  subroutine gather_zonal_array_2d_r4(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: local_array(:,:)
    real(4), intent(inout) :: array(:,:)

    integer ierr, i, k

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      array(1:size(local_array, 1),:) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r4(i,k), i - 1, 31, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_REAL, 0, 31, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_2d_r4

  subroutine gather_zonal_array_2d_r8(zonal_circle, local_array, array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: local_array(:,:)
    real(8), intent(inout) :: array(:,:)

    integer ierr, i, k

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      array(1:size(local_array, 1),:) = local_array
      do i = 2, zonal_circle%np
        call MPI_RECV(array, 1, zonal_circle%recv_type_r8(i,k), i - 1, 31, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
      end do
    else
      call MPI_SEND(local_array, size(local_array), MPI_DOUBLE, 0, 31, zonal_circle%comm, ierr)
    end if

  end subroutine gather_zonal_array_2d_r8

  subroutine scatter_zonal_array_1d_r4(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: array(:)
    real(4), intent(out) :: local_array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      local_array = array(1:size(local_array))
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r4(i,0), i - 1, 32, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_REAL, 0, 32, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine scatter_zonal_array_1d_r4

  subroutine scatter_zonal_array_1d_r8(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: array(:)
    real(8), intent(out) :: local_array(:)

    integer ierr, i

    if (zonal_circle%id == 0) then
      local_array = array(1:size(local_array))
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r8(i,0), i - 1, 32, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_DOUBLE, 0, 32, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine scatter_zonal_array_1d_r8

  subroutine scatter_zonal_array_2d_r4(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(4), intent(in) :: array(:,:)
    real(4), intent(out) :: local_array(:,:)

    integer ierr, i, k

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      local_array = array(1:size(local_array, 1),:)
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r4(i,k), i - 1, 33, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_REAL, 0, 33, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine scatter_zonal_array_2d_r4

  subroutine scatter_zonal_array_2d_r8(zonal_circle, array, local_array)

    type(zonal_circle_type), intent(in) :: zonal_circle
    real(8), intent(in) :: array(:,:)
    real(8), intent(out) :: local_array(:,:)

    integer ierr, i, k

    if (zonal_circle%id == 0) then
      k = merge(1, 2, size(local_array, 2) == global_mesh%num_full_lev)
      local_array = array(1:size(local_array, 1),:)
      do i = 2, zonal_circle%np
        call MPI_SEND(array, 1, zonal_circle%recv_type_r8(i,k), i - 1, 33, zonal_circle%comm, ierr)
      end do
    else
      call MPI_RECV(local_array, size(local_array), MPI_DOUBLE, 0, 33, zonal_circle%comm, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine scatter_zonal_array_2d_r8

end module parallel_mod
