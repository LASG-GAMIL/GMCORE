module parallel_mod

  use mpi
  use flogger
  use const_mod
  use mesh_mod
  use block_mod
  use process_mod
  use parallel_types_mod

  implicit none

  private

  public proc
  public process_init
  public process_stop
  public process_final
  public is_root_proc
  public async_type
  public async_v
  public async_u
  public async_gz
  public async_pv
  public async_pv_lon
  public async_pv_lat
  public async_ke
  public async_mf_lon_n
  public async_mf_lat_n
  public async_dpv_lon_t
  public async_dpv_lat_t
  public fill_zonal_halo
  public fill_halo
  public wait_halo
  public zero_halo
  public zonal_sum
  public global_sum
  public overlay_inner_halo
  public barrier

  interface fill_zonal_halo
    module procedure fill_zonal_halo_1d_r8
    module procedure fill_zonal_halo_1d_r8_async
    module procedure fill_zonal_halo_2d_r8
    module procedure fill_zonal_halo_2d_r8_async
  end interface fill_zonal_halo

  interface fill_halo
    module procedure fill_halo_2d_r8
    module procedure fill_halo_2d_r8_async
    module procedure fill_halo_3d_r8
    module procedure fill_halo_3d_r8_async
  end interface fill_halo

  interface zero_halo
    module procedure zero_halo_1d_r8
  end interface zero_halo

  interface zonal_sum
    module procedure zonal_sum_0d_r8
    module procedure zonal_sum_1d_r8
  end interface zonal_sum

  interface global_sum
    module procedure global_sum_0d_r8
  end interface global_sum

contains

  subroutine fill_zonal_halo_1d_r8(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer status(MPI_STATUS_SIZE), ierr
    integer i1, i2, i3, i4

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
      call MPI_SENDRECV(array(i1:i2), halo_width, MPI_DOUBLE, block%halo(east)%proc_id, 1, &
                        array(i3:i4), halo_width, MPI_DOUBLE, block%halo(west)%proc_id, 1, &
                        proc%comm, status, ierr)
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
      call MPI_SENDRECV(array(i1:i2), halo_width, MPI_DOUBLE, block%halo(west)%proc_id, 2, &
                        array(i3:i4), halo_width, MPI_DOUBLE, block%halo(east)%proc_id, 2, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_zonal_halo_1d_r8

  subroutine fill_zonal_halo_1d_r8_async(block, halo_width, array, async, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:)
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer ierr
    integer i1, i2, i3, i4

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
      call MPI_ISEND(array(i1:i2), halo_width, MPI_DOUBLE, block%halo(east)%proc_id, 1, proc%comm, async%west_send_req, ierr)
      call MPI_IRECV(array(i3:i4), halo_width, MPI_DOUBLE, block%halo(west)%proc_id, 1, proc%comm, async%west_recv_req, ierr)
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
      call MPI_ISEND(array(i1:i2), halo_width, MPI_DOUBLE, block%halo(west)%proc_id, 2, proc%comm, async%east_send_req, ierr)
      call MPI_IRECV(array(i3:i4), halo_width, MPI_DOUBLE, block%halo(east)%proc_id, 2, proc%comm, async%east_recv_req, ierr)
    end if

  end subroutine fill_zonal_halo_1d_r8_async

  subroutine fill_zonal_halo_2d_r8(block, halo_width, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer status(MPI_STATUS_SIZE), ierr
    integer nx, nz, i1, i2, i3, i4

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
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nz, MPI_DOUBLE, block%halo(east)%proc_id, 1, &
                        array(:,i3:i4), halo_width * nz, MPI_DOUBLE, block%halo(west)%proc_id, 1, &
                        proc%comm, status, ierr)
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
      call MPI_SENDRECV(array(:,i1:i2), halo_width * nz, MPI_DOUBLE, block%halo(west)%proc_id, 2, &
                        array(:,i3:i4), halo_width * nz, MPI_DOUBLE, block%halo(east)%proc_id, 2, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_zonal_halo_2d_r8

  subroutine fill_zonal_halo_2d_r8_async(block, halo_width, array, async, west_halo, east_halo)

    type(block_type), intent(in) :: block
    integer, intent(in) :: halo_width
    real(8), intent(inout)  :: array(:,:)
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer ierr
    integer nx, nz, i1, i2, i3, i4

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
      call MPI_ISEND(array(:,i1:i2), halo_width * nz, MPI_DOUBLE, block%halo(east)%proc_id, 1, proc%comm, async%west_send_req, ierr)
      call MPI_IRECV(array(:,i3:i4), halo_width * nz, MPI_DOUBLE, block%halo(west)%proc_id, 1, proc%comm, async%west_recv_req, ierr)
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
      call MPI_ISEND(array(:,i1:i2), halo_width * nz, MPI_DOUBLE, block%halo(west)%proc_id, 2, proc%comm, async%east_send_req, ierr)
      call MPI_IRECV(array(:,i3:i4), halo_width * nz, MPI_DOUBLE, block%halo(east)%proc_id, 2, proc%comm, async%east_recv_req, ierr)
    end if

  end subroutine fill_zonal_halo_2d_r8_async

  subroutine fill_halo_2d_r8(block, array, full_lon, full_lat, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer status(MPI_STATUS_SIZE), i, j, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_2d(i,j), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_2d(i,j), block%halo(west)%proc_id, 3, &
                        proc%comm, status, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_2d(i,j), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_2d(i,j), block%halo(east)%proc_id, 7, &
                        proc%comm, status, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11, &
                        proc%comm, status, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 15, &
                        array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 15, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_halo_2d_r8

  subroutine fill_halo_2d_r8_async(block, array, full_lon, full_lat, async, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    type(async_type), intent(inout) :: async
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: south_halo
    logical, intent(in), optional :: north_halo

    integer i, j, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)

    if (merge(west_halo, .true., present(west_halo))) then
      if (async%west_send_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%west_send_req, MPI_STATUS_IGNORE, ierr)
      if (async%west_recv_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%west_recv_req, MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(east)%send_type_2d(i,j), block%halo(east)%proc_id, 3, proc%comm, async%west_send_req, ierr)
      call MPI_IRECV(array, 1, block%halo(west)%recv_type_2d(i,j), block%halo(west)%proc_id, 3, proc%comm, async%west_recv_req, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      if (async%east_send_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%east_send_req, MPI_STATUS_IGNORE, ierr)
      if (async%east_recv_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%east_recv_req, MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(west)%send_type_2d(i,j), block%halo(west)%proc_id, 7, proc%comm, async%east_send_req, ierr)
      call MPI_IRECV(array, 1, block%halo(east)%recv_type_2d(i,j), block%halo(east)%proc_id, 7, proc%comm, async%east_recv_req, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      if (async%south_send_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%south_send_req, MPI_STATUS_IGNORE, ierr)
      if (async%south_recv_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%south_recv_req, MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(north)%send_type_2d(i,j), block%halo(north)%proc_id, 11, proc%comm, async%south_send_req, ierr)
      call MPI_IRECV(array, 1, block%halo(south)%recv_type_2d(i,j), block%halo(south)%proc_id, 11, proc%comm, async%south_recv_req, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      if (async%north_send_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%north_send_req, MPI_STATUS_IGNORE, ierr)
      if (async%north_recv_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%north_recv_req, MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(south)%send_type_2d(i,j), block%halo(south)%proc_id, 15, proc%comm, async%north_send_req, ierr)
      call MPI_IRECV(array, 1, block%halo(north)%recv_type_2d(i,j), block%halo(north)%proc_id, 15, proc%comm, async%north_recv_req, ierr)
    end if

  end subroutine fill_halo_2d_r8_async

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

    integer status(MPI_STATUS_SIZE), i, j, k, ierr

    i = merge(1, 2, full_lon)
    j = merge(1, 2, full_lat)
    k = merge(1, 2, merge(full_lev, .true., present(full_lev)))

    if (merge(west_halo, .true., present(west_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(east)%send_type_3d(i,j,k), block%halo(east)%proc_id, 3, &
                        array, 1, block%halo(west)%recv_type_3d(i,j,k), block%halo(west)%proc_id, 3, &
                        proc%comm, status, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(west)%send_type_3d(i,j,k), block%halo(west)%proc_id, 7, &
                        array, 1, block%halo(east)%recv_type_3d(i,j,k), block%halo(east)%proc_id, 7, &
                        proc%comm, status, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, &
                        array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, &
                        proc%comm, status, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      call MPI_SENDRECV(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15, &
                        array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15, &
                        proc%comm, status, ierr)
    end if

  end subroutine fill_halo_3d_r8

  subroutine fill_halo_3d_r8_async(block, array, full_lon, full_lat, async, full_lev, west_halo, east_halo, south_halo, north_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:,:)
    logical, intent(in) :: full_lon
    logical, intent(in) :: full_lat
    type(async_type), intent(inout) :: async
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
      if (async%west_send_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%west_send_req, MPI_STATUS_IGNORE, ierr)
      if (async%west_recv_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%west_recv_req, MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(east)%send_type_3d(i,j,k), block%halo(east)%proc_id, 3, proc%comm, async%west_send_req, ierr)
      call MPI_IRECV(array, 1, block%halo(west)%recv_type_3d(i,j,k), block%halo(west)%proc_id, 3, proc%comm, async%west_recv_req, ierr)
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      if (async%east_send_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%east_send_req, MPI_STATUS_IGNORE, ierr)
      if (async%east_recv_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%east_recv_req, MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(west)%send_type_3d(i,j,k), block%halo(west)%proc_id, 7, proc%comm, async%east_send_req, ierr)
      call MPI_IRECV(array, 1, block%halo(east)%recv_type_3d(i,j,k), block%halo(east)%proc_id, 7, proc%comm, async%east_recv_req, ierr)
    end if

    if (merge(south_halo, .true., present(south_halo))) then
      if (async%south_send_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%south_send_req, MPI_STATUS_IGNORE, ierr)
      if (async%south_recv_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%south_recv_req, MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(north)%send_type_3d(i,j,k), block%halo(north)%proc_id, 11, proc%comm, async%south_send_req, ierr)
      call MPI_IRECV(array, 1, block%halo(south)%recv_type_3d(i,j,k), block%halo(south)%proc_id, 11, proc%comm, async%south_recv_req, ierr)
    end if

    if (merge(north_halo, .true., present(north_halo))) then
      if (async%north_send_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%north_send_req, MPI_STATUS_IGNORE, ierr)
      if (async%north_recv_req /= MPI_REQUEST_NULL) call MPI_WAIT(async%north_recv_req, MPI_STATUS_IGNORE, ierr)
      call MPI_ISEND(array, 1, block%halo(south)%send_type_3d(i,j,k), block%halo(south)%proc_id, 15, proc%comm, async%north_send_req, ierr)
      call MPI_IRECV(array, 1, block%halo(north)%recv_type_3d(i,j,k), block%halo(north)%proc_id, 15, proc%comm, async%north_recv_req, ierr)
    end if

  end subroutine fill_halo_3d_r8_async

  subroutine wait_halo(async)

    type(async_type), intent(inout) :: async

    integer ierr

    if (async%west_send_req  /= MPI_REQUEST_NULL) then
      call MPI_WAIT(async%west_send_req , MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(async%west_recv_req , MPI_STATUS_IGNORE, ierr)
    end if
    if (async%east_send_req  /= MPI_REQUEST_NULL) then
      call MPI_WAIT(async%east_send_req , MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(async%east_recv_req , MPI_STATUS_IGNORE, ierr)
    end if
    if (async%south_send_req /= MPI_REQUEST_NULL) then
      call MPI_WAIT(async%south_send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(async%south_recv_req, MPI_STATUS_IGNORE, ierr)
    end if
    if (async%north_send_req /= MPI_REQUEST_NULL) then
      call MPI_WAIT(async%north_send_req, MPI_STATUS_IGNORE, ierr)
      call MPI_WAIT(async%north_recv_req, MPI_STATUS_IGNORE, ierr)
    end if

  end subroutine wait_halo

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

  subroutine overlay_inner_halo(block, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i1, i2, i3, i4
    integer status(MPI_STATUS_SIZE), ierr
    real(8) buffer(block%mesh%lon_halo_width)

    if (merge(west_halo, .false., present(west_halo))) then
      !   west halo |           |                       | east_halo
      !  ___________|___________|_______________________|___________
      ! |     |     |     |     |     |     |     |     |     |     |
      ! |     |     | i3  | i4  |     |     |     |     | i1  | i2  |
      ! |_____|_____|_____|_____|_____|_____|_____|_____|_____|_____|
      !                |                                   |
      !              1 + w                             n - w + 1
      i1 = size(array) - block%mesh%lon_halo_width + 1
      i2 = i1 + block%mesh%lon_halo_width - 1
      i3 =  1 + block%mesh%lon_halo_width
      i4 = i3 + block%mesh%lon_halo_width - 1
      call MPI_SENDRECV(array(i1:i2), block%mesh%lon_halo_width, MPI_DOUBLE, block%halo(east)%proc_id, 19, &
                        buffer      , block%mesh%lon_halo_width, MPI_DOUBLE, block%halo(west)%proc_id, 19, &
                        proc%comm, status, ierr)
      array(i3:i4) = array(i3:i4) + buffer
    end if

  end subroutine overlay_inner_halo

  subroutine zonal_sum_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine zonal_sum_0d_r8

  subroutine zonal_sum_1d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value(:)

    integer ierr
    real(8) res(size(value))

    call MPI_ALLREDUCE(value, res, size(value), MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine zonal_sum_1d_r8

  subroutine global_sum_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_0d_r8

  subroutine barrier()

    integer ierr

    call MPI_BARRIER(proc%comm, ierr)

  end subroutine barrier

end module parallel_mod
