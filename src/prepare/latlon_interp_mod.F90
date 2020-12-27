module latlon_interp_mod

  use flogger
  use const_mod
  use mesh_mod

  implicit none

  private

  public latlon_interp_bilinear_cell
  public latlon_interp_bilinear_lon_edge
  public latlon_interp_bilinear_lat_edge

contains

  subroutine latlon_interp_bilinear_cell(src_lon, src_lat, src_data, dst_mesh, dst_data, zero_pole, ierr)

    real(r8), intent(in) :: src_lon(:)
    real(r8), intent(in) :: src_lat(:)
    real(r8), intent(in) :: src_data(:,:)
    type(mesh_type), intent(in) :: dst_mesh
    real(r8), intent(out) :: dst_data(dst_mesh%full_lon_lb:dst_mesh%full_lon_ub, &
                                      dst_mesh%full_lat_lb:dst_mesh%full_lat_ub)
    logical, intent(in), optional :: zero_pole
    integer, intent(out), optional :: ierr

    integer nx1, ny1
    integer lb_x2, ub_x2, lb_y2, ub_y2
    integer i, ii, j, jj
    integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
    real(r8), allocatable :: xwgt1(:), xwgt2(:), ywgt1(:), ywgt2(:)
    real(r8) tmp1, tmp2
    logical is_found

    nx1 = size(src_lon)
    ny1 = size(src_lat)

    lb_x2 = dst_mesh%full_lon_ibeg
    ub_x2 = dst_mesh%full_lon_iend
    lb_y2 = dst_mesh%full_lat_ibeg
    ub_y2 = dst_mesh%full_lat_iend

    allocate(i1(lb_x2:ub_x2), i2(lb_x2:ub_x2))
    allocate(j1(lb_y2:ub_y2), j2(lb_y2:ub_y2))
    allocate(xwgt1(lb_x2:ub_x2), xwgt2(lb_x2:ub_x2))
    allocate(ywgt1(lb_y2:ub_y2), ywgt2(lb_y2:ub_y2))

    associate (x1 => src_lon, y1 => src_lat, &
               x2 => dst_mesh%full_lon_deg(lb_x2:ub_x2), &
               y2 => dst_mesh%full_lat_deg(lb_y2:ub_y2))
      do i = lb_x2, ub_x2
        is_found = .false.
        do ii = 1, nx1-1
          if (x2(i) >= x1(ii) .and. x2(i) < x1(ii+1)) then
            i1(i) = ii
            i2(i) = ii+1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i))
            is_found = .true.
            exit
          else if (x2(i) >= x1(nx1)) then
            i1(i) = nx1
            i2(i) = 1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i)) + 360.0d0
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Longitude mismatch!')
          end if
        end if
        xwgt1(i) = (tmp2 - x2(i)) / (tmp2 - tmp1)
        xwgt2(i) = (x2(i) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        is_found = .false.
        do jj = 1, ny1-1
          if ((y2(j) >= y1(jj) .and. y2(j) < y1(jj+1)) .or. &
            (jj == ny1 - 1 .and. (y2(j) <= y1(ny1) .or. abs(y2(j) - y1(ny1)) < 1.0e-6))) then
            j1(j) = jj
            j2(j) = jj + 1
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Latitude mismatch!')
          end if
        end if
        tmp1 = y1(j1(j))
        tmp2 = y1(j2(j))
        ywgt1(j) = (tmp2 - y2(j)) / (tmp2 - tmp1)
        ywgt2(j) = (y2(j) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        do i = lb_x2, ub_x2
          tmp1 = src_data(i1(i),j1(j)) * xwgt1(i) + src_data(i2(i),j1(j)) * xwgt2(i)
          tmp2 = src_data(i1(i),j2(j)) * xwgt1(i) + src_data(i2(i),j2(j)) * xwgt2(i)
          dst_data(i,j) = tmp1 * ywgt1(j) + tmp2 * ywgt2(j)
        end do
      end do
    end associate

    deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)

    ! Handle pole grids.
    if (dst_mesh%has_south_pole()) then
      if (merge(zero_pole, .false., present(zero_pole))) then
        dst_data(:,dst_mesh%full_lat_ibeg) = 0.0_r8
      else
        dst_data(:,dst_mesh%full_lat_ibeg) = sum(dst_data(:,dst_mesh%full_lat_ibeg)) / dst_mesh%num_full_lon
      end if
    end if
    if (dst_mesh%has_north_pole()) then
      if (merge(zero_pole, .false., present(zero_pole))) then
        dst_data(:,dst_mesh%full_lat_iend) = 0.0_r8
      else
        dst_data(:,dst_mesh%full_lat_iend) = sum(dst_data(:,dst_mesh%full_lat_iend)) / dst_mesh%num_full_lon
      end if
    end if

  end subroutine latlon_interp_bilinear_cell

  subroutine latlon_interp_bilinear_lon_edge(src_lon, src_lat, src_data, dst_mesh, dst_data, zero_pole, ierr)

    real(r8), intent(in) :: src_lon(:)
    real(r8), intent(in) :: src_lat(:)
    real(r8), intent(in) :: src_data(:,:)
    type(mesh_type), intent(in) :: dst_mesh
    real(r8), intent(out) :: dst_data(dst_mesh%half_lon_lb:dst_mesh%half_lon_ub, &
                                      dst_mesh%full_lat_lb:dst_mesh%full_lat_ub)
    logical, intent(in), optional :: zero_pole
    integer, intent(out), optional :: ierr

    integer nx1, ny1
    integer lb_x2, ub_x2, lb_y2, ub_y2
    integer i, ii, j, jj
    integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
    real(r8), allocatable :: xwgt1(:), xwgt2(:), ywgt1(:), ywgt2(:)
    real(r8) tmp1, tmp2
    logical is_found

    nx1 = size(src_lon)
    ny1 = size(src_lat)

    lb_x2 = dst_mesh%half_lon_ibeg
    ub_x2 = dst_mesh%half_lon_iend
    lb_y2 = dst_mesh%full_lat_ibeg
    ub_y2 = dst_mesh%full_lat_iend

    allocate(i1(lb_x2:ub_x2), i2(lb_x2:ub_x2))
    allocate(j1(lb_y2:ub_y2), j2(lb_y2:ub_y2))
    allocate(xwgt1(lb_x2:ub_x2), xwgt2(lb_x2:ub_x2))
    allocate(ywgt1(lb_y2:ub_y2), ywgt2(lb_y2:ub_y2))

    associate (x1 => src_lon, y1 => src_lat, &
               x2 => dst_mesh%half_lon_deg(lb_x2:ub_x2), &
               y2 => dst_mesh%full_lat_deg(lb_y2:ub_y2))
      do i = lb_x2, ub_x2
        is_found = .false.
        do ii = 1, nx1-1
          if (x2(i) >= x1(ii) .and. x2(i) < x1(ii+1)) then
            i1(i) = ii
            i2(i) = ii+1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i))
            is_found = .true.
            exit
          else if (x2(i) >= x1(nx1)) then
            i1(i) = nx1
            i2(i) = 1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i)) + 360.0d0
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Longitude mismatch!')
          end if
        end if
        xwgt1(i) = (tmp2 - x2(i)) / (tmp2 - tmp1)
        xwgt2(i) = (x2(i) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        is_found = .false.
        do jj = 1, ny1-1
          if ((y2(j) >= y1(jj) .and. y2(j) < y1(jj+1)) .or. &
            (jj == ny1 - 1 .and. (y2(j) <= y1(ny1) .or. abs(y2(j) - y1(ny1)) < 1.0e-6))) then
            j1(j) = jj
            j2(j) = jj + 1
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Latitude mismatch!')
          end if
        end if
        tmp1 = y1(j1(j))
        tmp2 = y1(j2(j))
        ywgt1(j) = (tmp2 - y2(j)) / (tmp2 - tmp1)
        ywgt2(j) = (y2(j) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        do i = lb_x2, ub_x2
          tmp1 = src_data(i1(i),j1(j)) * xwgt1(i) + src_data(i2(i),j1(j)) * xwgt2(i)
          tmp2 = src_data(i1(i),j2(j)) * xwgt1(i) + src_data(i2(i),j2(j)) * xwgt2(i)
          dst_data(i,j) = tmp1 * ywgt1(j) + tmp2 * ywgt2(j)
        end do
      end do
    end associate

    deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)

    ! Handle pole grids.
    if (dst_mesh%has_south_pole()) then
      if (merge(zero_pole, .false., present(zero_pole))) then
        dst_data(:,dst_mesh%full_lat_ibeg) = 0.0_r8
      else
        dst_data(:,dst_mesh%full_lat_ibeg) = sum(dst_data(:,dst_mesh%full_lat_ibeg)) / dst_mesh%num_full_lon
      end if
    end if
    if (dst_mesh%has_north_pole()) then
      if (merge(zero_pole, .false., present(zero_pole))) then
        dst_data(:,dst_mesh%full_lat_iend) = 0.0_r8
      else
        dst_data(:,dst_mesh%full_lat_iend) = sum(dst_data(:,dst_mesh%full_lat_iend)) / dst_mesh%num_full_lon
      end if
    end if

  end subroutine latlon_interp_bilinear_lon_edge

  subroutine latlon_interp_bilinear_lat_edge(src_lon, src_lat, src_data, dst_mesh, dst_data, zero_pole, ierr)

    real(r8), intent(in) :: src_lon(:)
    real(r8), intent(in) :: src_lat(:)
    real(r8), intent(in) :: src_data(:,:)
    type(mesh_type), intent(in) :: dst_mesh
    real(r8), intent(out) :: dst_data(dst_mesh%full_lon_lb:dst_mesh%full_lon_ub, &
                                      dst_mesh%half_lat_lb:dst_mesh%half_lat_ub)
    logical, intent(in), optional :: zero_pole
    integer, intent(out), optional :: ierr

    integer nx1, ny1
    integer lb_x2, ub_x2, lb_y2, ub_y2
    integer i, ii, j, jj
    integer, allocatable :: i1(:), j1(:), i2(:), j2(:)
    real(r8), allocatable :: xwgt1(:), xwgt2(:), ywgt1(:), ywgt2(:)
    real(r8) tmp1, tmp2
    logical is_found

    nx1 = size(src_lon)
    ny1 = size(src_lat)

    lb_x2 = dst_mesh%full_lon_ibeg
    ub_x2 = dst_mesh%full_lon_iend
    lb_y2 = dst_mesh%half_lat_ibeg
    ub_y2 = dst_mesh%half_lat_iend

    allocate(i1(lb_x2:ub_x2), i2(lb_x2:ub_x2))
    allocate(j1(lb_y2:ub_y2), j2(lb_y2:ub_y2))
    allocate(xwgt1(lb_x2:ub_x2), xwgt2(lb_x2:ub_x2))
    allocate(ywgt1(lb_y2:ub_y2), ywgt2(lb_y2:ub_y2))

    associate (x1 => src_lon, y1 => src_lat, &
               x2 => dst_mesh%full_lon_deg(lb_x2:ub_x2), &
               y2 => dst_mesh%half_lat_deg(lb_y2:ub_y2))
      do i = lb_x2, ub_x2
        is_found = .false.
        do ii = 1, nx1-1
          if (x2(i) >= x1(ii) .and. x2(i) < x1(ii+1)) then
            i1(i) = ii
            i2(i) = ii+1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i))
            is_found = .true.
            exit
          else if (x2(i) >= x1(nx1)) then
            i1(i) = nx1
            i2(i) = 1
            tmp1 = x1(i1(i))
            tmp2 = x1(i2(i)) + 360.0d0
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Longitude mismatch!')
          end if
        end if
        xwgt1(i) = (tmp2 - x2(i)) / (tmp2 - tmp1)
        xwgt2(i) = (x2(i) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        is_found = .false.
        do jj = 1, ny1-1
          if ((y2(j) >= y1(jj) .and. y2(j) < y1(jj+1)) .or. &
            (jj == ny1 - 1 .and. (y2(j) <= y1(ny1) .or. abs(y2(j) - y1(ny1)) < 1.0e-6))) then
            j1(j) = jj
            j2(j) = jj + 1
            is_found = .true.
            exit
          end if
        end do
        if (.not. is_found) then
          if (present(ierr)) then
            ierr = -1
            deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)
            return
          else
            call log_error('latlon_interp_bilinear: Latitude mismatch!')
          end if
        end if
        tmp1 = y1(j1(j))
        tmp2 = y1(j2(j))
        ywgt1(j) = (tmp2 - y2(j)) / (tmp2 - tmp1)
        ywgt2(j) = (y2(j) - tmp1) / (tmp2 - tmp1)
      end do

      do j = lb_y2, ub_y2
        do i = lb_x2, ub_x2
          tmp1 = src_data(i1(i),j1(j)) * xwgt1(i) + src_data(i2(i),j1(j)) * xwgt2(i)
          tmp2 = src_data(i1(i),j2(j)) * xwgt1(i) + src_data(i2(i),j2(j)) * xwgt2(i)
          dst_data(i,j) = tmp1 * ywgt1(j) + tmp2 * ywgt2(j)
        end do
      end do
    end associate

    deallocate(i1, i2, j1, j2, xwgt1, xwgt2, ywgt1, ywgt2)

    ! Handle pole grids.
    if (dst_mesh%has_south_pole()) then
      if (merge(zero_pole, .false., present(zero_pole))) then
        dst_data(:,dst_mesh%full_lat_ibeg) = 0.0_r8
      else
        dst_data(:,dst_mesh%full_lat_ibeg) = sum(dst_data(:,dst_mesh%full_lat_ibeg)) / dst_mesh%num_full_lon
      end if
    end if
    if (dst_mesh%has_north_pole()) then
      if (merge(zero_pole, .false., present(zero_pole))) then
        dst_data(:,dst_mesh%full_lat_iend) = 0.0_r8
      else
        dst_data(:,dst_mesh%full_lat_iend) = sum(dst_data(:,dst_mesh%full_lat_iend)) / dst_mesh%num_full_lon
      end if
    end if

  end subroutine latlon_interp_bilinear_lat_edge

end module latlon_interp_mod
