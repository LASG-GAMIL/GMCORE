module filter_mod

  use const_mod
  use block_mod
  use parallel_mod
  use filter_types_mod

  implicit none

  private

  public filter_type
  public filter_on_cell
  public filter_on_lon_edge
  public filter_on_lat_edge
  public filter_on_lev_edge

  interface filter_on_cell
    module procedure filter_on_cell_2d
    module procedure filter_on_cell_3d
  end interface filter_on_cell

contains

  subroutine filter_on_cell_2d(block, x, y)

    type(block_type), intent(in) :: block
    real(r8), intent(inout) :: x(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub)
    real(r8), intent(out), optional :: y(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                         block%mesh%full_lat_lb:block%mesh%full_lat_ub)

    real(r8) tmp(block%mesh%full_lon_lb:block%mesh%full_lon_ub)
    integer i, j, n, hn

    associate (mesh => block%mesh, filter => block%filter)
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      if (filter%ngrid_lon(j) > 1) then
        n  = filter%ngrid_lon(j)
        hn = (n - 1) / 2
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tmp(i) = sum(filter%wgt_lon(:n,j) * x(i-hn:i+hn,j))
        end do
        if (present(y)) then
          y(:,j) = tmp
        else
          x(:,j) = tmp
        end if
      else if (present(y)) then
        y(:,j) = x(:,j)
      end if
    end do
    end associate

  end subroutine filter_on_cell_2d

  subroutine filter_on_cell_3d(block, x, y)

    type(block_type), intent(in) :: block
    real(r8), intent(inout) :: x(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out), optional :: y(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                         block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                         block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    real(r8) tmp(block%mesh%full_lon_lb:block%mesh%full_lon_ub)
    integer i, j, k, n, hn

    associate (mesh => block%mesh, filter => block%filter)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        if (filter%ngrid_lon(j) > 1) then
          n  = filter%ngrid_lon(j)
          hn = (n - 1) / 2
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tmp(i) = sum(filter%wgt_lon(:n,j) * x(i-hn:i+hn,j,k))
          end do
          if (present(y)) then
            y(:,j,k) = tmp
          else
            x(:,j,k) = tmp
          end if
        else if (present(y)) then
          y(:,j,k) = x(:,j,k)
        end if
      end do
    end do
    end associate

  end subroutine filter_on_cell_3d

  subroutine filter_on_lon_edge(block, x, y)

    type(block_type), intent(in) :: block
    real(r8), intent(inout) :: x(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out), optional :: y(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                         block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                         block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    real(r8) tmp(block%mesh%half_lon_lb:block%mesh%half_lon_ub)
    integer i, j, k, n, hn

    associate (mesh => block%mesh, filter => block%filter)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        if (filter%ngrid_lon(j) > 1) then
          n  = filter%ngrid_lon(j)
          hn = (n - 1) / 2
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            tmp(i) = sum(filter%wgt_lon(:n,j) * x(i-hn:i+hn,j,k))
          end do
          if (present(y)) then
            y(:,j,k) = tmp
          else
            x(:,j,k) = tmp
          end if
        else if (present(y)) then
          y(:,j,k) = x(:,j,k)
        end if
      end do
    end do
    end associate

  end subroutine filter_on_lon_edge

  subroutine filter_on_lat_edge(block, x, y)

    type(block_type), intent(in) :: block
    real(r8), intent(inout) :: x(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out), optional :: y(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                         block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                         block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    real(r8) tmp(block%mesh%full_lon_lb:block%mesh%full_lon_ub)
    integer i, j, k, n, hn

    associate (mesh => block%mesh, filter => block%filter)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        if (filter%ngrid_lat(j) > 1) then
          n  = filter%ngrid_lat(j)
          hn = (n - 1) / 2
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tmp(i) = sum(filter%wgt_lat(:n,j) * x(i-hn:i+hn,j,k))
          end do
          if (present(y)) then
            y(:,j,k) = tmp
          else
            x(:,j,k) = tmp
          end if
        else if (present(y)) then
          y(:,j,k) = x(:,j,k)
        end if
      end do
    end do
    end associate

  end subroutine filter_on_lat_edge

  subroutine filter_on_lev_edge(block, x, y)

    type(block_type), intent(in) :: block
    real(r8), intent(inout) :: x(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8), intent(out), optional :: y(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                         block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                         block%mesh%half_lev_lb:block%mesh%half_lev_ub)

    real(r8) tmp(block%mesh%full_lon_lb:block%mesh%full_lon_ub)
    integer i, j, k, n, hn

    associate (mesh => block%mesh, filter => block%filter)
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        if (filter%ngrid_lon(j) > 1) then
          n  = filter%ngrid_lon(j)
          hn = (n - 1) / 2
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tmp(i) = sum(filter%wgt_lon(:n,j) * x(i-hn:i+hn,j,k))
          end do
          if (present(y)) then
            y(:,j,k) = tmp
          else
            x(:,j,k) = tmp
          end if
        else if (present(y)) then
          y(:,j,k) = x(:,j,k)
        end if
      end do
    end do
    end associate

  end subroutine filter_on_lev_edge

end module filter_mod
