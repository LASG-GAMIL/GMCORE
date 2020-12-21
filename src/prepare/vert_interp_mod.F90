module vert_interp_mod

  use flogger
  use const_mod
  use process_mod

  implicit none

  private

  public vert_interp_linear
  public vert_interp_log_linear

contains

  subroutine vert_interp_linear(x1, data1, x2, data2, allow_extrap, ierr)

    real(8), intent(in) :: x1(:), data1(:), x2(:)
    real(8), intent(out) :: data2(:)
    logical, intent(in) :: allow_extrap
    integer, intent(out), optional :: ierr

    integer nx1, nx2
    integer i, ii
    integer, allocatable :: i1(:), i2(:)
    real(8), allocatable :: wgt1(:), wgt2(:)
    real(8) tmp1, tmp2

    nx1 = size(x1)
    nx2 = size(x2)

    allocate(i1(nx2), i2(nx2))
    allocate(wgt1(nx2), wgt2(nx2))

    do i = 1, nx2
      do ii = 1, nx1 - 1
        if ((x2(i) >= x1(ii) .and. x2(i) < x1(ii+1))) then
          i1(i) = ii
          i2(i) = ii + 1
          exit
        else if (x2(i) < x1(1)) then
          if (allow_extrap) then
            i1(i) = 1
            i2(i) = 2
          else
            i1(i) = 1
            i2(i) = 1
          end if
          exit
        else if (x2(i) >= x1(nx1)) then
          if (allow_extrap) then
            i1(i) = nx1 - 1
            i2(i) = nx1
          else
            i1(i) = nx1
            i2(i) = nx1
          end if
          exit
        end if
      end do
      if (i1(i) == i2(i)) then
        wgt1(i) = 0.5d0
        wgt2(i) = 0.5d0
      else
        tmp1 = x1(i1(i))
        tmp2 = x1(i2(i))
        wgt1(i) = (tmp2 - x2(i)) / (tmp2 - tmp1)
        wgt2(i) = (x2(i) - tmp1) / (tmp2 - tmp1)
      end if
    end do

    do i = 1, nx2
      data2(i) = data1(i1(i)) * wgt1(i) + data1(i2(i)) * wgt2(i)
    end do

    deallocate(i1, i2, wgt1, wgt2)

  end subroutine vert_interp_linear

  subroutine vert_interp_log_linear(x1, data1, x2, data2, allow_extrap, ierr)

    real(8), intent(in) :: x1(:), data1(:), x2(:)
    real(8), intent(out) :: data2(:)
    logical, intent(in) :: allow_extrap
    integer, intent(out), optional :: ierr

    integer nx1, nx2
    integer i, ii
    integer, allocatable :: i1(:), i2(:)
    real(8), allocatable :: wgt1(:), wgt2(:)
    real(8) tmp1, tmp2

    if (any(x1 == 0.0) .or. any(x2 == 0.0)) then
      if (present(ierr)) then
        ierr = -1
        return
      else
        if (is_root_proc()) call log_error('vert_interp_log_linear: Input coordinate equals zero!')
      end if
    end if

    nx1 = size(x1)
    nx2 = size(x2)

    allocate(i1(nx2), i2(nx2))
    allocate(wgt1(nx2), wgt2(nx2))

    do i = 1, nx2
      do ii = 1, nx1 - 1
        if ((x2(i) >= x1(ii) .and. x2(i) < x1(ii+1))) then
          i1(i) = ii
          i2(i) = ii + 1
          exit
        else if (x2(i) < x1(1)) then
          if (allow_extrap) then
            i1(i) = 1
            i2(i) = 2
          else
            i1(i) = 1
            i2(i) = 1
          end if
          exit
        else if (x2(i) >= x1(nx1)) then
          if (allow_extrap) then
            i1(i) = nx1 - 1
            i2(i) = nx1
          else
            i1(i) = nx1
            i2(i) = nx1
          end if
          exit
        end if
      end do
      if (i1(i) == i2(i)) then
        wgt1(i) = 0.5d0
        wgt2(i) = 0.5d0
      else
        tmp1 = x1(i1(i))
        tmp2 = x1(i2(i))
        wgt1(i) = log(tmp2 / x2(i)) / log(tmp2 / tmp1)
        wgt2(i) = log(x2(i) / tmp1) / log(tmp2 / tmp1)
      end if
    end do

    do i = 1, nx2
      data2(i) = data1(i1(i)) * wgt1(i) + data1(i2(i)) * wgt2(i)
    end do

    deallocate(i1, i2, wgt1, wgt2)

  end subroutine vert_interp_log_linear

end module vert_interp_mod
