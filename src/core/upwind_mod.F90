module upwind_mod

  use const_mod

  implicit none

  private

  public upwind1
  public upwind3

contains

  pure real(r8) function upwind1(dir, wgt, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: f(2)
    real(r8), intent(in) :: wgt

    real(r8), parameter :: c11 =  0.5_r8
    real(r8), parameter :: c12 = -0.5_r8

    res = c11 * (f(2) + f(1)) + c12 * (f(2) - f(1)) * wgt * dir

  end function upwind1

  pure real(r8) function upwind3(dir, wgt, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: wgt
    real(r8), intent(in) :: f(4)

    real(r8), parameter :: c31 =  7.0_r8 / 12.0_r8
    real(r8), parameter :: c32 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c33 =  1.0_r8 / 12.0_r8

    res = c31 * (f(3) + f(2)) + c32 * (f(4) + f(1)) + c33 * (f(4) - f(1) - 3 * (f(3) - f(2))) * wgt * dir

  end function upwind3

end module upwind_mod
