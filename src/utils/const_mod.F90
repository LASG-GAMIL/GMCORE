module const_mod

  use, intrinsic :: ieee_arithmetic

  integer, parameter :: r8 = 8

  real(r8), parameter :: pi = atan(1.0_r8) * 4.0_r8
  real(r8), parameter :: pi2 = pi * 2
  real(r8), parameter :: pi05 = pi * 0.5_r8
  real(r8), parameter :: deg = 180.0_r8 / pi
  real(r8), parameter :: rad = pi / 180.0_r8
  real(r8), parameter :: omega = 2.0_r8 * pi / 86400.0_r8
  real(r8), parameter :: radius = 6.37122d6
  real(r8), parameter :: g = 9.80616_r8
  real(r8), parameter :: eps = epsilon(1.0_r8)
  real(r8), parameter :: inf = huge(1.0_r8)
  
  integer, parameter :: nest_ratio = 3

contains

  pure logical function is_inf(x) result(res)

    real(r8), intent(in) :: x

    res = x - 1 == x

  end function is_inf

end module const_mod
