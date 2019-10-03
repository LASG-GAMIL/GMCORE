module const_mod

  integer, parameter :: r8 = 8

  real(r8), parameter :: pi = atan(1.0d0) * 4.0d0
  real(r8), parameter :: pi2 = pi * 2
  real(r8), parameter :: pi05 = pi * 0.5_r8
  real(r8), parameter :: deg = 180.0d0 / pi
  real(r8), parameter :: rad = pi / 180.0d0
  real(r8), parameter :: omega = 2.0d0 * pi / 86400.0d0
  real(r8), parameter :: radius = 6.37122d6
  real(r8), parameter :: g = 9.80616d0
  real(r8), parameter :: eps = epsilon(1.0d0)
  real(r8), parameter :: inf = 1.0d33

end module const_mod
