module const_mod

  integer, parameter :: real_kind = 8

  real(real_kind), parameter :: pi = atan(1.0d0) * 4.0d0
  real(real_kind), parameter :: pi2 = pi * 2
  real(real_kind), parameter :: deg = 180.0d0 / pi
  real(real_kind), parameter :: rad = pi / 180.0d0
  real(real_kind), parameter :: omega = 2.0d0 * pi / 86400.0d0
  real(real_kind), parameter :: radius = 6.37122d6
  real(real_kind), parameter :: g = 9.80616d0
  real(real_kind), parameter :: eps = epsilon(1.0d0)

end module const_mod
