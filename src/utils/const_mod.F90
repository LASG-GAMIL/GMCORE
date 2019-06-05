module const_mod

  real, parameter :: pi = atan(1.0) * 4.0
  real, parameter :: pi2 = pi * 2
  real, parameter :: deg = 180.0 / pi
  real, parameter :: rad = pi / 180.0
  real, parameter :: omega = 2.0 * pi / 86400.0
  real, parameter :: radius = 6.37122e6
  real, parameter :: g = 9.80616
  real, parameter :: eps    = epsilon(1.0d0)

end module const_mod
