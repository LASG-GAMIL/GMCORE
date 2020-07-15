module formula_mod

  use const_mod

  implicit none

contains

  elemental real(r8) function potential_temperature(t, p) result(res)

    real(r8), intent(in) :: t
    real(r8), intent(in) :: p

    real(r8), parameter :: p0 = 1.0e5_r8
    real(r8), parameter :: k  = Rd / cp

    res = t * (p0 / p)**k

  end function potential_temperature

  elemental real(r8) function temperature(pt, p) result(res)

    real(r8), intent(in) :: pt
    real(r8), intent(in) :: p

    real(r8), parameter :: p0 = 1.0e5_r8
    real(r8), parameter :: k  = Rd / cp

    res = pt * (p / p0)**k

  end function temperature

end module formula_mod
