module formula_mod

  use const_mod

  implicit none

contains

  elemental real(r8) function potential_temperature(t, p) result(res)

    real(r8), intent(in) :: t
    real(r8), intent(in) :: p

    real(r8), parameter :: p0 = 1.0e5_r8

    res = t * (p0 / p)**Rd_o_cpd

  end function potential_temperature

  elemental real(r8) function temperature(pt, p) result(res)

    real(r8), intent(in) :: pt
    real(r8), intent(in) :: p

    real(r8), parameter :: p0 = 1.0e5_r8

    res = pt * (p / p0)**Rd_o_cpd

  end function temperature

  elemental real(r8) function dry_air_density(pt, p) result(res)

    real(r8), intent(in) :: pt
    real(r8), intent(in) :: p

    real(r8), parameter :: p0 = 1.0e5_r8

    res = p0 / Rd / pt * (p / p0)**cvd_o_cpd

  end function dry_air_density

end module formula_mod
