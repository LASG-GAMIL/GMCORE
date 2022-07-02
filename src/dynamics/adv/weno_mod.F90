module weno_mod

  use const_mod

  implicit none

  private

  public weno3
  public weno5

contains

  pure real(r8) function weno3(dir, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: f(4)

    real(r8), parameter :: c11 = -0.5_r8, c12 =  1.5_r8
    real(r8), parameter :: c21 =  0.5_r8, c22 =  0.5_r8
    real(r8), parameter :: weno_coef(2) = [1.0_r8 / 3.0_r8, 2.0_r8 / 3.0_r8]
    real(r8), parameter :: eps = 1.0e-16_r8

    real(r8) fs(2), beta(2), alpha(2)

    ! o-----o--x--o-----o
    ! 1     2     3     4
    select case (int(dir))
    case (1)
      fs(1) = c11 * f(1) + c12 * f(2); beta(1) = (f(2) - f(1))**2
      fs(2) = c21 * f(2) + c22 * f(3); beta(2) = (f(3) - f(2))**2
    case (-1)
      fs(1) = c11 * f(4) + c12 * f(3); beta(1) = (f(3) - f(4))**2
      fs(2) = c21 * f(3) + c22 * f(2); beta(2) = (f(2) - f(3))**2
    end select

    alpha = weno_coef / (eps + beta)**2
    res = dot_product(fs, alpha / sum(alpha))

  end function weno3

  pure real(r8) function weno5(dir, f) result(res)

    real(r8), intent(in) :: dir
    real(r8), intent(in) :: f(6)

    real(r8), parameter :: c11 =  1.0_r8 / 3.0_r8, c12 = -7.0_r8 / 6.0_r8, c13 = 11.0_r8 / 6.0_r8
    real(r8), parameter :: c21 = -1.0_r8 / 6.0_r8, c22 =  5.0_r8 / 6.0_r8, c23 =  1.0_r8 / 3.0_r8
    real(r8), parameter :: c31 =  1.0_r8 / 3.0_r8, c32 =  5.0_r8 / 6.0_r8, c33 = -1.0_r8 / 6.0_r8
    real(r8), parameter :: b1 = 13.0_r8 / 12.0_r8, b2 = 0.25_r8
    real(r8), parameter :: weno_coef(3) = [0.1_r8, 0.6_r8, 0.3_r8]
    real(r8), parameter :: eps = 1.0e-16_r8

    real(r8) fs(3), beta(3), alpha(3)

    ! o-----o-----o--x--o-----o-----o
    ! 1     2     3     4     5     6
    select case (int(dir))
    case (1)
      fs(1) = c11 * f(1) + c12 * f(2) + c13 * f(3)
      fs(2) = c21 * f(2) + c22 * f(3) + c23 * f(4)
      fs(3) = c31 * f(3) + c32 * f(4) + c33 * f(5)
      beta(1) = b1 * (f(1) - 2 * f(2) + f(3))**2 + b2 * (    f(1) - 4 * f(2) + 3 * f(3))**2
      beta(2) = b1 * (f(2) - 2 * f(3) + f(4))**2 + b2 * (    f(2)            -     f(4))**2
      beta(3) = b1 * (f(3) - 2 * f(4) + f(5))**2 + b2 * (3 * f(3) - 4 * f(4) +     f(5))**2
    case (-1)
      fs(1) = c11 * f(6) + c12 * f(5) + c13 * f(4)
      fs(2) = c21 * f(5) + c22 * f(4) + c23 * f(3)
      fs(3) = c31 * f(4) + c32 * f(3) + c33 * f(2)
      beta(1) = b1 * (f(6) - 2 * f(5) + f(4))**2 + b2 * (    f(6) - 4 * f(5) + 3 * f(4))**2
      beta(2) = b1 * (f(5) - 2 * f(4) + f(3))**2 + b2 * (    f(5)            -     f(3))**2
      beta(3) = b1 * (f(4) - 2 * f(3) + f(2))**2 + b2 * (3 * f(4) - 4 * f(3) +     f(2))**2
    end select

    alpha = weno_coef / (eps + beta)**2
    res = dot_product(fs, alpha / sum(alpha))

  end function weno5

end module weno_mod
