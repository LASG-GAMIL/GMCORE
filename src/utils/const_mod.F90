module const_mod

  use, intrinsic :: ieee_arithmetic
  use flogger

#ifdef REAL_KIND
  integer, parameter :: r8 = REAL_KIND
#else
  integer, parameter :: r8 = 8
#endif

  real(r8), parameter :: pi     = atan(1.0_r8) * 4.0_r8
  real(r8), parameter :: pi2    = pi * 2
  real(r8), parameter :: pi05   = pi * 0.5_r8
  real(r8), parameter :: deg    = 180.0_r8 / pi
  real(r8), parameter :: rad    = pi / 180.0_r8
  real(r8), parameter :: eps    = epsilon(1.0_r8)
  real(r8), parameter :: inf    = huge(1.0_r8)

  real(r8)            :: omega      ! s-1
  real(r8)            :: radius     ! m
  real(r8)            :: g          ! m2 s-2
  real(r8)            :: Rd         ! J kg-1 K-1
  real(r8)            :: Rv         ! J kg-1 K-1
  real(r8)            :: cp         ! J kg-1 K-1
  real(r8)            :: cv         ! J kg-1 K-1
  real(r8)            :: Rd_o_cp
  real(r8)            :: Rd_o_g
  real(r8)            :: cp_o_cv
  real(r8)            :: cv_o_cp
  real(r8)            :: lapse_rate ! K m-1

  integer, parameter :: nest_ratio = 3
  integer, parameter :: inf_i4 = 10000000

  ! Split schemes
  integer, parameter :: csp2          = 1

  integer, parameter :: no_wind_pass  = -1
  integer, parameter :: all_pass      = 0
  integer, parameter :: slow_pass     = 1
  integer, parameter :: fast_pass     = 2
  integer, parameter :: div_damp_pass = 3
  integer, parameter :: vor_damp_pass = 4
  integer, parameter :: nh_pass_1     = 5
  integer, parameter :: nh_pass_2     = 6

  integer, parameter :: west  = 1
  integer, parameter :: east  = 2
  integer, parameter :: south = 3
  integer, parameter :: north = 4

contains

  subroutine const_init(planet)

    character(*), intent(in) :: planet

    select case (planet)
    case ('earth')
      omega      = 2 * pi / 86400.0d0
      radius     = 6.37122d6
      g          = 9.80616d0
      Rd         = 287.04d0
      Rv         = 461.497d0
      cp         = 1004.0d0
      cv         = 717.0d0
      lapse_rate = 0.006d0
    case ('mars')
      omega      = 2 * pi / 88642.663d0
      radius     = 3.38992d6
      g          = 3.72d0
      Rd         = 191.84d0
      cp         = 735.0d0
      cv         = 543.16d0
      lapse_rate = 5.06d-3
    case default
      call log_error('Invalid planet!')
    end select

    Rd_o_g  = Rd / g
    Rd_o_cp = Rd / cp
    cp_o_cv = cp / cv
    cv_o_cp = cv / cp

  end subroutine const_init

  pure logical function is_inf(x) result(res)

    real(r8), intent(in) :: x

    res = x - 1 == x

  end function is_inf

end module const_mod
