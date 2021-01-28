module const_mod

  use, intrinsic :: ieee_arithmetic

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
  real(r8)            :: omega  = 2.0_r8 * pi / 86400.0_r8  ! s-1
  real(r8)            :: radius = 6.37122d6                 ! m
  real(r8), parameter :: g      = 9.80616_r8                ! m2 s-2
  real(r8), parameter :: eps    = epsilon(1.0_r8)
  real(r8), parameter :: inf    = huge(1.0_r8)

  real(r8), parameter :: Rd      = 287.04_r8                 ! J kg-1 K-1
  real(r8), parameter :: Rv      = 461.497_r8                ! J kg-1 K-1
  real(r8), parameter :: cp      = 1004.0_r8                 ! J kg-1 K-1
  real(r8), parameter :: cv      = 717.0_r8                  ! J kg-1 K-1
  real(r8), parameter :: Rd_o_cp = Rd / cp
  real(r8), parameter :: cp_o_cv = cp / cv
  real(r8), parameter :: cv_o_cp = cv / cp
  
  integer, parameter :: nest_ratio = 3
  integer, parameter :: inf_i4 = 10000000

  ! Split schemes
  integer, parameter :: csp2          = 1

  integer, parameter :: no_wind_pass  = -1
  integer, parameter :: nh_pass       = -2
  integer, parameter :: all_pass      = 0
  integer, parameter :: slow_pass     = 1
  integer, parameter :: fast_pass     = 2
  integer, parameter :: div_damp_pass = 3
  integer, parameter :: vor_damp_pass = 4

  integer, parameter :: west  = 1
  integer, parameter :: east  = 2
  integer, parameter :: south = 3
  integer, parameter :: north = 4

contains

  pure logical function is_inf(x) result(res)

    real(r8), intent(in) :: x

    res = x - 1 == x

  end function is_inf

end module const_mod
