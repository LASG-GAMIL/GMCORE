module refer_state_wrf_mod

  use const_mod
  use mesh_mod

  implicit none

  private

  public refer_state_wrf_set_phs

  real(r8), parameter :: p0          = 1.0e5_r8 ! Reference sea level pressure (Pa)
  real(r8), parameter :: t0          = 300.0_r8 ! Reference sea level temperature (270 to 300 K)
  real(r8), parameter :: a           = 50.0_r8  ! Temperature difference between pressure levels of p0 and p0/e (K)
  real(r8), parameter :: t_min       = 200.0_r8 ! Minimum temperature permitted (K)
  real(r8), parameter :: gamma_strat = -11.0_r8 ! Standard stratosphere lapse rate (K)
  real(r8), parameter :: p_strat     = 0.0_r8   ! Pressure at which stratospheric warming begins (Pa)

contains

  subroutine refer_state_wrf_set_phs(mesh, gzs, phs)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in ) :: gzs(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub)
    real(r8), intent(out) :: phs(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub)

    integer i, j

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        phs(i,j) = p0 * exp(-t0 / a + sqrt((t0 / a)**2 - 2 * gzs(i,j) / a / Rd))
      end do
    end do

  end subroutine refer_state_wrf_set_phs

end module refer_state_wrf_mod
