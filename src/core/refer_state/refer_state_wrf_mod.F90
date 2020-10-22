module refer_state_wrf_mod

  use const_mod
  use mesh_mod
  use static_mod
  use vert_coord_mod
  use refer_state_types_mod

  implicit none

  private

  public refer_state_wrf_init

  real(r8), parameter :: p0          = 1.0e5_r8 ! Reference sea level pressure (Pa)
  real(r8), parameter :: t0          = 300.0_r8 ! Reference sea level temperature (270 to 300 K)
  real(r8), parameter :: a           = 50.0_r8  ! Temperature difference between pressure levels of p0 and p0/e (K)
  real(r8), parameter :: t_min       = 200.0_r8 ! Minimum temperature permitted (K)
  real(r8), parameter :: gamma_strat = -11.0_r8 ! Standard stratosphere lapse rate (K)
  real(r8), parameter :: p_strat     = 0.0_r8   ! Pressure at which stratospheric warming begins (Pa)

contains

  subroutine refer_state_wrf_init(static, refer_state)

    type(static_type), intent(in) :: static
    type(refer_state_type), intent(out) :: refer_state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => static%mesh

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        refer_state%phs(i,j) = p0 * exp(-t0 / a + sqrt((t0 / a)**2 - 2 * static%gzs(i,j) / a / Rd))
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          refer_state%ph(i,j,k) = vert_coord_calc_ph(k, refer_state%phs(i,j))
          if (refer_state%ph(i,j,k) < p_strat) then
            refer_state%t(i,j,k) = t_min + gamma_strat * log(refer_state%ph(i,j,k) / p_strat)
          else
            refer_state%t(i,j,k) = max(t_min, t0 + a * log(refer_state%ph(i,j,k) / p0))
          end if
        end do
      end do
    end do

  end subroutine refer_state_wrf_init

end module refer_state_wrf_mod
