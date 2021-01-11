module rayleigh_damp_mod

  use const_mod
  use block_mod
  use parallel_mod

  implicit none

  private

  public rayleigh_damp_append_tend
  public rayleigh_damp_w

  real(r8), parameter :: p0 = 1.0e5_r8
  real(r8), parameter :: T0 = 290.0_r8
  real(r8), parameter :: z1 = 37.0e3_r8
  real(r8), parameter :: tau= 3 * 86400.0_r8
  real(r8), parameter :: h0 = 4.4e3_r8
  real(r8), parameter :: u0 = 0.0_r8
  real(r8), parameter :: v0 = 0.0_r8
  real(r8), parameter :: w0 = 0.0_r8
  real(r8), parameter :: H = Rd * T0 / g

contains

  subroutine rayleigh_damp_append_tend(block, state, tend)

    type(block_type), intent(in) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real(r8) p, z, kr
    integer i, j, k

    associate (mesh => block%mesh)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            p = 0.5_r8 * (state%ph(i,j,k) + state%ph(i+1,j,k))
            z = H * log(p0 / p)
            kr = (1 + tanh((z - z1) / h0)) / tau
            tend%du(i,j,k) = tend%du(i,j,k) - kr * (state%u(i,j,k) - u0)
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            p = 0.5_r8 * (state%ph(i,j,k) + state%ph(i,j+1,k))
#else
            p = 0.5_r8 * (state%ph(i,j,k) + state%ph(i,j-1,k))
#endif
            z = H * log(p0 / p)
            kr = (1 + tanh((z - z1) / h0)) / tau
            tend%dv(i,j,k) = tend%dv(i,j,k) - kr * (state%v(i,j,k) - v0)
          end do
        end do
      end do
    end associate

  end subroutine rayleigh_damp_append_tend

  subroutine rayleigh_damp_w(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    real(r8) p, z, kr
    integer i, j, k

    associate (mesh => block%mesh, w => state%w)
      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            p = 0.5_r8 * (state%ph(i,j,k) + state%ph(i+1,j,k))
            z = H * log(p0 / p)
            kr = (1 + tanh((z - z1) / h0)) / tau
            w(i,j,k) = w(i,j,k) - dt * kr * (w(i,j,k) - w0)
          end do
        end do
      end do
      call fill_halo(block, w, full_lon=.true., full_lat=.true., full_lev=.false.)
    end associate

  end subroutine rayleigh_damp_w

end module rayleigh_damp_mod
