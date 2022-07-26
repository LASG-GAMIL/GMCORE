module held_suarez_test_mod

  use flogger
  use const_mod, only: r8, Rd_o_cpd
  use formula_mod
  use vert_coord_mod
  use block_mod
  use parallel_mod
  use rossby_haurwitz_wave_3d_test_mod

  implicit none

  private

  public held_suarez_test_set_ic
  public held_suarez_test_apply_forcing

  real(r8), parameter :: sig_b    = 0.7_r8
  real(r8), parameter :: kf       = 1.0_r8 / 86400.0_r8   ! s-1
  real(r8), parameter :: ka       = 0.025_r8 / 86400.0_r8 ! s-1
  real(r8), parameter :: ks       = 0.25_r8 / 86400.0_r8  ! s-1
  real(r8), parameter :: dt_lat   = 60.0_r8               ! K
  real(r8), parameter :: dpt_lev  = 10.0_r8               ! K
  real(r8), parameter :: p0       = 1.0e5_r8              ! Pa

contains

  subroutine held_suarez_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(r8) random

    call rossby_haurwitz_wave_3d_test_set_ic(block)

    call random_seed()

    associate (mesh  => block%mesh, &
               phs   => block%state(1)%phs, &
               ph    => block%state(1)%ph , &
               t     => block%state(1)%t  , &
               pt    => block%state(1)%pt)
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        call random_number(random)
        phs(i,j) = phs(i,j) - (0.5_r8 + random) * mesh%full_cos_lat(j)**2
      end do
    end do
    call fill_halo(block, phs, full_lon=.true., full_lat=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          call random_number(random)
          pt(i,j,k) = potential_temperature(t(i,j,k), ph(i,j,k)) - (0.5_r8 + random) * mesh%full_cos_lat(j)**2
        end do
      end do
    end do
    call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine held_suarez_test_set_ic

  subroutine held_suarez_test_apply_forcing(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    real(r8) kv, kt, teq, p_p0
    integer i, j, k

    associate (mesh  => block%mesh , &
               u     => state%u_lon, &
               v     => state%v_lat, &
               ph    => state%ph   , &
               t     => state%t    , &
               pt    => state%pt   )
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      kv = kf * max(0.0_r8, (mesh%full_lev(k) - sig_b) / (1.0_r8 - sig_b))
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          u(i,j,k) = u(i,j,k) - dt * kv * u(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      kv = kf * max(0.0_r8, (mesh%full_lev(k) - sig_b) / (1.0_r8 - sig_b))
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          v(i,j,k) = v(i,j,k) - dt * kv * v(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        kt = ka + (ks - ka) * max(0.0_r8, (mesh%full_lev(k) - sig_b) / (1.0_r8 - sig_b)) * mesh%full_cos_lat(j)**4
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          p_p0 = ph(i,j,k) / p0
          teq = max(200.0_r8, (315.0_r8 - dt_lat * mesh%full_sin_lat(j)**2 - dpt_lev * log(p_p0) * mesh%full_cos_lat(j)**2) * p_p0**Rd_o_cpd)
          pt(i,j,k) = pt(i,j,k) - dt * kt * pt(i,j,k) * (1.0_r8 - teq / t(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine held_suarez_test_apply_forcing

end module held_suarez_test_mod
