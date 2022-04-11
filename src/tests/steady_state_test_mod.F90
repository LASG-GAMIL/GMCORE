module steady_state_test_mod

  use const_mod
  use formula_mod
  use parallel_mod
  use block_mod
  use vert_coord_mod

  implicit none

  private

  public steady_state_test_set_ic

  real(r8), parameter :: alpha = 0.0_r8
  real(r8), parameter :: u0    = 35        ! m s-1
  real(r8), parameter :: t0    = 288       ! K
  real(r8), parameter :: gamma = 0.005     ! K m-1
  real(r8), parameter :: dt    = 4.8e5     ! K
  real(r8), parameter :: eta0  = 0.252
  real(r8), parameter :: etat  = 0.2       ! Tropopause level

contains

  subroutine steady_state_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) etav, eta, tbar, gzbar, sin_lat, cos_lat
    integer i, j, k

    associate (mesh   => block%mesh           , &
               u      => block%state(1)%u_lon , &
               v      => block%state(1)%v_lat , &
               phs    => block%state(1)%phs   , &
               ph_lev => block%state(1)%ph_lev, &
               ph     => block%state(1)%ph    , &
               gz_lev => block%state(1)%gz_lev, &
               gz     => block%state(1)%gz    , &
               t      => block%state(1)%t     , &
               pt     => block%state(1)%pt    , &
               gzs    => block%static%gzs)
    phs = 1.0e5_r8
    v   = 0

    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, phs(i,j))
        end do
      end do
    end do
    call fill_halo(block, ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ph(i,j,k) = 0.5d0 * (ph_lev(i,j,k) + ph_lev(i,j,k+1))
        end do
      end do
    end do
    call fill_halo(block, ph, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi / 2d0
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          u(i,j,k) = u0 * cos(etav)**(1.5d0) * sin(2 * mesh%full_lat(j))**2
        end do
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        tbar = t0 * eta**(Rd * gamma / g)
      else
        tbar = t0 * eta**(Rd * gamma / g) + dt * (etat - eta)**5
      end if
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          t(i,j,k) = tbar + 3d0 / 4d0 * eta * pi * u0 / Rd * sin(etav) * sqrt(cos(etav)) * (         &
            (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * 2 * u0 * cos(etav)**1.5d0 + &
            (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega          &
          )
          pt(i,j,k) = potential_temperature(t(i,j,k), ph(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, t , full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      eta = mesh%half_lev(k)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g))
      else
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g)) - Rd * dt * (   &
            (log(eta / etat) + 137d0 / 60d0) * etat**5 - 5 * etat**4 * eta + &
            5 * etat**3 * eta**2 - 10d0 / 3d0 * etat**2 * eta**3 +           &
            5d0 / 4d0 * etat * eta**4 - 1d0 / 5d0 * eta**5                   &
          )
      end if
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gz_lev(i,j,k) = gzbar + u0 * cos(etav)**1.5d0 * (                                      &
            (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * u0 * cos(etav)**1.5d0 + &
            (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega      &
          )
        end do
      end do
    end do
    call fill_halo(block, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      eta = mesh%full_lev(k)
      etav = (eta - eta0) * pi / 2d0
      if (etat <= eta .and. eta <= 1) then
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g))
      else
        gzbar = t0 * g / gamma * (1 - eta**(Rd * gamma / g)) - Rd * dt * (   &
            (log(eta / etat) + 137d0 / 60d0) * etat**5 - 5 * etat**4 * eta + &
            5 * etat**3 * eta**2 - 10d0 / 3d0 * etat**2 * eta**3 +           &
            5d0 / 4d0 * etat * eta**4 - 1d0 / 5d0 * eta**5                   &
          )
      end if
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        sin_lat = mesh%full_sin_lat(j)
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gz(i,j,k) = gzbar + u0 * cos(etav)**1.5d0 * (                                          &
            (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * u0 * cos(etav)**1.5d0 + &
            (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega      &
          )
        end do
      end do
    end do
    call fill_halo(block, gz, full_lon=.true., full_lat=.true., full_lev=.true.)

    etav = (1 - eta0) * pi / 2
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        gzs(i,j) = u0 * cos(etav)**1.5d0 * (                                                   &
          (-2 * sin_lat**6 * (cos_lat**2 + 1d0 / 3d0) + 10d0 / 63d0) * u0 * cos(etav)**1.5d0 + &
          (8d0 / 5d0 * cos_lat**3 * (sin_lat**2 + 2d0 / 3d0) - pi / 4d0) * radius * omega      &
        )
      end do
    end do
    call fill_halo(block, gzs, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine steady_state_test_set_ic

end module steady_state_test_mod
