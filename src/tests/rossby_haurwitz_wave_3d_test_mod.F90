module rossby_haurwitz_wave_3d_test_mod

  use flogger
  use const_mod
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod

  implicit none

  private

  public rossby_haurwitz_wave_3d_test_set_ic

  integer , parameter :: n     = 4
  real(r8), parameter :: u0    = 50d0               ! m s-1
  real(r8)            :: gz0
  real(r8)            :: M
  real(r8), parameter :: t0    = 288d0              ! K
  real(r8), parameter :: gamma = 0.0065d0           ! K m-1
  real(r8), parameter :: pref  = 955d2              ! Pa

contains

  subroutine rossby_haurwitz_wave_3d_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) lon, cos_lat, sin_lat
    real(r8) a, a_lat, b_lat, c_lat, phi_p
    integer i, j, k

    gz0 = 8.0d3 * g
    M   = u0 / (n * radius)

    block%static%gzs(:,:) = 0d0

    a = radius

    associate (mesh   => block%mesh           , &
               u      => block%state(1)%u_lon , &
               v      => block%state(1)%v_lat , &
               phs    => block%state(1)%phs   , &
               pt     => block%state(1)%pt    , &
               t      => block%state(1)%t     , &
               ph_lev => block%state(1)%ph_lev, &
               ph     => block%state(1)%ph    , &
               gz_lev => block%state(1)%gz_lev, &
               gz     => block%state(1)%gz)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        cos_lat = mesh%full_cos_lat(j)
        sin_lat = mesh%full_sin_lat(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          lon = mesh%half_lon(i)
          u(i,j,k) = a * M * cos_lat + a * M * cos_lat**(n - 1) * cos(n * lon) * (n * sin_lat**2 - cos_lat**2)
        end do
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        cos_lat = mesh%half_cos_lat(j)
        sin_lat = mesh%half_sin_lat(j)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          lon = mesh%full_lon(i)
          v(i,j,k) = - a * M * n * cos_lat**(n-1) * sin_lat * sin(n * lon)
        end do
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        a_lat = M * (2 * omega + M) / 2 * cos_lat**2 + M**2 / 4d0 * cos_lat**(2 * n) * ((n + 1) * cos_lat**2 + (2 * n**2 - n - 2)) - n**2 * M**2 / 2 * cos_lat**(2 * (n - 1))
        b_lat = 2 * (omega + M) * M / (n + 1) / (n + 2) * cos_lat**n * ((n**2 + 2 * n + 2) - (n + 1)**2 * cos_lat**2)
        c_lat = M**2 / 4d0 * cos_lat**(2 * n) * ((n + 1) * cos_lat**2 - (n + 2))
        phi_p = a**2 * (a_lat + b_lat * cos(n * lon) + c_lat * cos(2 * n * lon))
        phs(i,j) = pref * (1 + gamma / g / t0 * phi_p)**(g / gamma / Rd)
      end do
    end do
    call fill_halo(block, phs, full_lon=.true., full_lat=.true.)

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
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          t (i,j,k) = t0 * (ph(i,j,k) / pref)**(gamma * Rd / g)
          pt(i,j,k) = potential_temperature(t(i,j,k), ph(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, t , full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gz(i,j,k) = g * t0 / gamma * (1 - (ph(i,j,k) / pref)**(gamma * Rd / g))
        end do
      end do
    end do
    call fill_halo(block, gz, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gz_lev(i,j,k) = g * t0 / gamma * (1 - (ph_lev(i,j,k) / pref)**(gamma * Rd / g))
        end do
      end do
    end do
    call fill_halo(block, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
    end associate

  end subroutine rossby_haurwitz_wave_3d_test_set_ic

end module rossby_haurwitz_wave_3d_test_mod
