module steady_state_pgf_test_mod

  use flogger
  use const_mod, only: r8, pi, Rd, g, omega
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod

  implicit none

  private

  public steady_state_pgf_test_set_params
  public steady_state_pgf_test_set_ic

  real(r8), parameter :: T0    = 300.d0      ! K
  real(r8), parameter :: h0    = 2000.d0     ! m
  real(r8), parameter :: p0    = 1.0e5       ! pa
  real(r8), parameter :: lonc  = 3.d0 * pi / 2
  real(r8), parameter :: latc  = 0.0
  real(r8), parameter :: Rm    = 3.d0 * pi / 4
  real(r8), parameter :: gamma = 0.0065d0
  real(r8), parameter :: osm   = pi / 16.d0

contains

  subroutine steady_state_pgf_test_set_params()

    omega = 0.0

  end subroutine steady_state_pgf_test_set_params

  subroutine steady_state_pgf_test_set_ic(block)

    type(block_type), intent(inout), target :: block
    real(r8) cos_lat, sin_lat, full_lon, r, height
    integer i, j, k

    associate (mesh   => block%mesh           , &
               u      => block%state(1)%u_lon , &
               v      => block%state(1)%v_lat , &
               phs    => block%state(1)%phs   , &
               ph_lev => block%state(1)%ph_lev, &
               ph     => block%state(1)%ph    , &
               t      => block%state(1)%t     , &
               pt     => block%state(1)%pt    , &
               gzs    => block%static%gzs)
    u = 0.0
    v = 0.0

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        full_lon = mesh%full_lon(i)
        r = acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(full_lon - lonc))
        if (r < Rm) gzs(i,j) = g * h0 / 2.d0 * (1.d0 + cos(pi * r / Rm)) * cos(pi * r / osm)**2
      end do
    end do
    call fill_halo(block, gzs, full_lon=.true., full_lat=.true.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        phs(i,j) = p0 * (1.d0 - gamma / T0 * gzs(i,j) / g)**(g / Rd / gamma) 
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
          t (i,j,k) = T0 * (ph(i,j,k) / p0)**(Rd * gamma / g)
          pt(i,j,k) = potential_temperature(t(i,j,k), ph(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, t , full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate
  
  end subroutine steady_state_pgf_test_set_ic

end module steady_state_pgf_test_mod
