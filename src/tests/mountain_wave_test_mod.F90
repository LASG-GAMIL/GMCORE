module mountain_wave_test_mod

  use flogger
  use namelist_mod
  use const_mod
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use operators_mod

  implicit none

  private

  public mountain_wave_test_set_ic

  real(r8), parameter :: T0   = 288.d0      ! K
  real(r8), parameter :: h0   = 2000.d0     ! m
  real(r8), parameter :: d    = 1.5e6 
  real(r8), parameter :: u0   = 20.d0       ! m s-1
  real(r8), parameter :: lonc = pi05
  real(r8), parameter :: latc = pi / 6.0
  real(r8), parameter :: kap  = 2.d0 / 7.d0
  real(r8), parameter :: psp  = 93000.d0    ! Pa
  real(r8), parameter :: N    = 0.0182      ! s-1

contains

  subroutine mountain_wave_test_set_ic(block)

    type(block_type), intent(inout), target :: block
    real(r8) cos_lat, sin_lat, full_lon, r
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
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          u(i,j,k) = u0 * cos_lat
        end do
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

    v = 0.0

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        full_lon = mesh%full_lon(i)
        r = radius * acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(full_lon - lonc))
        gzs(i,j) = g * h0 * exp(-(r / d)**2)
      end do
    end do
    call fill_halo(block, gzs, full_lon=.true., full_lat=.true.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        phs(i,j) = psp * exp(-0.5_r8 * radius * N**2 * u0 / g**2 / kap * (u0 / radius + 2.0_r8 * omega) * &
                   (sin_lat**2 - 1.0_r8) - N**2 / g**2 / kap * gzs(i,j))
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
          t (i,j,k) = 288.d0
          pt(i,j,k) = potential_temperature(t(i,j,k), ph(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, t , full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)

    if (nonhydrostatic) then
      call diag_gz_lev(block, block%state(1))
    end if
    end associate
  
  end subroutine mountain_wave_test_set_ic
  
end module mountain_wave_test_mod
