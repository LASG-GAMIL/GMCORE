module mountain_wave_test_mod

  use flogger
  use const_mod
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod

  implicit none

  private

  public mountain_wave_test_set_initial_condition

  real(r8), parameter :: T0   = 288.d0      ! K
  real(r8), parameter :: h0   = 2000.d0     ! m
  real(r8), parameter :: d    = 1.5e6 
  real(r8), parameter :: u0   = 20.d0       ! m s-1
  real(r8), parameter :: lonc = pi05
  real(r8), parameter :: latc = pi / 6.0
  real(r8), parameter :: kappa= 2.d0 / 7.d0
  real(r8), parameter :: psp  = 93000.d0    ! Pa
  real(r8), parameter :: N    = 0.0182      ! s-1

contains

  subroutine mountain_wave_test_set_initial_condition(block)

    type(block_type), intent(inout), target :: block
    real(r8) cos_lat, sin_lat, full_lon, r
    integer i, j, k
    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static

    mesh => block%mesh
    state => block%state(1)
    static => block%static
    
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        cos_lat = mesh%full_cos_lat(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%u(i,j,k) = u0 * cos_lat
        end do
      end do
    end do
    call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

    state%v = 0.0

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        full_lon = mesh%full_lon(i)
        r = radius * acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(full_lon - lonc))
        static%gzs(i,j) = g * h0 * exp(-(r / d)**2)
      end do
    end do
    call fill_halo(block, static%gzs, full_lon=.true., full_lat=.true.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%phs(i,j) = psp * exp(-0.5_r8 * radius * N**2 * u0 / g**2 / kappa * (u0 / radius + 2.0_r8 * omega) * &
                        (sin_lat**2 - 1.0_r8) - N**2 / g**2 / kappa * static%gzs(i,j))
      end do
    end do
    call fill_halo(block, state%phs, full_lon=.true., full_lat=.true.)

    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, state%phs(i,j))
        end do
      end do
    end do
    call fill_halo(block, state%ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%ph(i,j,k) = 0.5d0 * (state%ph_lev(i,j,k) + state%ph_lev(i,j,k+1))
        end do
      end do
    end do
    call fill_halo(block, state%ph, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%t(i,j,k) = 288.d0
          state%pt(i,j,k) = potential_temperature(state%t(i,j,k), state%ph(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, state%t, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
  
  end subroutine mountain_wave_test_set_initial_condition
  
end module mountain_wave_test_mod
