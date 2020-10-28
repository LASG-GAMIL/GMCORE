module steady_state_pgf_test_mod

  use flogger
  use const_mod
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod

  implicit none

  private

  public steady_state_pgf_test_set_initial_condition

  real(r8), parameter :: T0   = 300.d0      ! K
  real(r8), parameter :: h0   = 2000.d0     ! m
  real(r8), parameter :: p0   = 1.0e5       ! pa
  real(r8), parameter :: lonc = 3.d0 * pi / 2
  real(r8), parameter :: latc = 0.0
  real(r8), parameter :: Rm   = 3.d0 * pi / 4
  real(r8), parameter :: gamma= 0.0065d0
  real(r8), parameter :: osm  = pi / 16.d0   

contains
  subroutine steady_state_pgf_test_set_initial_condition(block)

    type(block_type), intent(inout), target :: block
    real(r8) cos_lat, sin_lat, full_lon, r, height
    integer i, j, k
    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static

    mesh => block%mesh
    state => block%state(1)
    static => block%static
    
    state%u = 0.0
    state%v = 0.0

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      sin_lat = mesh%full_sin_lat(j)
      cos_lat = mesh%full_cos_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        full_lon = mesh%full_lon(i)
        r = acos(sin(latc) * sin_lat + cos(latc) * cos_lat * cos(full_lon - lonc))
        if (r < Rm) static%gzs(i,j) = g * h0 / 2.d0 * (1.d0 + cos(pi * r / Rm)) * cos(pi * r / osm)**2
      end do
    end do
    call fill_halo(block, static%gzs, full_lon=.true., full_lat=.true.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%phs(i,j) = p0 * (1.d0 - gamma / T0 * static%gzs(i,j) / g)**(g / Rd / gamma) 
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
!          height = T0 / gamma * (1.d0 - (state%ph(i,j,k) / p0)**(Rd * gamma / g))
!          state%t(i,j,k) = T0 - gamma * height
        end do
      end do
    end do
    call fill_halo(block, state%ph, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%t(i,j,k) = T0 * (state%ph(i,j,k) / p0)**(Rd * gamma / g)
          state%pt(i,j,k) = potential_temperature(state%t(i,j,k), state%ph(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, state%t, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
  
  end subroutine steady_state_pgf_test_set_initial_condition

end module steady_state_pgf_test_mod
