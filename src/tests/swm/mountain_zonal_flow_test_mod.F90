module mountain_zonal_flow_test_mod

  use const_mod
  use parallel_mod
  use mesh_mod
  use state_mod
  use static_mod

  implicit none

  private

  public mountain_zonal_flow_test_set_initial_condition

  real, parameter :: alpha = 0.0
  real, parameter :: u0 = 20.0
  real, parameter :: gd0 = 5960.0 * g ! m2 s-2
  real, parameter :: lon0 = pi * 1.5
  real, parameter :: lat0 = pi / 6.0
  real, parameter :: ghs0 = 2000.0 * g
  real, parameter :: R = pi / 9.0

contains

  subroutine mountain_zonal_flow_test_set_initial_condition()

    real cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha, dlon, d
    integer i, j, k

    write(6, *) '[Notice]: Use mountain zonal flow initial condition.'

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)

#ifdef STAGGER_V_ON_POLE
    do j = static%mesh%half_lat_start_idx, static%mesh%half_lat_end_idx
#else
    do j = static%mesh%full_lat_start_idx, static%mesh%full_lat_end_idx
#endif
      do i = static%mesh%full_lon_start_idx, static%mesh%full_lon_end_idx
        dlon = abs(mesh%full_lon(i) - lon0)
        dlon = min(dlon, 2 * pi - dlon)
        d = min(R, sqrt(dlon**2 + (mesh%full_lat(j) - lat0)**2))
        static%ghs(i,j) = ghs0 * (1.0 - d / R)
      end do
    end do

    call parallel_fill_halo(static%ghs, all_halo=.true.)

#ifdef STAGGER_V_ON_POLE
    do j = state(1)%mesh%half_lat_start_idx, state(1)%mesh%half_lat_end_idx
#else
    do j = state(1)%mesh%full_lat_start_idx_no_pole, state(1)%mesh%full_lat_end_idx_no_pole
#endif
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = state(1)%mesh%half_lon_start_idx, state(1)%mesh%half_lon_end_idx
        cos_lon = mesh%half_cos_lon(i)
        state(1)%u(i,j) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
      end do
    end do

    call parallel_fill_halo(state(1)%u, all_halo=.true.)

#ifdef STAGGER_V_ON_POLE
    do j = state(1)%mesh%full_lat_start_idx_no_pole, state(1)%mesh%full_lat_end_idx_no_pole
#else
    do j = state(1)%mesh%half_lat_start_idx, state(1)%mesh%half_lat_end_idx
#endif
      do i = state(1)%mesh%full_lon_start_idx, state(1)%mesh%full_lon_end_idx
        sin_lon = mesh%full_cos_lon(i)
        state(1)%v(i,j) = - u0 * sin_lon * sin_alpha
      end do
    end do

    call parallel_fill_halo(state(1)%v, all_halo=.true.)

#ifdef STAGGER_V_ON_POLE
    do j = state(1)%mesh%half_lat_start_idx, state(1)%mesh%half_lat_end_idx
#else
    do j = state(1)%mesh%full_lat_start_idx, state(1)%mesh%full_lat_end_idx
#endif
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = state(1)%mesh%full_lon_start_idx, state(1)%mesh%full_lon_end_idx
        cos_lon = mesh%full_cos_lon(i)
        state(1)%gd(i,j) = gd0 - (radius * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2 - static%ghs(i,j)
      end do
    end do

    call parallel_fill_halo(state(1)%gd, all_halo=.true.)

  end subroutine mountain_zonal_flow_test_set_initial_condition

end module mountain_zonal_flow_test_mod