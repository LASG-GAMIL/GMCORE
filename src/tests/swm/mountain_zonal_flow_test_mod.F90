module mountain_zonal_flow_test_mod

  use flogger
  use const_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public mountain_zonal_flow_test_set_initial_condition

  real(r8), parameter :: alpha = 0.0
  real(r8), parameter :: u0 = 20.0
  real(r8), parameter :: gz0 = 5960.0 * g
  real(r8), parameter :: lon0 = pi * 1.5
  real(r8), parameter :: lat0 = pi / 6.0
  real(r8), parameter :: gzs0 = 2000.0 * g
  real(r8), parameter :: R = pi / 9.0

contains

  subroutine mountain_zonal_flow_test_set_initial_condition(block)

    type(block_type), intent(inout), target :: block

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha, dlon, d
    integer i, j, k
    type(mesh_type), pointer :: mesh

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)

    mesh => block%mesh

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        dlon = abs(mesh%full_lon(i) - lon0)
        dlon = min(dlon, 2 * pi - dlon)
        d = min(R, sqrt(dlon**2 + (mesh%full_lat(j) - lat0)**2))
        block%static%gzs(i,j) = gzs0 * (1.0 - d / R)
      end do
    end do
    call fill_halo(block, block%static%gzs, full_lon=.true., full_lat=.true.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        cos_lon = mesh%half_cos_lon(i)
        block%state(1)%u(i,j,1) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
      end do
    end do
    call fill_halo(block, block%state(1)%u, full_lon=.false., full_lat=.true.)

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        sin_lon = mesh%full_sin_lon(i)
        block%state(1)%v(i,j,1) = - u0 * sin_lon * sin_alpha
      end do
    end do
    call fill_halo(block, block%state(1)%v, full_lon=.true., full_lat=.false.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        cos_lon = mesh%full_cos_lon(i)
        block%state(1)%gz(i,j,1) = gz0 - (radius * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2
      end do
    end do
    call fill_halo(block, block%state(1)%gz, full_lon=.true., full_lat=.true.)

  end subroutine mountain_zonal_flow_test_set_initial_condition

end module mountain_zonal_flow_test_mod
