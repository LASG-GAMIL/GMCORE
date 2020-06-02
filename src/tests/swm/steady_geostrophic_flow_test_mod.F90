module steady_geostrophic_flow_test_mod

  use flogger
  use const_mod
  use mesh_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public steady_geostrophic_flow_test_set_initial_condition

  real(r8), parameter :: alpha = 0.0_r8
  real(r8), parameter :: u0 = 2.0_r8 * pi * radius / (12.0_r8 * 86400.0_r8)
  real(r8), parameter :: gz0 = 2.94e4_r8 ! m2 s-2

contains

  subroutine steady_geostrophic_flow_test_set_initial_condition(block)

    type(block_type), intent(inout), target :: block

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon, cos_alpha, sin_alpha
    integer i, j
    type(mesh_type), pointer :: mesh

    mesh => block%mesh

    cos_alpha = cos(alpha)
    sin_alpha = sin(alpha)
    block%static%gzs(:,:) = 0.0

    do j = block%state(1)%mesh%full_lat_ibeg, block%state(1)%mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = block%state(1)%mesh%half_lon_ibeg, block%state(1)%mesh%half_lon_iend
        cos_lon = mesh%half_cos_lon(i)
        block%state(1)%u(i,j,1) = u0 * (cos_lat * cos_alpha + cos_lon * sin_lat * sin_alpha)
      end do
    end do
    call fill_halo(block, block%state(1)%u(:,:,1), full_lon=.false., full_lat=.true.)

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        sin_lon = mesh%full_cos_lon(i)
        block%state(1)%v(i,j,1) = - u0 * sin_lon * sin_alpha
      end do
    end do
    call fill_halo(block, block%state(1)%v(:,:,1), full_lon=.true., full_lat=.false.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        cos_lon = mesh%full_cos_lon(i)
        block%state(1)%gz(i,j,1) = gz0 - (radius * omega * u0 + u0**2 * 0.5) * (sin_lat * cos_alpha - cos_lon * cos_lat * sin_alpha)**2
      end do
    end do
    call fill_halo(block, block%state(1)%gz(:,:,1), full_lon=.true., full_lat=.true.)

  end subroutine steady_geostrophic_flow_test_set_initial_condition

end module steady_geostrophic_flow_test_mod
