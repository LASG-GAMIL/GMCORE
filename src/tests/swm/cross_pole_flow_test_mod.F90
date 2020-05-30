module cross_pole_flow_test_mod
  
  use flogger
  use string
  use const_mod
  use mesh_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public cross_pole_flow_test_set_initial_condition

  real, parameter :: v0 = 20 ! m s-1
  real, parameter :: gz0 = 5.7684e4 ! m2 s-2

contains
  
  subroutine cross_pole_flow_test_set_initial_condition(block)

    type(block_type), intent(inout), target :: block

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon
    integer i, j
    type(mesh_type), pointer :: mesh

    mesh => block%mesh

    block%static%gzs(:,:) = 0.0

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        sin_lon = mesh%half_sin_lon(i)
        block%state(1)%u(i,j) = -v0 * sin_lon * sin_lat * (4.0 * cos_lat**2 - 1.0)
      end do
    end do
    call fill_halo(block, block%state(1)%u, full_lon=.false., full_lat=.true.)

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      sin_lat = mesh%half_sin_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        cos_lon = mesh%full_cos_lon(i)
        block%state(1)%v(i,j) = v0 * sin_lat**2 * cos_lon 
      end do 
    end do 
    call fill_halo(block, block%state(1)%v, full_lon=.true., full_lat=.false.)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        sin_lon = mesh%full_sin_lon(i)
        block%state(1)%gz(i,j) = gz0 + 2 * radius * omega * v0 * sin_lat**3 * cos_lat * sin_lon 
      end do 
    end do 
    call fill_halo(block, block%state(1)%gz, full_lon=.true., full_lat=.true.)

  end subroutine cross_pole_flow_test_set_initial_condition

end module cross_pole_flow_test_mod