module cross_pole_flow_test_mod
  
  use flogger
  use string
  use const_mod
  use parallel_mod
  use mesh_mod
  use state_mod
  use static_mod

  implicit none

  private

  public cross_pole_flow_test_set_initial_condition

  real, parameter :: v0 = 20 ! m s-1
  real, parameter :: gd0 = 5.7684e4 ! m2 s-2

contains
  
  subroutine cross_pole_flow_test_set_initial_condition(static, state)

    type(static_type), intent(inout) :: static
    type(state_type) , intent(inout) :: state

    real(r8) cos_lat, sin_lat, cos_lon, sin_lon
    integer i, j
    type(mesh_type), pointer :: mesh

    mesh => state%mesh

    static%ghs(:,:) = 0.0

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        sin_lon = mesh%half_sin_lon(i)
        state%u(i,j) = -v0 * sin_lon * sin_lat * (4.0 * cos_lat**2 - 1.0)
      end do
    end do
    call parallel_fill_halo(mesh, state%u)

    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      sin_lat = mesh%half_sin_lat(j)
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        cos_lon = mesh%full_cos_lon(i)
        state%v(i,j) = v0 * sin_lat**2 * cos_lon 
      end do 
    end do 
    call parallel_fill_halo(mesh, state%v)

    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        sin_lon = mesh%full_sin_lon(i)
        state%gd(i,j) = gd0 + 2 * radius * omega * v0 * sin_lat**3 * cos_lat * sin_lon 
      end do 
    end do 
    call parallel_fill_halo(mesh, state%gd)

  end subroutine cross_pole_flow_test_set_initial_condition

end module cross_pole_flow_test_mod