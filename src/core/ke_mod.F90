module ke_mod

  use const_mod
  use mesh_mod
  use state_mod
  use block_mod
  use process_mod
  use parallel_mod

  implicit none

  private

  public calc_ke_cell

contains

  subroutine calc_ke_cell(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j
    real(r8) pole

    mesh => state%mesh

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%ke(i,j) = (mesh%lon_edge_east_area(j  ) * state%u(i-1,j  )**2 + &
                         mesh%lon_edge_west_area(j  ) * state%u(i  ,j  )**2 + &
#ifdef V_POLE
                         mesh%lat_edge_north_area(j  ) * state%v(i  ,j  )**2 + &
                         mesh%lat_edge_south_area(j+1) * state%v(i  ,j+1)**2   &
#else
                         mesh%lat_edge_north_area(j-1) * state%v(i  ,j-1)**2 + &
                         mesh%lat_edge_south_area(j  ) * state%v(i  ,j  )**2   &
#endif
                        ) / mesh%cell_area(j)
      end do
    end do
#ifndef V_POLE
    ! Note: lat_edge_south_area and lat_edge_north_area at the Poles is the same as cell_area.
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      pole = 0.0d0
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        pole = pole + state%v(i,j)**2
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / mesh%num_full_lon * 0.5_r8
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%ke(i,j) = pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      pole = 0.0d0
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        pole = pole + state%v(i,j-1)**2
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / mesh%num_full_lon * 0.5_r8
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%ke(i,j) = pole
      end do
    end if
#endif
    call fill_halo(block, state%ke, full_lon=.true., full_lat=.true.)

  end subroutine calc_ke_cell

end module ke_mod