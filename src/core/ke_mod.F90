module ke_mod

  use const_mod
  use mesh_mod
  use allocator_mod
  use parallel_mod
  use state_mod
  use tend_mod

  implicit none

  private

  public calc_ke_on_cell

contains

  subroutine calc_ke_on_cell(state)

    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j
    real(r8) pole

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%ke(i,j) = (mesh%lon_edge_right_area(j  ) * state%u(i-1,j  )**2 + &
                         mesh%lon_edge_left_area (j  ) * state%u(i  ,j  )**2 + &
#ifdef STAGGER_V_ON_POLE
                         mesh%lat_edge_up_area   (j  ) * state%v(i  ,j  )**2 + &
                         mesh%lat_edge_down_area (j+1) * state%v(i  ,j+1)**2   &
#else
                         mesh%lat_edge_up_area   (j-1) * state%v(i  ,j-1)**2 + &
                         mesh%lat_edge_down_area (j  ) * state%v(i  ,j  )**2   &
#endif
                        ) / mesh%cell_area(j)
      end do
    end do
#ifndef STAGGER_V_ON_POLE
    ! Note: lat_edge_down_area and lat_edge_up_area at the Poles is the same as cell_area.
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_start_idx
      pole = 0.0d0
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%v(i,j)**2
      end do
      call parallel_zonal_sum(pole)
      pole = pole / mesh%num_full_lon * 0.5_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%ke(i,j) = pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_end_idx
      pole = 0.0d0
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%v(i,j-1)**2
      end do
      call parallel_zonal_sum(pole)
      pole = pole / mesh%num_full_lon * 0.5_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%ke(i,j) = pole
      end do
    end if
#endif
    call parallel_fill_halo(mesh, state%ke)

  end subroutine calc_ke_on_cell

end module ke_mod