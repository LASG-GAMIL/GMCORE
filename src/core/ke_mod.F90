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

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%ke(i,j) = (mesh%area_lon_west (j  ) * state%u(i-1,j  )**2 + &
                         mesh%area_lon_east (j  ) * state%u(i  ,j  )**2 + &
#ifdef V_POLE
                         mesh%area_lat_north(j  ) * state%v(i  ,j  )**2 + &
                         mesh%area_lat_south(j+1) * state%v(i  ,j+1)**2   &
#else
                         mesh%area_lat_north(j-1) * state%v(i  ,j-1)**2 + &
                         mesh%area_lat_south(j  ) * state%v(i  ,j  )**2   &
#endif
                        ) / mesh%area_cell(j)
      end do
    end do
!$OMP END PARALLEL DO
#ifndef V_POLE
    ! Note: area_lat_south and area_lat_north at the Poles is the same as area_cell.
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