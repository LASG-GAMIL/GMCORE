module ke_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use state_mod
  use block_mod
  use parallel_mod

  implicit none

  private

  public calc_ke

contains

  subroutine calc_ke(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j, k
    real(r8) ke_vtx(4), pole(state%mesh%num_full_lev)

    associate (mesh => block%mesh)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%ke(i,j,k) = (mesh%area_lon_west (j  ) * state%u(i-1,j  ,k)**2 + &
                               mesh%area_lon_east (j  ) * state%u(i  ,j  ,k)**2 + &
#ifdef V_POLE
                               mesh%area_lat_north(j  ) * state%v(i  ,j  ,k)**2 + &
                               mesh%area_lat_south(j+1) * state%v(i  ,j+1,k)**2   &
#else
                               mesh%area_lat_north(j-1) * state%v(i  ,j-1,k)**2 + &
                               mesh%area_lat_south(j  ) * state%v(i  ,j  ,k)**2   &
#endif
                              ) / mesh%area_cell(j)
          end do
        end do
      end do
#ifndef V_POLE
      ! Note: area_lat_south and area_lat_north at the Poles is the same as area_cell.
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        pole = 0.0d0
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pole(k) = pole(k) + state%v(i,j,k)**2
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole / mesh%num_full_lon * 0.5_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%ke(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        pole = 0.0d0
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pole(k) = pole(k) + state%v(i,j-1,k)**2
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole / mesh%num_full_lon * 0.5_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%ke(i,j,k) = pole(k)
          end do
        end do
      end if
#endif
#ifdef V_POLE
      call fill_halo(block, state%ke, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., north_halo=.false.)
#else
      call fill_halo(block, state%ke, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false.)
#endif
    end associate

  end subroutine calc_ke

end module ke_mod
