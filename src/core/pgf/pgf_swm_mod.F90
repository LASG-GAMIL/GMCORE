module pgf_swm_mod

  use const_mod
  use parallel_mod
  use block_mod
  use reduce_mod

  implicit none

contains

  subroutine pgf_swm_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

  end subroutine pgf_swm_prepare

  subroutine pgf_swm_run(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    type(mesh_type), pointer :: mesh
    integer i, j, k, move

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          tend%pgf_lon(:,j,k) = 0.0_r8
          do move = 1, block%reduced_mesh(j)%reduce_factor
            do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
              block%reduced_tend(j)%pgf_lon(i,k) = (                                            &
                block%reduced_state(j)%gz(k,i+1,0,move) - block%reduced_state(j)%gz(k,i,0,move) &
              ) / block%reduced_mesh(j)%de_lon(0)
            end do
            call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%pgf_lon(:,k), mesh, tend%pgf_lon(:,j,k))
          end do
          call overlay_inner_halo(block, tend%pgf_lon(:,j,k), west_halo=.true.)
        else
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            tend%pgf_lon(i,j,k) = (state%gz(i+1,j,k) - state%gz(i,j,k)) / mesh%de_lon(j)
          end do
        end if
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          tend%pgf_lat(i,j,k) = (state%gz(i,j  ,k) - state%gz(i,j-1,k)) / mesh%de_lat(j)
#else
          tend%pgf_lat(i,j,k) = (state%gz(i,j+1,k) - state%gz(i,j  ,k)) / mesh%de_lat(j)
#endif
        end do
      end do
    end do

  end subroutine pgf_swm_run

end module pgf_swm_mod
