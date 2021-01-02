module pgf_lin97_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use reduce_mod

  implicit none

contains

  subroutine pgf_lin97_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

  end subroutine pgf_lin97_prepare

  subroutine pgf_lin97_run(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    type(mesh_type), pointer :: mesh
    real(r8) dph1, dph2, dgz1, dgz2, dp1, dp2, dpdph
    integer i, j, k, move

    mesh => state%mesh

    !                    o
    !                   /|
    !                  / |
    !                 /  |
    !                /   |
    !   o-----------/------------o
    !   |          /|            |
    !   |         / |            |
    !   |        /  |            |
    !   |       /   |            |
    !   |      o    |            |
    !   o------|    -------------o
    !          |   /
    !          |  /
    !          | /
    !          |/
    !          o

    if (hydrostatic) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        !
        !   4             3
        ! i,j,k        i+1,j,k
        !   o-------------o
        !   |             !
        !   |             !
        !   |    i,j,k    !
        !   |             !
        !   |             !
        !   o-------------o
        ! i,j,k+1      i+1,j,k+1  --> east
        !   1             2
        !
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (block%reduced_mesh(j)%reduce_factor > 1) then
            tend%pgf_lon(:,j,k) = 0.0_r8
            do move = 1, block%reduced_mesh(j)%reduce_factor
              do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
                dph1 = block%reduced_state(j)%ph_lev(k+1,i+1,0,move)**Rd_o_cp - block%reduced_state(j)%ph_lev(k  ,i  ,0,move)**Rd_o_cp
                dph2 = block%reduced_state(j)%ph_lev(k+1,i  ,0,move)**Rd_o_cp - block%reduced_state(j)%ph_lev(k  ,i+1,0,move)**Rd_o_cp
                dgz1 = block%reduced_state(j)%gz_lev(k+1,i  ,0,move)          - block%reduced_state(j)%gz_lev(k  ,i+1,0,move)
                dgz2 = block%reduced_state(j)%gz_lev(k  ,i  ,0,move)          - block%reduced_state(j)%gz_lev(k+1,i+1,0,move)
                block%reduced_tend(j)%pgf_lon(i,k) = -(dph1 * dgz1 + dph2 * dgz2) / block%reduced_mesh(j)%de_lon(0) / (dph1 + dph2)
              end do
              call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%pgf_lon(:,k), mesh, tend%pgf_lon(:,j,k))
            end do
            call overlay_inner_halo(block, tend%pgf_lon(:,j,k), west_halo=.true.)
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              dph1 = state%ph_lev(i+1,j,k+1)**Rd_o_cp - state%ph_lev(i  ,j,k  )**Rd_o_cp ! 2 - 4
              dph2 = state%ph_lev(i  ,j,k+1)**Rd_o_cp - state%ph_lev(i+1,j,k  )**Rd_o_cp ! 1 - 3
              dgz1 = state%gz_lev(i  ,j,k+1)          - state%gz_lev(i+1,j,k  )          ! 1 - 3
              dgz2 = state%gz_lev(i  ,j,k  )          - state%gz_lev(i+1,j,k+1)          ! 4 - 2
              tend%pgf_lon(i,j,k) = -(dph1 * dgz1 + dph2 * dgz2) / mesh%de_lon(j) / (dph1 + dph2)
            end do
          end if
        end do
        !
        !   4             3
        ! i,j,k        i,j+1,k
        !   o-------------o
        !   |             !
        !   |             !
        !   |    i,j,k    !
        !   |             !
        !   |             !
        !   o-------------o
        ! i,j,k+1      i,j+1,k+1  --> north
        !   1             2
        !
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dph1 = state%ph_lev(i,j+1,k+1)**Rd_o_cp - state%ph_lev(i,j  ,k  )**Rd_o_cp ! 2 - 4
            dph2 = state%ph_lev(i,j  ,k+1)**Rd_o_cp - state%ph_lev(i,j+1,k  )**Rd_o_cp ! 1 - 3
            dgz1 = state%gz_lev(i,j  ,k+1)          - state%gz_lev(i,j+1,k  )          ! 1 - 3
            dgz2 = state%gz_lev(i,j  ,k  )          - state%gz_lev(i,j+1,k+1)          ! 4 - 2
            tend%pgf_lat(i,j,k) = -(dph1 * dgz1 + dph2 * dgz2) / mesh%de_lat(j) / (dph1 + dph2)
          end do
        end do
      end do
    else if (nonhydrostatic) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        !
        !   4             3
        ! i,j,k        i+1,j,k
        !   o-------------o
        !   |             !
        !   |             !
        !   |    i,j,k    !
        !   |             !
        !   |             !
        !   o-------------o
        ! i,j,k+1      i+1,j,k+1  --> east
        !   1             2
        !
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (block%reduced_mesh(j)%reduce_factor > 1) then
            tend%pgf_lon(:,j,k) = 0.0_r8
            do move = 1, block%reduced_mesh(j)%reduce_factor
              do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
                dph1 = block%reduced_state(j)%ph_lev(k+1,i+1,0,move)**Rd_o_cp - block%reduced_state(j)%ph_lev(k  ,i  ,0,move)**Rd_o_cp
                dph2 = block%reduced_state(j)%ph_lev(k+1,i  ,0,move)**Rd_o_cp - block%reduced_state(j)%ph_lev(k  ,i+1,0,move)**Rd_o_cp
                dgz1 = block%reduced_state(j)%gz_lev(k+1,i  ,0,move)          - block%reduced_state(j)%gz_lev(k  ,i+1,0,move)
                dgz2 = block%reduced_state(j)%gz_lev(k  ,i  ,0,move)          - block%reduced_state(j)%gz_lev(k+1,i+1,0,move)
                dp1  = block%reduced_state(j)%p_lev (k+1,i  ,0,move)          - block%reduced_state(j)%p_lev (k  ,i+1,0,move)
                dp2  = block%reduced_state(j)%p_lev (k  ,i  ,0,move)          - block%reduced_state(j)%p_lev (k+1,i+1,0,move)
                dpdph = (block%reduced_state(j)%p_lev_lon(k+1,i,0,move) - block%reduced_state(j)%p_lev_lon(k,i,0,move)) / &
                         block%reduced_state(j)%m_lon(k,i,0,move)
                block%reduced_tend(j)%pgf_lon(i,k) = -(                                     &
                  (dph1 * dgz1 + dph2 * dgz2) * dpdph +                                     &
                  (dph1 * dp1  + dph2 * dp2 ) / block%reduced_state(j)%rhod_lon(k,i,0,move) &
                ) / block%reduced_mesh(j)%de_lon(0) / (dph1 + dph2)
              end do
              call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%pgf_lon(:,k), mesh, tend%pgf_lon(:,j,k))
            end do
            call overlay_inner_halo(block, tend%pgf_lon(:,j,k), west_halo=.true.)
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              dph1 = state%ph_lev(i+1,j,k+1)**Rd_o_cp - state%ph_lev(i  ,j,k  )**Rd_o_cp ! 2 - 4
              dph2 = state%ph_lev(i  ,j,k+1)**Rd_o_cp - state%ph_lev(i+1,j,k  )**Rd_o_cp ! 1 - 3
              dgz1 = state%gz_lev(i  ,j,k+1)          - state%gz_lev(i+1,j,k  )          ! 1 - 3
              dgz2 = state%gz_lev(i  ,j,k  )          - state%gz_lev(i+1,j,k+1)          ! 4 - 2
              dp1  = state%p_lev (i  ,j,k+1)          - state%p_lev (i+1,j,k  )          ! 1 - 3
              dp2  = state%p_lev (i  ,j,k  )          - state%p_lev (i+1,j,k+1)          ! 4 - 2
              dpdph = (state%p_lev_lon(i,j,k+1) - state%p_lev_lon(i,j,k)) / state%m_lon(i,j,k)
              tend%pgf_lon(i,j,k) = -(                              &
                (dph1 * dgz1 + dph2 * dgz2) * dpdph +               &
                (dph1 * dp1  + dph2 * dp2 ) / state%rhod_lon(i,j,k) &
              ) / mesh%de_lon(j) / (dph1 + dph2)
            end do
          end if
        end do
        !
        !   4             3
        ! i,j,k        i,j+1,k
        !   o-------------o
        !   |             !
        !   |             !
        !   |    i,j,k    !
        !   |             !
        !   |             !
        !   o-------------o
        ! i,j,k+1      i,j+1,k+1  --> north
        !   1             2
        !
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dph1 = state%ph_lev(i,j+1,k+1)**Rd_o_cp - state%ph_lev(i,j  ,k  )**Rd_o_cp ! 2 - 4
            dph2 = state%ph_lev(i,j  ,k+1)**Rd_o_cp - state%ph_lev(i,j+1,k  )**Rd_o_cp ! 1 - 3
            dgz1 = state%gz_lev(i,j  ,k+1)          - state%gz_lev(i,j+1,k  )          ! 1 - 3
            dgz2 = state%gz_lev(i,j  ,k  )          - state%gz_lev(i,j+1,k+1)          ! 4 - 2
            dp1  = state%p_lev (i,j  ,k+1)          - state%p_lev (i,j+1,k  )          ! 1 - 3
            dp2  = state%p_lev (i,j  ,k  )          - state%p_lev (i,j+1,k+1)          ! 4 - 2
            dpdph = (state%p_lev_lat(i,j,k+1) - state%p_lev_lat(i,j,k)) / state%m_lat(i,j,k)
            tend%pgf_lat(i,j,k) = -(                              &
              (dph1 * dgz1 + dph2 * dgz2) * dpdph +               &
              (dph1 * dp1  + dph2 * dp2 ) / state%rhod_lat(i,j,k) &
            ) / mesh%de_lat(j) / (dph1 + dph2)
          end do
        end do
      end do
    end if

  end subroutine pgf_lin97_run

end module pgf_lin97_mod
