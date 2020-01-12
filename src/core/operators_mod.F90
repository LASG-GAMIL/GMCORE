module operators_mod

  use const_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod
  use block_mod
  use parallel_mod
  use namelist_mod
  use log_mod
  use pv_mod
  use reduce_mod

  implicit none

  private

  public operators_prepare
  public calc_qhu_qhv
  public calc_dkedlon_dkedlat
  public calc_dpedlon_dpedlat
  public calc_dmfdlon_dmfdlat

contains

  subroutine operators_prepare(blocks, itime)

    type(block_type), intent(inout) :: blocks(:)
    integer, intent(in) :: itime

    integer iblk

    do iblk = 1, size(blocks)
      call blocks(iblk)%state(itime)%update_all()
    end do

  end subroutine operators_prepare

  subroutine calc_qhu_qhv(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move

    mesh => state%mesh

    select case (pv_scheme)
    case (1)
      call calc_pv_on_edge_midpoint(state)
    case (3)
      call calc_pv_on_edge_apvm(state, dt)
    case default
      call log_error('Unknown PV scheme!')
    end select

#ifdef V_POLE
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%qhv(:,j) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
            block%reduced_tend(j)%qhv(i) = (                    &
              block%reduced_mesh(j)%full_tangent_wgt(1,0) * (   &
                block%reduced_state(j)%mf_lat_n(i  ,0,move) * ( &
                  block%reduced_state(j)%pv_lon(i  ,0,move) +   &
                  block%reduced_state(j)%pv_lat(i  ,0,move)     &
                ) +                                             &
                block%reduced_state(j)%mf_lat_n(i+1,0,move) * ( &
                  block%reduced_state(j)%pv_lon(i  ,0,move) +   &
                  block%reduced_state(j)%pv_lat(i+1,0,move)     &
                )                                               &
              ) +                                               &
              block%reduced_mesh(j)%full_tangent_wgt(2,0) * (   &
                block%reduced_state(j)%mf_lat_n(i  ,1,move) * ( &
                  block%reduced_state(j)%pv_lon(i  ,0,move) +   &
                  block%reduced_state(j)%pv_lat(i  ,1,move)     &
                ) +                                             &
                block%reduced_state(j)%mf_lat_n(i+1,1,move) * ( &
                  block%reduced_state(j)%pv_lon(i  ,0,move) +   &
                  block%reduced_state(j)%pv_lat(i+1,1,move)     &
                )                                               &
              )                                                 &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhv, mesh, tend%qhv(:,j))
        end do
        call overlay_inner_halo(mesh, tend%qhv(:,j), left_halo=.true.)
      else
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%qhv(i,j) = (                                                           &
            mesh%full_tangent_wgt(1,j) * (                                            &
              state%mf_lat_n(i  ,j  ) * (state%pv_lon(i,j) + state%pv_lat(i  ,j  )) + &
              state%mf_lat_n(i+1,j  ) * (state%pv_lon(i,j) + state%pv_lat(i+1,j  ))   &
            ) +                                                                       &
            mesh%full_tangent_wgt(2,j) * (                                            &
              state%mf_lat_n(i  ,j+1) * (state%pv_lon(i,j) + state%pv_lat(i  ,j+1)) + &
              state%mf_lat_n(i+1,j+1) * (state%pv_lon(i,j) + state%pv_lat(i+1,j+1))   &
            )                                                                         &
          ) * 0.5_r8
        end do
      end if
    end do
#else
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%qhv(:,j) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
            block%reduced_tend(j)%qhv(i) = (                     &
              block%reduced_mesh(j)%full_tangent_wgt(1,0) * (    &
                block%reduced_state(j)%mf_lat_n(i  ,-1,move) * ( &
                  block%reduced_state(j)%pv_lon(i  , 0,move) +   &
                  block%reduced_state(j)%pv_lat(i  ,-1,move)     &
                ) +                                              &
                block%reduced_state(j)%mf_lat_n(i+1,-1,move) * ( &
                  block%reduced_state(j)%pv_lon(i  , 0,move) +   &
                  block%reduced_state(j)%pv_lat(i+1,-1,move)     &
                )                                                &
              ) +                                                &
              block%reduced_mesh(j)%full_tangent_wgt(2,0) * (    &
                block%reduced_state(j)%mf_lat_n(i  , 0,move) * ( &
                  block%reduced_state(j)%pv_lon(i  , 0,move) +   &
                  block%reduced_state(j)%pv_lat(i  , 0,move)     &
                ) +                                              &
                block%reduced_state(j)%mf_lat_n(i+1, 0,move) * ( &
                  block%reduced_state(j)%pv_lon(i  , 0,move) +   &
                  block%reduced_state(j)%pv_lat(i+1, 0,move)     &
                )                                                &
              )                                                  &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhv, mesh, tend%qhv(:,j))
        end do
        call overlay_inner_halo(mesh, tend%qhv(:,j), left_halo=.true.)
      else
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%qhv(i,j) = (                                                           &
            mesh%full_tangent_wgt(1,j) * (                                            &
              state%mf_lat_n(i  ,j-1) * (state%pv_lon(i,j) + state%pv_lat(i  ,j-1)) + &
              state%mf_lat_n(i+1,j-1) * (state%pv_lon(i,j) + state%pv_lat(i+1,j-1))   &
            ) +                                                                       &
            mesh%full_tangent_wgt(2,j) * (                                            &
              state%mf_lat_n(i  ,j  ) * (state%pv_lon(i,j) + state%pv_lat(i  ,j  )) + &
              state%mf_lat_n(i+1,j  ) * (state%pv_lon(i,j) + state%pv_lat(i+1,j  ))   &
            )                                                                         &
          ) * 0.5_r8
        end do
      end if
    end do
#endif

#ifdef V_POLE
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      if (block%reduced_mesh(j-1)%reduce_factor > 0) then
        tend%qhu(:,j) = 0.0_r8
        do move = 1, block%reduced_mesh(j-1)%reduce_factor
          do i = block%reduced_mesh(j-1)%full_lon_ibeg, block%reduced_mesh(j-1)%full_lon_iend
            block%reduced_tend(j-1)%qhu(i) = (                    &
              block%reduced_mesh(j-1)%half_tangent_wgt(1,1) * (   &
                block%reduced_state(j-1)%mf_lon_n(i-1,0,move) * ( &
                  block%reduced_state(j-1)%pv_lat(i  ,1,move) +   &
                  block%reduced_state(j-1)%pv_lon(i-1,0,move)     &
                ) +                                         &
                block%reduced_state(j-1)%mf_lon_n(i  ,0,move) * ( &
                  block%reduced_state(j-1)%pv_lat(i  ,1,move) +   &
                  block%reduced_state(j-1)%pv_lon(i  ,0,move)     &
                )                                                 &
              )                                                   &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, block%reduced_mesh(j-1), block%reduced_tend(j-1)%qhu, mesh, tend%qhu(:,j))
        end do
        call overlay_inner_halo(mesh, tend%qhu(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%qhu(i,j) = ( &
            mesh%half_tangent_wgt(1,j) * (                                            &
              state%mf_lon_n(i-1,j-1) * (state%pv_lat(i,j) + state%pv_lon(i-1,j-1)) + &
              state%mf_lon_n(i  ,j-1) * (state%pv_lat(i,j) + state%pv_lon(i  ,j-1))   &
            )                                                                         &
          ) * 0.5_r8
        end do
      end if
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        call zero_halo(mesh, tend%qhu(:,j), right_halo=.true.)
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
            block%reduced_tend(j)%qhu(i) = (                    &
              block%reduced_mesh(j)%half_tangent_wgt(2,0) * (   &
                block%reduced_state(j)%mf_lon_n(i-1,0,move) * ( &
                  block%reduced_state(j)%pv_lat(i  ,0,move) +   &
                  block%reduced_state(j)%pv_lon(i-1,0,move)     &
                ) +                                       &
                block%reduced_state(j)%mf_lon_n(i  ,0,move) * ( &
                  block%reduced_state(j)%pv_lat(i  ,0,move) +   &
                  block%reduced_state(j)%pv_lon(i  ,0,move)     &
                )                                               &
              )                                                 &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhu, mesh, tend%qhu(:,j))
        end do
        call overlay_inner_halo(mesh, tend%qhu(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%qhu(i,j) = tend%qhu(i,j) + (                                           &
            mesh%half_tangent_wgt(2,j) * (                                            &
              state%mf_lon_n(i-1,j  ) * (state%pv_lat(i,j) + state%pv_lon(i-1,j  )) + &
              state%mf_lon_n(i  ,j  ) * (state%pv_lat(i,j) + state%pv_lon(i  ,j  ))   &
            )                                                                         &
          ) * 0.5_r8
        end do
      end if
    end do
#else
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%qhu(:,j) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
            block%reduced_tend(j)%qhu(i) = (                    &
              block%reduced_mesh(j)%half_tangent_wgt(1,0) * (   &
                block%reduced_state(j)%mf_lon_n(i-1,0,move) * ( &
                  block%reduced_state(j)%pv_lat(i  ,0,move) +   &
                  block%reduced_state(j)%pv_lon(i-1,0,move)     &
                ) +                                             &
                block%reduced_state(j)%mf_lon_n(i  ,0,move) * ( &
                  block%reduced_state(j)%pv_lat(i  ,0,move) +   &
                  block%reduced_state(j)%pv_lon(i  ,0,move)     &
                )                                               &
              )                                                 &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhu, mesh, tend%qhu(:,j))
        end do
        call overlay_inner_halo(mesh, tend%qhu(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%qhu(i,j) = (                                                         &
            mesh%half_tangent_wgt(1,j) * (                                          &
              state%mf_lon_n(i-1,j  ) * (state%pv_lat(i,j) + state%pv_lon(i-1,j)) + &
              state%mf_lon_n(i  ,j  ) * (state%pv_lat(i,j) + state%pv_lon(i  ,j))   &
            )                                                                       &
          ) * 0.5_r8
        end do
      end if
      if (block%reduced_mesh(j+1)%reduce_factor > 0) then
        call zero_halo(mesh, tend%qhu(:,j), right_halo=.true.)
        do move = 1, block%reduced_mesh(j+1)%reduce_factor
          do i = block%reduced_mesh(j+1)%full_lon_ibeg, block%reduced_mesh(j+1)%full_lon_iend
            block%reduced_tend(j+1)%qhu(i) = (                     &
              block%reduced_mesh(j+1)%half_tangent_wgt(2,-1) * (   &
                block%reduced_state(j+1)%mf_lon_n(i-1, 0,move) * ( &
                  block%reduced_state(j+1)%pv_lat(i  ,-1,move) +   &
                  block%reduced_state(j+1)%pv_lon(i-1, 0,move)     &
                ) +                                                &
                block%reduced_state(j+1)%mf_lon_n(i  , 0,move) * ( &
                  block%reduced_state(j+1)%pv_lat(i  ,-1,move) +   &
                  block%reduced_state(j+1)%pv_lon(i  , 0,move)     &
                )                                                  &
              )                                                    &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, block%reduced_mesh(j+1), block%reduced_tend(j+1)%qhu, mesh, tend%qhu(:,j))
        end do
        call overlay_inner_halo(mesh, tend%qhu(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%qhu(i,j) = tend%qhu(i,j) + (                                           &
            mesh%half_tangent_wgt(2,j) * (                                            &
              state%mf_lon_n(i-1,j+1) * (state%pv_lat(i,j) + state%pv_lon(i-1,j+1)) + &
              state%mf_lon_n(i  ,j+1) * (state%pv_lat(i,j) + state%pv_lon(i  ,j+1))   &
            )                                                                         &
          ) * 0.5_r8
        end do
      end if
    end do
#endif

  end subroutine calc_qhu_qhv

  subroutine calc_dkedlon_dkedlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move

    mesh => state%mesh

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%dkedlon(:,j) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
            block%reduced_tend(j)%dkedlon(i) = (                                    &
              block%reduced_state(j)%ke(i+1,0,move) - block%reduced_state(j)%ke(i,0,move) &
            ) / block%reduced_mesh(j)%de_lon(0)
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dkedlon, mesh, tend%dkedlon(:,j))
        end do
        call overlay_inner_halo(mesh, tend%dkedlon(:,j), left_halo=.true.)
      else
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%dkedlon(i,j) = (state%ke(i+1,j) - state%ke(i,j)) / mesh%de_lon(j)
        end do
      end if
    end do

    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        tend%dkedlat(i,j) = (state%ke(i,j) - state%ke(i,j-1)) / mesh%de_lat(j)
#else
        tend%dkedlat(i,j) = (state%ke(i,j+1) - state%ke(i,j)) / mesh%de_lat(j)
#endif
      end do
    end do

  end subroutine calc_dkedlon_dkedlat

  subroutine calc_dpedlon_dpedlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move

    mesh => state%mesh

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%dpedlon(:,j) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
            block%reduced_tend(j)%dpedlon(i) = (                                          &
              block%reduced_state (j)%gd (i+1,0,move) - block%reduced_state (j)%gd (i,0,move) + &
              block%reduced_static(j)%ghs(i+1,0,move) - block%reduced_static(j)%ghs(i,0,move)   &
            ) / block%reduced_mesh(j)%de_lon(0)
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dpedlon, mesh, tend%dpedlon(:,j))
        end do
        call overlay_inner_halo(mesh, tend%dpedlon(:,j), left_halo=.true.)
      else
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%dpedlon(i,j) = (                   &
                  state %gd (i+1,j) -       state %gd (i,j) + &
            block%static%ghs(i+1,j) - block%static%ghs(i,j)   &
          ) / mesh%de_lon(j)
        end do
      end if
    end do

    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        tend%dpedlat(i,j) = (                   &
                state %gd (i,j) -       state %gd (i,j-1) + &
          block%static%ghs(i,j) - block%static%ghs(i,j-1)   &
        ) / mesh%de_lat(j)
#else
        tend%dpedlat(i,j) = (                   &
                state %gd (i,j+1) - state %gd (i,j) + &
          block%static%ghs(i,j+1) - block%static%ghs(i,j)   &
        ) / mesh%de_lat(j)
#endif
      end do
    end do

  end subroutine calc_dpedlon_dpedlat

  subroutine calc_dmfdlon_dmfdlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move
    real(r8) pole

    mesh => state%mesh

    ! --------------------------------------------------------------------------
    !                       Zonal mass flux divergence
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%dmfdlon(:,j) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
            block%reduced_tend(j)%dmfdlon(i) = (                                                      &
              block%reduced_state(j)%mf_lon_n(i,0,move) - block%reduced_state(j)%mf_lon_n(i-1,0,move) &
            ) * block%reduced_mesh(j)%le_lon(0) / block%reduced_mesh(j)%cell_area(0)
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dmfdlon, mesh, tend%dmfdlon(:,j))
        end do
        call overlay_inner_halo(mesh, tend%dmfdlon(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dmfdlon(i,j) = (                         &
            state%mf_lon_n(i,j) - state%mf_lon_n(i-1,j) &
          ) * mesh%le_lon(j) / mesh%cell_area(j)
        end do
      end if
    end do

    ! --------------------------------------------------------------------------
    !                    Meridional mass flux divergence
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        tend%dmfdlat(i,j) = (                        &
          state%mf_lat_n(i,j+1) * mesh%le_lat(j+1) - &
          state%mf_lat_n(i,j  ) * mesh%le_lat(j  )   &
        ) / mesh%cell_area(j)
#else
        tend%dmfdlat(i,j) = (                        &
          state%mf_lat_n(i,j  ) * mesh%le_lat(j  ) - &
          state%mf_lat_n(i,j-1) * mesh%le_lat(j-1)   &
        ) / mesh%cell_area(j)
#endif
      end do
    end do

#ifndef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      pole = 0.0_r8
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        pole = pole + state%mf_lat_n(i,j)
      end do
      call zonal_sum(pole)
      pole = pole * mesh%le_lat(j) / mesh%num_full_lon / mesh%cell_area(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        tend%dmfdlat(i,j) = pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      pole = 0.0_r8
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        pole = pole - state%mf_lat_n(i,j-1)
      end do
      call zonal_sum(pole)
      pole = pole * mesh%le_lat(j-1) / mesh%num_full_lon / mesh%cell_area(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        tend%dmfdlat(i,j) = pole
      end do
    end if
#endif

  end subroutine calc_dmfdlon_dmfdlat

end module operators_mod
