module operators_mod

  use const_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod
  use parallel_mod
  use namelist_mod
  use log_mod
  use pv_mod
  use ke_mod
  use reduce_mod

  implicit none

  private

  public operators_prepare
  public calc_qhu_qhv
  public calc_dkedlon_dkedlat
  public calc_dpedlon_dpedlat
  public calc_dmfdlon_dmfdlat
  public calc_m_lon_m_lat       ! Shold be called when mass is updated.
  public calc_mf_lon_n_mf_lat_n ! Same as above.
  public calc_mf_lon_t_mf_lat_t ! Same as above.
  public calc_m_vtx             ! Same as above.

contains

  subroutine operators_prepare(state)

    type(state_type), intent(inout) :: state

    call calc_mf_lon_n_mf_lat_n(state)
    call calc_mf_lon_t_mf_lat_t(state)
    call calc_pv_on_vertex(state)
    call calc_ke_on_cell(state)

  end subroutine operators_prepare

  subroutine calc_m_lon_m_lat(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%m_lon(i,j) = (state%mesh%lon_edge_left_area (j) * state%gd(i,  j) + &
                            state%mesh%lon_edge_right_area(j) * state%gd(i+1,j)   &
                           ) / state%mesh%lon_edge_area(j) / g
      end do
    end do

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%m_lat(i,j) = (state%mesh%lat_edge_up_area  (j) * state%gd(i,j  ) + &
                            state%mesh%lat_edge_down_area(j) * state%gd(i,j-1)   &
                           ) / state%mesh%lat_edge_area(j) / g
#else
        state%m_lat(i,j) = (state%mesh%lat_edge_up_area  (j) * state%gd(i,j+1) + &
                            state%mesh%lat_edge_down_area(j) * state%gd(i,j  )   &
                           ) / state%mesh%lat_edge_area(j) / g
#endif
      end do
    end do

  end subroutine calc_m_lon_m_lat

  subroutine calc_m_vtx(state)

    type(state_type), intent(inout) :: state

    integer i, j
    real(r8) pole

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%m_vtx(i,j) = (                                                       &
          (state%gd(i,j-1) + state%gd(i+1,j-1)) * state%mesh%subcell_area(2,j-1) + &
          (state%gd(i,j  ) + state%gd(i+1,j  )) * state%mesh%subcell_area(1,j  )   &
        ) / state%mesh%vertex_area(j) / g
#else
        state%m_vtx(i,j) = (                                                       &
          (state%gd(i,j  ) + state%gd(i+1,j  )) * state%mesh%subcell_area(2,j  ) + &
          (state%gd(i,j+1) + state%gd(i+1,j+1)) * state%mesh%subcell_area(1,j+1)   &
        ) / state%mesh%vertex_area(j) / g
#endif
      end do
    end do
#ifdef STAGGER_V_ON_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%half_lat_start_idx
      pole = 0.0_r8
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%gd(i,j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / g
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%m_vtx(i,j) = pole
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%half_lat_end_idx
      pole = 0.0_r8
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%gd(i,j-1)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / g
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%m_vtx(i,j) = pole
      end do
    end if
#endif

  end subroutine calc_m_vtx

  subroutine calc_mf_lon_n_mf_lat_n(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%mf_lon_n(i,j) = state%m_lon(i,j) * state%u(i,j)
      end do
    end do
    call parallel_fill_halo(state%mesh, state%mf_lon_n)

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%mf_lat_n(i,j) = state%m_lat(i,j) * state%v(i,j)
      end do
    end do
    call parallel_fill_halo(state%mesh, state%mf_lat_n)

  end subroutine calc_mf_lon_n_mf_lat_n

  subroutine calc_mf_lon_t_mf_lat_t(state)

    type(state_type), intent(inout), target :: state

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%mf_lon_t(i,j) = state%mesh%full_tangent_wgt(1,j) * (state%mf_lat_n(i,j  ) + state%mf_lat_n(i+1,j  )) + &
                              state%mesh%full_tangent_wgt(2,j) * (state%mf_lat_n(i,j+1) + state%mf_lat_n(i+1,j+1))
#else
        state%mf_lon_t(i,j) = state%mesh%full_tangent_wgt(1,j) * (state%mf_lat_n(i,j-1) + state%mf_lat_n(i+1,j-1)) + &
                              state%mesh%full_tangent_wgt(2,j) * (state%mf_lat_n(i,j  ) + state%mf_lat_n(i+1,j  ))
#endif
      end do
    end do
    call parallel_fill_halo(state%mesh, state%mf_lon_t)

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%mf_lat_t(i,j) = state%mesh%half_tangent_wgt(1,j) * (state%mf_lon_n(i-1,j-1) + state%mf_lon_n(i,j-1)) + &
                              state%mesh%half_tangent_wgt(2,j) * (state%mf_lon_n(i-1,j  ) + state%mf_lon_n(i,j  ))
#else
        state%mf_lat_t(i,j) = state%mesh%half_tangent_wgt(1,j) * (state%mf_lon_n(i-1,j  ) + state%mf_lon_n(i,j  )) + &
                              state%mesh%half_tangent_wgt(2,j) * (state%mf_lon_n(i-1,j+1) + state%mf_lon_n(i,j+1))
#endif
      end do
    end do
    call parallel_fill_halo(state%mesh, state%mf_lat_t)

  end subroutine calc_mf_lon_t_mf_lat_t

  subroutine calc_qhu_qhv(state, tend, dt)

    type(state_type), intent(inout) :: state
    type(tend_type ), intent(inout) :: tend
    real(r8)        , intent(in   ) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move

    mesh => state%mesh

    select case (pv_scheme)
    case (1)
      call calc_pv_on_edge_midpoint(state)
    case (2)
      call calc_pv_on_edge_upwind(state)
    case (3)
      call calc_pv_on_edge_apvm(state, dt)
    case (4)
      call calc_pv_on_edge_scale_aware_apvm(state)
    case default
      call log_error('Unknown PV scheme!')
    end select

#ifdef STAGGER_V_ON_POLE
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      if (reduced_mesh(j)%reduce_factor > 0) then
        tend%qhv(:,j) = 0.0_r8
        do move = 1, reduced_mesh(j)%reduce_factor
          do i = reduced_mesh(j)%half_lon_start_idx, reduced_mesh(j)%half_lon_end_idx
            reduced_tend(j)%qhv(i) = (                    &
              reduced_mesh(j)%full_tangent_wgt(1,0) * (   &
                reduced_state(j)%mf_lat_n(i  ,0,move) * ( &
                  reduced_state(j)%pv_lon(i  ,0,move) +   &
                  reduced_state(j)%pv_lat(i  ,0,move)     &
                ) +                                       &
                reduced_state(j)%mf_lat_n(i+1,0,move) * ( &
                  reduced_state(j)%pv_lon(i  ,0,move) +   &
                  reduced_state(j)%pv_lat(i+1,0,move)     &
                )                                         &
              ) +                                         &
              reduced_mesh(j)%full_tangent_wgt(2,0) * (   &
                reduced_state(j)%mf_lat_n(i  ,1,move) * ( &
                  reduced_state(j)%pv_lon(i  ,0,move) +   &
                  reduced_state(j)%pv_lat(i  ,1,move)     &
                ) +                                       &
                reduced_state(j)%mf_lat_n(i+1,1,move) * ( &
                  reduced_state(j)%pv_lon(i  ,0,move) +   &
                  reduced_state(j)%pv_lat(i+1,1,move)     &
                )                                         &
              )                                           &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%qhv, mesh, tend%qhv(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%qhv(:,j), left_halo=.true.)
      else
        do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
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
    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      if (reduced_mesh(j)%reduce_factor > 0) then
        tend%qhv(:,j) = 0.0_r8
        do move = 1, reduced_mesh(j)%reduce_factor
          do i = reduced_mesh(j)%half_lon_start_idx, reduced_mesh(j)%half_lon_end_idx
            reduced_tend(j)%qhv(i) = (                     &
              reduced_mesh(j)%full_tangent_wgt(1,0) * (    &
                reduced_state(j)%mf_lat_n(i  ,-1,move) * ( &
                  reduced_state(j)%pv_lon(i  , 0,move) +   &
                  reduced_state(j)%pv_lat(i  ,-1,move)     &
                ) +                                        &
                reduced_state(j)%mf_lat_n(i+1,-1,move) * ( &
                  reduced_state(j)%pv_lon(i  , 0,move) +   &
                  reduced_state(j)%pv_lat(i+1,-1,move)     &
                )                                          &
              ) +                                          &
              reduced_mesh(j)%full_tangent_wgt(2,0) * (    &
                reduced_state(j)%mf_lat_n(i  , 0,move) * ( &
                  reduced_state(j)%pv_lon(i  , 0,move) +   &
                  reduced_state(j)%pv_lat(i  , 0,move)     &
                ) +                                        &
                reduced_state(j)%mf_lat_n(i+1, 0,move) * ( &
                  reduced_state(j)%pv_lon(i  , 0,move) +   &
                  reduced_state(j)%pv_lat(i+1, 0,move)     &
                )                                          &
              )                                            &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%qhv, mesh, tend%qhv(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%qhv(:,j), left_halo=.true.)
      else
        do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
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

#ifdef STAGGER_V_ON_POLE
    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      if (reduced_mesh(j-1)%reduce_factor > 0) then
        tend%qhu(:,j) = 0.0_r8
        do move = 1, reduced_mesh(j-1)%reduce_factor
          do i = reduced_mesh(j-1)%full_lon_start_idx, reduced_mesh(j-1)%full_lon_end_idx
            reduced_tend(j-1)%qhu(i) = (                    &
              reduced_mesh(j-1)%half_tangent_wgt(1,1) * (   &
                reduced_state(j-1)%mf_lon_n(i-1,0,move) * ( &
                  reduced_state(j-1)%pv_lat(i  ,1,move) +   &
                  reduced_state(j-1)%pv_lon(i-1,0,move)     &
                ) +                                         &
                reduced_state(j-1)%mf_lon_n(i  ,0,move) * ( &
                  reduced_state(j-1)%pv_lat(i  ,1,move) +   &
                  reduced_state(j-1)%pv_lon(i  ,0,move)     &
                )                                           &
              )                                             &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, reduced_mesh(j-1), reduced_tend(j-1)%qhu, mesh, tend%qhu(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%qhu(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%qhu(i,j) = ( &
            mesh%half_tangent_wgt(1,j) * (                                            &
              state%mf_lon_n(i-1,j-1) * (state%pv_lat(i,j) + state%pv_lon(i-1,j-1)) + &
              state%mf_lon_n(i  ,j-1) * (state%pv_lat(i,j) + state%pv_lon(i  ,j-1))   &
            )                                                                         &
          ) * 0.5_r8
        end do
      end if
      if (reduced_mesh(j)%reduce_factor > 0) then
        call parallel_zero_halo(mesh, tend%qhu(:,j), right_halo=.true.)
        do move = 1, reduced_mesh(j)%reduce_factor
          do i = reduced_mesh(j)%full_lon_start_idx, reduced_mesh(j)%full_lon_end_idx
            reduced_tend(j)%qhu(i) = (                    &
              reduced_mesh(j)%half_tangent_wgt(2,0) * (   &
                reduced_state(j)%mf_lon_n(i-1,0,move) * ( &
                  reduced_state(j)%pv_lat(i  ,0,move) +   &
                  reduced_state(j)%pv_lon(i-1,0,move)     &
                ) +                                       &
                reduced_state(j)%mf_lon_n(i  ,0,move) * ( &
                  reduced_state(j)%pv_lat(i  ,0,move) +   &
                  reduced_state(j)%pv_lon(i  ,0,move)     &
                )                                         &
              )                                           &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%qhu, mesh, tend%qhu(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%qhu(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
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
    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      if (reduced_mesh(j)%reduce_factor > 0) then
        tend%qhu(:,j) = 0.0_r8
        do move = 1, reduced_mesh(j)%reduce_factor
          do i = reduced_mesh(j)%full_lon_start_idx, reduced_mesh(j)%full_lon_end_idx
            reduced_tend(j)%qhu(i) = (                    &
              reduced_mesh(j)%half_tangent_wgt(1,0) * (   &
                reduced_state(j)%mf_lon_n(i-1,0,move) * ( &
                  reduced_state(j)%pv_lat(i  ,0,move) +   &
                  reduced_state(j)%pv_lon(i-1,0,move)     &
                ) +                                       &
                reduced_state(j)%mf_lon_n(i  ,0,move) * ( &
                  reduced_state(j)%pv_lat(i  ,0,move) +   &
                  reduced_state(j)%pv_lon(i  ,0,move)     &
                )                                         &
              )                                           &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%qhu, mesh, tend%qhu(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%qhu(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%qhu(i,j) = (                                                         &
            mesh%half_tangent_wgt(1,j) * (                                          &
              state%mf_lon_n(i-1,j  ) * (state%pv_lat(i,j) + state%pv_lon(i-1,j)) + &
              state%mf_lon_n(i  ,j  ) * (state%pv_lat(i,j) + state%pv_lon(i  ,j))   &
            )                                                                       &
          ) * 0.5_r8
        end do
      end if
      if (reduced_mesh(j+1)%reduce_factor > 0) then
        call parallel_zero_halo(mesh, tend%qhu(:,j), right_halo=.true.)
        do move = 1, reduced_mesh(j+1)%reduce_factor
          do i = reduced_mesh(j+1)%full_lon_start_idx, reduced_mesh(j+1)%full_lon_end_idx
            reduced_tend(j+1)%qhu(i) = (                     &
              reduced_mesh(j+1)%half_tangent_wgt(2,-1) * (   &
                reduced_state(j+1)%mf_lon_n(i-1, 0,move) * ( &
                  reduced_state(j+1)%pv_lat(i  ,-1,move) +   &
                  reduced_state(j+1)%pv_lon(i-1, 0,move)     &
                ) +                                          &
                reduced_state(j+1)%mf_lon_n(i  , 0,move) * ( &
                  reduced_state(j+1)%pv_lat(i  ,-1,move) +   &
                  reduced_state(j+1)%pv_lon(i  , 0,move)     &
                )                                            &
              )                                              &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, reduced_mesh(j+1), reduced_tend(j+1)%qhu, mesh, tend%qhu(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%qhu(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
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

  subroutine calc_dkedlon_dkedlat(state, tend, dt)

    type(state_type ), intent(inout) :: state
    type(tend_type  ), intent(inout) :: tend
    real(r8)         , intent(in   ) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      if (reduced_mesh(j)%reduce_factor > 0) then
        tend%dkedlon(:,j) = 0.0_r8
        do move = 1, reduced_mesh(j)%reduce_factor
          do i = reduced_mesh(j)%full_lon_start_idx, reduced_mesh(j)%full_lon_end_idx
            reduced_tend(j)%dkedlon(i) = (                                    &
              reduced_state(j)%ke(i+1,0,move) - reduced_state(j)%ke(i,0,move) &
            ) / reduced_mesh(j)%de_lon(0)
          end do
          call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%dkedlon, mesh, tend%dkedlon(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%dkedlon(:,j), left_halo=.true.)
      else
        do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
          tend%dkedlon(i,j) = (state%ke(i+1,j) - state%ke(i,j)) / mesh%de_lon(j)
        end do
      end if
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%dkedlat(i,j) = (state%ke(i,j) - state%ke(i,j-1)) / mesh%de_lat(j)
#else
        tend%dkedlat(i,j) = (state%ke(i,j+1) - state%ke(i,j)) / mesh%de_lat(j)
#endif
      end do
    end do

  end subroutine calc_dkedlon_dkedlat

  subroutine calc_dpedlon_dpedlat(static, state, tend, dt)

    type(static_type), intent(in   ) :: static
    type(state_type ), intent(inout) :: state
    type(tend_type  ), intent(inout) :: tend
    real(r8)         , intent(in   ) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move

    mesh => state%mesh

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      if (reduced_mesh(j)%reduce_factor > 0) then
        tend%dpedlon(:,j) = 0.0_r8
        do move = 1, reduced_mesh(j)%reduce_factor
          do i = reduced_mesh(j)%full_lon_start_idx, reduced_mesh(j)%full_lon_end_idx
            reduced_tend(j)%dpedlon(i) = (                                          &
              reduced_state (j)%gd (i+1,0,move) - reduced_state (j)%gd (i,0,move) + &
              reduced_static(j)%ghs(i+1,0,move) - reduced_static(j)%ghs(i,0,move)   &
            ) / reduced_mesh(j)%de_lon(0)
          end do
          call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%dpedlon, mesh, tend%dpedlon(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%dpedlon(:,j), left_halo=.true.)
      else
        do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
          tend%dpedlon(i,j) = (                   &
            state %gd (i+1,j) - state %gd (i,j) + &
            static%ghs(i+1,j) - static%ghs(i,j)   &
          ) / mesh%de_lon(j)
        end do
      end if
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%dpedlat(i,j) = (                   &
          state %gd (i,j) - state %gd (i,j-1) + &
          static%ghs(i,j) - static%ghs(i,j-1)   &
        ) / mesh%de_lat(j)
#else
        tend%dpedlat(i,j) = (                   &
          state %gd (i,j+1) - state %gd (i,j) + &
          static%ghs(i,j+1) - static%ghs(i,j)   &
        ) / mesh%de_lat(j)
#endif
      end do
    end do

  end subroutine calc_dpedlon_dpedlat

  subroutine calc_dmfdlon_dmfdlat(state, tend, dt)

    type(state_type), intent(in   ) :: state
    type(tend_type) , intent(inout) :: tend
    real(r8)        , intent(in   ) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move
    real(r8) pole

    mesh => state%mesh

    ! --------------------------------------------------------------------------
    !                       Zonal mass flux divergence
    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      if (reduced_mesh(j)%reduce_factor > 0) then
        tend%dmfdlon(:,j) = 0.0_r8
        do move = 1, reduced_mesh(j)%reduce_factor
          do i = reduced_mesh(j)%full_lon_start_idx, reduced_mesh(j)%full_lon_end_idx
            reduced_tend(j)%dmfdlon(i) = (                                                &
              reduced_state(j)%mf_lon_n(i,0,move) - reduced_state(j)%mf_lon_n(i-1,0,move) &
            ) * reduced_mesh(j)%le_lon(0) / reduced_mesh(j)%cell_area(0)
          end do
          call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%dmfdlon, mesh, tend%dmfdlon(:,j))
        end do
        call parallel_overlay_inner_halo(mesh, tend%dmfdlon(:,j), left_halo=.true.)
      else
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%dmfdlon(i,j) = (                         &
            state%mf_lon_n(i,j) - state%mf_lon_n(i-1,j) &
          ) * mesh%le_lon(j) / mesh%cell_area(j)
        end do
      end if
    end do

    ! --------------------------------------------------------------------------
    !                    Meridional mass flux divergence
    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
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

#ifndef STAGGER_V_ON_POLE
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_start_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole + state%mf_lat_n(i,j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole * mesh%le_lat(j) / mesh%num_full_lon / mesh%cell_area(j)
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%dmfdlat(i,j) = pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_end_idx
      pole = 0.0_r8
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        pole = pole - state%mf_lat_n(i,j-1)
      end do
      call parallel_zonal_sum(pole)
      pole = pole * mesh%le_lat(j-1) / mesh%num_full_lon / mesh%cell_area(j)
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        tend%dmfdlat(i,j) = pole
      end do
    end if
#endif

  end subroutine calc_dmfdlon_dmfdlat

end module operators_mod
