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

  implicit none

  private

  public operators_prepare
  public nonlinear_coriolis_operator
  public energy_gradient_operator
  public mass_flux_divergence_operator

contains

  subroutine operators_prepare(state)

    type(state_type), intent(inout) :: state

    call calc_normal_mass_flux(state)
    call calc_tangent_mass_flux(state)
    call calc_pv_on_vertex(state)
    call calc_ke_on_cell(state)

  end subroutine operators_prepare

  subroutine calc_normal_mass_flux(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%mass_lon(i,j) = (state%mesh%lon_edge_left_area (j) * state%hd(i,  j) + &
                               state%mesh%lon_edge_right_area(j) * state%hd(i+1,j)   &
                              ) / state%mesh%lon_edge_area(j)
        state%mass_flux_lon_n(i,j) = state%mass_lon(i,j) * state%u(i,j)
      end do
    end do

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%mass_lat(i,j) = (state%mesh%lat_edge_up_area  (j) * state%hd(i,j  ) + &
                               state%mesh%lat_edge_down_area(j) * state%hd(i,j-1)   &
                              ) / state%mesh%lat_edge_area(j)
#else
        state%mass_lat(i,j) = (state%mesh%lat_edge_up_area  (j) * state%hd(i,j+1) + &
                               state%mesh%lat_edge_down_area(j) * state%hd(i,j  )   &
                              ) / state%mesh%lat_edge_area(j)
#endif
        state%mass_flux_lat_n(i,j) = state%mass_lat(i,j) * state%v(i,j)
      end do
    end do

    call parallel_fill_halo(state%mesh, state%mass_flux_lon_n, all_halo=.true.)
    call parallel_fill_halo(state%mesh, state%mass_flux_lat_n, all_halo=.true.)

  end subroutine calc_normal_mass_flux

  subroutine calc_tangent_mass_flux(state)

    type(state_type), intent(inout), target :: state

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%mass_flux_lat_t(i,j) = state%mesh%full_tangent_wgt(1,j) * (state%mass_flux_lat_n(i,j  ) + state%mass_flux_lat_n(i+1,j  )) + &
                                     state%mesh%full_tangent_wgt(2,j) * (state%mass_flux_lat_n(i,j+1) + state%mass_flux_lat_n(i+1,j+1))
#else
        state%mass_flux_lat_t(i,j) = state%mesh%full_tangent_wgt(1,j) * (state%mass_flux_lat_n(i,j-1) + state%mass_flux_lat_n(i+1,j-1)) + &
                                     state%mesh%full_tangent_wgt(2,j) * (state%mass_flux_lat_n(i,j  ) + state%mass_flux_lat_n(i+1,j  ))
#endif
      end do
    end do

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%mass_flux_lon_t(i,j) = state%mesh%half_tangent_wgt(1,j) * (state%mass_flux_lon_n(i-1,j-1) + state%mass_flux_lon_n(i,j-1)) + &
                                     state%mesh%half_tangent_wgt(2,j) * (state%mass_flux_lon_n(i-1,j  ) + state%mass_flux_lon_n(i,j  ))
#else
        state%mass_flux_lon_t(i,j) = state%mesh%half_tangent_wgt(1,j) * (state%mass_flux_lon_n(i-1,j  ) + state%mass_flux_lon_n(i,j  )) + &
                                     state%mesh%half_tangent_wgt(2,j) * (state%mass_flux_lon_n(i-1,j+1) + state%mass_flux_lon_n(i,j+1))
#endif
      end do
    end do

    call parallel_fill_halo(state%mesh, state%mass_flux_lat_t, all_halo=.true.)
    call parallel_fill_halo(state%mesh, state%mass_flux_lon_t, all_halo=.true.)

  end subroutine calc_tangent_mass_flux

  subroutine nonlinear_coriolis_operator(state, tend)

    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    select case (pv_scheme)
    case (1)
      call calc_pv_on_edge_midpoint(state)
    case (2)
      call calc_pv_on_edge_upwind(state)
    case default
      call log_error('Unknown PV scheme!')
    end select

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%qhv(i,j) = (state%mesh%full_tangent_wgt(1,j) * (state%mass_flux_lat_n(i  ,j  ) * (state%pv_lon(i,j) + state%pv_lat(i  ,j  ))  + &
                                                             state%mass_flux_lat_n(i+1,j  ) * (state%pv_lon(i,j) + state%pv_lat(i+1,j  ))) + &
                         state%mesh%full_tangent_wgt(2,j) * (state%mass_flux_lat_n(i  ,j+1) * (state%pv_lon(i,j) + state%pv_lat(i  ,j+1))  + &
                                                             state%mass_flux_lat_n(i+1,j+1) * (state%pv_lon(i,j) + state%pv_lat(i+1,j+1)))   &
                        ) * 0.5d0
#else
        tend%qhv(i,j) = (state%mesh%full_tangent_wgt(1,j) * (state%mass_flux_lat_n(i  ,j-1) * (state%pv_lon(i,j) + state%pv_lat(i  ,j-1))  + &
                                                             state%mass_flux_lat_n(i+1,j-1) * (state%pv_lon(i,j) + state%pv_lat(i+1,j-1))) + &
                         state%mesh%full_tangent_wgt(2,j) * (state%mass_flux_lat_n(i  ,j  ) * (state%pv_lon(i,j) + state%pv_lat(i  ,j  ))  + &
                                                             state%mass_flux_lat_n(i+1,j  ) * (state%pv_lon(i,j) + state%pv_lat(i+1,j  )))   &
                        ) * 0.5d0
#endif
      end do
    end do

    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%qhu(i,j) = (state%mesh%half_tangent_wgt(1,j) * (state%mass_flux_lon_n(i-1,j-1) * (state%pv_lat(i,j) + state%pv_lon(i-1,j-1))  + &
                                                             state%mass_flux_lon_n(i  ,j-1) * (state%pv_lat(i,j) + state%pv_lon(i  ,j-1))) + &
                         state%mesh%half_tangent_wgt(2,j) * (state%mass_flux_lon_n(i-1,j  ) * (state%pv_lat(i,j) + state%pv_lon(i-1,j  ))  + &
                                                             state%mass_flux_lon_n(i  ,j  ) * (state%pv_lat(i,j) + state%pv_lon(i  ,j  )))   &
                        ) * 0.5d0
#else
        tend%qhu(i,j) = (state%mesh%half_tangent_wgt(1,j) * (state%mass_flux_lon_n(i-1,j  ) * (state%pv_lat(i,j) + state%pv_lon(i-1,j  ))  + &
                                                             state%mass_flux_lon_n(i  ,j  ) * (state%pv_lat(i,j) + state%pv_lon(i  ,j  ))) + &
                         state%mesh%half_tangent_wgt(2,j) * (state%mass_flux_lon_n(i-1,j+1) * (state%pv_lat(i,j) + state%pv_lon(i-1,j+1))  + &
                                                             state%mass_flux_lon_n(i  ,j+1) * (state%pv_lat(i,j) + state%pv_lon(i  ,j+1)))   &
                        ) * 0.5d0
#endif
      end do
    end do

  end subroutine nonlinear_coriolis_operator

  subroutine energy_gradient_operator(static, state, tend)

    type(static_type), intent(in   ) :: static
    type(state_type ), intent(inout) :: state
    type(tend_type  ), intent(inout) :: tend

    integer i, j

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        tend%dEdlon(i,j) = (state%ke_cell(i+1,j) - state%ke_cell(i,j)  +   &
                            (state%hd    (i+1,j) - state%hd     (i,j)  +   &
                             static%hs   (i+1,j) - static%hs    (i,j)) * g &
                           ) / state%mesh%cell_lon_distance(j)
      end do
    end do

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        tend%dEdlat(i,j) = (state%ke_cell(i,j) - state%ke_cell(i,j-1)  +   &
                            (state%hd    (i,j) - state%hd     (i,j-1)  +   &
                             static%hs   (i,j) - static%hs    (i,j-1)) * g &
#else
        tend%dEdlat(i,j) = (state%ke_cell(i,j+1) - state%ke_cell(i,j)  +   &
                            (state%hd    (i,j+1) - state%hd     (i,j)  +   &
                             static%hs   (i,j+1) - static%hs    (i,j)) * g &
#endif
                           ) / state%mesh%cell_lat_distance(j)
      end do
    end do

  end subroutine energy_gradient_operator

  subroutine mass_flux_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j
    real(real_kind) pole

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        tend%div_mass_flux(i,j) = ((state%mass_flux_lon_n(i,j  ) - state%mass_flux_lon_n(i-1,j)) * mesh%vertex_lat_distance(j) + &
#ifdef STAGGER_V_ON_POLE
                                   (state%mass_flux_lat_n(i,j+1) * mesh%vertex_lon_distance(j+1) -                               &
                                    state%mass_flux_lat_n(i,j  ) * mesh%vertex_lon_distance(j  ))                                &
#else
                                   (state%mass_flux_lat_n(i,j  ) * mesh%vertex_lon_distance(j  ) -                               &
                                    state%mass_flux_lat_n(i,j-1) * mesh%vertex_lon_distance(j-1))                                &
#endif
                                  ) / mesh%cell_area(j)
      end do
    end do
#ifndef STAGGER_V_ON_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%full_lat_start_idx
      pole = 0.0d0
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%mass_flux_lat_n(i,j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole * mesh%vertex_lon_distance(j) / state%mesh%num_full_lon / mesh%cell_area(j)
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        tend%div_mass_flux(i,j) = pole
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%full_lat_end_idx
      pole = 0.0d0
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole - state%mass_flux_lat_n(i,j-1)
      end do
      call parallel_zonal_sum(pole)
      pole = pole * mesh%vertex_lon_distance(j-1) / state%mesh%num_full_lon / mesh%cell_area(j)
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        tend%div_mass_flux(i,j) = pole
      end do
    end if
#endif

  end subroutine mass_flux_divergence_operator

end module operators_mod