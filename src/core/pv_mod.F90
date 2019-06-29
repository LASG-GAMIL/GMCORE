module pv_mod

  use const_mod
  use mesh_mod
  use allocator_mod
  use parallel_mod
  use state_mod
  use tend_mod

  implicit none

  private

  public calc_pv_on_vertex
  public calc_pv_on_edge_midpoint
  public calc_pv_on_edge_upwind

contains

  subroutine calc_pv_on_vertex(state)

    type(state_type), intent(inout) :: state

    integer i, j
    real(real_kind) pole

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%mass_vertex(i,j) = ((state%gd(i,j-1) + state%gd(i+1,j-1)) * state%mesh%subcell_area(2,j-1) + &
                                  (state%gd(i,j  ) + state%gd(i+1,j  )) * state%mesh%subcell_area(1,j  )   &
                                 ) / state%mesh%vertex_area(j) / g
#else
        state%mass_vertex(i,j) = ((state%gd(i,j  ) + state%gd(i+1,j  )) * state%mesh%subcell_area(2,j  ) + &
                                  (state%gd(i,j+1) + state%gd(i+1,j+1)) * state%mesh%subcell_area(1,j+1)   &
                                 ) / state%mesh%vertex_area(j) / g
#endif
      end do
    end do
#ifdef STAGGER_V_ON_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%half_lat_start_idx
      pole = 0.0d0
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%gd(i,j)
      end do
      state%mass_vertex(:,j) = pole / state%mesh%num_half_lon / g
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%half_lat_end_idx
      pole = 0.0d0
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        pole = pole + state%gd(i,j-1)
      end do
      state%mass_vertex(:,j) = pole / state%mesh%num_half_lon / g
    end if
#endif

    do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%pv(i,j) = ((state%u(i  ,j-1) * state%mesh%cell_lon_distance(j-1) - &
                          state%u(i  ,j  ) * state%mesh%cell_lon_distance(j  ) + &
                          state%v(i+1,j  ) * state%mesh%cell_lat_distance(j  ) - &
                          state%v(i  ,j  ) * state%mesh%cell_lat_distance(j  )   &
                         ) / state%mesh%vertex_area(j) + state%mesh%half_f(j)    &
                        ) / state%mass_vertex(i,j)
#else
        state%pv(i,j) = ((state%u(i  ,j  ) * state%mesh%cell_lon_distance(j  ) - &
                          state%u(i  ,j+1) * state%mesh%cell_lon_distance(j+1) + &
                          state%v(i+1,j  ) * state%mesh%cell_lat_distance(j  ) - &
                          state%v(i  ,j  ) * state%mesh%cell_lat_distance(j  )   &
                         ) / state%mesh%vertex_area(j) + state%mesh%half_f(j)    &
                        ) / state%mass_vertex(i,j)
#endif
      end do
    end do
#ifdef STAGGER_V_ON_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%half_lat_start_idx
      pole = 0.0d0
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        pole = pole - state%u(i,j) * state%mesh%cell_lon_distance(j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / state%mesh%vertex_area(j)
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%pv(i,j) = (pole + state%mesh%half_f(j)) / state%mass_vertex(i,j)
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%half_lat_end_idx
      pole = 0.0d0
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        pole = pole + state%u(i,j-1) * state%mesh%cell_lon_distance(j-1)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / state%mesh%vertex_area(j)
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%pv(i,j) = (pole + state%mesh%half_f(j)) / state%mass_vertex(i,j)
      end do
    end if
#else
    ! Special treatment of vorticity around Poles
    if (state%mesh%has_south_pole()) then
      j = state%mesh%half_lat_start_idx
      pole = 0.0d0
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        pole = pole - state%u(i,j+1) * state%mesh%cell_lon_distance(j+1)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / state%mesh%vertex_area(j)
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%pv(i,j) = (pole + state%mesh%half_f(j)) / state%mass_vertex(i,j)
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%half_lat_end_idx
      pole = 0.0d0
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        pole = pole + state%u(i,j) * state%mesh%cell_lon_distance(j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / state%mesh%num_half_lon / state%mesh%vertex_area(j)
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%pv(i,j) = (pole + state%mesh%half_f(j)) / state%mass_vertex(i,j)
      end do
    end if
#endif
    call parallel_fill_halo(state%mesh, state%pv, all_halo=.true.)

  end subroutine calc_pv_on_vertex

  subroutine calc_pv_on_edge_midpoint(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%pv_lat(i,j) = 0.5 * (state%pv(i-1,j) + state%pv(i,j))
      end do 
    end do 
      
    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%pv_lon(i,j) = 0.5 * (state%pv(i,j) + state%pv(i,j+1))
#else
        state%pv_lon(i,j) = 0.5 * (state%pv(i,j) + state%pv(i,j-1))
#endif
      end do 
    end do 

    call parallel_fill_halo(state%mesh, state%pv_lon, all_halo = .true.)
    call parallel_fill_halo(state%mesh, state%pv_lat, all_halo = .true.)

  end subroutine calc_pv_on_edge_midpoint

  subroutine calc_pv_on_edge_upwind(state)

    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%pv_lat(i,j) = 0.5 * (state%pv(i,j) + state%pv(i-1,j)) - state%mesh%half_upwind_beta(j) * &
                            0.5 * (state%pv(i,j) - state%pv(i-1,j)) * &
                            sign(1.0_real_kind, state%mass_flux_lon_t(i,j))  
      end do 
    end do 

    do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole 
     do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
#ifdef STAGGER_V_ON_POLE
        state%pv_lon(i,j) = 0.5 * (state%pv(i,j+1) + state%pv(i,j)) - state%mesh%full_upwind_beta(j) * &
                            0.5 * (state%pv(i,j+1) - state%pv(i,j)) * &
                            sign(1.0_real_kind, state%mass_flux_lat_t(i,j))
#else
        state%pv_lon(i,j) = 0.5 * (state%pv(i,j) + state%pv(i,j-1)) - state%mesh%full_upwind_beta(j) * &
                            0.5 * (state%pv(i,j) - state%pv(i,j-1)) * &
                            sign(1.0_real_kind, state%mass_flux_lat_t(i,j))
#endif
      end do 
    end do

    call parallel_fill_halo(state%mesh, state%pv_lon, all_halo = .true.)
    call parallel_fill_halo(state%mesh, state%pv_lat, all_halo = .true.)

  end subroutine calc_pv_on_edge_upwind

end module pv_mod