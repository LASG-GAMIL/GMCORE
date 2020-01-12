module pv_mod

  use const_mod
  use mesh_mod
  use allocator_mod
  use namelist_mod
  use parallel_mod
  use state_mod
  use tend_mod

  implicit none

  private

  public calc_pv_on_edge_midpoint
  public calc_pv_on_edge_apvm

contains

  subroutine calc_dpv_on_edge(state)

    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    ! Tangent pv difference
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%dpv_lat_t(i,j) = state%pv(i,j) - state%pv(i-1,j)
      end do
    end do
    call fill_halo(mesh, state%dpv_lat_t)

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
        state%dpv_lon_t(i,j) = state%pv(i,j+1) - state%pv(i,j)
#else
        state%dpv_lon_t(i,j) = state%pv(i,j) - state%pv(i,j-1)
#endif
      end do
    end do
    call fill_halo(mesh, state%dpv_lon_t)

    ! Normal pv difference
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        state%dpv_lat_n(i,j) = 0.25_r8 * (state%dpv_lon_t(i-1,j-1) + state%dpv_lon_t(i,j-1) + &
                                          state%dpv_lon_t(i-1,j  ) + state%dpv_lon_t(i,j  ))
#else
        state%dpv_lat_n(i,j) = 0.25_r8 * (state%dpv_lon_t(i-1,j  ) + state%dpv_lon_t(i,j  ) + &
                                          state%dpv_lon_t(i-1,j+1) + state%dpv_lon_t(i,j+1))
#endif
      end do
    end do

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
        state%dpv_lon_n(i,j) = 0.25_r8 * (state%dpv_lat_t(i,j  ) + state%dpv_lat_t(i+1,j  ) + &
                                          state%dpv_lat_t(i,j+1) + state%dpv_lat_t(i+1,j+1))
#else
        state%dpv_lon_n(i,j) = 0.25_r8 * (state%dpv_lat_t(i,j-1) + state%dpv_lat_t(i+1,j-1) + &
                                          state%dpv_lat_t(i,j  ) + state%dpv_lat_t(i+1,j  ))
#endif
      end do
    end do

  end subroutine calc_dpv_on_edge

  subroutine calc_pv_on_edge_midpoint(state)

    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%pv_lat(i,j) = 0.5_r8 * (state%pv(i-1,j) + state%pv(i,j))
      end do 
    end do 
    call fill_halo(mesh, state%pv_lon)

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i,j+1))
#else
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i,j-1))
#endif
      end do 
    end do 
    call fill_halo(mesh, state%pv_lat)

  end subroutine calc_pv_on_edge_midpoint

  subroutine calc_pv_on_edge_apvm(state, dt)

    type(state_type), intent(inout) :: state
    real(r8)        , intent(in   ) :: dt

    type(mesh_type), pointer :: mesh
    real(r8) u, v, le, de
    integer i, j

    call calc_dpv_on_edge(state)

    mesh => state%mesh

    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      le = mesh%le_lat(j)
      de = mesh%de_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        u = state%mf_lat_t(i,j) / state%m_lat(i,j)
        v = state%v(i,j)
        state%pv_lat(i,j) = 0.5_r8 * (state%pv(i,j) + state%pv(i-1,j)) - &
          0.5_r8 * (u * state%dpv_lat_t(i,j) / le + v * state%dpv_lat_n(i,j) / de) * dt
      end do
    end do
#ifdef V_POLE
    state%pv_lat(:,mesh%half_lat_ibeg) = state%pv(:,mesh%half_lat_ibeg)
    state%pv_lat(:,mesh%half_lat_iend  ) = state%pv(:,mesh%half_lat_iend  )
#endif
    call fill_halo(mesh, state%pv_lat)

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      le = mesh%le_lon(j)
      de = mesh%de_lon(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        u = state%u(i,j)
        v = state%mf_lon_t(i,j) / state%m_lon(i,j)
#ifdef V_POLE
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j+1) + state%pv(i,j)) - &
          0.5_r8 * (u * state%dpv_lon_n(i,j) / de + v * state%dpv_lon_t(i,j) / le) * dt
#else
        state%pv_lon(i,j) = 0.5_r8 * (state%pv(i,j-1) + state%pv(i,j)) - &
          0.5_r8 * (u * state%dpv_lon_n(i,j) / de + v * state%dpv_lon_t(i,j) / le) * dt
#endif
      end do
    end do
    call fill_halo(mesh, state%pv_lon)

  end subroutine calc_pv_on_edge_apvm

end module pv_mod