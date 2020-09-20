module pv_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use state_mod
  use block_mod
  use parallel_mod

  implicit none

  private

  public calc_pv_vtx
  public calc_pv_edge_midpoint
  public calc_pv_edge_apvm

contains

  subroutine calc_vor_vtx(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) pole(state%mesh%num_full_lev)
    integer i, j, k

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%vor(i,j,k) = (                                                                &
            state%u(i  ,j-1,k) * mesh%de_lon(j-1) - state%u(i  ,j  ,k) * mesh%de_lon(j  ) + &
            state%v(i+1,j  ,k) * mesh%de_lat(j  ) - state%v(i  ,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
#else
          state%vor(i,j,k) = (                                                                 &
            state%u(i  ,j  ,k) * mesh%de_lon(j  ) - state%u(i  ,j+1,k) * mesh%de_lon(j+1) + &
            state%v(i+1,j  ,k) * mesh%de_lat(j  ) - state%v(i  ,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
#endif
        end do
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      pole = 0.0_r8
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pole(k) = pole(k) - state%u(i,j,k) * mesh%de_lon(j)
        end do
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / mesh%num_half_lon / mesh%area_vtx(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%vor(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      pole = 0.0_r8
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pole(k) = pole(k) + state%u(i,j-1,k) * mesh%de_lon(j-1)
        end do
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / mesh%num_half_lon / mesh%area_vtx(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%vor(i,j,k) = pole(k)
        end do
      end do
    end if
#else
    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            pole(k) = pole(k) - state%u(i,j+1,k) * mesh%de_lon(j+1)
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole / global_mesh%num_half_lon / mesh%area_vtx(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            state%vor(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            pole(k) = pole(k) + state%u(i,j,k) * mesh%de_lon(j)
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole / global_mesh%num_half_lon / mesh%area_vtx(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            state%vor(i,j,k) = pole(k)
          end do
        end do
      end if
    end if
#endif
    call fill_halo(block, state%vor, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine calc_vor_vtx

  subroutine calc_pv_vtx(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) pole(state%mesh%num_full_lev)
    integer i, j, k

    call calc_vor_vtx(block, state)

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%pv(i,j,k) = (state%vor(i,j,k) + mesh%half_f(j)) / state%m_vtx(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, state%pv, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine calc_pv_vtx

  subroutine calc_dpv_edge(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    ! Tangent pv difference
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%dpv_lat_t(i,j,k) = state%pv(i,j,k) - state%pv(i-1,j,k)
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%dpv_lat_t, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%dpv_lat_t, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)
#endif

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%dpv_lon_t(i,j,k) = state%pv(i,j+1,k) - state%pv(i,j  ,k)
#else
          state%dpv_lon_t(i,j,k) = state%pv(i,j  ,k) - state%pv(i,j-1,k)
#endif
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%dpv_lon_t, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%dpv_lon_t, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)
#endif

    ! Normal pv difference
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          state%dpv_lat_n(i,j,k) = 0.25_r8 * (state%dpv_lon_t(i-1,j-1,k) + state%dpv_lon_t(i,j-1,k) + &
                                              state%dpv_lon_t(i-1,j  ,k) + state%dpv_lon_t(i,j  ,k))
#else
          state%dpv_lat_n(i,j,k) = 0.25_r8 * (state%dpv_lon_t(i-1,j  ,k) + state%dpv_lon_t(i,j  ,k) + &
                                              state%dpv_lon_t(i-1,j+1,k) + state%dpv_lon_t(i,j+1,k))
#endif
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%dpv_lon_n(i,j,k) = 0.25_r8 * (state%dpv_lat_t(i,j  ,k) + state%dpv_lat_t(i+1,j  ,k) + &
                                              state%dpv_lat_t(i,j+1,k) + state%dpv_lat_t(i+1,j+1,k))
#else
          state%dpv_lon_n(i,j,k) = 0.25_r8 * (state%dpv_lat_t(i,j-1,k) + state%dpv_lat_t(i+1,j-1,k) + &
                                              state%dpv_lat_t(i,j  ,k) + state%dpv_lat_t(i+1,j  ,k))
#endif
        end do
      end do
    end do

  end subroutine calc_dpv_edge

  subroutine calc_pv_edge_midpoint(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%pv_lat(i,j,k) = 0.5_r8 * (state%pv(i-1,j,k) + state%pv(i,j,k))
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)
#endif

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%pv_lon(i,j,k) = 0.5_r8 * (state%pv(i,j,k) + state%pv(i,j+1,k))
#else
          state%pv_lon(i,j,k) = 0.5_r8 * (state%pv(i,j,k) + state%pv(i,j-1,k))
#endif
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)
#endif

  end subroutine calc_pv_edge_midpoint

  subroutine calc_pv_edge_apvm(block, state, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    real(r8) u, v
    integer i, j, k

    call calc_dpv_edge(block, state)

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          u = state%mf_lat_t(i,j,k) / state%m_lat(i,j,k)
          v = state%v(i,j,k)
          state%pv_lat(i,j,k) = 0.5_r8 * (                &
            state%pv(i,j,k) + state%pv(i-1,j,k)           &
          ) - 0.5_r8 * (                                  &
            u * state%dpv_lat_t(i,j,k) / mesh%le_lat(j) + &
            v * state%dpv_lat_n(i,j,k) / mesh%de_lat(j)   &
          ) * dt
        end do
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) state%pv_lat(:,mesh%half_lat_ibeg,:) = state%pv(:,mesh%half_lat_ibeg,:)
    if (mesh%has_north_pole()) state%pv_lat(:,mesh%half_lat_iend,:) = state%pv(:,mesh%half_lat_iend,:)
#endif
#ifdef V_POLE
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)
#endif

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          u = state%u(i,j,k)
          v = state%mf_lon_t(i,j,k) / state%m_lon(i,j,k)
#ifdef V_POLE
          state%pv_lon(i,j,k) = 0.5_r8 * (                &
            state%pv(i,j+1,k) + state%pv(i,j,k)           &
          ) - 0.5_r8 * (                                  &
            u * state%dpv_lon_n(i,j,k) / mesh%de_lon(j) + &
            v * state%dpv_lon_t(i,j,k) / mesh%le_lon(j)   &
          ) * dt
#else
          state%pv_lon(i,j,k) = 0.5_r8 * (                &
            state%pv(i,j-1,k) + state%pv(i,j,k)           &
          ) - 0.5_r8 * (                                  &
            u * state%dpv_lon_n(i,j,k) / mesh%de_lon(j) + &
            v * state%dpv_lon_t(i,j,k) / mesh%le_lon(j)   &
          ) * dt
#endif
        end do
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)
#endif

  end subroutine calc_pv_edge_apvm

end module pv_mod
