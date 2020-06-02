module pv_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use state_mod
  use block_mod
  use parallel_mod
  use damp_mod

  implicit none

  private

  public calc_pv_vtx
  public calc_pv_edge_midpoint
  public calc_pv_edge_apvm

contains

  subroutine calc_pv_vtx(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) pole
    integer i, j

    mesh => state%mesh

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
        state%pv(i,j,1) = (                                                                 &
          (                                                                                 &
            state%u(i  ,j-1,1) * mesh%de_lon(j-1) - state%u(i  ,j  ,1) * mesh%de_lon(j  ) + &
            state%v(i+1,j  ,1) * mesh%de_lat(j  ) - state%v(i  ,j  ,1) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j) + mesh%half_f(j)                                             &
        ) / state%m_vtx(i,j,1)
#else
        state%pv(i,j,1) = (                                                                 &
          (                                                                                 &
            state%u(i  ,j  ,1) * mesh%de_lon(j  ) - state%u(i  ,j+1,1) * mesh%de_lon(j+1) + &
            state%v(i+1,j  ,1) * mesh%de_lat(j  ) - state%v(i  ,j  ,1) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j) + mesh%half_f(j)                                             &
        ) / state%m_vtx(i,j,1)
#endif
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      pole = 0.0_r8
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        pole = pole - state%u(i,j,1) * mesh%de_lon(j)
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / mesh%num_half_lon / mesh%area_vtx(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        state%pv(i,j) = (pole + mesh%half_f(j)) / state%m_vtx(i,j)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      pole = 0.0_r8
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        pole = pole + state%u(i,j-1,1) * mesh%de_lon(j-1)
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / mesh%num_half_lon / mesh%area_vtx(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        state%pv(i,j) = (pole + mesh%half_f(j)) / state%m_vtx(i,j)
      end do
    end if
#else
    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        pole = 0.0_r8
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pole = pole - state%u(i,j+1,1) * mesh%de_lon(j+1)
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole / global_mesh%num_half_lon / mesh%area_vtx(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%pv(i,j,1) = (pole + mesh%half_f(j)) / state%m_vtx(i,j,1)
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        pole = 0.0_r8
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pole = pole + state%u(i,j,1) * mesh%de_lon(j)
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole / global_mesh%num_half_lon / mesh%area_vtx(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%pv(i,j,1) = (pole + mesh%half_f(j)) / state%m_vtx(i,j,1)
        end do
      end if
    end if
#endif
#ifdef V_POLE
    call fill_halo(block, state%pv, full_lon=.false., full_lat=.false., async=state%async(async_pv), south_halo=.false.)
#else
    call fill_halo(block, state%pv, full_lon=.false., full_lat=.false., async=state%async(async_pv), north_halo=.false.)
#endif

  end subroutine calc_pv_vtx

  subroutine calc_dpv_edge(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

    call wait_halo(state%async(async_pv))

    ! Tangent pv difference
!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%dpv_lat_t(i,j,1) = state%pv(i,j,1) - state%pv(i-1,j,1)
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    call fill_halo(block, state%dpv_lat_t, full_lon=.true., full_lat=.false., west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%dpv_lat_t, full_lon=.true., full_lat=.false., west_halo=.false., north_halo=.false.)
#endif

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
        state%dpv_lon_t(i,j,1) = state%pv(i,j+1,1) - state%pv(i,j  ,1)
#else
        state%dpv_lon_t(i,j,1) = state%pv(i,j  ,1) - state%pv(i,j-1,1)
#endif
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    call fill_halo(block, state%dpv_lon_t, full_lon=.false., full_lat=.true., east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%dpv_lon_t, full_lon=.false., full_lat=.true., east_halo=.false., south_halo=.false.)
#endif

    ! Normal pv difference
!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        state%dpv_lat_n(i,j,1) = 0.25_r8 * (state%dpv_lon_t(i-1,j-1,1) + state%dpv_lon_t(i,j-1,1) + &
                                            state%dpv_lon_t(i-1,j  ,1) + state%dpv_lon_t(i,j  ,1))
#else
        state%dpv_lat_n(i,j,1) = 0.25_r8 * (state%dpv_lon_t(i-1,j  ,1) + state%dpv_lon_t(i,j  ,1) + &
                                            state%dpv_lon_t(i-1,j+1,1) + state%dpv_lon_t(i,j+1,1))
#endif
      end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
        state%dpv_lon_n(i,j,1) = 0.25_r8 * (state%dpv_lat_t(i,j  ,1) + state%dpv_lat_t(i+1,j  ,1) + &
                                            state%dpv_lat_t(i,j+1,1) + state%dpv_lat_t(i+1,j+1,1))
#else
        state%dpv_lon_n(i,j,1) = 0.25_r8 * (state%dpv_lat_t(i,j-1,1) + state%dpv_lat_t(i+1,j-1,1) + &
                                            state%dpv_lat_t(i,j  ,1) + state%dpv_lat_t(i+1,j  ,1))
#endif
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine calc_dpv_edge

  subroutine calc_pv_edge_midpoint(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => state%mesh

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        state%pv_lat(i,j,1) = 0.5_r8 * (state%pv(i-1,j,1) + state%pv(i,j,1))
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., async=state%async(async_pv_lat), west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., async=state%async(async_pv_lat), west_halo=.false., north_halo=.false.)
#endif

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
        state%pv_lon(i,j,1) = 0.5_r8 * (state%pv(i,j,1) + state%pv(i,j+1,1))
#else
        state%pv_lon(i,j,1) = 0.5_r8 * (state%pv(i,j,1) + state%pv(i,j-1,1))
#endif
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., async=state%async(async_pv_lon), east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., async=state%async(async_pv_lon), east_halo=.false., south_halo=.false.)
#endif

  end subroutine calc_pv_edge_midpoint

  subroutine calc_pv_edge_apvm(block, state, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    real(r8) u, v
    integer i, j

    call calc_dpv_edge(block, state)

    mesh => state%mesh

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        u = state%mf_lat_t(i,j,1) / state%m_lat(i,j,1)
        v = state%v(i,j,1)
        state%pv_lat(i,j,1) = 0.5_r8 * (                &
          state%pv(i,j,1) + state%pv(i-1,j,1)           &
        ) - 0.5_r8 * (                                  &
          u * state%dpv_lat_t(i,j,1) / mesh%le_lat(j) + &
          v * state%dpv_lat_n(i,j,1) / mesh%de_lat(j)   &
        ) * dt
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    if (mesh%has_south_pole()) state%pv_lat(:,mesh%half_lat_ibeg,1) = state%pv(:,mesh%half_lat_ibeg,1)
    if (mesh%has_north_pole()) state%pv_lat(:,mesh%half_lat_iend,1) = state%pv(:,mesh%half_lat_iend,1)
#endif
#ifdef V_POLE
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., async=state%async(async_pv_lat), west_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., async=state%async(async_pv_lat), west_halo=.false., north_halo=.false.)
#endif

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        u = state%u(i,j,1)
        v = state%mf_lon_t(i,j,1) / state%m_lon(i,j,1)
#ifdef V_POLE
        state%pv_lon(i,j,1) = 0.5_r8 * (                &
          state%pv(i,j+1,1) + state%pv(i,j,1)           &
        ) - 0.5_r8 * (                                  &
          u * state%dpv_lon_n(i,j,1) / mesh%de_lon(j) + &
          v * state%dpv_lon_t(i,j,1) / mesh%le_lon(j)   &
        ) * dt
#else
        state%pv_lon(i,j,1) = 0.5_r8 * (                &
          state%pv(i,j-1,1) + state%pv(i,j,1)           &
        ) - 0.5_r8 * (                                  &
          u * state%dpv_lon_n(i,j,1) / mesh%de_lon(j) + &
          v * state%dpv_lon_t(i,j,1) / mesh%le_lon(j)   &
        ) * dt
#endif
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., async=state%async(async_pv_lon), east_halo=.false., north_halo=.false.)
#else
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., async=state%async(async_pv_lon), east_halo=.false., south_halo=.false.)
#endif

  end subroutine calc_pv_edge_apvm

end module pv_mod
