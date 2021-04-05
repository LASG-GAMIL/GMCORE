module pv_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use state_mod
  use block_mod
  use parallel_mod
  use upwind_mod
  use weno_mod

  implicit none

  private

  public calc_vor
  public diag_pv
  public interp_pv_midpoint
  public interp_pv_upwind
  public interp_pv_weno
  public interp_pv_apvm

contains

  subroutine calc_vor(block, state, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    real(r8) work(state%mesh%half_lon_ibeg:state%mesh%half_lon_iend,state%mesh%num_full_lev)
    real(r8) pole(state%mesh%num_full_lev)
    integer i, j, k

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%vor(i,j,k) = (                                                                 &
            state%u(i  ,j  ,k) * mesh%de_lon(j  ) - state%u(i  ,j+1,k) * mesh%de_lon(j+1) + &
            state%v(i+1,j  ,k) * mesh%de_lat(j  ) - state%v(i  ,j  ,k) * mesh%de_lat(j  )   &
          ) / mesh%area_vtx(j)
        end do
      end do
    end do
    if (pv_pole_stokes) then
      ! Special treatment of vorticity around Poles
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            work(i,k) = - state%u(i,j+1,k) * mesh%de_lon(j+1)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole / global_mesh%num_half_lon / mesh%area_vtx(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            state%vor(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            work(i,k) = state%u(i,j,k) * mesh%de_lon(j)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole / global_mesh%num_half_lon / mesh%area_vtx(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            state%vor(i,j,k) = pole(k)
          end do
        end do
      end if
    end if
    call fill_halo(block, state%vor, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine calc_vor

  subroutine diag_pv(block, state, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    real(r8) pole(state%mesh%num_full_lev)
    integer i, j, k

    call calc_vor(block, state, dt)

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%pv(i,j,k) = (state%vor(i,j,k) + mesh%half_f(j)) / state%m_vtx(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, state%pv, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine diag_pv

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
    call fill_halo(block, state%dpv_lat_t, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., north_halo=.false.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%dpv_lon_t(i,j,k) = state%pv(i,j  ,k) - state%pv(i,j-1,k)
        end do
      end do
    end do
    call fill_halo(block, state%dpv_lon_t, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false., south_halo=.false.)

    ! Normal pv difference
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%dpv_lat_n(i,j,k) = 0.25_r8 * (state%dpv_lon_t(i-1,j  ,k) + state%dpv_lon_t(i,j  ,k) + &
                                              state%dpv_lon_t(i-1,j+1,k) + state%dpv_lon_t(i,j+1,k))
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%dpv_lon_n(i,j,k) = 0.25_r8 * (state%dpv_lat_t(i,j-1,k) + state%dpv_lat_t(i+1,j-1,k) + &
                                              state%dpv_lat_t(i,j  ,k) + state%dpv_lat_t(i+1,j  ,k))
        end do
      end do
    end do

  end subroutine calc_dpv_edge

  subroutine interp_pv_midpoint(block, state)

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
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%pv_lon(i,j,k) = 0.5_r8 * (state%pv(i,j,k) + state%pv(i,j-1,k))
        end do
      end do
    end do
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false.)

  end subroutine interp_pv_midpoint

  subroutine interp_pv_upwind(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    real(r8) ut, vt, b
    integer i, j, k

    associate (mesh     => block%mesh    , &
               un       => state%u       , & ! in
               vn       => state%v       , & ! in
               m_lon    => state%m_lon   , & ! in
               m_lat    => state%m_lat   , & ! in
               mf_lon_t => state%mf_lon_t, & ! in
               mf_lat_t => state%mf_lat_t, & ! in
               pv       => state%pv      , & ! in
               pv_lon   => state%pv_lon  , & ! out
               pv_lat   => state%pv_lat)     ! out
      select case (upwind_order_pv)
      case (3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            if (mesh%is_full_lat_next_to_pole(j)) then
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                vt = mf_lon_t(i,j,k) / m_lon(i,j,k)
                b  = abs(vt) / (sqrt(un(i,j,k)**2 + vt**2) + eps)
                pv_lon(i,j,k) = b * upwind1(sign(1.0_r8, vt), upwind_wgt_pv, pv(i,j-1:j,k)) + &
                                (1 - b) * 0.5_r8 * (pv(i,j-1,k) + pv(i,j,k))
              end do
            else
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                vt = mf_lon_t(i,j,k) / m_lon(i,j,k)
                b  = abs(vt) / (sqrt(un(i,j,k)**2 + vt**2) + eps)
                pv_lon(i,j,k) = b * upwind3(sign(1.0_r8, vt), upwind_wgt_pv, pv(i,j-2:j+1,k)) + &
                                (1 - b) * 0.5_r8 * (pv(i,j-1,k) + pv(i,j,k))
              end do
            end if
          end do
        end do
        call fill_halo(block, pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false.)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              ut = mf_lat_t(i,j,k) / m_lat(i,j,k)
              b  = abs(ut) / (sqrt(ut**2 + vn(i,j,k)**2) + eps)
              pv_lat(i,j,k) = b * upwind3(sign(1.0_r8, ut), upwind_wgt_pv, pv(i-2:i+1,j,k)) + &
                              (1 - b) * 0.5_r8 * (pv(i-1,j,k) + pv(i,j,k))
            end do
          end do
        end do
        call fill_halo(block, pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)
      end select
    end associate

  end subroutine interp_pv_upwind

  subroutine interp_pv_weno(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j, k

    associate (mesh   => block%mesh    , &
               u      => state%mf_lat_t, &
               v      => state%mf_lon_t, &
               pv     => state%pv      , &
               pv_lon => state%pv_lon  , &
               pv_lat => state%pv_lat)
      select case (weno_order_pv)
      case (3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            if (mesh%is_full_lat_next_to_pole(j)) then
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                pv_lon(i,j,k) = upwind1(sign(1.0_r8, v(i,j,k)), upwind_wgt_pv, pv(i,j-1:j,k))
              end do
            else
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                pv_lon(i,j,k) = weno3(sign(1.0_r8, v(i,j,k)), pv(i,j-2:j+1,k))
              end do
            end if
          end do
        end do
        call fill_halo(block, pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false.)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              pv_lat(i,j,k) = weno3(sign(1.0_r8, u(i,j,k)), pv(i-2:i+1,j,k))
            end do
          end do
        end do
        call fill_halo(block, pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)
      end select
    end associate

  end subroutine interp_pv_weno

  subroutine interp_pv_apvm(block, state, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: dt

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
    call fill_halo(block, state%pv_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          u = state%u(i,j,k)
          v = state%mf_lon_t(i,j,k) / state%m_lon(i,j,k)
          state%pv_lon(i,j,k) = 0.5_r8 * (                &
            state%pv(i,j-1,k) + state%pv(i,j,k)           &
          ) - 0.5_r8 * (                                  &
            u * state%dpv_lon_n(i,j,k) / mesh%de_lon(j) + &
            v * state%dpv_lon_t(i,j,k) / mesh%le_lon(j)   &
          ) * dt
        end do
      end do
    end do
    call fill_halo(block, state%pv_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false.)

  end subroutine interp_pv_apvm

end module pv_mod
