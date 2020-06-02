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
  use ke_mod
  use reduce_mod

  implicit none

  private

  public operators_prepare
  public calc_m_lon_m_lat
  public calc_m_vtx
  public calc_mf_lon_n_mf_lat_n
  public calc_mf_lon_t_mf_lat_t
  public calc_qhu_qhv
  public calc_dkedlon_dkedlat
  public calc_dpedlon_dpedlat
  public calc_dmfdlon_dmfdlat

  interface operators_prepare
    module procedure operators_prepare_1
    module procedure operators_prepare_2
  end interface operators_prepare

contains

  subroutine operators_prepare_1(blocks, itime)

    type(block_type), intent(inout) :: blocks(:)
    integer, intent(in) :: itime

    integer iblk

    do iblk = 1, size(blocks)
      call calc_m_lon_m_lat(blocks(iblk), blocks(iblk)%state(itime))
      call calc_m_vtx(blocks(iblk), blocks(iblk)%state(itime))
      call calc_mf_lon_n_mf_lat_n(blocks(iblk), blocks(iblk)%state(itime))
      call calc_mf_lon_t_mf_lat_t(blocks(iblk), blocks(iblk)%state(itime))
      call calc_pv_vtx(blocks(iblk), blocks(iblk)%state(itime))
      call calc_ke_cell(blocks(iblk), blocks(iblk)%state(itime))
    end do

  end subroutine operators_prepare_1

  subroutine operators_prepare_2(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call calc_m_lon_m_lat(block, state)
    call calc_m_vtx(block, state)
    call calc_mf_lon_n_mf_lat_n(block, state)
    call calc_mf_lon_t_mf_lat_t(block, state)
    call calc_pv_vtx(block, state)
    call calc_ke_cell(block, state)

  end subroutine operators_prepare_2

  subroutine calc_m_lon_m_lat(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j

    do j = state%mesh%full_lat_lb, state%mesh%full_lat_ub
      do i = state%mesh%full_lon_lb, state%mesh%full_lon_ub
        state%m(i,j,1) = (state%gz(i,j,1) - block%static%gzs(i,j)) / g
      end do
    end do

!$OMP PARALLEL DO COLLAPSE(2)
    do j = state%mesh%full_lat_ibeg_no_pole, state%mesh%full_lat_iend_no_pole
      do i = state%mesh%half_lon_ibeg, state%mesh%half_lon_iend
        state%m_lon(i,j,1) = (state%mesh%area_lon_west(j) * state%m(i,  j,1) + &
                              state%mesh%area_lon_east(j) * state%m(i+1,j,1)   &
                             ) / state%mesh%area_lon(j)
      end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2)
    do j = state%mesh%half_lat_ibeg_no_pole, state%mesh%half_lat_iend_no_pole
      do i = state%mesh%full_lon_ibeg, state%mesh%full_lon_iend
#ifdef V_POLE
        state%m_lat(i,j,1) = (state%mesh%area_lat_north(j) * state%m(i,j  ,1) + &
                              state%mesh%area_lat_south(j) * state%m(i,j-1,1)   &
                             ) / state%mesh%area_lat(j)
#else
        state%m_lat(i,j,1) = (state%mesh%area_lat_north(j) * state%m(i,j+1,1) + &
                              state%mesh%area_lat_south(j) * state%m(i,j  ,1)   &
                             ) / state%mesh%area_lat(j)
#endif
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine calc_m_lon_m_lat

  subroutine calc_m_vtx(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j
    real(r8) pole

!$OMP PARALLEL DO COLLAPSE(2)
    do j = state%mesh%half_lat_ibeg_no_pole, state%mesh%half_lat_iend_no_pole
      do i = state%mesh%half_lon_ibeg, state%mesh%half_lon_iend
#ifdef V_POLE
        state%m_vtx(i,j,1) = (                                                       &
          (state%m(i,j-1,1) + state%m(i+1,j-1,1)) * state%mesh%area_subcell(2,j-1) + &
          (state%m(i,j  ,1) + state%m(i+1,j  ,1)) * state%mesh%area_subcell(1,j  )   &
        ) / state%mesh%area_vtx(j)
#else
        state%m_vtx(i,j,1) = (                                                       &
          (state%m(i,j  ,1) + state%m(i+1,j  ,1)) * state%mesh%area_subcell(2,j  ) + &
          (state%m(i,j+1,1) + state%m(i+1,j+1,1)) * state%mesh%area_subcell(1,j+1)   &
        ) / state%mesh%area_vtx(j)
#endif
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    if (state%mesh%has_south_pole()) then
      j = state%mesh%half_lat_ibeg
      pole = 0.0_r8
      do i = state%mesh%full_lon_ibeg, state%mesh%full_lon_iend
        pole = pole + state%m(i,j,1)
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / state%mesh%num_half_lon
      do i = state%mesh%half_lon_ibeg, state%mesh%half_lon_iend
        state%m_vtx(i,j,1) = pole
      end do
    end if
    if (state%mesh%has_north_pole()) then
      j = state%mesh%half_lat_iend
      pole = 0.0_r8
      do i = state%mesh%full_lon_ibeg, state%mesh%full_lon_iend
        pole = pole + state%m(i,j-1,1)
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole / state%mesh%num_half_lon
      do i = state%mesh%half_lon_ibeg, state%mesh%half_lon_iend
        state%m_vtx(i,j,1) = pole
      end do
    end if
#endif

  end subroutine calc_m_vtx

  subroutine calc_mf_lon_n_mf_lat_n(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j

!$OMP PARALLEL DO COLLAPSE(2)
    do j = state%mesh%full_lat_ibeg_no_pole, state%mesh%full_lat_iend_no_pole
      do i = state%mesh%half_lon_ibeg, state%mesh%half_lon_iend
        state%mf_lon_n(i,j,1) = state%m_lon(i,j,1) * state%u(i,j,1)
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    call fill_halo(block, state%mf_lon_n, full_lon=.false., full_lat=.true., async=state%async(async_mf_lon_n), north_halo=.false.)
#else
    call fill_halo(block, state%mf_lon_n, full_lon=.false., full_lat=.true., async=state%async(async_mf_lon_n), south_halo=.false.)
#endif

!$OMP PARALLEL DO COLLAPSE(2)
    do j = state%mesh%half_lat_ibeg_no_pole, state%mesh%half_lat_iend_no_pole
      do i = state%mesh%full_lon_ibeg, state%mesh%full_lon_iend
        state%mf_lat_n(i,j,1) = state%m_lat(i,j,1) * state%v(i,j,1)
      end do
    end do
!$OMP END PARALLEL DO
#ifdef V_POLE
    call fill_halo(block, state%mf_lat_n, full_lon=.true., full_lat=.false., async=state%async(async_mf_lat_n), south_halo=.false.)
#else
    call fill_halo(block, state%mf_lat_n, full_lon=.true., full_lat=.false., async=state%async(async_mf_lat_n), north_halo=.false.)
#endif

  end subroutine calc_mf_lon_n_mf_lat_n

  subroutine calc_mf_lon_t_mf_lat_t(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j

    call wait_halo(state%async(async_mf_lon_n))
!$OMP PARALLEL DO COLLAPSE(2)
    do j = state%mesh%half_lat_ibeg_no_pole, state%mesh%half_lat_iend_no_pole
      do i = state%mesh%full_lon_ibeg, state%mesh%full_lon_iend
#ifdef V_POLE
        state%mf_lat_t(i,j,1) = state%mesh%half_tangent_wgt(1,j) * (state%mf_lon_n(i-1,j-1,1) + state%mf_lon_n(i,j-1,1)) + &
                                state%mesh%half_tangent_wgt(2,j) * (state%mf_lon_n(i-1,j  ,1) + state%mf_lon_n(i,j  ,1))
#else
        state%mf_lat_t(i,j,1) = state%mesh%half_tangent_wgt(1,j) * (state%mf_lon_n(i-1,j  ,1) + state%mf_lon_n(i,j  ,1)) + &
                                state%mesh%half_tangent_wgt(2,j) * (state%mf_lon_n(i-1,j+1,1) + state%mf_lon_n(i,j+1,1))
#endif
      end do
    end do
!$OMP END PARALLEL DO

    call wait_halo(state%async(async_mf_lat_n))
!$OMP PARALLEL DO COLLAPSE(2)
    do j = state%mesh%full_lat_ibeg_no_pole, state%mesh%full_lat_iend_no_pole
      do i = state%mesh%half_lon_ibeg, state%mesh%half_lon_iend
#ifdef V_POLE
        state%mf_lon_t(i,j,1) = state%mesh%full_tangent_wgt(1,j) * (state%mf_lat_n(i,j  ,1) + state%mf_lat_n(i+1,j  ,1)) + &
                                state%mesh%full_tangent_wgt(2,j) * (state%mf_lat_n(i,j+1,1) + state%mf_lat_n(i+1,j+1,1))
#else
        state%mf_lon_t(i,j,1) = state%mesh%full_tangent_wgt(1,j) * (state%mf_lat_n(i,j-1,1) + state%mf_lat_n(i+1,j-1,1)) + &
                                state%mesh%full_tangent_wgt(2,j) * (state%mf_lat_n(i,j  ,1) + state%mf_lat_n(i+1,j  ,1))
#endif
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine calc_mf_lon_t_mf_lat_t

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
      call calc_pv_edge_midpoint(block, state)
    case (3)
      call calc_pv_edge_apvm(block, state, dt)
    case default
      call log_error('Unknown PV scheme!')
    end select

    call wait_halo(state%async(async_pv_lon))
!$OMP PARALLEL DO
#ifdef V_POLE
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      if (block%reduced_mesh(j-1)%reduce_factor > 0) then
        tend%qhu(:,j,1) = 0.0_r8
        do move = 1, block%reduced_mesh(j-1)%reduce_factor
          call wait_halo(block%reduced_state(j-1)%async(async_pv_lon,0,move))
          do i = block%reduced_mesh(j-1)%full_lon_ibeg, block%reduced_mesh(j-1)%full_lon_iend
            block%reduced_tend(j-1)%qhu(i) = (                    &
              block%reduced_mesh(j-1)%half_tangent_wgt(1,1) * (   &
                block%reduced_state(j-1)%mf_lon_n(i-1,0,move) * ( &
                  block%reduced_state(j-1)%pv_lat(i  ,1,move) +   &
                  block%reduced_state(j-1)%pv_lon(i-1,0,move)     &
                ) +                                               &
                block%reduced_state(j-1)%mf_lon_n(i  ,0,move) * ( &
                  block%reduced_state(j-1)%pv_lat(i  ,1,move) +   &
                  block%reduced_state(j-1)%pv_lon(i  ,0,move)     &
                )                                                 &
              )                                                   &
            ) * 0.5_r8
          end do
          call reduce_append_array(move, block%reduced_mesh(j-1), block%reduced_tend(j-1)%qhu, mesh, tend%qhu(:,j,1))
        end do
        call overlay_inner_halo(block, tend%qhu(:,j,1), west_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          if (coriolis_scheme == 1) then
            tend%qhu(i,j) = (                                                           &
              mesh%half_tangent_wgt(1,j) * (                                            &
                state%mf_lon_n(i-1,j-1) * (state%pv_lat(i,j) + state%pv_lon(i-1,j-1)) + &
                state%mf_lon_n(i  ,j-1) * (state%pv_lat(i,j) + state%pv_lon(i  ,j-1))   &
              )                                                                         &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            tend%qhu(i,j) = state%mf_lat_t(i,j) * state%pv_lat(i,j)
          end if 
        end do
      end if
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        call zero_halo(block, tend%qhu(:,j,1), east_halo=.true.)
        do move = 1, block%reduced_mesh(j)%reduce_factor
          call wait_halo(block%reduced_state(j)%async(async_pv_lon,0,move))
          do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
            block%reduced_tend(j)%qhu(i) = (                    &
              block%reduced_mesh(j)%half_tangent_wgt(2,0) * (   &
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
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhu, mesh, tend%qhu(:,j,1))
        end do
        call overlay_inner_halo(block, tend%qhu(:,j,1), west_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          if (coriolis_scheme == 1) then
            tend%qhu(i,j,1) = tend%qhu(i,j) + (                                         &
              mesh%half_tangent_wgt(2,j) * (                                            &
                state%mf_lon_n(i-1,j  ) * (state%pv_lat(i,j) + state%pv_lon(i-1,j  )) + &
                state%mf_lon_n(i  ,j  ) * (state%pv_lat(i,j) + state%pv_lon(i  ,j  ))   &
              )                                                                         &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            tend%qhu(i,j,1) = state%mf_lat_t(i,j) * state%pv_lat(i,j)
          end if 
        end do
      end if
    end do
#else
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%qhu(:,j,1) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          call wait_halo(block%reduced_state(j)%async(async_pv_lon,0,move))
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
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhu, mesh, tend%qhu(:,j,1))
        end do
        call overlay_inner_halo(block, tend%qhu(:,j,1), west_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          if (coriolis_scheme == 1) then
            tend%qhu(i,j,1) = (                                                             &
              mesh%half_tangent_wgt(1,j) * (                                                &
                state%mf_lon_n(i-1,j  ,1) * (state%pv_lat(i,j,1) + state%pv_lon(i-1,j,1)) + &
                state%mf_lon_n(i  ,j  ,1) * (state%pv_lat(i,j,1) + state%pv_lon(i  ,j,1))   &
              )                                                                             &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            tend%qhu(i,j,1) = state%mf_lat_t(i,j,1) * state%pv_lat(i,j,1)
          end if 
        end do
      end if
      if (block%reduced_mesh(j+1)%reduce_factor > 0) then
        call zero_halo(block, tend%qhu(:,j,1), east_halo=.true.)
        do move = 1, block%reduced_mesh(j+1)%reduce_factor
          call wait_halo(block%reduced_state(j+1)%async(async_pv_lon,0,move))
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
          call reduce_append_array(move, block%reduced_mesh(j+1), block%reduced_tend(j+1)%qhu, mesh, tend%qhu(:,j,1))
        end do
        call overlay_inner_halo(block, tend%qhu(:,j,1), west_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          if (coriolis_scheme == 1) then
            tend%qhu(i,j,1) = tend%qhu(i,j,1) + (                                             &
              mesh%half_tangent_wgt(2,j) * (                                                  &
                state%mf_lon_n(i-1,j+1,1) * (state%pv_lat(i,j,1) + state%pv_lon(i-1,j+1,1)) + &
                state%mf_lon_n(i  ,j+1,1) * (state%pv_lat(i,j,1) + state%pv_lon(i  ,j+1,1))   &
              )                                                                               &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then 
            tend%qhu(i,j,1) = state%mf_lat_t(i,j,1) * state%pv_lat(i,j,1)
          end if 
        end do
      end if
    end do
#endif
!$OMP END PARALLEL DO

    call wait_halo(state%async(async_pv_lat))
!$OMP PARALLEL DO
#ifdef V_POLE
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%qhv(:,j,1) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          call wait_halo(block%reduced_state(j)%async(async_pv_lat,0,move))
          call wait_halo(block%reduced_state(j)%async(async_pv_lat,1,move))
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
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhv, mesh, tend%qhv(:,j,1))
        end do
        call overlay_inner_halo(block, tend%qhv(:,j,1), west_halo=.true.)
      else
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          if (coriolis_scheme == 1) then
            tend%qhv(i,j,1) = (                                                         &
              mesh%full_tangent_wgt(1,j) * (                                            &
                state%mf_lat_n(i  ,j  ) * (state%pv_lon(i,j) + state%pv_lat(i  ,j  )) + &
                state%mf_lat_n(i+1,j  ) * (state%pv_lon(i,j) + state%pv_lat(i+1,j  ))   &
              ) +                                                                       &
              mesh%full_tangent_wgt(2,j) * (                                            &
                state%mf_lat_n(i  ,j+1) * (state%pv_lon(i,j) + state%pv_lat(i  ,j+1)) + &
                state%mf_lat_n(i+1,j+1) * (state%pv_lon(i,j) + state%pv_lat(i+1,j+1))   &
              )                                                                         &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            tend%qhv(i,j,1) = state%mf_lon_t(i,j) * state%pv_lon(i,j)
          end if
        end do
      end if
    end do
#else
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%qhv(:,j,1) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          call wait_halo(block%reduced_state(j)%async(async_pv_lat,-1,move))
          call wait_halo(block%reduced_state(j)%async(async_pv_lat, 0,move))
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
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhv, mesh, tend%qhv(:,j,1))
        end do
        call overlay_inner_halo(block, tend%qhv(:,j,1), west_halo=.true.)
      else
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          if (coriolis_scheme == 1) then
            tend%qhv(i,j,1) = (                                                           &
              mesh%full_tangent_wgt(1,j) * (                                              &
                state%mf_lat_n(i  ,j-1,1) * (state%pv_lon(i,j,1) + state%pv_lat(i  ,j-1,1)) + &
                state%mf_lat_n(i+1,j-1,1) * (state%pv_lon(i,j,1) + state%pv_lat(i+1,j-1,1))   &
              ) +                                                                         &
              mesh%full_tangent_wgt(2,j) * (                                              &
                state%mf_lat_n(i  ,j  ,1) * (state%pv_lon(i,j,1) + state%pv_lat(i  ,j  ,1)) + &
                state%mf_lat_n(i+1,j  ,1) * (state%pv_lon(i,j,1) + state%pv_lat(i+1,j  ,1))   &
              )                                                                           &
            ) * 0.5_r8
          else if (coriolis_scheme == 2) then
            tend%qhv(i,j,1) = state%mf_lon_t(i,j,1) * state%pv_lon(i,j,1)
          end if 
        end do
      end if
    end do
#endif
!$OMP END PARALLEL DO

  end subroutine calc_qhu_qhv

  subroutine calc_dkedlon_dkedlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move

    mesh => state%mesh

    call wait_halo(state%async(async_ke))
!$OMP PARALLEL DO
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%dkedlon(:,j,1) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          call wait_halo(block%reduced_state(j)%async(async_ke,0,move))
          do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
            block%reduced_tend(j)%dkedlon(i) = (                                    &
              block%reduced_state(j)%ke(i+1,0,move) - block%reduced_state(j)%ke(i,0,move) &
            ) / block%reduced_mesh(j)%de_lon(0)
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dkedlon, mesh, tend%dkedlon(:,j,1))
        end do
        call overlay_inner_halo(block, tend%dkedlon(:,j,1), west_halo=.true.)
      else
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%dkedlon(i,j,1) = (state%ke(i+1,j,1) - state%ke(i,j,1)) / mesh%de_lon(j)
        end do
      end if
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        tend%dkedlat(i,j,1) = (state%ke(i,j  ,1) - state%ke(i,j-1,1)) / mesh%de_lat(j)
#else
        tend%dkedlat(i,j,1) = (state%ke(i,j+1,1) - state%ke(i,j  ,1)) / mesh%de_lat(j)
#endif
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine calc_dkedlon_dkedlat

  subroutine calc_dpedlon_dpedlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, move

    mesh => state%mesh

    call wait_halo(state%async(async_gz))
!$OMP PARALLEL DO
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%dpedlon(:,j,1) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
            block%reduced_tend(j)%dpedlon(i) = (                                          &
              block%reduced_state(j)%gz(i+1,0,move) - block%reduced_state(j)%gz(i,0,move) &
            ) / block%reduced_mesh(j)%de_lon(0)
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dpedlon, mesh, tend%dpedlon(:,j,1))
        end do
        call overlay_inner_halo(block, tend%dpedlon(:,j,1), west_halo=.true.)
      else
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%dpedlon(i,j,1) = (state%gz(i+1,j,1) - state%gz(i,j,1)) / mesh%de_lon(j)
        end do
      end if
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        tend%dpedlat(i,j,1) = (state%gz(i,j  ,1) - state%gz(i,j-1,1)) / mesh%de_lat(j)
#else
        tend%dpedlat(i,j,1) = (state%gz(i,j+1,1) - state%gz(i,j  ,1)) / mesh%de_lat(j)
#endif
      end do
    end do
!$OMP END PARALLEL DO

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
!$OMP PARALLEL DO
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        tend%dmfdlon(:,j,1) = 0.0_r8
        do move = 1, block%reduced_mesh(j)%reduce_factor
          do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
            block%reduced_tend(j)%dmfdlon(i) = (                                                      &
              block%reduced_state(j)%mf_lon_n(i,0,move) - block%reduced_state(j)%mf_lon_n(i-1,0,move) &
            ) * block%reduced_mesh(j)%le_lon(0) / block%reduced_mesh(j)%area_cell(0)
          end do
          call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dmfdlon, mesh, tend%dmfdlon(:,j,1))
        end do
        call overlay_inner_halo(block, tend%dmfdlon(:,j,1), west_halo=.true.)
      else
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dmfdlon(i,j,1) = (                           &
            state%mf_lon_n(i,j,1) - state%mf_lon_n(i-1,j,1) &
          ) * mesh%le_lon(j) / mesh%area_cell(j)
        end do
      end if
    end do
!$OMP END PARALLEL DO

    ! --------------------------------------------------------------------------
    !                    Meridional mass flux divergence
!$OMP PARALLEL DO COLLAPSE(2)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
        tend%dmfdlat(i,j,1) = (                        &
          state%mf_lat_n(i,j+1,1) * mesh%le_lat(j+1) - &
          state%mf_lat_n(i,j  ,1) * mesh%le_lat(j  )   &
        ) / mesh%area_cell(j)
#else
        tend%dmfdlat(i,j,1) = (                        &
          state%mf_lat_n(i,j  ,1) * mesh%le_lat(j  ) - &
          state%mf_lat_n(i,j-1,1) * mesh%le_lat(j-1)   &
        ) / mesh%area_cell(j)
#endif
      end do
    end do
!$OMP END PARALLEL DO

#ifndef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      pole = 0.0_r8
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        pole = pole + state%mf_lat_n(i,j,1)
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole * mesh%le_lat(j) / mesh%num_full_lon / mesh%area_cell(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        tend%dmfdlat(i,j,1) = pole
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      pole = 0.0_r8
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        pole = pole - state%mf_lat_n(i,j-1,1)
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole * mesh%le_lat(j-1) / mesh%num_full_lon / mesh%area_cell(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        tend%dmfdlat(i,j,1) = pole
      end do
    end if
#endif

  end subroutine calc_dmfdlon_dmfdlat

end module operators_mod
