module operators_mod

  use const_mod
  use vert_coord_mod
  use block_mod
  use parallel_mod
  use formula_mod
  use namelist_mod
  use log_mod
  use pgf_mod
  use pv_mod
  use ke_mod
  use pgf_mod
  use interp_mod
  use reduce_mod
  use zonal_damp_mod

  implicit none

  private

  public operators_prepare
  public calc_ph_lev_ph
  public calc_gz_lev_gz
  public calc_t
  public calc_wp
  public calc_wedphdlev
  public calc_div
  public calc_m
  public calc_m_lon_m_lat
  public calc_m_vtx
  public calc_mf_lon_n_mf_lat_n
  public calc_mf_lon_t_mf_lat_t
  public calc_qhu_qhv
  public calc_dkedlon_dkedlat
  public calc_dmfdlon_dmfdlat
  public calc_dptfdlon_dptfdlat
  public calc_dptfdlev
  public calc_dphs
  public calc_wedudlev_wedvdlev

  interface operators_prepare
    module procedure operators_prepare_1
    module procedure operators_prepare_2
  end interface operators_prepare

contains

  subroutine operators_prepare_1(blocks, itime, dt)

    type(block_type), intent(inout) :: blocks(:)
    integer, intent(in) :: itime
    real(r8), intent(in) :: dt

    integer iblk

    do iblk = 1, size(blocks)
      call calc_ph_lev_ph           (blocks(iblk), blocks(iblk)%state(itime))
      call calc_m                   (blocks(iblk), blocks(iblk)%state(itime))
      call calc_t                   (blocks(iblk), blocks(iblk)%state(itime))
      call pgf_prepare              (blocks(iblk), blocks(iblk)%state(itime))
      call calc_m_lon_m_lat         (blocks(iblk), blocks(iblk)%state(itime))
      call calc_m_vtx               (blocks(iblk), blocks(iblk)%state(itime))
      call calc_mf_lon_n_mf_lat_n   (blocks(iblk), blocks(iblk)%state(itime))
      call calc_mf_lon_t_mf_lat_t   (blocks(iblk), blocks(iblk)%state(itime))
      call calc_pv_vtx              (blocks(iblk), blocks(iblk)%state(itime))
      call calc_pv_edge             (blocks(iblk), blocks(iblk)%state(itime), dt)
      call calc_ke_cell             (blocks(iblk), blocks(iblk)%state(itime))
      call calc_gz_lev_gz           (blocks(iblk), blocks(iblk)%state(itime))
      call calc_pt_lon_pt_lat_pt_lev(blocks(iblk), blocks(iblk)%state(itime))
      call calc_div                 (blocks(iblk), blocks(iblk)%state(itime))

      call reduce_run(blocks(iblk), blocks(iblk)%state(itime), dt, all_pass)
    end do

  end subroutine operators_prepare_1

  subroutine operators_prepare_2(block, state, dt, pass)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass

    call calc_ph_lev_ph           (block, state)
    call calc_m                   (block, state)
    call calc_t                   (block, state)
    call pgf_prepare              (block, state)
    call calc_m_lon_m_lat         (block, state)
    call calc_mf_lon_n_mf_lat_n   (block, state)
    call calc_mf_lon_t_mf_lat_t   (block, state)
    call calc_ke_cell             (block, state)
    call calc_gz_lev_gz           (block, state)
    call calc_pt_lon_pt_lat_pt_lev(block, state)
    if (pass == all_pass .or. pass == slow_pass) then
      call calc_m_vtx             (block, state)
      call calc_pv_vtx            (block, state)
      call calc_pv_edge           (block, state, dt)
      call calc_div               (block, state)
    end if

    call reduce_run(block, state, dt, pass)

  end subroutine operators_prepare_2

  subroutine calc_ph_lev_ph(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    if (baroclinic) then
      mesh => state%mesh

      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, state%phs(i,j))
          end do
        end do
      end do
      call fill_halo(block, state%ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%ph(i,j,k) = 0.5_r8 * (state%ph_lev(i,j,k) + state%ph_lev(i,j,k+1))
          end do
        end do
      end do
      call fill_halo(block, state%ph, full_lon=.true., full_lat=.true., full_lev=.true.)
    end if

  end subroutine calc_ph_lev_ph

  subroutine calc_t(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    if (baroclinic) then
      mesh => state%mesh

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%t(i,j,k) = temperature(state%pt(i,j,k), state%ph(i,j,k))
          end do
        end do
      end do
      call fill_halo(block, state%t, full_lon=.true., full_lat=.true., full_lev=.true.)
    end if

  end subroutine calc_t

  subroutine calc_wp(block, state, tend, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, l
    real(r8) mf

    if (baroclinic) then
      mesh => state%mesh

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            mf = 0.5_r8 * (tend%dmfdlon(i,j,k) + tend%dmfdlat(i,j,k))
            do l = 1, k - 1
              mf = mf + tend%dmfdlon(i,j,l) + tend%dmfdlat(i,j,l)
            end do
#ifdef V_POLE
            state%wp(i,j,k) = - mf + 0.5_r8 * (                                                     &
              state%u(i-1,j  ,k) * (state%ph(i  ,j  ,k) - state%ph(i-1,j  ,k)) / mesh%de_lon(j  ) + &
              state%u(i  ,j  ,k) * (state%ph(i+1,j  ,k) - state%ph(i  ,j  ,k)) / mesh%de_lon(j  ) + &
              state%v(i  ,j  ,k) * (state%ph(i  ,j  ,k) - state%ph(i  ,j-1,k)) / mesh%de_lat(j  ) + &
              state%v(i  ,j+1,k) * (state%ph(i  ,j+1,k) - state%ph(i  ,j  ,k)) / mesh%de_lat(j+1)   &
            )
#else
            state%wp(i,j,k) = - mf + 0.5_r8 * (                                                     &
              state%u(i-1,j  ,k) * (state%ph(i  ,j  ,k) - state%ph(i-1,j  ,k)) / mesh%de_lon(j  ) + &
              state%u(i  ,j  ,k) * (state%ph(i+1,j  ,k) - state%ph(i  ,j  ,k)) / mesh%de_lon(j  ) + &
              state%v(i  ,j-1,k) * (state%ph(i  ,j  ,k) - state%ph(i  ,j-1,k)) / mesh%de_lat(j-1) + &
              state%v(i  ,j  ,k) * (state%ph(i  ,j+1,k) - state%ph(i  ,j  ,k)) / mesh%de_lat(j  )   &
            )
#endif
          end do
        end do
      end do
    end if

  end subroutine calc_wp

  subroutine calc_wedphdlev(block, state, tend, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(in) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, l
    real(r8) mf

    if (baroclinic) then
      mesh => state%mesh

      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            mf = 0.0_r8
            do l = 1, k - 1
              mf = mf + tend%dmfdlon(i,j,l) + tend%dmfdlat(i,j,l)
            end do
            state%wedphdlev(i,j,k) = - vert_coord_calc_dphdt_lev(k, tend%dphs(i,j)) - mf
          end do
        end do
      end do
      ! Set vertical boundary conditions.
      state%wedphdlev(:,:,mesh%half_lev_ibeg) = 0.0_r8
      state%wedphdlev(:,:,mesh%half_lev_iend) = 0.0_r8
#ifdef V_POLE
      call fill_halo(block, state%wedphdlev, full_lon=.true., full_lat=.true., full_lev=.false., west_halo=.false., north_halo=.false.)
#else
      call fill_halo(block, state%wedphdlev, full_lon=.true., full_lat=.true., full_lev=.false., west_halo=.false., south_halo=.false.)
#endif

      call interp_cell_to_edge_on_half_level(mesh, state%wedphdlev, state%wedphdlev_lon, state%wedphdlev_lat)
    end if

  end subroutine calc_wedphdlev

  subroutine calc_div(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) pole(global_mesh%num_full_lev)
    integer i, j, k

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          state%div(i,j,k) = (                                                          &
            (state%u(i,j  ,k) * mesh%le_lon(j  ) - state%u(i-1,j,k) * mesh%le_lon(j)) + &
            (state%v(i,j+1,k) * mesh%le_lat(j+1) - state%v(i  ,j,k) * mesh%le_lat(j))   &
          ) / mesh%area_cell(j)
#else
          state%div(i,j,k) = (                                                          &
            (state%u(i,j,k) * mesh%le_lon(j) - state%u(i-1,  j,k) * mesh%le_lon(j  )) + &
            (state%v(i,j,k) * mesh%le_lat(j) - state%v(i  ,j-1,k) * mesh%le_lat(j-1))   &
          ) / mesh%area_cell(j)
#endif
        end do
      end do
    end do
#ifndef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      pole = 0.0_r8
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          pole(k) = pole(k) + state%v(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%div(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      pole = 0.0_r8
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          pole(k) = pole(k) - state%v(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%div(i,j,k) = pole(k)
        end do
      end do
    end if
#endif
    if (div_damp_order == 4) then
      call fill_halo(block, state%div, full_lon=.true., full_lat=.true., full_lev=.true.)
    else
#ifdef V_POLE
      call fill_halo(block, state%div, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., north_halo=.false.)
#else
      call fill_halo(block, state%div, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false.)
#endif
    end if

    if (div_damp_order == 4) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%div2(i,j,k) = (                                                                     &
              state%div(i+1,j,k) - 2 * state%div(i,j,k) + state%div(i-1,j,k)                          &
            ) / mesh%de_lon(j)**2 + (                                                                 &
              (state%div(i,j+1,k) - state%div(i,j  ,k)) * mesh%half_cos_lat(j  ) / mesh%de_lat(j  ) - &
              (state%div(i,j  ,k) - state%div(i,j-1,k)) * mesh%half_cos_lat(j-1) / mesh%de_lat(j-1)   &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
          end do
        end do
      end do
      call fill_halo(block, state%div2, full_lon=.true., full_lat=.true., full_lev=.true.)
    end if

  end subroutine calc_div

  subroutine calc_gz_lev_gz(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k, l
    real(r8) dgz

    if (baroclinic .and. hydrostatic) then
      mesh => state%mesh

      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dgz = 0.0_r8
            do l = k, mesh%num_full_lev
              dgz = dgz + Rd * state%t(i,j,l) * log(state%ph_lev(i,j,l+1) / state%ph_lev(i,j,l))
            end do
            state%gz_lev(i,j,k) = block%static%gzs(i,j) + dgz
          end do
        end do
      end do
      call fill_halo(block, state%gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ! state%gz(i,j,k) = 0.5_r8 * (state%gz_lev(i,j,k) + state%gz_lev(i,j,k+1)) ! Simplified version
            state%gz(i,j,k) = state%gz_lev(i,j,k+1) + state%ak(i,j,k) * Rd * state%t(i,j,k) ! Simmons and Burridge (1981)
          end do
        end do
      end do
      call fill_halo(block, state%gz, full_lon=.true., full_lat=.true., full_lev=.true.)
    end if

  end subroutine calc_gz_lev_gz

  subroutine calc_m(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    if (baroclinic) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%m(i,j,k) = state%ph_lev(i,j,k+1) - state%ph_lev(i,j,k)
#ifndef NDEBUG
            if (state%m(i,j,k) <= 0) then
              print *, i, j, k
              print *, 'phs    =', state%phs(i,j)
              call process_stop(1)
            end if
#endif
          end do
        end do
      end do
    else
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%m(i,j,1) = (state%gz(i,j,1) - block%static%gzs(i,j)) / g
        end do
      end do
    end if
    call fill_halo(block, state%m, full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine calc_m

  subroutine calc_m_lon_m_lat(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh

    mesh => state%mesh

    call interp_cell_to_edge_on_full_level(mesh, state%m, state%m_lon, state%m_lat)
    call fill_halo(block, state%m_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
    call fill_halo(block, state%m_lat, full_lon=.true., full_lat=.false., full_lev=.true.)

  end subroutine calc_m_lon_m_lat

  subroutine calc_pt_lon_pt_lat_pt_lev(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh

    if (baroclinic) then
      mesh => state%mesh

      call interp_cell_to_edge_on_full_level(mesh, state%pt, state%pt_lon, state%pt_lat, reversed_area=.true.)
      call fill_halo(block, state%pt_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
#ifdef V_POLE
      call fill_halo(block, state%pt_lat, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false.)
#else
      call fill_halo(block, state%pt_lat, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)
#endif
      call interp_full_level_to_half_level_on_cell(mesh, state%pt, state%pt_lev)
    end if

  end subroutine calc_pt_lon_pt_lat_pt_lev

  subroutine calc_m_vtx(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k
    real(r8) pole(state%mesh%num_full_lev)

    mesh => state%mesh

    call interp_cell_to_vertex_on_full_level(mesh, state%m, state%m_vtx)

  end subroutine calc_m_vtx

  subroutine calc_mf_lon_n_mf_lat_n(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%mf_lon_n(i,j,k) = state%m_lon(i,j,k) * state%u(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, state%mf_lon_n, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%mf_lat_n(i,j,k) = state%m_lat(i,j,k) * state%v(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, state%mf_lat_n, full_lon=.true., full_lat=.false., full_lev=.true.)

  end subroutine calc_mf_lon_n_mf_lat_n

  subroutine calc_mf_lon_t_mf_lat_t(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          state%mf_lat_t(i,j,k) = mesh%half_tangent_wgt(1,j) * (state%mf_lon_n(i-1,j-1,k) + state%mf_lon_n(i,j-1,k)) + &
                                  mesh%half_tangent_wgt(2,j) * (state%mf_lon_n(i-1,j  ,k) + state%mf_lon_n(i,j  ,k))
#else
          state%mf_lat_t(i,j,k) = mesh%half_tangent_wgt(1,j) * (state%mf_lon_n(i-1,j  ,k) + state%mf_lon_n(i,j  ,k)) + &
                                  mesh%half_tangent_wgt(2,j) * (state%mf_lon_n(i-1,j+1,k) + state%mf_lon_n(i,j+1,k))
#endif
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%mf_lon_t(i,j,k) = mesh%full_tangent_wgt(1,j) * (state%mf_lat_n(i,j  ,k) + state%mf_lat_n(i+1,j  ,k)) + &
                                  mesh%full_tangent_wgt(2,j) * (state%mf_lat_n(i,j+1,k) + state%mf_lat_n(i+1,j+1,k))
#else
          state%mf_lon_t(i,j,k) = mesh%full_tangent_wgt(1,j) * (state%mf_lat_n(i,j-1,k) + state%mf_lat_n(i+1,j-1,k)) + &
                                  mesh%full_tangent_wgt(2,j) * (state%mf_lat_n(i,j  ,k) + state%mf_lat_n(i+1,j  ,k))
#endif
        end do
      end do
    end do

  end subroutine calc_mf_lon_t_mf_lat_t

  subroutine calc_pv_edge(block, state, dt)
    
    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state 
    real(r8), intent(in) :: dt

    select case (pv_scheme)
    case (1)
      call calc_pv_edge_midpoint(block, state)
    case (3)
      call calc_pv_edge_apvm(block, state, dt)
    case default
      call log_error('Unknown PV scheme!')
    end select

  end subroutine calc_pv_edge

  subroutine calc_qhu_qhv(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, move

    mesh => state%mesh

    call calc_pv_edge(block, state, dt)

#ifdef V_POLE
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        if (block%reduced_mesh(j-1)%reduce_factor > 1) then
          tend%qhu(:,j,k) = 0.0_r8
          do move = 1, block%reduced_mesh(j-1)%reduce_factor
            do i = block%reduced_mesh(j-1)%full_lon_ibeg, block%reduced_mesh(j-1)%full_lon_iend
              block%reduced_tend(j-1)%qhu(i,k) = (                    &
                block%reduced_mesh(j-1)%half_tangent_wgt(1,1) * (     &
                  block%reduced_state(j-1)%mf_lon_n(k,i-1,0,move) * ( &
                    block%reduced_state(j-1)%pv_lat(k,i  ,1,move) +   &
                    block%reduced_state(j-1)%pv_lon(k,i-1,0,move)     &
                  ) +                                                 &
                  block%reduced_state(j-1)%mf_lon_n(k,i  ,0,move) * ( &
                    block%reduced_state(j-1)%pv_lat(k,i  ,1,move) +   &
                    block%reduced_state(j-1)%pv_lon(k,i  ,0,move)     &
                  )                                                   &
                )                                                     &
              ) * 0.5_r8
            end do
            call reduce_append_array(move, block%reduced_mesh(j-1), block%reduced_tend(j-1)%qhu(:,k), mesh, tend%qhu(:,j,k))
          end do
          call overlay_inner_halo(block, tend%qhu(:,j,k), west_halo=.true.)
        else
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            if (coriolis_scheme == 1) then
              tend%qhu(i,j,k) = (                                                               &
                mesh%half_tangent_wgt(1,j) * (                                                  &
                  state%mf_lon_n(i-1,j-1,k) * (state%pv_lat(i,j,k) + state%pv_lon(i-1,j-1,k)) + &
                  state%mf_lon_n(i  ,j-1,k) * (state%pv_lat(i,j,k) + state%pv_lon(i  ,j-1,k))   &
                )                                                                               &
              ) * 0.5_r8
            else if (coriolis_scheme == 2) then
              tend%qhu(i,j,k) = state%mf_lat_t(i,j,k) * state%pv_lat(i,j,k)
            end if
          end do
        end if
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          call zero_halo(block, tend%qhu(:,j,k), east_halo=.true.)
          do move = 1, block%reduced_mesh(j)%reduce_factor
            do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
              block%reduced_tend(j)%qhu(i,k) = (                    &
                block%reduced_mesh(j)%half_tangent_wgt(2,0) * (     &
                  block%reduced_state(j)%mf_lon_n(k,i-1,0,move) * ( &
                    block%reduced_state(j)%pv_lat(k,i  ,0,move) +   &
                    block%reduced_state(j)%pv_lon(k,i-1,0,move)     &
                  ) +                                               &
                  block%reduced_state(j)%mf_lon_n(k,i  ,0,move) * ( &
                    block%reduced_state(j)%pv_lat(k,i  ,0,move) +   &
                    block%reduced_state(j)%pv_lon(k,i  ,0,move)     &
                  )                                                 &
                )                                                   &
              ) * 0.5_r8
            end do
            call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhu(:,k), mesh, tend%qhu(:,j,k))
          end do
          call overlay_inner_halo(block, tend%qhu(:,j,k), west_halo=.true.)
        else
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            if (coriolis_scheme == 1) then
              tend%qhu(i,j,k) = tend%qhu(i,j,k) + (                                             &
                mesh%half_tangent_wgt(2,j) * (                                                  &
                  state%mf_lon_n(i-1,j  ,k) * (state%pv_lat(i,j,k) + state%pv_lon(i-1,j  ,k)) + &
                  state%mf_lon_n(i  ,j  ,k) * (state%pv_lat(i,j,k) + state%pv_lon(i  ,j  ,k))   &
                )                                                                               &
              ) * 0.5_r8
            else if (coriolis_scheme == 2) then
              tend%qhu(i,j,k) = state%mf_lat_t(i,j,k) * state%pv_lat(i,j,k)
            end if
          end do
        end if
      end do
    end do
#else
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          tend%qhu(:,j,k) = 0.0_r8
          do move = 1, block%reduced_mesh(j)%reduce_factor
            do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
              block%reduced_tend(j)%qhu(i,k) = (                    &
                block%reduced_mesh(j)%half_tangent_wgt(1,0) * (     &
                  block%reduced_state(j)%mf_lon_n(k,i-1,0,move) * ( &
                    block%reduced_state(j)%pv_lat(k,i  ,0,move) +   &
                    block%reduced_state(j)%pv_lon(k,i-1,0,move)     &
                  ) +                                               &
                  block%reduced_state(j)%mf_lon_n(k,i  ,0,move) * ( &
                    block%reduced_state(j)%pv_lat(k,i  ,0,move) +   &
                    block%reduced_state(j)%pv_lon(k,i  ,0,move)     &
                  )                                                 &
                )                                                   &
              ) * 0.5_r8
            end do
            call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhu(:,k), mesh, tend%qhu(:,j,k))
          end do
          call overlay_inner_halo(block, tend%qhu(:,j,k), west_halo=.true.)
        else
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            if (coriolis_scheme == 1) then
              tend%qhu(i,j,k) = (                                                             &
                mesh%half_tangent_wgt(1,j) * (                                                &
                  state%mf_lon_n(i-1,j  ,k) * (state%pv_lat(i,j,k) + state%pv_lon(i-1,j,k)) + &
                  state%mf_lon_n(i  ,j  ,k) * (state%pv_lat(i,j,k) + state%pv_lon(i  ,j,k))   &
                )                                                                             &
              ) * 0.5_r8
            else if (coriolis_scheme == 2) then
              tend%qhu(i,j,k) = state%mf_lat_t(i,j,k) * state%pv_lat(i,j,k)
            end if
          end do
        end if
        if (block%reduced_mesh(j+1)%reduce_factor > 1) then
          call zero_halo(block, tend%qhu(:,j,k), east_halo=.true.)
          do move = 1, block%reduced_mesh(j+1)%reduce_factor
            do i = block%reduced_mesh(j+1)%full_lon_ibeg, block%reduced_mesh(j+1)%full_lon_iend
              block%reduced_tend(j+1)%qhu(i,k) = (                     &
                block%reduced_mesh(j+1)%half_tangent_wgt(2,-1) * (     &
                  block%reduced_state(j+1)%mf_lon_n(k,i-1, 0,move) * ( &
                    block%reduced_state(j+1)%pv_lat(k,i  ,-1,move) +   &
                    block%reduced_state(j+1)%pv_lon(k,i-1, 0,move)     &
                  ) +                                                  &
                  block%reduced_state(j+1)%mf_lon_n(k,i  , 0,move) * ( &
                    block%reduced_state(j+1)%pv_lat(k,i  ,-1,move) +   &
                    block%reduced_state(j+1)%pv_lon(k,i  , 0,move)     &
                  )                                                    &
                )                                                      &
              ) * 0.5_r8
            end do
            call reduce_append_array(move, block%reduced_mesh(j+1), block%reduced_tend(j+1)%qhu(:,k), mesh, tend%qhu(:,j,k))
          end do
          call overlay_inner_halo(block, tend%qhu(:,j,k), west_halo=.true.)
        else
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            if (coriolis_scheme == 1) then
              tend%qhu(i,j,k) = tend%qhu(i,j,k) + (                                             &
                mesh%half_tangent_wgt(2,j) * (                                                  &
                  state%mf_lon_n(i-1,j+1,k) * (state%pv_lat(i,j,k) + state%pv_lon(i-1,j+1,k)) + &
                  state%mf_lon_n(i  ,j+1,k) * (state%pv_lat(i,j,k) + state%pv_lon(i  ,j+1,k))   &
                )                                                                               &
              ) * 0.5_r8
            else if (coriolis_scheme == 2) then
              tend%qhu(i,j,k) = state%mf_lat_t(i,j,k) * state%pv_lat(i,j,k)
            end if
          end do
        end if
      end do
    end do
#endif

#ifdef V_POLE
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          tend%qhv(:,j,k) = 0.0_r8
          do move = 1, block%reduced_mesh(j)%reduce_factor
            do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
              block%reduced_tend(j)%qhv(i,k) = (                    &
                block%reduced_mesh(j)%full_tangent_wgt(1,0) * (     &
                  block%reduced_state(j)%mf_lat_n(k,i  ,0,move) * ( &
                    block%reduced_state(j)%pv_lon(k,i  ,0,move) +   &
                    block%reduced_state(j)%pv_lat(k,i  ,0,move)     &
                  ) +                                               &
                  block%reduced_state(j)%mf_lat_n(k,i+1,0,move) * ( &
                    block%reduced_state(j)%pv_lon(k,i  ,0,move) +   &
                    block%reduced_state(j)%pv_lat(k,i+1,0,move)     &
                  )                                                 &
                ) +                                                 &
                block%reduced_mesh(j)%full_tangent_wgt(2,0) * (     &
                  block%reduced_state(j)%mf_lat_n(k,i  ,1,move) * ( &
                    block%reduced_state(j)%pv_lon(k,i  ,0,move) +   &
                    block%reduced_state(j)%pv_lat(k,i  ,1,move)     &
                  ) +                                               &
                  block%reduced_state(j)%mf_lat_n(k,i+1,1,move) * ( &
                    block%reduced_state(j)%pv_lon(k,i  ,0,move) +   &
                    block%reduced_state(j)%pv_lat(k,i+1,1,move)     &
                  )                                                 &
                )                                                   &
              ) * 0.5_r8
            end do
            call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhv(:,k), mesh, tend%qhv(:,j,k))
          end do
          call overlay_inner_halo(block, tend%qhv(:,j,k), west_halo=.true.)
        else
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            if (coriolis_scheme == 1) then
              tend%qhv(i,j,k) = (                                                               &
                mesh%full_tangent_wgt(1,j) * (                                                  &
                  state%mf_lat_n(i  ,j  ,k) * (state%pv_lon(i,j,k) + state%pv_lat(i  ,j  ,k)) + &
                  state%mf_lat_n(i+1,j  ,k) * (state%pv_lon(i,j,k) + state%pv_lat(i+1,j  ,k))   &
                ) +                                                                             &
                mesh%full_tangent_wgt(2,j) * (                                                  &
                  state%mf_lat_n(i  ,j+1,k) * (state%pv_lon(i,j,k) + state%pv_lat(i  ,j+1,k)) + &
                  state%mf_lat_n(i+1,j+1,k) * (state%pv_lon(i,j,k) + state%pv_lat(i+1,j+1,k))   &
                )                                                                               &
              ) * 0.5_r8
            else if (coriolis_scheme == 2) then
              tend%qhv(i,j,k) = state%mf_lon_t(i,j,k) * state%pv_lon(i,j,k)
            end if
          end do
        end if
      end do
    end do
#else
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          tend%qhv(:,j,k) = 0.0_r8
          do move = 1, block%reduced_mesh(j)%reduce_factor
            do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
              block%reduced_tend(j)%qhv(i,k) = (                     &
                block%reduced_mesh(j)%full_tangent_wgt(1,0) * (      &
                  block%reduced_state(j)%mf_lat_n(k,i  ,-1,move) * ( &
                    block%reduced_state(j)%pv_lon(k,i  , 0,move) +   &
                    block%reduced_state(j)%pv_lat(k,i  ,-1,move)     &
                  ) +                                                &
                  block%reduced_state(j)%mf_lat_n(k,i+1,-1,move) * ( &
                    block%reduced_state(j)%pv_lon(k,i  , 0,move) +   &
                    block%reduced_state(j)%pv_lat(k,i+1,-1,move)     &
                  )                                                  &
                ) +                                                  &
                block%reduced_mesh(j)%full_tangent_wgt(2,0) * (      &
                  block%reduced_state(j)%mf_lat_n(k,i  , 0,move) * ( &
                    block%reduced_state(j)%pv_lon(k,i  , 0,move) +   &
                    block%reduced_state(j)%pv_lat(k,i  , 0,move)     &
                  ) +                                                &
                  block%reduced_state(j)%mf_lat_n(k,i+1, 0,move) * ( &
                    block%reduced_state(j)%pv_lon(k,i  , 0,move) +   &
                    block%reduced_state(j)%pv_lat(k,i+1, 0,move)     &
                  )                                                  &
                )                                                    &
              ) * 0.5_r8
            end do
            call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%qhv(:,k), mesh, tend%qhv(:,j,k))
          end do
          call overlay_inner_halo(block, tend%qhv(:,j,k), west_halo=.true.)
        else
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            if (coriolis_scheme == 1) then
              tend%qhv(i,j,k) = (                                                           &
                mesh%full_tangent_wgt(1,j) * (                                              &
                  state%mf_lat_n(i  ,j-1,k) * (state%pv_lon(i,j,k) + state%pv_lat(i  ,j-1,k)) + &
                  state%mf_lat_n(i+1,j-1,k) * (state%pv_lon(i,j,k) + state%pv_lat(i+1,j-1,k))   &
                ) +                                                                         &
                mesh%full_tangent_wgt(2,j) * (                                              &
                  state%mf_lat_n(i  ,j  ,k) * (state%pv_lon(i,j,k) + state%pv_lat(i  ,j  ,k)) + &
                  state%mf_lat_n(i+1,j  ,k) * (state%pv_lon(i,j,k) + state%pv_lat(i+1,j  ,k))   &
                )                                                                           &
              ) * 0.5_r8
            else if (coriolis_scheme == 2) then
              tend%qhv(i,j,k) = state%mf_lon_t(i,j,k) * state%pv_lon(i,j,k)
            end if
          end do
        end if
      end do
    end do
#endif

  end subroutine calc_qhu_qhv

  subroutine calc_dkedlon_dkedlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, move, cyc

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (.false. .and. block%reduced_mesh(j)%reduce_factor > 1) then
          tend%dkedlon(:,j,k) = 0.0_r8
          do move = 1, block%reduced_mesh(j)%reduce_factor
            do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
              block%reduced_tend(j)%dkedlon(i,k) = (                                            &
                block%reduced_state(j)%ke(k,i+1,0,move) - block%reduced_state(j)%ke(k,i,0,move) &
              ) / block%reduced_mesh(j)%de_lon(0)
            end do
            call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dkedlon(:,k), mesh, tend%dkedlon(:,j,k))
          end do
          call overlay_inner_halo(block, tend%dkedlon(:,j,k), west_halo=.true.)
        else
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            tend%dkedlon(i,j,k) = (state%ke(i+1,j,k) - state%ke(i,j,k)) / mesh%de_lon(j)
          end do
          if (abs(mesh%full_lat_deg(j)) > 85.0_r8) then
            call fill_zonal_halo(block, mesh%lon_halo_width, tend%dkedlon(:,j,k))
            do cyc = 1, 5
              call zonal_damp_1d(block, 4, dt, mesh%half_lon_lb, mesh%half_lon_ub, mesh%lon_halo_width, tend%dkedlon(:,j,k))
            end do
          end if
        end if
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          tend%dkedlat(i,j,k) = (state%ke(i,j  ,k) - state%ke(i,j-1,k)) / mesh%de_lat(j)
#else
          tend%dkedlat(i,j,k) = (state%ke(i,j+1,k) - state%ke(i,j  ,k)) / mesh%de_lat(j)
#endif
        end do
        if (abs(mesh%half_lat_deg(j)) > 85.0_r8) then
          call fill_zonal_halo(block, mesh%lon_halo_width, tend%dkedlat(:,j,k))
          do cyc = 1, 5
            call zonal_damp_1d(block, 4, dt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, tend%dkedlat(:,j,k))
          end do
        end if
      end do
    end do

  end subroutine calc_dkedlon_dkedlat

  subroutine calc_dmfdlon_dmfdlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, move
    real(r8) pole(state%mesh%num_full_lev)

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          tend%dmfdlon(:,j,k) = 0.0_r8
          do move = 1, block%reduced_mesh(j)%reduce_factor
            do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
              block%reduced_tend(j)%dmfdlon(i,k) = (                                                        &
                block%reduced_state(j)%mf_lon_n(k,i,0,move) - block%reduced_state(j)%mf_lon_n(k,i-1,0,move) &
              ) * block%reduced_mesh(j)%le_lon(0) / block%reduced_mesh(j)%area_cell(0)
            end do
            call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dmfdlon(:,k), mesh, tend%dmfdlon(:,j,k))
          end do
          call overlay_inner_halo(block, tend%dmfdlon(:,j,k), west_halo=.true.)
        else
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%dmfdlon(i,j,k) = (                           &
              state%mf_lon_n(i,j,k) - state%mf_lon_n(i-1,j,k) &
            ) * mesh%le_lon(j) / mesh%area_cell(j)
          end do
        end if
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          tend%dmfdlat(i,j,k) = (                        &
            state%mf_lat_n(i,j+1,k) * mesh%le_lat(j+1) - &
            state%mf_lat_n(i,j  ,k) * mesh%le_lat(j  )   &
          ) / mesh%area_cell(j)
#else
          tend%dmfdlat(i,j,k) = (                        &
            state%mf_lat_n(i,j  ,k) * mesh%le_lat(j  ) - &
            state%mf_lat_n(i,j-1,k) * mesh%le_lat(j-1)   &
          ) / mesh%area_cell(j)
#endif
        end do
      end do
    end do

#ifndef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      pole = 0.0_r8
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          pole(k) = pole(k) + state%mf_lat_n(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dmfdlat(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      pole = 0.0_r8
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          pole(k) = pole(k) - state%mf_lat_n(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_comm, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dmfdlat(i,j,k) = pole(k)
        end do
      end do
    end if
#endif

  end subroutine calc_dmfdlon_dmfdlat

  subroutine calc_dptfdlon_dptfdlat(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, move
    real(r8) pole(state%mesh%num_full_lev)

    if (baroclinic) then
      mesh => state%mesh

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (block%reduced_mesh(j)%reduce_factor > 1) then
            tend%dptfdlon(:,j,k) = 0.0_r8
            do move = 1, block%reduced_mesh(j)%reduce_factor
              do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
                block%reduced_tend(j)%dptfdlon(i,k) = (          &
                  block%reduced_state(j)%ptf_lon(k,i  ,0,move) - &
                  block%reduced_state(j)%ptf_lon(k,i-1,0,move)   &
                ) * block%reduced_mesh(j)%le_lon(0) / block%reduced_mesh(j)%area_cell(0)
              end do
              call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%dptfdlon(:,k), mesh, tend%dptfdlon(:,j,k))
            end do
            call overlay_inner_halo(block, tend%dptfdlon(:,j,k), west_halo=.true.)
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dptfdlon(i,j,k) = (                            &
                state%mf_lon_n(i  ,j,k) * state%pt_lon(i  ,j,k) - &
                state%mf_lon_n(i-1,j,k) * state%pt_lon(i-1,j,k)   &
              ) * mesh%le_lon(j) / mesh%area_cell(j)
            end do
          end if
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            tend%dptfdlat(i,j,k) = (                                               &
              state%mf_lat_n(i,j+1,k) * state%pt_lat(i,j+1,k) * mesh%le_lat(j+1) - &
              state%mf_lat_n(i,j  ,k) * state%pt_lat(i,j  ,k) * mesh%le_lat(j  )   &
            ) / mesh%area_cell(j)
#else
            tend%dptfdlat(i,j,k) = (                                               &
              state%mf_lat_n(i,j  ,k) * state%pt_lat(i,j  ,k) * mesh%le_lat(j  ) - &
              state%mf_lat_n(i,j-1,k) * state%pt_lat(i,j-1,k) * mesh%le_lat(j-1)   &
            ) / mesh%area_cell(j)
#endif
          end do
        end do
      end do

#ifndef V_POLE
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pole(k) = pole(k) + state%mf_lat_n(i,j,k) * state%pt_lat(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%dptfdlat(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pole(k) = pole(k) - state%mf_lat_n(i,j-1,k) * state%pt_lat(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%dptfdlat(i,j,k) = pole(k)
          end do
        end do
      end if
#endif
    end if

  end subroutine calc_dptfdlon_dptfdlat

  subroutine calc_dptfdlev(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, move
    real(r8) pole(state%mesh%num_full_lev)

    if (baroclinic) then
      mesh => state%mesh

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%dptfdlev(i,j,k) = state%wedphdlev(i,j,k+1) * state%pt_lev(i,j,k+1) - &
                                   state%wedphdlev(i,j,k  ) * state%pt_lev(i,j,k  )
          end do
        end do
      end do
    end if

  end subroutine calc_dptfdlev

  subroutine calc_dphs(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    tend%dphs = 0.0_r8
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dphs(i,j) = tend%dphs(i,j) - tend%dmfdlon(i,j,k) - tend%dmfdlat(i,j,k)
        end do
      end do
    end do

  end subroutine calc_dphs

  subroutine calc_wedudlev_wedvdlev(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k

    if (baroclinic) then
      mesh => state%mesh

      do k = mesh%full_lev_ibeg + 1, mesh%full_lev_iend - 1
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            tend%wedudlev(i,j,k) = (                                                   &
                state%wedphdlev_lon(i,j,k+1) * (state%u(i,j,k+1) - state%u(i,j,k  )) + &
                state%wedphdlev_lon(i,j,k  ) * (state%u(i,j,k  ) - state%u(i,j,k-1))   &
              ) / state%m_lon(i,j,k) / 2.0_r8
          end do
        end do
      end do
      k = mesh%full_lev_ibeg
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%wedudlev(i,j,k) = (state%wedphdlev_lon(i,j,k+1) * &
            (state%u(i,j,k+1) - state%u(i,j,k  ))) / state%m_lon(i,j,k) / 2.0_r8
        end do
      end do
      k = mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%wedudlev(i,j,k) = (state%wedphdlev_lon(i,j,k  ) * &
            (state%u(i,j,k  ) - state%u(i,j,k-1))) / state%m_lon(i,j,k) / 2.0_r8
        end do
      end do

      do k = mesh%full_lev_ibeg + 1, mesh%full_lev_iend - 1
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%wedvdlev(i,j,k) = (                                                 &
              state%wedphdlev_lat(i,j,k+1) * (state%v(i,j,k+1) - state%v(i,j,k  )) + &
              state%wedphdlev_lat(i,j,k  ) * (state%v(i,j,k  ) - state%v(i,j,k-1))   &
            ) / state%m_lat(i,j,k) / 2.0_r8
          end do
        end do
      end do
      k = mesh%full_lev_ibeg
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%wedvdlev(i,j,k) = (state%wedphdlev_lat(i,j,k+1) * &
            (state%v(i,j,k+1) - state%v(i,j,k  ))) / state%m_lat(i,j,k) / 2.0_r8
        end do
      end do
      k = mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%wedvdlev(i,j,k) = (state%wedphdlev_lat(i,j,k  ) * &
            (state%v(i,j,k  ) - state%v(i,j,k-1))) / state%m_lat(i,j,k) / 2.0_r8
        end do
      end do
    end if

  end subroutine calc_wedudlev_wedvdlev

end module operators_mod
