module reduce_mod

  ! Here we only work on full latitudes.
  !
  ! ===: full lat
  ! ---: half lat
  !
  !    V on Poles           U on Poles
  !
  ! ---q---v---q---  2    ===u===h===u===  1
  !
  ! ===u===h===u===  1    ---q---v---q---  0
  !
  ! ---q---v---q---  1    ===u===h===u===  0
  !
  ! ===u===h===u===  0    ---q---v---q--- -1
  !
  ! ---q---v---q---  0    ===u===h===u=== -1    South Pole
  !
  ! =============== -1    ---------------       virtual boundary

  use flogger
  use string
  use const_mod
  use namelist_mod
  use mesh_mod
  use static_mod
  use state_mod
  use parallel_mod
  use block_mod
  use reduced_types_mod

  implicit none

  private

  public reduce_init
  public reduce_run
  public reduce_append_array

  interface
    subroutine reduce_sub_interface(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)
      import block_type, mesh_type, state_type, reduced_mesh_type, reduced_static_type, reduced_state_type, r8
      integer, intent(in) :: j
      integer, intent(in) :: buf_j
      integer, intent(in) :: move
      type(block_type), intent(in) :: block
      type(mesh_type), intent(in) :: raw_mesh
      type(state_type), intent(inout) :: raw_state
      type(reduced_mesh_type), intent(in) :: reduced_mesh
      type(reduced_static_type), intent(in) :: reduced_static
      type(reduced_state_type), intent(inout) :: reduced_state
      real(r8), intent(in) :: dt
    end subroutine reduce_sub_interface
  end interface

contains

  subroutine reduce_init(blocks)

    type(block_type), intent(inout) :: blocks(:)

    integer iblk, itime, j, full_j

    do iblk = 1, size(blocks)
      allocate(blocks(iblk)%reduced_mesh  (blocks(iblk)%mesh%full_lat_lb    :blocks(iblk)%mesh%full_lat_ub    ))
      allocate(blocks(iblk)%reduced_static(blocks(iblk)%mesh%full_lat_ibeg-1:blocks(iblk)%mesh%full_lat_iend+1))
      allocate(blocks(iblk)%reduced_state (blocks(iblk)%mesh%full_lat_ibeg-1:blocks(iblk)%mesh%full_lat_iend+1))
      allocate(blocks(iblk)%reduced_tend  (blocks(iblk)%mesh%full_lat_ibeg-1:blocks(iblk)%mesh%full_lat_iend+1))
      do j = 1, size(reduce_factors)
        if (reduce_factors(j) == 0) cycle
        if (mod(global_mesh%num_full_lon, reduce_factors(j)) /= 0) then
          call log_error('Zonal reduce factor ' // to_string(reduce_factors(j)) // &
            ' cannot divide zonal grid number ' // to_string(global_mesh%num_full_lon) // '!')
        end if
        ! South Pole
#ifdef V_POLE
        full_j = j
#else
        full_j = j + 1
#endif
        if (full_j >= blocks(iblk)%mesh%full_lat_lb .and. full_j <= blocks(iblk)%mesh%full_lat_ub) then
          call reduce_mesh(reduce_factors(j), full_j, blocks(iblk)%mesh, blocks(iblk)%reduced_mesh(full_j))
          blocks(iblk)%reduced_mesh(full_j)%damp_order = damp_orders(j)
        end if
        ! North Pole
#ifdef V_POLE
        full_j = global_mesh%full_lat_iend - j + 1
#else
        full_j = global_mesh%full_lat_iend - j
#endif
        if (full_j >= blocks(iblk)%mesh%full_lat_lb .and. full_j <= blocks(iblk)%mesh%full_lat_ub) then
          call reduce_mesh(reduce_factors(j), full_j, blocks(iblk)%mesh, blocks(iblk)%reduced_mesh(full_j))
          blocks(iblk)%reduced_mesh(full_j)%damp_order = damp_orders(j)
        end if
      end do
      do j = blocks(iblk)%mesh%full_lat_ibeg - 1, blocks(iblk)%mesh%full_lat_iend + 1
        if (blocks(iblk)%reduced_mesh(j)%reduce_factor > 0) then
          call allocate_reduced_static(blocks(iblk)%reduced_mesh(j), blocks(iblk)%reduced_static(j))
          call reduce_static(j, blocks(iblk), blocks(iblk)%mesh, blocks(iblk)%static, blocks(iblk)%reduced_mesh(j), blocks(iblk)%reduced_static(j))
          call allocate_reduced_state(blocks(iblk)%reduced_mesh(j), blocks(iblk)%reduced_state(j))
          call allocate_reduced_tend(blocks(iblk)%reduced_mesh(j), blocks(iblk)%reduced_tend(j))
        end if
      end do
    end do

  end subroutine reduce_init

  subroutine reduce_run(block, state, dt, pass)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass

    integer j

    ! Extend loop range by 1 is for Coriolis forces.
    do j = block%mesh%full_lat_ibeg - 1, block%mesh%full_lat_iend + 1
      if (block%reduced_mesh(j)%reduce_factor > 0) then
        call reduce_state(j, block, block%mesh, state, block%reduced_mesh(j), block%reduced_static(j), block%reduced_state(j), dt, pass)
      end if
    end do

  end subroutine reduce_run

  subroutine reduce_mesh(reduce_factor, j, raw_mesh, reduced_mesh)

    integer, intent(in) :: reduce_factor
    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(reduced_mesh_type), intent(inout) :: reduced_mesh

    real(r8) x(3), y(3), z(3)
    integer i, buf_j

    ! Check if decomposition is OK for reduce.
    if (mod(raw_mesh%num_full_lon / reduce_factor, proc%cart_dims(1)) /= 0) then
      call log_error('Parallel zonal decomposition cannot divide reduced zonal grids!')
    end if

    reduced_mesh%reduce_factor = reduce_factor
    reduced_mesh%halo_width    = 1
    reduced_mesh%num_full_lon  = raw_mesh%num_full_lon / reduce_factor
    reduced_mesh%num_half_lon  = raw_mesh%num_half_lon / reduce_factor
    reduced_mesh%full_lon_ibeg = raw_mesh%full_lon_ibeg / reduce_factor + 1
    reduced_mesh%full_lon_iend = reduced_mesh%full_lon_ibeg + reduced_mesh%num_full_lon - 1
    reduced_mesh%half_lon_ibeg = reduced_mesh%full_lon_ibeg
    reduced_mesh%half_lon_iend = reduced_mesh%half_lon_ibeg + reduced_mesh%num_half_lon - 1
    reduced_mesh%full_lon_lb   = reduced_mesh%full_lon_ibeg - reduced_mesh%halo_width
    reduced_mesh%full_lon_ub   = reduced_mesh%full_lon_iend + reduced_mesh%halo_width
    reduced_mesh%half_lon_lb   = reduced_mesh%half_lon_ibeg - reduced_mesh%halo_width
    reduced_mesh%half_lon_ub   = reduced_mesh%half_lon_iend + reduced_mesh%halo_width

    reduced_mesh%num_full_lev  = raw_mesh%num_full_lev
    reduced_mesh%num_half_lev  = raw_mesh%num_half_lev
    reduced_mesh%full_lev_ibeg = raw_mesh%full_lev_ibeg
    reduced_mesh%full_lev_iend = raw_mesh%full_lev_iend
    reduced_mesh%half_lev_ibeg = raw_mesh%half_lev_ibeg
    reduced_mesh%half_lev_iend = raw_mesh%half_lev_iend
    reduced_mesh%full_lev_lb   = raw_mesh%full_lev_lb
    reduced_mesh%full_lev_ub   = raw_mesh%full_lev_ub
    reduced_mesh%half_lev_lb   = raw_mesh%half_lev_lb
    reduced_mesh%half_lev_ub   = raw_mesh%half_lev_ub

    do buf_j = lbound(reduced_mesh%full_lat, 1), ubound(reduced_mesh%full_lat, 1)
      if (j + buf_j >= raw_mesh%full_lat_lb .and. j + buf_j <= raw_mesh%full_lat_ub) then
        reduced_mesh%full_lat(buf_j) = raw_mesh%full_lat(j+buf_j)
      end if
    end do
    do buf_j = lbound(reduced_mesh%half_lat, 1), ubound(reduced_mesh%half_lat, 1)
      if (j + buf_j >= raw_mesh%half_lat_lb .and. j + buf_j <= raw_mesh%half_lat_ub) then
        reduced_mesh%half_lat(buf_j) = raw_mesh%half_lat(j+buf_j)
      end if
    end do
    do buf_j = lbound(reduced_mesh%half_f, 1), ubound(reduced_mesh%half_f, 1)
      if (j + buf_j >= raw_mesh%half_lat_lb .and. j + buf_j <= raw_mesh%half_lat_ub) then
        reduced_mesh%half_f(buf_j) = raw_mesh%half_f(j+buf_j)
      end if
    end do

    ! Cell area
    do buf_j = lbound(reduced_mesh%area_cell, 1), ubound(reduced_mesh%area_cell, 1)
      if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
        reduced_mesh%area_cell(buf_j) = raw_mesh%area_cell(j+buf_j) * reduce_factor
      end if
    end do
    do buf_j = lbound(reduced_mesh%area_subcell, 2), ubound(reduced_mesh%area_subcell, 2)
      if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
        reduced_mesh%area_subcell(1,buf_j) = raw_mesh%area_subcell(1,j+buf_j) * reduce_factor
        reduced_mesh%area_subcell(2,buf_j) = raw_mesh%area_subcell(2,j+buf_j) * reduce_factor
      end if
    end do
    do buf_j = lbound(reduced_mesh%area_lon, 1), ubound(reduced_mesh%area_lon, 1)
      if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
        reduced_mesh%area_lon_west(buf_j) = raw_mesh%area_lon_west(j+buf_j) * reduce_factor
        reduced_mesh%area_lon_east(buf_j) = raw_mesh%area_lon_east(j+buf_j) * reduce_factor
        reduced_mesh%area_lon     (buf_j) = raw_mesh%area_lon     (j+buf_j) * reduce_factor
      end if
    end do
    ! Vertex area
    do buf_j = lbound(reduced_mesh%area_vtx, 1), ubound(reduced_mesh%area_vtx, 1)
      if (raw_mesh%is_inside_with_halo_half_lat(j+buf_j)) then
        reduced_mesh%area_vtx(buf_j) = raw_mesh%area_vtx(j+buf_j) * reduce_factor
      end if
    end do
    do buf_j = lbound(reduced_mesh%area_lat, 1), ubound(reduced_mesh%area_lat, 1)
      if (raw_mesh%is_inside_with_halo_half_lat(j+buf_j)) then
        reduced_mesh%area_lat_north(buf_j) = raw_mesh%area_lat_north(j+buf_j) * reduce_factor
        reduced_mesh%area_lat_south(buf_j) = raw_mesh%area_lat_south(j+buf_j) * reduce_factor
        reduced_mesh%area_lat      (buf_j) = raw_mesh%area_lat      (j+buf_j) * reduce_factor
      end if
    end do
    ! Edge lengths and cell distances
    do buf_j = lbound(reduced_mesh%le_lat, 1), ubound(reduced_mesh%le_lat, 1)
      if (raw_mesh%is_inside_with_halo_half_lat(j+buf_j)) then
        reduced_mesh%le_lat(buf_j) = raw_mesh%le_lat(j+buf_j) * reduce_factor
      end if
    end do
    do buf_j = lbound(reduced_mesh%de_lat, 1), ubound(reduced_mesh%de_lat, 1)
      if (raw_mesh%is_inside_with_halo_half_lat(j+buf_j)) then
        reduced_mesh%de_lat(buf_j) = raw_mesh%de_lat(j+buf_j)
      end if
    end do
    do buf_j = lbound(reduced_mesh%le_lon, 1), ubound(reduced_mesh%le_lon, 1)
      if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
        reduced_mesh%le_lon(buf_j) = raw_mesh%le_lon(j+buf_j)
      end if
    end do
    do buf_j = lbound(reduced_mesh%de_lon, 1), ubound(reduced_mesh%de_lon, 1)
      if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
        reduced_mesh%de_lon(buf_j) = raw_mesh%de_lon(j+buf_j) * reduce_factor
      end if
    end do

#ifdef V_POLE
    do buf_j = lbound(reduced_mesh%full_tangent_wgt, 2), ubound(reduced_mesh%full_tangent_wgt, 2)
      if (reduced_mesh%le_lat(buf_j  ) /= 0 .and. reduced_mesh%de_lon(buf_j) /= 0) then
        reduced_mesh%full_tangent_wgt(1,buf_j) = reduced_mesh%le_lat(buf_j  ) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lat(buf_j+1) /= 0 .and. reduced_mesh%de_lon(buf_j) /= 0) then
        reduced_mesh%full_tangent_wgt(2,buf_j) = reduced_mesh%le_lat(buf_j+1) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
    end do
#else
    do buf_j = lbound(reduced_mesh%full_tangent_wgt, 2), ubound(reduced_mesh%full_tangent_wgt, 2)
      if (reduced_mesh%le_lat(buf_j-1) /= 0 .and. reduced_mesh%de_lon(buf_j) /= 0) then
        reduced_mesh%full_tangent_wgt(1,buf_j) = reduced_mesh%le_lat(buf_j-1) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lat(buf_j  ) /= 0 .and. reduced_mesh%de_lon(buf_j) /= 0) then
        reduced_mesh%full_tangent_wgt(2,buf_j) = reduced_mesh%le_lat(buf_j  ) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
    end do
#endif
#ifdef V_POLE
    do buf_j = lbound(reduced_mesh%half_tangent_wgt, 2), ubound(reduced_mesh%half_tangent_wgt, 2)
      if (raw_mesh%is_south_pole(j+buf_j+1)) cycle
      if (reduced_mesh%le_lon(buf_j-1) /= 0 .and. reduced_mesh%de_lat(buf_j) /= 0) then
        reduced_mesh%half_tangent_wgt(1,buf_j) = reduced_mesh%le_lon(buf_j-1) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lon(buf_j  ) /= 0 .and. reduced_mesh%de_lat(buf_j) /= 0) then
        reduced_mesh%half_tangent_wgt(2,buf_j) = reduced_mesh%le_lon(buf_j  ) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
    end do
#else
    do buf_j = lbound(reduced_mesh%half_tangent_wgt, 2), ubound(reduced_mesh%half_tangent_wgt, 2)
      if (reduced_mesh%le_lon(buf_j  ) /= 0 .and. reduced_mesh%de_lat(buf_j) /= 0) then
        reduced_mesh%half_tangent_wgt(1,buf_j) = reduced_mesh%le_lon(buf_j  ) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lon(buf_j+1) /= 0 .and. reduced_mesh%de_lat(buf_j) /= 0) then
        reduced_mesh%half_tangent_wgt(2,buf_j) = reduced_mesh%le_lon(buf_j+1) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
    end do
#endif

  end subroutine reduce_mesh

  subroutine allocate_reduced_state(reduced_mesh, reduced_state)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state

#ifdef V_POLE
    allocate(reduced_state%mf_lon_n (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_n (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lon_t (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_t (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%gz       (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m        (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lon    (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lat    (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%u        (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%v        (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv       (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lon   (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lat   (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_t(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_n(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_t(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_n(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%ke       (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:0,reduced_mesh%reduce_factor))
#else
    allocate(reduced_state%mf_lon_n (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_n (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lon_t (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_t (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%gz       (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m        (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lon    (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lat    (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%u        (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%v        (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv       (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lon   (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lat   (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_t(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_n(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_t(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_n(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%ke       (reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:0,reduced_mesh%reduce_factor))
#endif
    allocate(reduced_state%async    (11,-2:2,reduced_mesh%reduce_factor))

  end subroutine allocate_reduced_state

  subroutine reduce_state(j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt, pass)

    integer, intent(in) :: j
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass

    call apply_reduce(lbound(reduced_state%mf_lon_n, 3), ubound(reduced_state%mf_lon_n, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_mf_lon_n, dt)
    call apply_reduce(lbound(reduced_state%mf_lat_n, 3), ubound(reduced_state%mf_lat_n, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_mf_lat_n, dt)
    call apply_reduce(lbound(reduced_state%mf_lon_t, 3), ubound(reduced_state%mf_lon_t, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_mf_lon_t, dt)
    call apply_reduce(lbound(reduced_state%mf_lat_t, 3), ubound(reduced_state%mf_lat_t, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_mf_lat_t, dt)
    call apply_reduce(lbound(reduced_state%gz      , 3), ubound(reduced_state%gz      , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_gz      , dt)
    if (pass == all_pass .or. pass == slow_pass) then
      call apply_reduce(lbound(reduced_state%m_lon    , 3), ubound(reduced_state%m_lon    , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_m_lon      , dt)
      call apply_reduce(lbound(reduced_state%m_lat    , 3), ubound(reduced_state%m_lat    , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_m_lat      , dt)
      call apply_reduce(lbound(reduced_state%u        , 3), ubound(reduced_state%u        , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_u          , dt)
      call apply_reduce(lbound(reduced_state%v        , 3), ubound(reduced_state%v        , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_v          , dt)
      call apply_reduce(lbound(reduced_state%pv       , 3), ubound(reduced_state%pv       , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_pv         , dt)
      call apply_reduce(lbound(reduced_state%dpv_lon_t, 3), ubound(reduced_state%dpv_lon_t, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_dpv_lon_t  , dt)
      call apply_reduce(lbound(reduced_state%dpv_lat_n, 3), ubound(reduced_state%dpv_lat_n, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_dpv_lat_n  , dt)
      call apply_reduce(lbound(reduced_state%dpv_lat_t, 3), ubound(reduced_state%dpv_lat_t, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_dpv_lat_t  , dt)
      call apply_reduce(lbound(reduced_state%dpv_lon_n, 3), ubound(reduced_state%dpv_lon_n, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_dpv_lon_n  , dt)
      call apply_reduce(lbound(reduced_state%pv_lon   , 3), ubound(reduced_state%pv_lon   , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_pv_lon_apvm, dt)
      call apply_reduce(lbound(reduced_state%pv_lat   , 3), ubound(reduced_state%pv_lat   , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_pv_lat_apvm, dt)
    end if
    if (pass == all_pass .or. pass == fast_pass) then
      call apply_reduce(lbound(reduced_state%ke       , 3), ubound(reduced_state%ke       , 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_ke         , dt)
    end if

  end subroutine reduce_state

  subroutine reduce_gzs(j, buf_j, move, block, raw_mesh, raw_static, reduced_mesh, reduced_static)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(static_type), intent(in) :: raw_static
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    integer raw_i, i

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        reduced_static%gzs(i,buf_j,move) = sum(raw_static%gzs(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      reduced_static%gzs(:,buf_j,move) = reduced_static%gzs(:,buf_j,move) / reduced_mesh%reduce_factor
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_static%gzs(:,buf_j,move))
    end if

  end subroutine reduce_gzs

  subroutine apply_reduce(buf_lb, buf_ub, j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, reduce_sub, dt)

    integer, intent(in) :: buf_lb
    integer, intent(in) :: buf_ub
    integer, intent(in) :: j
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    procedure(reduce_sub_interface) reduce_sub
    real(r8), intent(in) :: dt

    integer buf_j, move

    do move = 1, reduced_mesh%reduce_factor
      do buf_j = buf_lb, buf_ub
        call reduce_sub(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)
      end do
    end do

  end subroutine apply_reduce

  subroutine reduce_u(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%area_lon(buf_j) == 0) return
    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%u(k,i,buf_j,move) = reduced_state%mf_lon_n(k,i,buf_j,move) / reduced_state%m_lon(k,i,buf_j,move)
      end do
    end do

  end subroutine reduce_u

  subroutine reduce_v(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%area_lat(buf_j) == 0) return
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%v(k,i,buf_j,move) = reduced_state%mf_lat_n(k,i,buf_j,move) / reduced_state%m_lat(k,i,buf_j,move)
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%v(:,:,buf_j,move), west_halo=.false.)

  end subroutine reduce_v

  subroutine reduce_gz(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%gz(k,i,buf_j,move) = sum(raw_state%gz(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k))
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      reduced_state%gz(:,:,buf_j,move) = reduced_state%gz(:,:,buf_j,move) / reduced_mesh%reduce_factor
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%gz(:,:,buf_j,move), west_halo=.false.)
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%m(k,:,buf_j,move) = (reduced_state%gz(k,:,buf_j,move) - reduced_static%gzs(:,buf_j,move)) / g
      end do
    end if

  end subroutine reduce_gz

  subroutine reduce_pv(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    real(r8) m_vtx
    integer raw_i, i, k

#ifdef V_POLE
    if (raw_mesh%is_outside_pole_half_lat(j+buf_j)) then
      return
    else if (raw_mesh%is_south_pole(j+buf_j)) then
      reduced_state%pv(:,:,buf_j,move) = raw_state%pv(raw_mesh%full_lon_ibeg,raw_mesh%half_lat_ibeg,:)
    else if (raw_mesh%is_north_pole(j+buf_j)) then
      reduced_state%pv(:,:,buf_j,move) = raw_state%pv(raw_mesh%full_lon_ibeg,raw_mesh%half_lat_iend,:)
    else
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          m_vtx = (                                                                                                            &
            (reduced_state%m(k,i,buf_j-1,move) + reduced_state%m(k,i+1,buf_j-1,move)) * reduced_mesh%area_subcell(2,buf_j-1) + &
            (reduced_state%m(k,i,buf_j  ,move) + reduced_state%m(k,i+1,buf_j  ,move)) * reduced_mesh%area_subcell(1,buf_j  )   &
          ) / reduced_mesh%area_vtx(buf_j)
          reduced_state%pv(k,i,buf_j,move) = (                                     &
            (                                                                      &
              reduced_state%u(k,i  ,buf_j-1,move) * reduced_mesh%de_lon(buf_j-1) - &
              reduced_state%u(k,i  ,buf_j  ,move) * reduced_mesh%de_lon(buf_j  ) + &
              reduced_state%v(k,i+1,buf_j  ,move) * reduced_mesh%de_lat(buf_j  ) - &
              reduced_state%v(k,i  ,buf_j  ,move) * reduced_mesh%de_lat(buf_j  )   &
            ) / reduced_mesh%area_vtx(buf_j) + reduced_mesh%half_f(buf_j)          &
          ) / m_vtx
        end do
      end do
    end if
#else
    if (raw_mesh%is_outside_pole_half_lat(j+buf_j)) then
      return
    else if (raw_mesh%is_south_pole(j+buf_j) .or. raw_mesh%is_north_pole(j+buf_j+1)) then
      raw_i = raw_mesh%half_lon_ibeg + move - 1
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%pv(k,i,buf_j,move) = sum(raw_state%pv(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k))
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      reduced_state%pv(:,:,buf_j,move) = reduced_state%pv(:,:,buf_j,move) / reduced_mesh%reduce_factor
    else
      if (reduced_mesh%area_vtx(buf_j) == 0) return
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          m_vtx = (                                                                                                            &
            (reduced_state%m(k,i,buf_j  ,move) + reduced_state%m(k,i+1,buf_j  ,move)) * reduced_mesh%area_subcell(2,buf_j  ) + &
            (reduced_state%m(k,i,buf_j+1,move) + reduced_state%m(k,i+1,buf_j+1,move)) * reduced_mesh%area_subcell(1,buf_j+1)   &
          ) / reduced_mesh%area_vtx(buf_j)
          reduced_state%pv(k,i,buf_j,move) = (                                     &
            (                                                                      &
              reduced_state%u(k,i  ,buf_j  ,move) * reduced_mesh%de_lon(buf_j  ) - &
              reduced_state%u(k,i  ,buf_j+1,move) * reduced_mesh%de_lon(buf_j+1) + &
              reduced_state%v(k,i+1,buf_j  ,move) * reduced_mesh%de_lat(buf_j  ) - &
              reduced_state%v(k,i  ,buf_j  ,move) * reduced_mesh%de_lat(buf_j  )   &
            ) / reduced_mesh%area_vtx(buf_j) + reduced_mesh%half_f(buf_j)          &
          ) / m_vtx
        end do
      end do
    end if
#endif
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%pv(:,:,buf_j,move))

  end subroutine reduce_pv

  subroutine reduce_m_lon(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%area_lon(buf_j) == 0) return
    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%m_lon(k,i,buf_j,move) = (                                   &
          reduced_mesh%area_lon_west(buf_j) * reduced_state%m(k,i  ,buf_j,move) + &
          reduced_mesh%area_lon_east(buf_j) * reduced_state%m(k,i+1,buf_j,move)   &
        ) / reduced_mesh%area_lon(buf_j)
      end do
    end do

  end subroutine reduce_m_lon

  subroutine reduce_m_lat(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%area_lat(buf_j) == 0) return
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
#ifdef V_POLE
        reduced_state%m_lat(k,i,buf_j,move) = (                                    &
          reduced_mesh%area_lat_north(buf_j) * reduced_state%m(k,i,buf_j  ,move) + &
          reduced_mesh%area_lat_south(buf_j) * reduced_state%m(k,i,buf_j-1,move)   &
        ) / reduced_mesh%area_lat(buf_j)
#else
        reduced_state%m_lat(k,i,buf_j,move) = (                                    &
          reduced_mesh%area_lat_north(buf_j) * reduced_state%m(k,i,buf_j+1,move) + &
          reduced_mesh%area_lat_south(buf_j) * reduced_state%m(k,i,buf_j  ,move)   &
        ) / reduced_mesh%area_lat(buf_j)
#endif
      end do
    end do

  end subroutine reduce_m_lat

  subroutine reduce_mf_lon_n(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%mf_lon_n(k,i,buf_j,move) = sum(raw_state%mf_lon_n(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k))
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      reduced_state%mf_lon_n(:,:,buf_j,move) = reduced_state%mf_lon_n(:,:,buf_j,move) / reduced_mesh%reduce_factor
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%mf_lon_n(:,:,buf_j,move), east_halo=.false.)
    end if

  end subroutine reduce_mf_lon_n

  subroutine reduce_mf_lat_n(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_half_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%mf_lat_n(k,i,buf_j,move) = sum(raw_state%mf_lat_n(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k))
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      reduced_state%mf_lat_n(:,:,buf_j,move) = reduced_state%mf_lat_n(:,:,buf_j,move) / reduced_mesh%reduce_factor
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%mf_lat_n(:,:,buf_j,move), west_halo=.false.)
    end if

  end subroutine reduce_mf_lat_n

  subroutine reduce_mf_lon_t(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
#ifdef V_POLE
        reduced_state%mf_lon_t(k,i,buf_j,move) =                                                                                             &
          reduced_mesh%full_tangent_wgt(1,buf_j) * (reduced_state%mf_lat_n(k,i,buf_j  ,move) + reduced_state%mf_lat_n(k,i+1,buf_j  ,move)) + &
          reduced_mesh%full_tangent_wgt(2,buf_j) * (reduced_state%mf_lat_n(k,i,buf_j+1,move) + reduced_state%mf_lat_n(k,i+1,buf_j+1,move))
#else
        reduced_state%mf_lon_t(k,i,buf_j,move) =                                                                                             &
          reduced_mesh%full_tangent_wgt(1,buf_j) * (reduced_state%mf_lat_n(k,i,buf_j-1,move) + reduced_state%mf_lat_n(k,i+1,buf_j-1,move)) + &
          reduced_mesh%full_tangent_wgt(2,buf_j) * (reduced_state%mf_lat_n(k,i,buf_j  ,move) + reduced_state%mf_lat_n(k,i+1,buf_j  ,move))
#endif
      end do
    end do

  end subroutine reduce_mf_lon_t

  subroutine reduce_mf_lat_t(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    if (is_inf(reduced_mesh%half_tangent_wgt(1,buf_j)) .or. is_inf(reduced_mesh%half_tangent_wgt(2,buf_j))) return
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
#ifdef V_POLE
        reduced_state%mf_lat_t(k,i,buf_j,move) =                                                                                             &
          reduced_mesh%half_tangent_wgt(1,buf_j) * (reduced_state%mf_lon_n(k,i-1,buf_j-1,move) + reduced_state%mf_lon_n(k,i,buf_j-1,move)) + &
          reduced_mesh%half_tangent_wgt(2,buf_j) * (reduced_state%mf_lon_n(k,i-1,buf_j  ,move) + reduced_state%mf_lon_n(k,i,buf_j  ,move))
#else
        reduced_state%mf_lat_t(k,i,buf_j,move) =                                                                                             &
          reduced_mesh%half_tangent_wgt(1,buf_j) * (reduced_state%mf_lon_n(k,i-1,buf_j  ,move) + reduced_state%mf_lon_n(k,i,buf_j  ,move)) + &
          reduced_mesh%half_tangent_wgt(2,buf_j) * (reduced_state%mf_lon_n(k,i-1,buf_j+1,move) + reduced_state%mf_lon_n(k,i,buf_j+1,move))
#endif
      end do
    end do

  end subroutine reduce_mf_lat_t

  subroutine reduce_dpv_lon_t(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
#ifdef V_POLE
        reduced_state%dpv_lon_t(k,i,buf_j,move) = reduced_state%pv(k,i,buf_j+1,move) - reduced_state%pv(k,i,buf_j  ,move)
#else
        reduced_state%dpv_lon_t(k,i,buf_j,move) = reduced_state%pv(k,i,buf_j  ,move) - reduced_state%pv(k,i,buf_j-1,move)
#endif
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%dpv_lon_t(:,:,buf_j,move))

  end subroutine reduce_dpv_lon_t

  subroutine reduce_dpv_lat_t(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%dpv_lat_t(k,i,buf_j,move) = reduced_state%pv(k,i,buf_j,move) - reduced_state%pv(k,i-1,buf_j,move)
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%dpv_lat_t(:,:,buf_j,move))

  end subroutine reduce_dpv_lat_t

  subroutine reduce_dpv_lon_n(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
#ifdef V_POLE
        reduced_state%dpv_lon_n(k,i,buf_j,move) = 0.25_r8 * ( &
          reduced_state%dpv_lat_t(k,i  ,buf_j  ,move) +       &
          reduced_state%dpv_lat_t(k,i+1,buf_j  ,move) +       &
          reduced_state%dpv_lat_t(k,i  ,buf_j+1,move) +       &
          reduced_state%dpv_lat_t(k,i+1,buf_j+1,move)         &
        )
#else
        reduced_state%dpv_lon_n(k,i,buf_j,move) = 0.25_r8 * ( &
          reduced_state%dpv_lat_t(k,i  ,buf_j-1,move) +       &
          reduced_state%dpv_lat_t(k,i+1,buf_j-1,move) +       &
          reduced_state%dpv_lat_t(k,i  ,buf_j  ,move) +       &
          reduced_state%dpv_lat_t(k,i+1,buf_j  ,move)         &
        )
#endif
      end do
    end do

  end subroutine reduce_dpv_lon_n

  subroutine reduce_dpv_lat_n(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i, k

    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
#ifdef V_POLE
        reduced_state%dpv_lat_n(k,i,buf_j,move) = 0.25_r8 * ( &
          reduced_state%dpv_lon_t(k,i-1,buf_j-1,move) +       &
          reduced_state%dpv_lon_t(k,i  ,buf_j-1,move) +       &
          reduced_state%dpv_lon_t(k,i-1,buf_j  ,move) +       &
          reduced_state%dpv_lon_t(k,i  ,buf_j  ,move)         &
        )
#else
        reduced_state%dpv_lat_n(k,i,buf_j,move) = 0.25_r8 * ( &
          reduced_state%dpv_lon_t(k,i-1,buf_j  ,move) +       &
          reduced_state%dpv_lon_t(k,i  ,buf_j  ,move) +       &
          reduced_state%dpv_lon_t(k,i-1,buf_j+1,move) +       &
          reduced_state%dpv_lon_t(k,i  ,buf_j+1,move)         &
        )
#endif
      end do
    end do

  end subroutine reduce_dpv_lat_n

  subroutine reduce_pv_lon_apvm(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    real(r8) u, v, le, de
    integer i, k

    le = reduced_mesh%le_lon(buf_j)
    de = reduced_mesh%de_lon(buf_j)
    if (le == 0 .or. de == 0) return
    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        u = reduced_state%u(k,i,buf_j,move)
        v = reduced_state%mf_lon_t(k,i,buf_j,move) / reduced_state%m_lon(k,i,buf_j,move)
#ifdef V_POLE
        reduced_state%pv_lon(k,i,buf_j,move) = 0.5_r8 * (    &
          reduced_state%pv(k,i,buf_j+1,move) +               &
          reduced_state%pv(k,i,buf_j  ,move)                 &
        ) - 0.5_r8 * (                                       &
          u * reduced_state%dpv_lon_n(k,i,buf_j,move) / de + &
          v * reduced_state%dpv_lon_t(k,i,buf_j,move) / le   &
        ) * dt
#else
        reduced_state%pv_lon(k,i,buf_j,move) = 0.5_r8 * (    &
          reduced_state%pv(k,i,buf_j-1,move) +               &
          reduced_state%pv(k,i,buf_j  ,move)                 &
        ) - 0.5_r8 * (                                       &
          u * reduced_state%dpv_lon_n(k,i,buf_j,move) / de + &
          v * reduced_state%dpv_lon_t(k,i,buf_j,move) / le   &
        ) * dt
#endif
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%pv_lon(:,:,buf_j,move), &
                         reduced_state%async(async_pv_lon,buf_j,move), east_halo=.false.)

  end subroutine reduce_pv_lon_apvm

  subroutine reduce_pv_lat_apvm(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    real(r8) u, v, le, de
    integer i, k

    le = reduced_mesh%le_lat(buf_j)
    de = reduced_mesh%de_lat(buf_j)
    if (le == 0 .or. de == 0) return
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        u = reduced_state%mf_lat_t(k,i,buf_j,move) / reduced_state%m_lat(k,i,buf_j,move)
        v = reduced_state%v(k,i,buf_j,move)
        reduced_state%pv_lat(k,i,buf_j,move) = 0.5_r8 * (    &
          reduced_state%pv(k,i  ,buf_j,move) +               &
          reduced_state%pv(k,i-1,buf_j,move)                 &
        ) - 0.5_r8 * (                                       &
          u * reduced_state%dpv_lat_t(k,i,buf_j,move) / le + &
          v * reduced_state%dpv_lat_n(k,i,buf_j,move) / de   &
        ) * dt
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%pv_lat(:,:,buf_j,move), &
                         reduced_state%async(async_pv_lat,buf_j,move), west_halo=.false.)

  end subroutine reduce_pv_lat_apvm

  subroutine reduce_ke(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i, k

    raw_i = raw_mesh%full_lon_ibeg + move - 1
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%ke(k,i,buf_j,move) = sum(raw_state%ke(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k))
      end do
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%ke(:,:,buf_j,move) = reduced_state%ke(:,:,buf_j,move) / reduced_mesh%reduce_factor
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%ke(:,:,buf_j,move), &
                         reduced_state%async(async_ke,buf_j,move), west_halo=.false.)

  end subroutine reduce_ke

  subroutine reduce_append_array(move, reduced_mesh, reduced_array, raw_mesh, raw_array)

    integer, intent(in) :: move
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    real(r8), intent(in) :: reduced_array(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub)
    type(mesh_type), intent(in) :: raw_mesh
    real(r8), intent(inout) :: raw_array(raw_mesh%full_lon_lb:raw_mesh%full_lon_ub)

    integer raw_i, i

    raw_i = raw_mesh%full_lon_ibeg + move - 1
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1) = raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1) + reduced_array(i) / reduced_mesh%reduce_factor
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do

  end subroutine reduce_append_array

  subroutine allocate_reduced_static(reduced_mesh, reduced_static)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    allocate(reduced_static%gzs(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, -2:2,reduced_mesh%reduce_factor))

  end subroutine allocate_reduced_static

  subroutine allocate_reduced_tend(reduced_mesh, reduced_tend)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_tend_type), intent(inout) :: reduced_tend

    allocate(reduced_tend%qhu    (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub))
    allocate(reduced_tend%qhv    (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub))
    allocate(reduced_tend%dmfdlon(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub))
    allocate(reduced_tend%dpedlon(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub))
    allocate(reduced_tend%dkedlon(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub))

  end subroutine allocate_reduced_tend

  subroutine reduce_static(j, block, raw_mesh, raw_static, reduced_mesh, reduced_static)

    integer, intent(in) :: j
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(static_type), intent(in) :: raw_static
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    integer buf_j, move

    do move = 1, reduced_mesh%reduce_factor
      do buf_j = lbound(reduced_static%gzs, 2), ubound(reduced_static%gzs, 2)
        call reduce_gzs(j, buf_j, move, block, raw_mesh, raw_static, reduced_mesh, reduced_static)
      end do
    end do

  end subroutine reduce_static

end module reduce_mod
