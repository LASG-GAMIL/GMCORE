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
  use formula_mod
  use mesh_mod
  use static_mod
  use state_mod
  use parallel_mod
  use block_mod
  use reduced_types_mod
  use interp_mod, only: upwind_pole_wgt
  use upwind_mod
  use weno_mod

  implicit none

  private

  public reduce_init
  public reduce_run
  public reduce_append_array
  public reduced_mesh_type
  public reduced_tend_type
  public reduced_state_type

  interface
    subroutine reduce_sub_interface(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)
      import block_type, mesh_type, state_type, reduced_mesh_type, reduced_static_type, reduced_state_type
      integer, intent(in) :: j
      integer, intent(in) :: buf_j
      integer, intent(in) :: move
      type(block_type), intent(in) :: block
      type(mesh_type), intent(in) :: raw_mesh
      type(state_type), intent(inout) :: raw_state
      type(reduced_mesh_type), intent(in) :: reduced_mesh
      type(reduced_static_type), intent(in) :: reduced_static
      type(reduced_state_type), intent(inout) :: reduced_state
      real(8), intent(in) :: dt
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
          call log_error('Zonal reduce factor ' // to_str(reduce_factors(j)) // &
            ' cannot divide zonal grid number ' // to_str(global_mesh%num_full_lon) // '!')
        end if
        ! South Pole
#ifdef V_POLE
        full_j = j
#else
        full_j = j + 1
#endif
        if (full_j >= blocks(iblk)%mesh%full_lat_lb .and. full_j <= blocks(iblk)%mesh%full_lat_ub) then
          call reduce_mesh(reduce_factors(j), full_j, blocks(iblk)%mesh, blocks(iblk)%reduced_mesh(full_j))
        end if
        ! North Pole
#ifdef V_POLE
        full_j = global_mesh%full_lat_iend - j + 1
#else
        full_j = global_mesh%full_lat_iend - j
#endif
        if (full_j >= blocks(iblk)%mesh%full_lat_lb .and. full_j <= blocks(iblk)%mesh%full_lat_ub) then
          call reduce_mesh(reduce_factors(j), full_j, blocks(iblk)%mesh, blocks(iblk)%reduced_mesh(full_j))
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
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    integer j

    select case (pass)
    case (nh_pass_1)
      do j = block%mesh%full_lat_ibeg, block%mesh%full_lat_iend
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          call reduce_nh_state_1(j, block, block%mesh, state, block%reduced_mesh(j), block%reduced_static(j), block%reduced_state(j), dt)
        end if
      end do
    case (nh_pass_2)
      do j = block%mesh%full_lat_ibeg, block%mesh%full_lat_iend
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          call reduce_nh_state_2(j, block, block%mesh, state, block%reduced_mesh(j), block%reduced_static(j), block%reduced_state(j), dt)
        end if
      end do
    case default
      ! Extend loop range by 1 is for Coriolis forces. FIXME: Revise it.
      do j = block%mesh%full_lat_ibeg - 1, block%mesh%full_lat_iend + 1
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          call reduce_state(j, block, block%mesh, state, block%reduced_mesh(j), block%reduced_static(j), block%reduced_state(j), dt, pass)
        end if
      end do
    end select

  end subroutine reduce_run

  subroutine reduce_mesh(reduce_factor, j, raw_mesh, reduced_mesh)

    integer, intent(in) :: reduce_factor
    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(reduced_mesh_type), intent(inout) :: reduced_mesh

    real(r8) r0
    integer i, buf_j

    ! Check if decomposition is OK for reduce.
    if (mod(raw_mesh%num_full_lon, reduce_factor) /= 0) then
      if (is_root_proc()) call log_error('Parallel zonal decomposition cannot divide reduced zonal grids!')
    end if

    reduced_mesh%reduce_factor = reduce_factor

    allocate(reduced_mesh%weights(reduce_factor))
    reduced_mesh%weights = 1.0_r8 / reduce_factor

    reduced_mesh%halo_width    = merge(2, 1, pv_scheme == 2 .or. pv_scheme == 3)
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

    integer is, ie, ks, ke, i1, i2, i3

#define full_lev_full_lon_dims(x, buf_bound) reduced_state%x(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,buf_bound,reduced_mesh%reduce_factor)
#define full_lev_half_lon_dims(x, buf_bound) reduced_state%x(reduced_mesh%full_lev_lb:reduced_mesh%full_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,buf_bound,reduced_mesh%reduce_factor)
#define half_lev_full_lon_dims(x, buf_bound) reduced_state%x(reduced_mesh%half_lev_lb:reduced_mesh%half_lev_ub,reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,buf_bound,reduced_mesh%reduce_factor)
#define half_lev_half_lon_dims(x, buf_bound) reduced_state%x(reduced_mesh%half_lev_lb:reduced_mesh%half_lev_ub,reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,buf_bound,reduced_mesh%reduce_factor)

#ifdef V_POLE
    call log_error('Prepare to remove V_POLE codes.')
#else
    if (baroclinic) then
      allocate(full_lev_half_lon_dims(mf_lon_n      , -2:2))
      allocate(full_lev_full_lon_dims(mf_lat_n      , -2:1))
      allocate(full_lev_full_lon_dims(m             , -2:2))
      allocate(full_lev_half_lon_dims(m_lon         , -2:2))
      allocate(full_lev_full_lon_dims(m_lat         , -2:1))
      allocate(full_lev_full_lon_dims(ke            ,  0:0))
      allocate(full_lev_full_lon_dims(gz            ,  0:0))
      allocate(full_lev_full_lon_dims(pt            ,  0:0))
      allocate(full_lev_half_lon_dims(pt_lon        ,  0:0))
      allocate(full_lev_half_lon_dims(ptf_lon       ,  0:0))
      allocate(full_lev_full_lon_dims(ph            ,  0:0))
      allocate(half_lev_full_lon_dims(ph_lev        ,  0:0))
      allocate(full_lev_full_lon_dims(t             ,  0:0))
      if (pgf_scheme == 'sb81') then
        allocate(full_lev_full_lon_dims(ak          ,  0:0))
        allocate(full_lev_half_lon_dims(t_lnpop_lon ,  0:0))
        allocate(full_lev_half_lon_dims(ak_t_lon    ,  0:0))
      end if
      if (pgf_scheme == 'lin97' .or. nonhydrostatic) then
        allocate(half_lev_full_lon_dims(gz_lev      ,  0:0))
      end if
      if (nonhydrostatic) then
        allocate(half_lev_full_lon_dims(m_lev       ,  0:0))
        allocate(half_lev_half_lon_dims(mf_lev_lon_n,  0:0))
        allocate(half_lev_full_lon_dims(w_lev       ,  0:0))
        allocate(half_lev_half_lon_dims(w_lev_lon   ,  0:0))
        allocate(half_lev_half_lon_dims(p_lev       ,  0:0))
        allocate(half_lev_half_lon_dims(p_lev_lon   ,  0:0))
        allocate(half_lev_half_lon_dims(gz_lev_lon  ,  0:0))
        allocate(full_lev_half_lon_dims(rhod_lon    ,  0:0))
      end if
    else ! barotropic
      allocate(full_lev_half_lon_dims(mf_lon_n      , -2:2))
      allocate(full_lev_full_lon_dims(mf_lat_n      , -2:1))
      allocate(full_lev_full_lon_dims(m             , -2:2))
      allocate(full_lev_half_lon_dims(m_lon         , -2:2))
      allocate(full_lev_full_lon_dims(m_lat         , -2:1))
      allocate(full_lev_full_lon_dims(ke            ,  0:0))
      allocate(full_lev_full_lon_dims(gz            , -2:2))
    end if

    ! Potential vorticity
    allocate(full_lev_half_lon_dims(pv_lon          ,  0:0))
    allocate(full_lev_full_lon_dims(pv_lat          , -1:0))
    if (.not. reduce_pv_directly .or. pv_scheme == 4) then
      allocate(full_lev_half_lon_dims(u             , -2:2))
      allocate(full_lev_full_lon_dims(v             , -2:1))
    end if
    select case (pv_scheme)
    case (1)
      allocate(full_lev_half_lon_dims(pv            , -1:0))
    case (2, 3)
      allocate(full_lev_half_lon_dims(pv            , -2:1))
      allocate(full_lev_half_lon_dims(mf_lon_t      ,  0:0))
      allocate(full_lev_full_lon_dims(mf_lat_t      , -1:0))
    case (4)
      allocate(full_lev_half_lon_dims(pv            , -2:1))
      allocate(full_lev_half_lon_dims(mf_lon_t      ,  0:0))
      allocate(full_lev_full_lon_dims(mf_lat_t      , -1:0))
      allocate(full_lev_half_lon_dims(dpv_lon_n     ,  0:0))
      allocate(full_lev_half_lon_dims(dpv_lat_t     , -1:0))
      allocate(full_lev_full_lon_dims(dpv_lon_t     , -1:1))
      allocate(full_lev_full_lon_dims(dpv_lat_n     , -1:0))
    end select
#endif
    allocate(reduced_state%async(11,-2:2,reduced_mesh%reduce_factor))

    ! Initialize async objects.
    do i1 = 1, size(reduced_state%async, 1)
      do i2 = lbound(reduced_state%async, 2), ubound(reduced_state%async, 2)
        do i3 = 1, size(reduced_state%async, 3)
          call reduced_state%async(i1,i2,i3)%init(size(proc%ngb))
        end do
      end do
    end do

  end subroutine allocate_reduced_state

  subroutine reduce_state(j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt, pass)

    integer, intent(in) :: j
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

#define reduce_args(x, sub) lbound(reduced_state%x, 3), ubound(reduced_state%x, 3), j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, sub, dt

    if (.not. reduce_pv_directly .or. pv_scheme == 4 .or. pgf_scheme == 'sb81') then
      call apply_reduce(reduce_args(m, reduce_m))
    end if
    if (pass /= no_wind_pass) then
      ! Normal mass flux
      call apply_reduce(reduce_args(mf_lon_n, reduce_mf_lon_n))
      call apply_reduce(reduce_args(mf_lat_n, reduce_mf_lat_n))
      ! Kinetic energy
      if (do_reduce_ke) then
        call apply_reduce(reduce_args(ke    , reduce_ke      ))
      end if
      ! Potential vorticity
      if (pass == all_pass .or. pass == slow_pass) then
        if (.not. reduce_pv_directly .or. pv_scheme == 4) then
          call apply_reduce(reduce_args(m_lon      , reduce_m_lon      ))
          call apply_reduce(reduce_args(m_lat      , reduce_m_lat      ))
          call apply_reduce(reduce_args(u          , reduce_u          ))
          call apply_reduce(reduce_args(v          , reduce_v          ))
        end if
        call apply_reduce(reduce_args(pv           , reduce_pv         ))
        select case (pv_scheme)
        case (1)
          call apply_reduce(reduce_args(pv_lon     , reduce_pv_lon_midpoint))
          call apply_reduce(reduce_args(pv_lat     , reduce_pv_lat_midpoint))
        case (2)
          call apply_reduce(reduce_args(mf_lon_t   , reduce_mf_lon_t     ))
          call apply_reduce(reduce_args(mf_lat_t   , reduce_mf_lat_t     ))
          call apply_reduce(reduce_args(pv_lon     , reduce_pv_lon_upwind))
          call apply_reduce(reduce_args(pv_lat     , reduce_pv_lat_upwind))
        case (3)
          call apply_reduce(reduce_args(mf_lon_t   , reduce_mf_lon_t   ))
          call apply_reduce(reduce_args(mf_lat_t   , reduce_mf_lat_t   ))
          call apply_reduce(reduce_args(pv_lon     , reduce_pv_lon_weno))
          call apply_reduce(reduce_args(pv_lat     , reduce_pv_lat_weno))
        case (4)
          call apply_reduce(reduce_args(mf_lon_t   , reduce_mf_lon_t   ))
          call apply_reduce(reduce_args(mf_lat_t   , reduce_mf_lat_t   ))
          call apply_reduce(reduce_args(dpv_lon_t  , reduce_dpv_lon_t  ))
          call apply_reduce(reduce_args(dpv_lat_n  , reduce_dpv_lat_n  ))
          call apply_reduce(reduce_args(dpv_lat_t  , reduce_dpv_lat_t  ))
          call apply_reduce(reduce_args(dpv_lon_n  , reduce_dpv_lon_n  ))
          call apply_reduce(reduce_args(pv_lon     , reduce_pv_lon_apvm))
          call apply_reduce(reduce_args(pv_lat     , reduce_pv_lat_apvm))
        end select
      end if
    end if
    ! Potential temperature flux
    if (baroclinic) then
      call apply_reduce(reduce_args(pt             , reduce_pt          ))
      call apply_reduce(reduce_args(ptf_lon        , reduce_ptf_lon     ))
    end if
    ! Horizontal pressure gradient force
    if (hydrostatic) then ! Nonhydrostatic PGF terms is in reduce_nh_state_2.
      call apply_reduce(reduce_args(t              , reduce_t           ))
      call apply_reduce(reduce_args(gz             , reduce_gz          ))
      select case (pgf_scheme)
      case ('sb81')
        call apply_reduce(reduce_args(t_lnpop_lon  , reduce_t_lnpop_lon ))
        call apply_reduce(reduce_args(ak_t_lon     , reduce_ak_t_lon    ))
      case ('lin97')
        call apply_reduce(reduce_args(ph_lev       , reduce_ph_lev      ))
        call apply_reduce(reduce_args(gz_lev       , reduce_gz_lev      ))
      end select
    else if (.not. baroclinic) then
      call apply_reduce(reduce_args(gz             , reduce_gz          ))
    end if

  end subroutine reduce_state

  subroutine reduce_nh_state_1(j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    ! Advection of w and gz
    call apply_reduce(reduce_args(m_lev        , reduce_m_lev       ))
    call apply_reduce(reduce_args(mf_lev_lon_n , reduce_mf_lev_lon_n))
    call apply_reduce(reduce_args(gz_lev       , reduce_gz_lev      )) ! TODO: This should be called at the first time.
    call apply_reduce(reduce_args(gz_lev_lon   , reduce_gz_lev_lon  ))
    call apply_reduce(reduce_args(w_lev        , reduce_w_lev       ))
    call apply_reduce(reduce_args(w_lev_lon    , reduce_w_lev_lon   ))

  end subroutine reduce_nh_state_1

  subroutine reduce_nh_state_2(j, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    ! Horizontal pressure gradient force
    call apply_reduce(reduce_args(m_lon        , reduce_m_lon       ))
    call apply_reduce(reduce_args(ph_lev       , reduce_ph_lev      ))
    call apply_reduce(reduce_args(gz_lev       , reduce_gz_lev      ))
    call apply_reduce(reduce_args(p_lev        , reduce_p_lev       ))
    call apply_reduce(reduce_args(p_lev_lon    , reduce_p_lev_lon   ))
    call apply_reduce(reduce_args(rhod_lon     , reduce_rhod_lon    ))

  end subroutine reduce_nh_state_2

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
        reduced_static%gzs(i,buf_j,move) = sum(raw_static%gzs(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j) * reduced_mesh%weights)
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
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
    real(8), intent(in) :: dt

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
    real(8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%area_lon(buf_j) == 0) return
    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%u(k,i,buf_j,move) = reduced_state%mf_lon_n(k,i,buf_j,move) / reduced_state%m_lon(k,i,buf_j,move)
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%u(:,:,buf_j,move))

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
    real(8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%area_lat(buf_j) == 0) return
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%v(k,i,buf_j,move) = reduced_state%mf_lat_n(k,i,buf_j,move) / reduced_state%m_lat(k,i,buf_j,move)
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%v(:,:,buf_j,move))

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
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%gz(k,i,buf_j,move) = sum(raw_state%gz(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%gz(:,:,buf_j,move), west_halo=.false.)
    end if

  end subroutine reduce_gz

  subroutine reduce_gz_lev(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
          reduced_state%gz_lev(k,i,buf_j,move) = sum(raw_state%gz_lev(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%gz_lev(:,:,buf_j,move), west_halo=.false.)
    end if

  end subroutine reduce_gz_lev

  subroutine reduce_gz_lev_lon(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%half_lon_ibeg + move - 1
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
          reduced_state%gz_lev_lon(k,i,buf_j,move) = sum(raw_state%gz_lev_lon(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%gz_lev_lon(:,:,buf_j,move))
    end if

  end subroutine reduce_gz_lev_lon

  subroutine reduce_w_lev(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
          reduced_state%w_lev(k,i,buf_j,move) = sum(raw_state%w_lev(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
    end if

  end subroutine reduce_w_lev

  subroutine reduce_w_lev_lon(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%half_lon_ibeg + move - 1
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
          reduced_state%w_lev_lon(k,i,buf_j,move) = sum(raw_state%w_lev_lon(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%w_lev_lon(:,:,buf_j,move), east_halo=.false.)
    end if

  end subroutine reduce_w_lev_lon

  subroutine reduce_p_lev(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
          reduced_state%p_lev(k,i,buf_j,move) = sum(raw_state%p_lev(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%p_lev(:,:,buf_j,move), west_halo=.false.)
    end if

  end subroutine reduce_p_lev

  subroutine reduce_p_lev_lon(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%half_lon_ibeg + move - 1
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
          reduced_state%p_lev_lon(k,i,buf_j,move) = sum(raw_state%p_lev_lon(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
    end if

  end subroutine reduce_p_lev_lon

  subroutine reduce_rhod_lon(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%half_lon_ibeg + move - 1
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%rhod_lon(k,i,buf_j,move) = sum(raw_state%rhod_lon(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
    end if

  end subroutine reduce_rhod_lon

  subroutine reduce_m(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%m(k,i,buf_j,move) = sum(raw_state%m(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%m(:,:,buf_j,move), west_halo=.false.)
    end if

  end subroutine reduce_m

  subroutine reduce_m_lev(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
          reduced_state%m_lev(k,i,buf_j,move) = sum(raw_state%m_lev(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
    end if

  end subroutine reduce_m_lev

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
    real(8), intent(in) :: dt

    real(r8) m_vtx
    integer raw_i, i, k

    if (reduce_pv_directly) then
      ! Test directly reduce PV to save computation.
      if (.not. raw_mesh%is_outside_pole_half_lat(j+buf_j)) then
        raw_i = raw_mesh%half_lon_ibeg + move - 1
        do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
          do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
            reduced_state%pv(k,i,buf_j,move) = sum(raw_state%pv(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
          end do
          raw_i = raw_i + reduced_mesh%reduce_factor
        end do
      end if
    else
#ifdef V_POLE
#else
      if (raw_mesh%is_outside_pole_half_lat(j+buf_j)) then
        return
      else if (raw_mesh%is_south_pole(j+buf_j) .or. raw_mesh%is_north_pole(j+buf_j+1)) then
        raw_i = raw_mesh%half_lon_ibeg + move - 1
        do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
          do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
            reduced_state%pv (k,i,buf_j,move) = sum(raw_state%pv (raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
          end do
          raw_i = raw_i + reduced_mesh%reduce_factor
        end do
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
    end if
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
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (reduced_mesh%area_lon(buf_j) == 0) return
    raw_i = raw_mesh%half_lon_ibeg + move - 1
    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%m_lon(k,i,buf_j,move) = sum(raw_state%m_lon(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
      end do
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%m_lon(:,:,buf_j,move), east_halo=.false.)

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
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (reduced_mesh%area_lat(buf_j) == 0) return
    raw_i = raw_mesh%full_lon_ibeg + move - 1
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%m_lat(k,i,buf_j,move) = sum(raw_state%m_lat(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
      end do
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%m_lat(:,:,buf_j,move), west_halo=.false.)

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
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%half_lon_ibeg + move - 1
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%mf_lon_n(k,i,buf_j,move) = sum(raw_state%mf_lon_n(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%mf_lon_n(:,:,buf_j,move), east_halo=.false.)
    end if

  end subroutine reduce_mf_lon_n

  subroutine reduce_mf_lev_lon_n(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_full_lat(j+buf_j)) then
      raw_i = raw_mesh%half_lon_ibeg + move - 1
      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
          reduced_state%mf_lev_lon_n(k,i,buf_j,move) = sum(raw_state%mf_lev_lon_n(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%mf_lev_lon_n(:,:,buf_j,move), east_halo=.false.)
    end if

  end subroutine reduce_mf_lev_lon_n

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
    real(8), intent(in) :: dt

    integer raw_i, i, k

    if (raw_mesh%is_inside_with_halo_half_lat(j+buf_j)) then
      raw_i = raw_mesh%full_lon_ibeg + move - 1
      do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          reduced_state%mf_lat_n(k,i,buf_j,move) = sum(raw_state%mf_lat_n(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
        end do
        raw_i = raw_i + reduced_mesh%reduce_factor
      end do
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
    real(8), intent(in) :: dt

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
    real(8), intent(in) :: dt

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
    real(8), intent(in) :: dt

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
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%dpv_lon_t(:,:,buf_j,move), east_halo=.false.)

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
    real(8), intent(in) :: dt

    integer i, k

    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%dpv_lat_t(k,i,buf_j,move) = reduced_state%pv(k,i,buf_j,move) - reduced_state%pv(k,i-1,buf_j,move)
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%dpv_lat_t(:,:,buf_j,move), west_halo=.false.)

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
    real(8), intent(in) :: dt

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
    real(8), intent(in) :: dt

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

  subroutine reduce_pv_lon_midpoint(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%le_lon(buf_j) == 0 .or. reduced_mesh%de_lon(buf_j) == 0) return
    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
#ifdef V_POLE
        reduced_state%pv_lon(k,i,buf_j,move) = 0.5_r8 * (    &
          reduced_state%pv(k,i,buf_j+1,move) +               &
          reduced_state%pv(k,i,buf_j  ,move)                 &
        )
#else
        reduced_state%pv_lon(k,i,buf_j,move) = 0.5_r8 * (    &
          reduced_state%pv(k,i,buf_j-1,move) +               &
          reduced_state%pv(k,i,buf_j  ,move)                 &
        )
#endif
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%pv_lon(:,:,buf_j,move), east_halo=.false.)

  end subroutine reduce_pv_lon_midpoint

  subroutine reduce_pv_lon_upwind(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%le_lon(buf_j) == 0 .or. reduced_mesh%de_lon(buf_j) == 0) return
    associate (v => reduced_state%mf_lon_t, pv => reduced_state%pv, pv_lon => reduced_state%pv_lon)
      select case (upwind_order_pv)
      case (3)
        if (raw_mesh%is_full_lat_next_to_pole(j+buf_j)) then
          do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
            do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
              pv_lon(k,i,buf_j,move) = upwind1(sign(1.0_r8, v(k,i,buf_j,move)), upwind_wgt_pv, pv(k,i,buf_j-1:buf_j,move))
            end do
          end do
        else
          do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
            do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
              pv_lon(k,i,buf_j,move) = upwind3(sign(1.0_r8, v(k,i,buf_j,move)), upwind_wgt_pv, pv(k,i,buf_j-2:buf_j+1,move))
            end do
          end do
        end if
      end select
      call fill_zonal_halo(block, reduced_mesh%halo_width, pv_lon(:,:,buf_j,move), east_halo=.false.)
    end associate

  end subroutine reduce_pv_lon_upwind

  subroutine reduce_pv_lon_weno(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%le_lon(buf_j) == 0 .or. reduced_mesh%de_lon(buf_j) == 0) return
    associate (v => reduced_state%mf_lon_t, pv => reduced_state%pv, pv_lon => reduced_state%pv_lon)
      select case (weno_order_pv)
      case (3)
        if (raw_mesh%is_full_lat_next_to_pole(j+buf_j)) then
          do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
            do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
              pv_lon(k,i,buf_j,move) = upwind1(sign(1.0_r8, v(k,i,buf_j,move)), upwind_wgt_pv, pv(k,i,buf_j-1:buf_j,move))
            end do
          end do
        else
          do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
            do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
              pv_lon(k,i,buf_j,move) = weno3(sign(1.0_r8, v(k,i,buf_j,move)), pv(k,i,buf_j-2:buf_j+1,move))
            end do
          end do
        end if
      end select
      call fill_zonal_halo(block, reduced_mesh%halo_width, pv_lon(:,:,buf_j,move), east_halo=.false.)
    end associate

  end subroutine reduce_pv_lon_weno

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
    real(8), intent(in) :: dt

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
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%pv_lon(:,:,buf_j,move), east_halo=.false.)

  end subroutine reduce_pv_lon_apvm

  subroutine reduce_pv_lat_midpoint(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%le_lat(buf_j) == 0 .or. reduced_mesh%de_lat(buf_j) == 0) return
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%pv_lat(k,i,buf_j,move) = 0.5_r8 * (    &
          reduced_state%pv(k,i  ,buf_j,move) +               &
          reduced_state%pv(k,i-1,buf_j,move)                 &
        )
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%pv_lat(:,:,buf_j,move), west_halo=.false.)

  end subroutine reduce_pv_lat_midpoint

  subroutine reduce_pv_lat_upwind(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%le_lat(buf_j) == 0 .or. reduced_mesh%de_lat(buf_j) == 0) return
    associate (u => reduced_state%mf_lat_t, pv => reduced_state%pv, pv_lat => reduced_state%pv_lat)
      select case (upwind_order_pv)
      case (3)
        do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
          do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
            pv_lat(k,i,buf_j,move) = upwind3(sign(1.0_r8, u(k,i,buf_j,move)), upwind_wgt_pv, pv(k,i-2:i+1,buf_j,move))
          end do
        end do
      end select
      call fill_zonal_halo(block, reduced_mesh%halo_width, pv_lat(:,:,buf_j,move), west_halo=.false.)
    end associate

  end subroutine reduce_pv_lat_upwind

  subroutine reduce_pv_lat_weno(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer i, k

    if (reduced_mesh%le_lat(buf_j) == 0 .or. reduced_mesh%de_lat(buf_j) == 0) return
    associate (u => reduced_state%mf_lat_t, pv => reduced_state%pv, pv_lat => reduced_state%pv_lat)
      select case (weno_order_pv)
      case (3)
        do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
          do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
            pv_lat(k,i,buf_j,move) = weno3(sign(1.0_r8, u(k,i,buf_j,move)), pv(k,i-2:i+1,buf_j,move))
          end do
        end do
      end select
      call fill_zonal_halo(block, reduced_mesh%halo_width, pv_lat(:,:,buf_j,move), west_halo=.false.)
    end associate

  end subroutine reduce_pv_lat_weno

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
    real(8), intent(in) :: dt

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
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%pv_lat(:,:,buf_j,move), west_halo=.false.)

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
    real(8), intent(in) :: dt

    integer raw_i, i, k

    raw_i = raw_mesh%full_lon_ibeg + move - 1
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%ke(k,i,buf_j,move) = sum(raw_state%ke(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
      end do
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%ke(:,:,buf_j,move), west_halo=.false.)

  end subroutine reduce_ke

  subroutine reduce_pt(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    raw_i = raw_mesh%full_lon_ibeg + move - 1
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%pt(k,i,buf_j,move) = sum(raw_state%pt(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * reduced_mesh%weights)
      end do
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%pt(:,:,buf_j,move))

  end subroutine reduce_pt

  ! FIXME: This variable is not necessary!
  subroutine reduce_ptf_lon(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    real(r8) beta
    integer i, k

    if (reduced_mesh%area_lon(buf_j) == 0) return
    associate (pt       => reduced_state%pt      , &
               pt_lon   => reduced_state%pt_lon  , &
               mf_lon_n => reduced_state%mf_lon_n, &
               ptf_lon  => reduced_state%ptf_lon)
      if (upwind_order == -1 .and. weno_order == -1) then
        do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
          do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
            pt_lon(k,i,buf_j,move) = 0.5_r8 * (pt(k,i,buf_j,move) + pt(k,i+1,buf_j,move))
          end do
        end do
      else
        select case (upwind_order)
        case (3)
          beta = upwind_wgt_pt * upwind_pole_wgt(j+buf_j)
          do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
            do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
              pt_lon(k,i,buf_j,move) = upwind3(sign(1.0_r8, mf_lon_n(k,i,buf_j,move)), beta, pt(k,i-1:i+2,buf_j,move))
            end do
          end do
        end select
        select case (weno_order)
        case (3)
          do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
            do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
              pt_lon(k,i,buf_j,move) = weno3(sign(1.0_r8, mf_lon_n(k,i,buf_j,move)), pt(k,i-1:i+2,buf_j,move))
            end do
          end do
        end select
      end if

      do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
        do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
          ptf_lon(k,i,buf_j,move) = mf_lon_n(k,i,buf_j,move) * pt_lon(k,i,buf_j,move)
        end do
      end do
      call fill_zonal_halo(block, reduced_mesh%halo_width, ptf_lon(:,:,buf_j,move), east_halo=.false.)
    end associate

  end subroutine reduce_ptf_lon

  subroutine reduce_t(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer i, k

    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%t(k,i,buf_j,move) = temperature(reduced_state%pt(k,i,buf_j,move), reduced_state%ph(k,i,buf_j,move))
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%t(:,:,buf_j,move), west_halo=.false.)

  end subroutine reduce_t

  subroutine reduce_t_lnpop_lon(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer i, k

    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%t_lnpop_lon(k,i,buf_j,move) = (                                         &
          reduced_mesh%area_lon_west(buf_j) * reduced_state%t(k,i  ,buf_j,move) * log(        &
            reduced_state%ph_lev(k+1,i  ,buf_j,move) / reduced_state%ph_lev(k,i  ,buf_j,move) &
          ) +                                                                                 &
          reduced_mesh%area_lon_east(buf_j) * reduced_state%t(k,i+1,buf_j,move) * log(        &
            reduced_state%ph_lev(k+1,i+1,buf_j,move) / reduced_state%ph_lev(k,i+1,buf_j,move) &
          )                                                                                   &
        ) / reduced_mesh%area_lon(buf_j)
      end do
    end do

  end subroutine reduce_t_lnpop_lon

  subroutine reduce_ak_t_lon(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k
    real(r8), parameter :: ln2 = log(2.0_r8)

    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        if (k == 1) then
          reduced_state%ak(k,i,buf_j,move) = ln2
        else
          reduced_state%ak(k,i,buf_j,move) = 1.0_r8 - reduced_state%ph_lev(k,i,buf_j,move) / reduced_state%m(k,i,buf_j,move) * &
            log(reduced_state%ph_lev(k+1,i,buf_j,move) / reduced_state%ph_lev(k,i,buf_j,move))
        end if
      end do
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%ak(:,:,buf_j,move), west_halo=.false.)

    do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%ak_t_lon(k,i,buf_j,move) = (                                                                     &
          reduced_mesh%area_lon_west(buf_j) * reduced_state%ak(k,i  ,buf_j,move) * reduced_state%t(k,i  ,buf_j,move) + &
          reduced_mesh%area_lon_east(buf_j) * reduced_state%ak(k,i+1,buf_j,move) * reduced_state%t(k,i+1,buf_j,move)   &
        ) / reduced_mesh%area_lon(buf_j)
      end do
    end do

  end subroutine reduce_ak_t_lon

  subroutine reduce_ph_lev(j, buf_j, move, block, raw_mesh, raw_state, reduced_mesh, reduced_static, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(block_type), intent(in) :: block
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(in) :: reduced_static
    type(reduced_state_type), intent(inout) :: reduced_state
    real(8), intent(in) :: dt

    integer raw_i, i, k

    raw_i = raw_mesh%full_lon_ibeg + move - 1
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%half_lev_ibeg, reduced_mesh%half_lev_iend
        reduced_state%ph_lev(k,i,buf_j,move) = sum(                              &
          raw_state%ph_lev(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j,k) * &
          reduced_mesh%weights                                                   &
        )
      end do
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    call fill_zonal_halo(block, reduced_mesh%halo_width, reduced_state%ph_lev(:,:,buf_j,move), west_halo=.false.)

    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      do k = reduced_mesh%full_lev_ibeg, reduced_mesh%full_lev_iend
        reduced_state%ph(k,i,buf_j,move) = 0.5_r8 * (reduced_state%ph_lev(k,i,buf_j,move) + reduced_state%ph_lev(k+1,i,buf_j,move))
      end do
    end do

  end subroutine reduce_ph_lev

  subroutine reduce_append_array(move, reduced_mesh, reduced_array, raw_mesh, raw_array)

    integer, intent(in) :: move
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    real(r8), intent(in) :: reduced_array(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub)
    type(mesh_type), intent(in) :: raw_mesh
    real(r8), intent(inout) :: raw_array(raw_mesh%full_lon_lb:raw_mesh%full_lon_ub)

    integer raw_i, i

    raw_i = raw_mesh%full_lon_ibeg + move - 1
    do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
      raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1) = raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1) + reduced_array(i) * reduced_mesh%weights
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

    integer is, ie, ks, ke

    is = reduced_mesh%full_lon_lb; ie = reduced_mesh%full_lon_ub
    ks = reduced_mesh%full_lev_lb; ke = reduced_mesh%full_lev_ub

    allocate(reduced_tend%qhu     (is:ie,ks:ke))
    allocate(reduced_tend%fu      (is:ie,ks:ke))
    allocate(reduced_tend%dmfdlon (is:ie,ks:ke))
    allocate(reduced_tend%dptfdlon(is:ie,ks:ke))

    is = reduced_mesh%half_lon_lb; ie = reduced_mesh%half_lon_ub
    ks = reduced_mesh%full_lev_lb; ke = reduced_mesh%full_lev_ub

    allocate(reduced_tend%qhv     (is:ie,ks:ke))
    allocate(reduced_tend%fv      (is:ie,ks:ke))
    allocate(reduced_tend%pgf_lon (is:ie,ks:ke))
    allocate(reduced_tend%dkedlon (is:ie,ks:ke))

    is = reduced_mesh%half_lon_lb; ie = reduced_mesh%half_lon_ub
    ks = reduced_mesh%half_lev_lb; ke = reduced_mesh%half_lev_ub

    allocate(reduced_tend%adv_gz_lon(is:ie,ks:ke))
    allocate(reduced_tend%adv_w_lon (is:ie,ks:ke))

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
