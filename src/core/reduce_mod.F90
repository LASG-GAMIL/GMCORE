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
  use sphere_geometry_mod
  use mesh_mod
  use static_mod
  use state_mod
  use parallel_mod

  implicit none

  private

  public reduced_mesh_type
  public reduced_state_type
  public reduce_init
  public reduce_run
  public reduce_append_array
  public reduce_replace_pv
  public reduce_final

  public reduced_mesh
  public reduced_static
  public reduced_state
  public reduced_tend

  type reduced_mesh_type
    integer :: reduce_factor = 0
    integer :: damp_order = 0
    integer halo_width
    integer num_full_lon
    integer num_half_lon
    integer full_lon_start_idx
    integer full_lon_end_idx
    integer half_lon_start_idx
    integer half_lon_end_idx
    integer full_lon_lb
    integer full_lon_ub
    integer half_lon_lb
    integer half_lon_ub
#ifdef V_POLE
    real(r8), dimension(  -1:1) :: full_lat            = inf
    real(r8), dimension(  -1:1) :: half_lat            = inf
    real(r8), dimension(  -1:1) :: cell_area           = 0
    real(r8), dimension(2,-2:2) :: subcell_area        = 0
    real(r8), dimension(  -2:2) :: lon_edge_area       = 0
    real(r8), dimension(  -2:2) :: lon_edge_left_area  = 0
    real(r8), dimension(  -2:2) :: lon_edge_right_area = 0
    real(r8), dimension(  -1:2) :: vertex_area         = 0
    real(r8), dimension(  -1:2) :: lat_edge_area       = 0
    real(r8), dimension(  -1:2) :: lat_edge_up_area    = 0
    real(r8), dimension(  -1:2) :: lat_edge_down_area  = 0
    real(r8), dimension(  -1:1) :: le_lon              = inf
    real(r8), dimension(  -2:2) :: de_lon              = inf
    real(r8), dimension(  -1:2) :: le_lat              = inf
    real(r8), dimension(  -1:2) :: de_lat              = inf
    real(r8), dimension(2,-1:1) :: full_tangent_wgt    = inf
    real(r8), dimension(2, 0:1) :: half_tangent_wgt    = inf
    real(r8), dimension(  -1:2) :: half_f              = inf
#else
    real(r8), dimension(  -1:1) :: full_lat            = inf
    real(r8), dimension(  -1:1) :: half_lat            = inf
    real(r8), dimension(  -1:1) :: cell_area           = 0
    real(r8), dimension(2,-2:2) :: subcell_area        = 0
    real(r8), dimension(  -2:2) :: lon_edge_area       = 0
    real(r8), dimension(  -2:2) :: lon_edge_left_area  = 0
    real(r8), dimension(  -2:2) :: lon_edge_right_area = 0
    real(r8), dimension(  -2:1) :: vertex_area         = 0
    real(r8), dimension(  -2:1) :: lat_edge_area       = 0
    real(r8), dimension(  -2:1) :: lat_edge_up_area    = 0
    real(r8), dimension(  -2:1) :: lat_edge_down_area  = 0
    real(r8), dimension(  -1:1) :: le_lon              = inf
    real(r8), dimension(  -2:2) :: de_lon              = inf
    real(r8), dimension(  -2:1) :: le_lat              = inf
    real(r8), dimension(  -2:1) :: de_lat              = inf
    real(r8), dimension(2,-1:1) :: full_tangent_wgt    = inf
    real(r8), dimension(2,-1:0) :: half_tangent_wgt    = inf
    real(r8), dimension(  -2:1) :: half_f              = inf
#endif
  end type reduced_mesh_type

  type(reduced_mesh_type), allocatable :: reduced_mesh(:)

  type reduced_static_type
    real(r8), allocatable, dimension(:,:,:) :: ghs
  contains
    final :: reduced_static_final
  end type reduced_static_type

  type(reduced_static_type), allocatable :: reduced_static(:)

  type reduced_state_type
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:) :: gd
    real(r8), allocatable, dimension(:,:,:) :: pv
    real(r8), allocatable, dimension(:,:,:) :: pv_lon
    real(r8), allocatable, dimension(:,:,:) :: pv_lat
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_t
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_t
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_n
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_n
    real(r8), allocatable, dimension(:,:,:) :: m_lon
    real(r8), allocatable, dimension(:,:,:) :: m_lat
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_n
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_t
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_n
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_t
    real(r8), allocatable, dimension(:,:,:) :: ke
  contains
    final :: reduced_state_final
  end type reduced_state_type

  type(reduced_state_type), allocatable :: reduced_state(:)

  type reduced_tend_type
    real(r8), allocatable, dimension(:) :: qhv
    real(r8), allocatable, dimension(:) :: qhu
    real(r8), allocatable, dimension(:) :: dmfdlon
    real(r8), allocatable, dimension(:) :: dpedlon
    real(r8), allocatable, dimension(:) :: dkedlon
  contains
    final :: reduced_tend_final
  end type reduced_tend_type

  type(reduced_tend_type), allocatable :: reduced_tend(:)

  interface
    subroutine reduce_sub_interface(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)
      import mesh_type, state_type, reduced_mesh_type, reduced_state_type, r8
      integer, intent(in) :: j
      integer, intent(in) :: buf_j
      integer, intent(in) :: move
      type(mesh_type), intent(in) :: raw_mesh
      type(state_type), intent(in) :: raw_state
      type(reduced_mesh_type), intent(in) :: reduced_mesh
      type(reduced_state_type), intent(inout) :: reduced_state
      real(r8), intent(in) :: dt
    end subroutine reduce_sub_interface
  end interface

contains

  subroutine reduce_init()

    integer j, full_j

    allocate(reduced_mesh  (global_mesh%full_lat_lb:global_mesh%full_lat_ub))
    allocate(reduced_static(global_mesh%full_lat_start_idx:global_mesh%full_lat_end_idx))
    allocate(reduced_state (global_mesh%full_lat_start_idx:global_mesh%full_lat_end_idx))
    allocate(reduced_tend  (global_mesh%full_lat_start_idx:global_mesh%full_lat_end_idx))

    do j = 1, size(reduce_factors)
      if (reduce_factors(j) == 0) exit
      if (mod(global_mesh%num_full_lon, reduce_factors(j)) /= 0) then
        call log_error('Zonal reduce factor ' // to_string(reduce_factors(j)) // ' cannot divide zonal grid number ' // to_string(global_mesh%num_full_lon) // '!')
      end if
      if (global_mesh%has_south_pole()) then
#ifdef V_POLE
        full_j = global_mesh%full_lat_start_idx+j-1
#else
        full_j = global_mesh%full_lat_start_idx+j
#endif
        call reduce_mesh(reduce_factors(j), full_j, global_mesh, reduced_mesh(full_j))
        reduced_mesh(full_j)%damp_order = damp_orders(j)
      end if
      if (global_mesh%has_north_pole()) then
#ifdef V_POLE
        full_j = global_mesh%full_lat_end_idx-j+1
#else
        full_j = global_mesh%full_lat_end_idx-j
#endif
        call reduce_mesh(reduce_factors(j), full_j, global_mesh, reduced_mesh(full_j))
        reduced_mesh(full_j)%damp_order = damp_orders(j)
      end if
    end do

    do j = global_mesh%full_lat_start_idx, global_mesh%full_lat_end_idx
      if (reduced_mesh(j)%reduce_factor > 0) then
        call allocate_reduced_static(reduced_mesh(j), reduced_static(j))
        call reduce_static(j, global_mesh, static, reduced_mesh(j), reduced_static(j))
        call allocate_reduced_state(reduced_mesh(j), reduced_state(j))
        call allocate_reduced_tend(reduced_mesh(j), reduced_tend(j))
      end if
    end do

  end subroutine reduce_init

  subroutine reduce_run(state, dt)

    type(state_type), intent(in) :: state
    real(r8), intent(in) :: dt

    integer j

    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      if (reduced_mesh(j)%reduce_factor > 0) then
        call reduce_state(j, state%mesh, state, reduced_mesh(j), reduced_state(j), dt)
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

    reduced_mesh%reduce_factor      = reduce_factor
    reduced_mesh%halo_width         = raw_mesh%halo_width
    reduced_mesh%num_full_lon       = raw_mesh%num_full_lon / reduce_factor
    reduced_mesh%num_half_lon       = raw_mesh%num_half_lon / reduce_factor
    reduced_mesh%full_lon_start_idx = raw_mesh%full_lon_start_idx                                 ! FIXME: This is wrong in parallel.
    reduced_mesh%full_lon_end_idx   = raw_mesh%full_lon_start_idx + reduced_mesh%num_full_lon - 1 ! FIXME: This is wrong in parallel.
    reduced_mesh%half_lon_start_idx = raw_mesh%half_lon_start_idx                                 ! FIXME: This is wrong in parallel.
    reduced_mesh%half_lon_end_idx   = raw_mesh%half_lon_start_idx + reduced_mesh%num_half_lon - 1 ! FIXME: This is wrong in parallel.
    reduced_mesh%full_lon_lb        = reduced_mesh%full_lon_start_idx - raw_mesh%halo_width
    reduced_mesh%full_lon_ub        = reduced_mesh%full_lon_end_idx + raw_mesh%halo_width
    reduced_mesh%half_lon_lb        = reduced_mesh%half_lon_start_idx - raw_mesh%halo_width
    reduced_mesh%half_lon_ub        = reduced_mesh%half_lon_end_idx + raw_mesh%halo_width

    reduced_mesh%full_lat = raw_mesh%full_lat(j+lbound(reduced_mesh%full_lat, 1):j+ubound(reduced_mesh%full_lat, 1))
    reduced_mesh%half_lat = raw_mesh%half_lat(j+lbound(reduced_mesh%half_lat, 1):j+ubound(reduced_mesh%half_lat, 1))
    reduced_mesh%half_f   = raw_mesh%half_f  (j+lbound(reduced_mesh%half_f  , 1):j+ubound(reduced_mesh%half_f  , 1))

    ! Cell area
    do buf_j = lbound(reduced_mesh%cell_area, 1), ubound(reduced_mesh%cell_area, 1)
      if (.not. is_inf(raw_mesh%cell_area(j+buf_j))) then
        reduced_mesh%cell_area(buf_j) = raw_mesh%cell_area(j+buf_j) * reduce_factor
      end if
    end do
    do buf_j = lbound(reduced_mesh%subcell_area, 2), ubound(reduced_mesh%subcell_area, 2)
      if (raw_mesh%is_outside_full_lat(j+buf_j)) cycle
      reduced_mesh%subcell_area(1,buf_j) = raw_mesh%subcell_area(1,j+buf_j) * reduce_factor
      reduced_mesh%subcell_area(2,buf_j) = raw_mesh%subcell_area(2,j+buf_j) * reduce_factor
    end do
    do buf_j = lbound(reduced_mesh%lon_edge_area, 1), ubound(reduced_mesh%lon_edge_area, 1)
      if (raw_mesh%is_outside_full_lat(j+buf_j)) cycle
      reduced_mesh%lon_edge_left_area (buf_j) = raw_mesh%lon_edge_left_area (j+buf_j) * reduce_factor
      reduced_mesh%lon_edge_right_area(buf_j) = raw_mesh%lon_edge_right_area(j+buf_j) * reduce_factor
      reduced_mesh%lon_edge_area      (buf_j) = raw_mesh%lon_edge_area      (j+buf_j) * reduce_factor
    end do
    ! Vertex area
    do buf_j = lbound(reduced_mesh%vertex_area, 1), ubound(reduced_mesh%vertex_area, 1)
      if (raw_mesh%is_outside_half_lat(j+buf_j)) cycle
      reduced_mesh%vertex_area(buf_j) = raw_mesh%vertex_area(j+buf_j) * reduce_factor
    end do
    do buf_j = lbound(reduced_mesh%lat_edge_area, 1), ubound(reduced_mesh%lat_edge_area, 1)
      if (raw_mesh%is_outside_half_lat(j+buf_j)) cycle
      reduced_mesh%lat_edge_up_area  (buf_j) = raw_mesh%lat_edge_up_area  (j+buf_j) * reduce_factor
      reduced_mesh%lat_edge_down_area(buf_j) = raw_mesh%lat_edge_down_area(j+buf_j) * reduce_factor
      reduced_mesh%lat_edge_area     (buf_j) = raw_mesh%lat_edge_area     (j+buf_j) * reduce_factor
    end do
    ! Edge lengths and cell distances
    do buf_j = lbound(reduced_mesh%le_lat, 1), ubound(reduced_mesh%le_lat, 1)
      if (raw_mesh%is_outside_half_lat(j+buf_j)) cycle
      reduced_mesh%le_lat(buf_j) = raw_mesh%le_lat(j+buf_j) * reduce_factor
    end do
    do buf_j = lbound(reduced_mesh%de_lat, 1), ubound(reduced_mesh%de_lat, 1)
      if (raw_mesh%is_outside_half_lat(j+buf_j)) cycle
      reduced_mesh%de_lat(buf_j) = raw_mesh%de_lat(j+buf_j)
    end do
    do buf_j = lbound(reduced_mesh%le_lon, 1), ubound(reduced_mesh%le_lon, 1)
      if (raw_mesh%is_outside_full_lat(j+buf_j) .or. is_inf(raw_mesh%le_lon(j+buf_j))) cycle
      reduced_mesh%le_lon(buf_j) = raw_mesh%le_lon(j+buf_j)
    end do
    do buf_j = lbound(reduced_mesh%de_lon, 1), ubound(reduced_mesh%de_lon, 1)
      if (raw_mesh%is_outside_full_lat(j+buf_j) .or. is_inf(raw_mesh%de_lon(j+buf_j))) cycle
      reduced_mesh%de_lon(buf_j) = raw_mesh%de_lon(j+buf_j) * reduce_factor
    end do

#ifdef V_POLE
    do buf_j = lbound(reduced_mesh%full_tangent_wgt, 2), ubound(reduced_mesh%full_tangent_wgt, 2)
      if (.not. is_inf(reduced_mesh%le_lat(buf_j  )) .and. .not. is_inf(reduced_mesh%de_lon(buf_j))) then
        reduced_mesh%full_tangent_wgt(1,buf_j) = reduced_mesh%le_lat(buf_j  ) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
      if (.not. is_inf(reduced_mesh%le_lat(buf_j+1)) .and. .not. is_inf(reduced_mesh%de_lon(buf_j))) then
        reduced_mesh%full_tangent_wgt(2,buf_j) = reduced_mesh%le_lat(buf_j+1) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
    end do
#else
    do buf_j = lbound(reduced_mesh%full_tangent_wgt, 2), ubound(reduced_mesh%full_tangent_wgt, 2)
      if (reduced_mesh%le_lat(buf_j-1) /= inf .and. reduced_mesh%de_lon(buf_j) /= inf) then
        reduced_mesh%full_tangent_wgt(1,buf_j) = reduced_mesh%le_lat(buf_j-1) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lat(buf_j  ) /= inf .and. reduced_mesh%de_lon(buf_j) /= inf) then
        reduced_mesh%full_tangent_wgt(2,buf_j) = reduced_mesh%le_lat(buf_j  ) / reduced_mesh%de_lon(buf_j) * 0.25_r8
      end if
    end do
#endif
#ifdef V_POLE
    do buf_j = lbound(reduced_mesh%half_tangent_wgt, 2), ubound(reduced_mesh%half_tangent_wgt, 2)
      if (raw_mesh%is_south_pole(j+buf_j+1)) cycle
      if (.not. is_inf(reduced_mesh%le_lon(buf_j-1)) .and. .not. is_inf(reduced_mesh%de_lat(buf_j))) then
        reduced_mesh%half_tangent_wgt(1,buf_j) = reduced_mesh%le_lon(buf_j-1) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
      if (.not. is_inf(reduced_mesh%le_lon(buf_j  )) .and. .not. is_inf(reduced_mesh%de_lat(buf_j))) then
        reduced_mesh%half_tangent_wgt(2,buf_j) = reduced_mesh%le_lon(buf_j  ) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
    end do
#else
    do buf_j = lbound(reduced_mesh%half_tangent_wgt, 2), ubound(reduced_mesh%half_tangent_wgt, 2)
      if (reduced_mesh%le_lon(buf_j  ) /= inf .and. reduced_mesh%de_lat(buf_j) /= inf) then
        reduced_mesh%half_tangent_wgt(1,buf_j) = reduced_mesh%le_lon(buf_j  ) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
      if (reduced_mesh%le_lon(buf_j+1) /= inf .and. reduced_mesh%de_lat(buf_j) /= inf) then
        reduced_mesh%half_tangent_wgt(2,buf_j) = reduced_mesh%le_lon(buf_j+1) / reduced_mesh%de_lat(buf_j) * 0.25_r8
      end if
    end do
#endif

  end subroutine reduce_mesh

  subroutine allocate_reduced_state(reduced_mesh, reduced_state)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state

#ifdef V_POLE
    allocate(reduced_state%mf_lon_n (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_n (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lon_t (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_t (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%gd       (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lon    (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lat    (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%u        (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%v        (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv       (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lon   (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lat   (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_t(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_n(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_t(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_n(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%ke       (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:0,reduced_mesh%reduce_factor))
#else
    allocate(reduced_state%mf_lon_n (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_n (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lon_t (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%mf_lat_t (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%gd       (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lon    (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%m_lat    (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%u        (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:2,reduced_mesh%reduce_factor))
    allocate(reduced_state%v        (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-2:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv       (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-2:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lon   (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%pv_lat   (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_t(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:1,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lon_n(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub, 0:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_t(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub,-1:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%dpv_lat_n(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub,-1:0,reduced_mesh%reduce_factor))
    allocate(reduced_state%ke       (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:0,reduced_mesh%reduce_factor))
#endif

  end subroutine allocate_reduced_state

  subroutine reduce_state(j, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    reduced_state%pv       (:,:,:) = 0
    reduced_state%pv_lon   (:,:,:) = 0
    reduced_state%pv_lat   (:,:,:) = 0
    reduced_state%dpv_lon_t(:,:,:) = inf
    reduced_state%dpv_lat_t(:,:,:) = inf
    reduced_state%dpv_lon_n(:,:,:) = inf
    reduced_state%dpv_lat_n(:,:,:) = inf
    reduced_state%m_lon    (:,:,:) = inf
    reduced_state%m_lat    (:,:,:) = inf
    reduced_state%mf_lon_n (:,:,:) = inf
    reduced_state%mf_lon_t (:,:,:) = inf
    reduced_state%mf_lat_n (:,:,:) = inf
    reduced_state%mf_lat_t (:,:,:) = inf

    call apply_reduce(lbound(reduced_state%mf_lon_n , 2), ubound(reduced_state%mf_lon_n , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_mf_lon_n   , dt)
    call apply_reduce(lbound(reduced_state%mf_lat_n , 2), ubound(reduced_state%mf_lat_n , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_mf_lat_n   , dt)
    call apply_reduce(lbound(reduced_state%mf_lon_t , 2), ubound(reduced_state%mf_lon_t , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_mf_lon_t   , dt)
    call apply_reduce(lbound(reduced_state%mf_lat_t , 2), ubound(reduced_state%mf_lat_t , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_mf_lat_t   , dt)
    call apply_reduce(lbound(reduced_state%gd       , 2), ubound(reduced_state%gd       , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_gd         , dt)
    call apply_reduce(lbound(reduced_state%m_lon    , 2), ubound(reduced_state%m_lon    , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_m_lon      , dt)
    call apply_reduce(lbound(reduced_state%m_lat    , 2), ubound(reduced_state%m_lat    , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_m_lat      , dt)
    call apply_reduce(lbound(reduced_state%u        , 2), ubound(reduced_state%u        , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_u          , dt)
    call apply_reduce(lbound(reduced_state%v        , 2), ubound(reduced_state%v        , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_v          , dt)
    call apply_reduce(lbound(reduced_state%pv       , 2), ubound(reduced_state%pv       , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_pv         , dt)
    call apply_reduce(lbound(reduced_state%dpv_lon_t, 2), ubound(reduced_state%dpv_lon_t, 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_dpv_lon_t  , dt)
    call apply_reduce(lbound(reduced_state%dpv_lat_n, 2), ubound(reduced_state%dpv_lat_n, 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_dpv_lat_n  , dt)
    call apply_reduce(lbound(reduced_state%dpv_lat_t, 2), ubound(reduced_state%dpv_lat_t, 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_dpv_lat_t  , dt)
    call apply_reduce(lbound(reduced_state%dpv_lon_n, 2), ubound(reduced_state%dpv_lon_n, 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_dpv_lon_n  , dt)
    call apply_reduce(lbound(reduced_state%pv_lon   , 2), ubound(reduced_state%pv_lon   , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_pv_lon_apvm, dt)
    call apply_reduce(lbound(reduced_state%pv_lat   , 2), ubound(reduced_state%pv_lat   , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_pv_lat_apvm, dt)
    call apply_reduce(lbound(reduced_state%ke       , 2), ubound(reduced_state%ke       , 2), j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_ke         , dt)

  end subroutine reduce_state

  subroutine reduce_ghs(j, buf_j, move, raw_mesh, raw_static, reduced_mesh, reduced_static)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(static_type), intent(in) :: raw_static
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    integer raw_i, i

    raw_i = raw_mesh%full_lon_start_idx + move - 1
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_static%ghs(i,buf_j,move) = sum(raw_static%ghs(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_static%ghs(:,buf_j,move) = reduced_static%ghs(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_static%ghs(:,buf_j,move))

  end subroutine reduce_ghs

  subroutine apply_reduce(buf_lb, buf_ub, j, raw_mesh, raw_state, reduced_mesh, reduced_state, reduce_sub, dt)

    integer, intent(in) :: buf_lb
    integer, intent(in) :: buf_ub
    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    procedure(reduce_sub_interface) reduce_sub
    real(r8), intent(in) :: dt

    integer buf_j, move

    do move = 1, reduced_mesh%reduce_factor
      do buf_j = buf_lb, buf_ub
        call reduce_sub(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)
      end do
    end do

  end subroutine apply_reduce

  subroutine reduce_u(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_state%u(i,buf_j,move) = reduced_state%mf_lon_n(i,buf_j,move) / reduced_state%m_lon(i,buf_j,move)
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%u(:,buf_j,move))

  end subroutine reduce_u

  subroutine reduce_v(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%v(i,buf_j,move) = reduced_state%mf_lat_n(i,buf_j,move) / reduced_state%m_lat(i,buf_j,move)
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%v(:,buf_j,move))

  end subroutine reduce_v

  subroutine reduce_gd(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i

    if (raw_mesh%is_outside_full_lat(j+buf_j)) return
    raw_i = raw_mesh%full_lon_start_idx + move - 1
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%gd(i,buf_j,move) = sum(raw_state%gd(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%gd(:,buf_j,move) = reduced_state%gd(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%gd(:,buf_j,move))

  end subroutine reduce_gd

  subroutine reduce_pv(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    real(r8) m_vtx, pole, sign
    integer i, u_j

#ifdef V_POLE
    if (raw_mesh%is_outside_half_lat(j+buf_j)) then
      return
    else if (raw_mesh%is_south_pole(j+buf_j)) then
      reduced_state%pv(:,buf_j,move) = raw_state%pv(raw_mesh%full_lon_start_idx,raw_mesh%half_lat_start_idx)
    else if (raw_mesh%is_north_pole(j+buf_j)) then
      reduced_state%pv(:,buf_j,move) = raw_state%pv(raw_mesh%full_lon_start_idx,raw_mesh%half_lat_end_idx)
    else
      do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
        m_vtx = (                                                                                                          &
          (reduced_state%gd(i,buf_j-1,move) + reduced_state%gd(i+1,buf_j-1,move)) * reduced_mesh%subcell_area(2,buf_j-1) + &
          (reduced_state%gd(i,buf_j  ,move) + reduced_state%gd(i+1,buf_j  ,move)) * reduced_mesh%subcell_area(1,buf_j  )   &
        ) / reduced_mesh%vertex_area(buf_j) / g
        reduced_state%pv(i,buf_j,move) = (                                     &
          (                                                                    &
            reduced_state%u(i  ,buf_j-1,move) * reduced_mesh%de_lon(buf_j-1) - &
            reduced_state%u(i  ,buf_j  ,move) * reduced_mesh%de_lon(buf_j  ) + &
            reduced_state%v(i+1,buf_j  ,move) * reduced_mesh%de_lat(buf_j  ) - &
            reduced_state%v(i  ,buf_j  ,move) * reduced_mesh%de_lat(buf_j  )   &
          ) / reduced_mesh%vertex_area(buf_j) + reduced_mesh%half_f(buf_j)     &
        ) / m_vtx
      end do
    end if
#else
    if (raw_mesh%is_outside_half_lat(j+buf_j)) then
      return
    else if (raw_mesh%is_south_pole(j+buf_j) .or. raw_mesh%is_north_pole(j+buf_j+1)) then
      sign = merge(-1.0_r8, 1.0_r8, raw_mesh%full_lat(j) < 0.0_r8)
      u_j  = merge(buf_j+1, buf_j , raw_mesh%full_lat(j) < 0.0_r8)
      pole = 0.0_r8
      do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
        pole = pole + sign * reduced_state%u(i,u_j,move) * reduced_mesh%de_lon(u_j)
      end do
      call parallel_zonal_sum(pole)
      pole = pole / reduced_mesh%num_half_lon / reduced_mesh%vertex_area(buf_j)
      do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
        m_vtx = (                                                                                                          &
          (reduced_state%gd(i,buf_j  ,move) + reduced_state%gd(i+1,buf_j  ,move)) * reduced_mesh%subcell_area(2,buf_j  ) + &
          (reduced_state%gd(i,buf_j+1,move) + reduced_state%gd(i+1,buf_j+1,move)) * reduced_mesh%subcell_area(1,buf_j+1)   &
        ) / reduced_mesh%vertex_area(buf_j) / g
        reduced_state%pv(i,buf_j,move) = (pole + reduced_mesh%half_f(buf_j)) / m_vtx
      end do
    else
      do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
        m_vtx = (                                                                                                          &
          (reduced_state%gd(i,buf_j  ,move) + reduced_state%gd(i+1,buf_j  ,move)) * reduced_mesh%subcell_area(2,buf_j  ) + &
          (reduced_state%gd(i,buf_j+1,move) + reduced_state%gd(i+1,buf_j+1,move)) * reduced_mesh%subcell_area(1,buf_j+1)   &
        ) / reduced_mesh%vertex_area(buf_j) / g
        reduced_state%pv(i,buf_j,move) = (                                     &
          (                                                                    &
            reduced_state%u(i  ,buf_j  ,move) * reduced_mesh%de_lon(buf_j  ) - &
            reduced_state%u(i  ,buf_j+1,move) * reduced_mesh%de_lon(buf_j+1) + &
            reduced_state%v(i+1,buf_j  ,move) * reduced_mesh%de_lat(buf_j  ) - &
            reduced_state%v(i  ,buf_j  ,move) * reduced_mesh%de_lat(buf_j  )   &
          ) / reduced_mesh%vertex_area(buf_j) + reduced_mesh%half_f(buf_j)     &
        ) / m_vtx
      end do
    end if
#endif
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%pv(:,buf_j,move))

  end subroutine reduce_pv

  subroutine reduce_m_lon(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    if (raw_mesh%is_outside_full_lat(j+buf_j)) return
    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_state%m_lon(i,buf_j,move) = (                                          &
        reduced_mesh%lon_edge_left_area (buf_j) * reduced_state%gd(i  ,buf_j,move) + &
        reduced_mesh%lon_edge_right_area(buf_j) * reduced_state%gd(i+1,buf_j,move)   &
      ) / reduced_mesh%lon_edge_area(buf_j) / g
    end do

  end subroutine reduce_m_lon

  subroutine reduce_m_lat(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    if (reduced_mesh%lat_edge_area(buf_j) == 0) return
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
#ifdef V_POLE
      reduced_state%m_lat(i,buf_j,move) = (                                         &
        reduced_mesh%lat_edge_up_area  (buf_j) * reduced_state%gd(i,buf_j  ,move) + &
        reduced_mesh%lat_edge_down_area(buf_j) * reduced_state%gd(i,buf_j-1,move)   &
      ) / reduced_mesh%lat_edge_area(buf_j) / g
#else
      reduced_state%m_lat(i,buf_j,move) = (                                         &
        reduced_mesh%lat_edge_up_area  (buf_j) * reduced_state%gd(i,buf_j+1,move) + &
        reduced_mesh%lat_edge_down_area(buf_j) * reduced_state%gd(i,buf_j  ,move)   &
      ) / reduced_mesh%lat_edge_area(buf_j) / g
#endif
    end do

  end subroutine reduce_m_lat

  subroutine reduce_mf_lon_n(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i

    if (raw_mesh%is_outside_full_lat(j+buf_j)) return
    raw_i = raw_mesh%full_lon_start_idx + move - 1
    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      reduced_state%mf_lon_n(i,buf_j,move) = sum(raw_state%mf_lon_n(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%mf_lon_n(:,buf_j,move) = reduced_state%mf_lon_n(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lon_n(:,buf_j,move))

  end subroutine reduce_mf_lon_n

  subroutine reduce_mf_lat_n(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i

    raw_i = move
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%mf_lat_n(i,buf_j,move) = sum(raw_state%mf_lat_n(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%mf_lat_n(:,buf_j,move) = reduced_state%mf_lat_n(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lat_n(:,buf_j,move))

  end subroutine reduce_mf_lat_n

  subroutine reduce_mf_lon_t(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
#ifdef V_POLE
      reduced_state%mf_lon_t(i,buf_j,move) =                                                                                           &
        reduced_mesh%full_tangent_wgt(1,buf_j) * (reduced_state%mf_lat_n(i,buf_j  ,move) + reduced_state%mf_lat_n(i+1,buf_j  ,move)) + &
        reduced_mesh%full_tangent_wgt(2,buf_j) * (reduced_state%mf_lat_n(i,buf_j+1,move) + reduced_state%mf_lat_n(i+1,buf_j+1,move))
#else
      reduced_state%mf_lon_t(i,buf_j,move) =                                                                                           &
        reduced_mesh%full_tangent_wgt(1,buf_j) * (reduced_state%mf_lat_n(i,buf_j-1,move) + reduced_state%mf_lat_n(i+1,buf_j-1,move)) + &
        reduced_mesh%full_tangent_wgt(2,buf_j) * (reduced_state%mf_lat_n(i,buf_j  ,move) + reduced_state%mf_lat_n(i+1,buf_j  ,move))
#endif
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lon_t(:,buf_j,move))

  end subroutine reduce_mf_lon_t

  subroutine reduce_mf_lat_t(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    if (is_inf(reduced_mesh%half_tangent_wgt(1,buf_j)) .or. is_inf(reduced_mesh%half_tangent_wgt(2,buf_j))) return
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
#ifdef V_POLE
      reduced_state%mf_lat_t(i,buf_j,move) =                                                                                           &
        reduced_mesh%half_tangent_wgt(1,buf_j) * (reduced_state%mf_lon_n(i-1,buf_j-1,move) + reduced_state%mf_lon_n(i,buf_j-1,move)) + &
        reduced_mesh%half_tangent_wgt(2,buf_j) * (reduced_state%mf_lon_n(i-1,buf_j  ,move) + reduced_state%mf_lon_n(i,buf_j  ,move))
#else
      reduced_state%mf_lat_t(i,buf_j,move) =                                                                                           &
        reduced_mesh%half_tangent_wgt(1,buf_j) * (reduced_state%mf_lon_n(i-1,buf_j  ,move) + reduced_state%mf_lon_n(i,buf_j  ,move)) + &
        reduced_mesh%half_tangent_wgt(2,buf_j) * (reduced_state%mf_lon_n(i-1,buf_j+1,move) + reduced_state%mf_lon_n(i,buf_j+1,move))
#endif
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%mf_lat_t(:,buf_j,move))

  end subroutine reduce_mf_lat_t

  subroutine reduce_dpv_lon_t(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
#ifdef V_POLE
      reduced_state%dpv_lon_t(i,buf_j,move) = reduced_state%pv(i,buf_j+1,move) - reduced_state%pv(i,buf_j  ,move)
#else
      reduced_state%dpv_lon_t(i,buf_j,move) = reduced_state%pv(i,buf_j  ,move) - reduced_state%pv(i,buf_j-1,move)
#endif
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%dpv_lon_t(:,buf_j,move))

  end subroutine reduce_dpv_lon_t

  subroutine reduce_dpv_lat_t(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%dpv_lat_t(i,buf_j,move) = reduced_state%pv(i+1,buf_j,move) - reduced_state%pv(i,buf_j,move)
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%dpv_lat_t(:,buf_j,move))

  end subroutine reduce_dpv_lat_t

  subroutine reduce_dpv_lon_n(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
#ifdef V_POLE
      reduced_state%dpv_lon_n(i,buf_j,move) = 0.25_r8 * ( &
        reduced_state%dpv_lat_t(i  ,buf_j  ,move) +       &
        reduced_state%dpv_lat_t(i+1,buf_j  ,move) +       &
        reduced_state%dpv_lat_t(i  ,buf_j+1,move) +       &
        reduced_state%dpv_lat_t(i+1,buf_j+1,move)         &
      )
#else
      reduced_state%dpv_lon_n(i,buf_j,move) = 0.25_r8 * ( &
        reduced_state%dpv_lat_t(i  ,buf_j-1,move) +       &
        reduced_state%dpv_lat_t(i+1,buf_j-1,move) +       &
        reduced_state%dpv_lat_t(i  ,buf_j  ,move) +       &
        reduced_state%dpv_lat_t(i+1,buf_j  ,move)         &
      )
#endif
    end do

  end subroutine reduce_dpv_lon_n

  subroutine reduce_dpv_lat_n(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer i

    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
#ifdef V_POLE
      reduced_state%dpv_lat_n(i,buf_j,move) = 0.25_r8 * ( &
        reduced_state%dpv_lon_t(i-1,buf_j-1,move) +       &
        reduced_state%dpv_lon_t(i  ,buf_j-1,move) +       &
        reduced_state%dpv_lon_t(i-1,buf_j  ,move) +       &
        reduced_state%dpv_lon_t(i  ,buf_j  ,move)         &
      )
#else
      reduced_state%dpv_lat_n(i,buf_j,move) = 0.25_r8 * ( &
        reduced_state%dpv_lon_t(i-1,buf_j  ,move) +       &
        reduced_state%dpv_lon_t(i  ,buf_j  ,move) +       &
        reduced_state%dpv_lon_t(i-1,buf_j+1,move) +       &
        reduced_state%dpv_lon_t(i  ,buf_j+1,move)         &
      )
#endif
    end do

  end subroutine reduce_dpv_lat_n

  subroutine reduce_pv_lon_apvm(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    real(r8) u, v, le, de
    integer i

    le = reduced_mesh%le_lon(buf_j)
    de = reduced_mesh%de_lon(buf_j)
    if (le == inf .or. de == inf) return
    do i = reduced_mesh%half_lon_start_idx, reduced_mesh%half_lon_end_idx
      u = reduced_state%u(i,buf_j,move)
      v = reduced_state%mf_lon_t(i,buf_j,move) / reduced_state%m_lon(i,buf_j,move)
#ifdef V_POLE
      reduced_state%pv_lon(i,buf_j,move) = 0.5_r8 * (    &
        reduced_state%pv(i,buf_j+1,move) +               &
        reduced_state%pv(i,buf_j  ,move)                 &
      ) - 0.5_r8 * (                                     &
        u * reduced_state%dpv_lon_n(i,buf_j,move) / de + &
        v * reduced_state%dpv_lon_t(i,buf_j,move) / le   &
      ) * dt
#else
      reduced_state%pv_lon(i,buf_j,move) = 0.5_r8 * (    &
        reduced_state%pv(i,buf_j-1,move) +               &
        reduced_state%pv(i,buf_j  ,move)                 &
      ) - 0.5_r8 * (                                     &
        u * reduced_state%dpv_lon_n(i,buf_j,move) / de + &
        v * reduced_state%dpv_lon_t(i,buf_j,move) / le   &
      ) * dt
#endif
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%pv_lon(:,buf_j,move))

  end subroutine reduce_pv_lon_apvm

  subroutine reduce_pv_lat_apvm(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    real(r8) u, v, le, de
    integer i

    le = reduced_mesh%le_lat(buf_j)
    de = reduced_mesh%de_lat(buf_j)
    if (le == inf .or. de == inf) return
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      u = reduced_state%mf_lat_t(i,buf_j,move) / reduced_state%m_lat(i,buf_j,move)
      v = reduced_state%v(i,buf_j,move)
      reduced_state%pv_lat(i,buf_j,move) = 0.5_r8 * (    &
        reduced_state%pv(i  ,buf_j,move) +               &
        reduced_state%pv(i-1,buf_j,move)                 &
      ) - 0.5_r8 * (                                     &
        u * reduced_state%dpv_lat_t(i,buf_j,move) / le + &
        v * reduced_state%dpv_lat_n(i,buf_j,move) / de   &
      ) * dt
    end do
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%pv_lat(:,buf_j,move))

  end subroutine reduce_pv_lat_apvm

  subroutine reduce_ke(j, buf_j, move, raw_mesh, raw_state, reduced_mesh, reduced_state, dt)

    integer, intent(in) :: j
    integer, intent(in) :: buf_j
    integer, intent(in) :: move
    type(mesh_type), intent(in) :: raw_mesh
    type(state_type), intent(in) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_state_type), intent(inout) :: reduced_state
    real(r8), intent(in) :: dt

    integer raw_i, i

    raw_i = raw_mesh%full_lon_start_idx + move - 1
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      reduced_state%ke(i,buf_j,move) = sum(raw_state%ke(raw_i:raw_i+reduced_mesh%reduce_factor-1,j+buf_j))
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do
    reduced_state%ke(:,buf_j,move) = reduced_state%ke(:,buf_j,move) / reduced_mesh%reduce_factor
    call parallel_fill_halo(reduced_mesh%halo_width, reduced_state%ke(:,buf_j,move))

  end subroutine reduce_ke

  subroutine reduce_append_array(move, reduced_mesh, reduced_array, raw_mesh, raw_array)

    integer, intent(in) :: move
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    real(r8), intent(in) :: reduced_array(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub)
    type(mesh_type), intent(in) :: raw_mesh
    real(r8), intent(inout) :: raw_array(raw_mesh%full_lon_lb:raw_mesh%full_lon_ub)

    integer i, raw_i

    raw_i = raw_mesh%full_lon_start_idx + move - 1
    do i = reduced_mesh%full_lon_start_idx, reduced_mesh%full_lon_end_idx
      raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1) = raw_array(raw_i:raw_i+reduced_mesh%reduce_factor-1) + reduced_array(i) / reduced_mesh%reduce_factor
      raw_i = raw_i + reduced_mesh%reduce_factor
    end do

  end subroutine reduce_append_array

  subroutine allocate_reduced_static(reduced_mesh, reduced_static)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    allocate(reduced_static%ghs(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub, 0:0,reduced_mesh%reduce_factor))

  end subroutine allocate_reduced_static

  subroutine allocate_reduced_tend(reduced_mesh, reduced_tend)

    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_tend_type), intent(inout) :: reduced_tend

    allocate(reduced_tend%qhu    (reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub))
    allocate(reduced_tend%qhv    (reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub))
    allocate(reduced_tend%dmfdlon(reduced_mesh%full_lon_lb:reduced_mesh%full_lon_ub))
    allocate(reduced_tend%dpedlon(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub))
    allocate(reduced_tend%dkedlon(reduced_mesh%half_lon_lb:reduced_mesh%half_lon_ub))

  end subroutine allocate_reduced_tend

  subroutine reduce_static(j, raw_mesh, raw_static, reduced_mesh, reduced_static)

    integer, intent(in) :: j
    type(mesh_type), intent(in) :: raw_mesh
    type(static_type), intent(in) :: raw_static
    type(reduced_mesh_type), intent(in) :: reduced_mesh
    type(reduced_static_type), intent(inout) :: reduced_static

    integer buf_j, move

    do move = 1, reduced_mesh%reduce_factor
      do buf_j = lbound(reduced_static%ghs, 2), ubound(reduced_static%ghs, 2)
        call reduce_ghs(j, buf_j, move, raw_mesh, raw_static, reduced_mesh, reduced_static)
      end do
    end do

  end subroutine reduce_static

  subroutine reduce_replace_pv(reduced_mesh, reduced_state, raw_state)

    type(state_type), intent(inout) :: raw_state
    type(reduced_mesh_type), intent(in) :: reduced_mesh(raw_state%mesh%full_lat_lb:raw_state%mesh%full_lat_ub)
    type(reduced_state_type), intent(in) :: reduced_state(raw_state%mesh%full_lat_start_idx:raw_state%mesh%full_lat_end_idx)

    integer j, move

    do j = raw_state%mesh%half_lat_start_idx_no_pole, raw_state%mesh%half_lat_end_idx_no_pole
      if (reduced_mesh(j)%reduce_factor > 0) then
        raw_state%pv(:,j) = 0.0
        do move = 1, reduced_mesh(j)%reduce_factor
          call reduce_append_array(move, reduced_mesh(j), reduced_state(j)%pv(:,0,move), raw_state%mesh, raw_state%pv(:,j))
        end do
        call parallel_overlay_inner_halo(raw_state%mesh, raw_state%pv(:,j), left_halo=.true.)
      end if
    end do

  end subroutine reduce_replace_pv

  subroutine reduced_static_final(this)

    type(reduced_static_type), intent(inout) :: this

    if (allocated(this%ghs)) deallocate(this%ghs)

  end subroutine reduced_static_final

  subroutine reduced_state_final(this)

    type(reduced_state_type), intent(inout) :: this

    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%gd       )) deallocate(this%gd       )
    if (allocated(this%pv       )) deallocate(this%pv       )
    if (allocated(this%pv_lon   )) deallocate(this%pv_lon   )
    if (allocated(this%pv_lat   )) deallocate(this%pv_lat   )
    if (allocated(this%dpv_lon_n)) deallocate(this%dpv_lon_n)
    if (allocated(this%dpv_lat_n)) deallocate(this%dpv_lat_n)
    if (allocated(this%dpv_lon_t)) deallocate(this%dpv_lon_t)
    if (allocated(this%dpv_lat_t)) deallocate(this%dpv_lat_t)
    if (allocated(this%m_lon    )) deallocate(this%m_lon    )
    if (allocated(this%m_lat    )) deallocate(this%m_lat    )
    if (allocated(this%mf_lon_n )) deallocate(this%mf_lon_n )
    if (allocated(this%mf_lon_t )) deallocate(this%mf_lon_t )
    if (allocated(this%mf_lat_n )) deallocate(this%mf_lat_n )
    if (allocated(this%mf_lat_t )) deallocate(this%mf_lat_t )
    if (allocated(this%ke       )) deallocate(this%ke       )

  end subroutine reduced_state_final

  subroutine reduced_tend_final(this)

    type(reduced_tend_type), intent(inout) :: this

    if (allocated(this%qhv    )) deallocate(this%qhv    )
    if (allocated(this%qhu    )) deallocate(this%qhu    )
    if (allocated(this%dmfdlon)) deallocate(this%dmfdlon)
    if (allocated(this%dpedlon)) deallocate(this%dpedlon)
    if (allocated(this%dkedlon)) deallocate(this%dkedlon)

  end subroutine reduced_tend_final

  subroutine reduce_final()

    if (allocated(reduced_mesh  )) deallocate(reduced_mesh  )
    if (allocated(reduced_static)) deallocate(reduced_static)
    if (allocated(reduced_state )) deallocate(reduced_state )
    if (allocated(reduced_tend  )) deallocate(reduced_tend  )

  end subroutine reduce_final

end module reduce_mod
