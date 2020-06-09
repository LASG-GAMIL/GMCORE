module reduced_types_mod

  use const_mod
  use parallel_types_mod

  implicit none

  type reduced_mesh_type
    integer :: reduce_factor = 0
    integer :: damp_order = 0
    integer halo_width
    integer num_full_lon
    integer num_half_lon
    integer full_lon_ibeg
    integer full_lon_iend
    integer half_lon_ibeg
    integer half_lon_iend
    integer full_lon_lb
    integer full_lon_ub
    integer half_lon_lb
    integer half_lon_ub
    integer num_full_lev
    integer num_half_lev
    integer full_lev_ibeg
    integer full_lev_iend
    integer half_lev_ibeg
    integer half_lev_iend
    integer full_lev_lb
    integer full_lev_ub
    integer half_lev_lb
    integer half_lev_ub
#ifdef V_POLE
    real(r8), dimension(  -1:1) :: full_lat         = inf
    real(r8), dimension(  -1:1) :: half_lat         = inf
    real(r8), dimension(  -1:1) :: area_cell        = 0
    real(r8), dimension(2,-2:2) :: area_subcell     = 0
    real(r8), dimension(  -2:2) :: area_lon         = 0
    real(r8), dimension(  -2:2) :: area_lon_west    = 0
    real(r8), dimension(  -2:2) :: area_lon_east    = 0
    real(r8), dimension(  -1:2) :: area_vtx         = 0
    real(r8), dimension(  -1:2) :: area_lat         = 0
    real(r8), dimension(  -1:2) :: area_lat_north   = 0
    real(r8), dimension(  -1:2) :: area_lat_south   = 0
    real(r8), dimension(  -1:1) :: le_lon           = 0
    real(r8), dimension(  -2:2) :: de_lon           = 0
    real(r8), dimension(  -1:2) :: le_lat           = 0
    real(r8), dimension(  -1:2) :: de_lat           = 0
    real(r8), dimension(2,-1:1) :: full_tangent_wgt = inf
    real(r8), dimension(2, 0:1) :: half_tangent_wgt = inf
    real(r8), dimension(  -1:2) :: half_f           = inf
#else
    real(r8), dimension(  -1:1) :: full_lat         = inf
    real(r8), dimension(  -1:1) :: half_lat         = inf
    real(r8), dimension(  -1:1) :: area_cell        = 0
    real(r8), dimension(2,-2:2) :: area_subcell     = 0
    real(r8), dimension(  -2:2) :: area_lon         = 0
    real(r8), dimension(  -2:2) :: area_lon_west    = 0
    real(r8), dimension(  -2:2) :: area_lon_east    = 0
    real(r8), dimension(  -2:1) :: area_vtx         = 0
    real(r8), dimension(  -2:1) :: area_lat         = 0
    real(r8), dimension(  -2:1) :: area_lat_north   = 0
    real(r8), dimension(  -2:1) :: area_lat_south   = 0
    real(r8), dimension(  -1:1) :: le_lon           = inf
    real(r8), dimension(  -2:2) :: de_lon           = inf
    real(r8), dimension(  -2:1) :: le_lat           = inf
    real(r8), dimension(  -2:1) :: de_lat           = inf
    real(r8), dimension(2,-1:1) :: full_tangent_wgt = inf
    real(r8), dimension(2,-1:0) :: half_tangent_wgt = inf
    real(r8), dimension(  -2:1) :: half_f           = inf
#endif
  end type reduced_mesh_type

  type reduced_static_type
    real(r8), allocatable, dimension(:,:,:) :: gzs
  contains
    final :: reduced_static_final
  end type reduced_static_type

  type reduced_state_type
    real(r8), allocatable, dimension(:,:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:,:) :: gz
    real(r8), allocatable, dimension(:,:,:,:) :: m
    real(r8), allocatable, dimension(:,:,:,:) :: m_lon
    real(r8), allocatable, dimension(:,:,:,:) :: m_lat
    real(r8), allocatable, dimension(:,:,:,:) :: mf_lon_n
    real(r8), allocatable, dimension(:,:,:,:) :: mf_lon_t
    real(r8), allocatable, dimension(:,:,:,:) :: mf_lat_n
    real(r8), allocatable, dimension(:,:,:,:) :: mf_lat_t
    real(r8), allocatable, dimension(:,:,:,:) :: pv
    real(r8), allocatable, dimension(:,:,:,:) :: pv_lon
    real(r8), allocatable, dimension(:,:,:,:) :: pv_lat
    real(r8), allocatable, dimension(:,:,:,:) :: dpv_lon_t
    real(r8), allocatable, dimension(:,:,:,:) :: dpv_lat_t
    real(r8), allocatable, dimension(:,:,:,:) :: dpv_lon_n
    real(r8), allocatable, dimension(:,:,:,:) :: dpv_lat_n
    real(r8), allocatable, dimension(:,:,:,:) :: ke
    type(async_type), allocatable :: async(:,:,:)
  contains
    final :: reduced_state_final
  end type reduced_state_type

  type reduced_tend_type
    real(r8), allocatable, dimension(:,:) :: qhv
    real(r8), allocatable, dimension(:,:) :: qhu
    real(r8), allocatable, dimension(:,:) :: dmfdlon
    real(r8), allocatable, dimension(:,:) :: dpedlon
    real(r8), allocatable, dimension(:,:) :: dkedlon
  contains
    final :: reduced_tend_final
  end type reduced_tend_type

contains

  subroutine reduced_static_final(this)

    type(reduced_static_type), intent(inout) :: this

    if (allocated(this%gzs)) deallocate(this%gzs)

  end subroutine reduced_static_final

  subroutine reduced_state_final(this)

    type(reduced_state_type), intent(inout) :: this

    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%gz       )) deallocate(this%gz       )
    if (allocated(this%m        )) deallocate(this%m        )
    if (allocated(this%m_lon    )) deallocate(this%m_lon    )
    if (allocated(this%m_lat    )) deallocate(this%m_lat    )
    if (allocated(this%mf_lon_n )) deallocate(this%mf_lon_n )
    if (allocated(this%mf_lon_t )) deallocate(this%mf_lon_t )
    if (allocated(this%mf_lat_n )) deallocate(this%mf_lat_n )
    if (allocated(this%mf_lat_t )) deallocate(this%mf_lat_t )
    if (allocated(this%pv       )) deallocate(this%pv       )
    if (allocated(this%pv_lon   )) deallocate(this%pv_lon   )
    if (allocated(this%pv_lat   )) deallocate(this%pv_lat   )
    if (allocated(this%dpv_lon_n)) deallocate(this%dpv_lon_n)
    if (allocated(this%dpv_lat_n)) deallocate(this%dpv_lat_n)
    if (allocated(this%dpv_lon_t)) deallocate(this%dpv_lon_t)
    if (allocated(this%dpv_lat_t)) deallocate(this%dpv_lat_t)
    if (allocated(this%ke       )) deallocate(this%ke       )
    if (allocated(this%async    )) deallocate(this%async    )

  end subroutine reduced_state_final

  subroutine reduced_tend_final(this)

    type(reduced_tend_type), intent(inout) :: this

    if (allocated(this%qhv    )) deallocate(this%qhv    )
    if (allocated(this%qhu    )) deallocate(this%qhu    )
    if (allocated(this%dmfdlon)) deallocate(this%dmfdlon)
    if (allocated(this%dpedlon)) deallocate(this%dpedlon)
    if (allocated(this%dkedlon)) deallocate(this%dkedlon)

  end subroutine reduced_tend_final

end module reduced_types_mod
