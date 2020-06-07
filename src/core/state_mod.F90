module state_mod

  use const_mod
  use mesh_mod
  use allocator_mod
  use parallel_types_mod

  implicit none

  private

  public state_type

  type state_type
    type(mesh_type), pointer :: mesh => null()
    ! For nesting
    integer :: id = 0
    type(state_type), pointer :: parent => null()
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:) :: w          ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: gz
    real(r8), allocatable, dimension(:,:,:) :: gzi        ! Geopotential height on half levels
    real(r8), allocatable, dimension(:,:,:) :: m
    real(r8), allocatable, dimension(:,:,:) :: m_vtx
    real(r8), allocatable, dimension(:,:,:) :: m_lon
    real(r8), allocatable, dimension(:,:,:) :: m_lat
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_n
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_n
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_t
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_t
    real(r8), allocatable, dimension(:,:,:) :: pv         ! Potential vorticity
    real(r8), allocatable, dimension(:,:,:) :: pv_lon
    real(r8), allocatable, dimension(:,:,:) :: pv_lat
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_t
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_n
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_t
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_n
    real(r8), allocatable, dimension(:,:,:) :: ke         ! Kinetic energy
    real(r8), allocatable, dimension(:,:,:) :: detadt     ! Vertical coordinate speed
    real(r8), allocatable, dimension(:,:,:) :: tp         ! Potential temperature
    real(r8), allocatable, dimension(:,:,:) :: t          ! Temperature
    real(r8), allocatable, dimension(:,:,:) :: p          ! Pressure on full levels
    real(r8), allocatable, dimension(:,:,:) :: pi         ! Pressure on half levels
    real(r8), allocatable, dimension(:,:,:) :: ph         ! Hydrostatic pressure on full levels
    real(r8), allocatable, dimension(:,:,:) :: phi        ! Hydrostatic pressure on half levels
    real(r8), allocatable, dimension(:,:  ) :: phs        ! Surface hydrostatic pressure
    real(r8) tm
    real(r8) te
    real(r8) tpe
    real(r8) tav
    type(async_type) async(11)
  contains
    procedure :: init => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

contains

  subroutine state_init(this, mesh)

    class(state_type), intent(inout)         :: this
    type(mesh_type  ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%u        , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v        , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%gz       , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m        , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_vtx    , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_lon    , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_lat    , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lon_n , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lon_t , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lat_n , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lat_t , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv       , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv_lon   , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv_lat   , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpv_lon_t, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpv_lon_n, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpv_lat_t, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpv_lat_n, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ke       , full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%u        )) deallocate(this%u        )
    if (allocated(this%v        )) deallocate(this%v        )
    if (allocated(this%gz       )) deallocate(this%gz       )
    if (allocated(this%m        )) deallocate(this%m        )
    if (allocated(this%m_vtx    )) deallocate(this%m_vtx    )
    if (allocated(this%m_lon    )) deallocate(this%m_lon    )
    if (allocated(this%m_lat    )) deallocate(this%m_lat    )
    if (allocated(this%mf_lon_n )) deallocate(this%mf_lon_n )
    if (allocated(this%mf_lat_n )) deallocate(this%mf_lat_n )
    if (allocated(this%mf_lat_t )) deallocate(this%mf_lat_t )
    if (allocated(this%mf_lon_t )) deallocate(this%mf_lon_t )
    if (allocated(this%pv       )) deallocate(this%pv       )
    if (allocated(this%pv_lon   )) deallocate(this%pv_lon   )
    if (allocated(this%pv_lat   )) deallocate(this%pv_lat   )
    if (allocated(this%dpv_lon_t)) deallocate(this%dpv_lon_t)
    if (allocated(this%dpv_lon_n)) deallocate(this%dpv_lon_n)
    if (allocated(this%dpv_lat_t)) deallocate(this%dpv_lat_t)
    if (allocated(this%dpv_lat_n)) deallocate(this%dpv_lat_n)
    if (allocated(this%ke       )) deallocate(this%ke       )

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
