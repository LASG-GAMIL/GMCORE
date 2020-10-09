module tend_mod

  use const_mod
  use mesh_mod
  use namelist_mod
  use allocator_mod

  implicit none

  private

  public tend_type

  type tend_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:,:) :: du
    real(r8), allocatable, dimension(:,:,:) :: dv
    real(r8), allocatable, dimension(:,:,:) :: dgz
    real(r8), allocatable, dimension(:,:,:) :: dpt
    real(r8), allocatable, dimension(:,:  ) :: dphs
    logical :: updated_du   = .false.
    logical :: updated_dv   = .false.
    logical :: updated_dgz  = .false.
    logical :: updated_dpt  = .false.
    logical :: updated_dphs = .false.
    logical :: copy_gz  = .false.
    logical :: copy_pt  = .false.
    logical :: copy_phs = .false.
    ! Individual tendencies
    real(r8), allocatable, dimension(:,:,:) :: qhv
    real(r8), allocatable, dimension(:,:,:) :: qhu
    real(r8), allocatable, dimension(:,:,:) :: dkedlon
    real(r8), allocatable, dimension(:,:,:) :: dkedlat
    real(r8), allocatable, dimension(:,:,:) :: dmfdlon
    real(r8), allocatable, dimension(:,:,:) :: dmfdlat
    real(r8), allocatable, dimension(:,:,:) :: dptfdlon ! Zonal potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: dptfdlat ! Meridional potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: dptfdlev ! Vertical potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: pgf_lon
    real(r8), allocatable, dimension(:,:,:) :: pgf_lat
    real(r8), allocatable, dimension(:,:,:) :: wedudlev
    real(r8), allocatable, dimension(:,:,:) :: wedvdlev
    real(r8), allocatable, dimension(:,:,:) :: dvordlon
    real(r8), allocatable, dimension(:,:,:) :: dvordlat
  contains
    procedure :: init => tend_init
    procedure :: reset_flags => tend_reset_flags
    procedure :: clear => tend_clear
    final :: tend_final
  end type tend_type

contains

  subroutine tend_init(this, mesh)

    class(tend_type), intent(inout)         :: this
    type(mesh_type ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%du      , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dv      , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dgz     , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpt     , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dphs    , full_lon=.true., full_lat=.true.                 )
    call allocate_array(mesh, this%qhv     , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%qhu     , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dkedlon , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dkedlat , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmfdlon , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmfdlat , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlon, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlat, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlev, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pgf_lon , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pgf_lat , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%wedudlev, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%wedvdlev, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dvordlon, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dvordlat, half_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine tend_init

  subroutine tend_reset_flags(this)

    class(tend_type), intent(inout) :: this

    this%updated_du   = .false.
    this%updated_dv   = .false.
    this%updated_dgz  = .false.
    this%updated_dpt  = .false.
    this%updated_dphs = .false.

    this%copy_gz  = .false.
    this%copy_pt  = .false.
    this%copy_phs = .false.

  end subroutine tend_reset_flags

  subroutine tend_clear(this)

    class(tend_type), intent(inout) :: this

    if (allocated(this%du      )) deallocate(this%du      )
    if (allocated(this%dv      )) deallocate(this%dv      )
    if (allocated(this%dgz     )) deallocate(this%dgz     )
    if (allocated(this%dpt     )) deallocate(this%dpt     )
    if (allocated(this%dphs    )) deallocate(this%dphs    )
    if (allocated(this%qhv     )) deallocate(this%qhv     )
    if (allocated(this%qhu     )) deallocate(this%qhu     )
    if (allocated(this%dkedlon )) deallocate(this%dkedlon )
    if (allocated(this%dkedlat )) deallocate(this%dkedlat )
    if (allocated(this%dmfdlon )) deallocate(this%dmfdlon )
    if (allocated(this%dmfdlat )) deallocate(this%dmfdlat )
    if (allocated(this%dptfdlon)) deallocate(this%dptfdlon)
    if (allocated(this%dptfdlat)) deallocate(this%dptfdlat)
    if (allocated(this%dptfdlev)) deallocate(this%dptfdlev)
    if (allocated(this%pgf_lon )) deallocate(this%pgf_lon )
    if (allocated(this%pgf_lat )) deallocate(this%pgf_lat )
    if (allocated(this%wedudlev)) deallocate(this%wedudlev)
    if (allocated(this%wedvdlev)) deallocate(this%wedvdlev)
    if (allocated(this%dvordlon)) deallocate(this%dvordlon)
    if (allocated(this%dvordlat)) deallocate(this%dvordlat)

  end subroutine tend_clear

  subroutine tend_final(this)

    type(tend_type), intent(inout) :: this

    call this%clear()

  end subroutine tend_final

end module tend_mod
