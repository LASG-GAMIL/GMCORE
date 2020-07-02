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
    real(r8), allocatable, dimension(:,:,:) :: dphs
    ! Individual tendencies
    real(r8), allocatable, dimension(:,:,:) :: qhv
    real(r8), allocatable, dimension(:,:,:) :: qhu
    real(r8), allocatable, dimension(:,:,:) :: dpedlon
    real(r8), allocatable, dimension(:,:,:) :: dkedlon
    real(r8), allocatable, dimension(:,:,:) :: dpedlat
    real(r8), allocatable, dimension(:,:,:) :: dkedlat
    real(r8), allocatable, dimension(:,:,:) :: dmfdlon
    real(r8), allocatable, dimension(:,:,:) :: dmfdlat
    real(r8), allocatable, dimension(:,:,:) :: dptfdlon ! Zonal potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: dptfdlat ! Meridional potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: dptfdlev ! Vertical potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: dpdlon
    real(r8), allocatable, dimension(:,:,:) :: dpdlat
    real(r8), allocatable, dimension(:,:,:) :: wedudlev
    real(r8), allocatable, dimension(:,:,:) :: wedvdlev
  contains
    procedure :: init => tend_init
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
    call allocate_array(mesh, this%dpedlon , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dkedlon , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpedlat , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dkedlat , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmfdlon , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmfdlat , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlon, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlon, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlon, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpdlon  , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpdlat  , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%wedudlev, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%wedvdlev, full_lon=.true., half_lat=.true., full_lev=.true.)

  end subroutine tend_init

  subroutine tend_clear(this)

    class(tend_type), intent(inout) :: this

    if (allocated(this%du      )) deallocate(this%du      )
    if (allocated(this%dv      )) deallocate(this%dv      )
    if (allocated(this%dgz     )) deallocate(this%dgz     )
    if (allocated(this%dpt     )) deallocate(this%dpt     )
    if (allocated(this%dphs    )) deallocate(this%dphs    )
    if (allocated(this%qhv     )) deallocate(this%qhv     )
    if (allocated(this%qhu     )) deallocate(this%qhu     )
    if (allocated(this%dpedlon )) deallocate(this%dpedlon )
    if (allocated(this%dpedlon )) deallocate(this%dkedlon )
    if (allocated(this%dkedlat )) deallocate(this%dpedlat )
    if (allocated(this%dkedlat )) deallocate(this%dkedlat )
    if (allocated(this%dmfdlon )) deallocate(this%dmfdlon )
    if (allocated(this%dmfdlat )) deallocate(this%dmfdlat )
    if (allocated(this%dptfdlon)) deallocate(this%dptfdlon)
    if (allocated(this%dptfdlat)) deallocate(this%dptfdlat)
    if (allocated(this%dptfdlev)) deallocate(this%dptfdlev)
    if (allocated(this%dpdlon  )) deallocate(this%dpdlon  )
    if (allocated(this%dpdlat  )) deallocate(this%dpdlat  )
    if (allocated(this%wedudlev)) deallocate(this%wedudlev)
    if (allocated(this%wedvdlev)) deallocate(this%wedvdlev)

  end subroutine tend_clear

  subroutine tend_final(this)

    type(tend_type), intent(inout) :: this

    call this%clear()

  end subroutine tend_final

end module tend_mod
