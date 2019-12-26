module tend_mod

  use const_mod
  use mesh_mod
  use namelist_mod
  use allocator_mod

  implicit none

  private

  public tend_type
  public tends
  public tend_init_root

  type tend_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:) :: du
    real(r8), allocatable, dimension(:,:) :: dv
    real(r8), allocatable, dimension(:,:) :: dgd
    ! Individual tendencies
    real(r8), allocatable, dimension(:,:) :: qhv
    real(r8), allocatable, dimension(:,:) :: qhu
    real(r8), allocatable, dimension(:,:) :: dpedlon
    real(r8), allocatable, dimension(:,:) :: dkedlon
    real(r8), allocatable, dimension(:,:) :: dpedlat
    real(r8), allocatable, dimension(:,:) :: dkedlat
    real(r8), allocatable, dimension(:,:) :: dmfdlon
    real(r8), allocatable, dimension(:,:) :: dmfdlat
  contains
    procedure :: init => tend_init
    procedure :: clear => tend_clear
    final :: tend_final
  end type tend_type

  type(tend_type), allocatable :: tends(:)

contains

  subroutine tend_init_root()

    integer i

    if (.not. allocated(tends)) then
      select case (trim(time_scheme))
      case ('pc2', 'rk2')
        allocate(tends(3))
      case ('rk3')
        allocate(tends(4))
      case ('rk4')
        allocate(tends(5))
      end select
      do i = lbound(tends, 1), ubound(tends, 1)
        call tends(i)%init(global_mesh)
      end do
    end if

  end subroutine tend_init_root

  subroutine tend_init(this, mesh)

    class(tend_type), intent(inout)         :: this
    type(mesh_type ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%du     , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dv     , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dgd    , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%qhu    , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%qhv    , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dpedlon, half_lon=.true., full_lon=.true.)
    call allocate_array(mesh, this%dkedlon, half_lon=.true., full_lon=.true.)
    call allocate_array(mesh, this%dpedlat, full_lon=.true., half_lon=.true.)
    call allocate_array(mesh, this%dkedlat, full_lon=.true., half_lon=.true.)
    call allocate_array(mesh, this%dmfdlon, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dmfdlat, full_lon=.true., full_lat=.true.)

  end subroutine tend_init

  subroutine tend_clear(this)

    class(tend_type), intent(inout) :: this

    if (allocated(this%du     )) deallocate(this%du     )
    if (allocated(this%dv     )) deallocate(this%dv     )
    if (allocated(this%dgd    )) deallocate(this%dgd    )
    if (allocated(this%qhu    )) deallocate(this%qhu    )
    if (allocated(this%qhv    )) deallocate(this%qhv    )
    if (allocated(this%dpedlon)) deallocate(this%dpedlon)
    if (allocated(this%dpedlon)) deallocate(this%dkedlon)
    if (allocated(this%dkedlat)) deallocate(this%dpedlat)
    if (allocated(this%dkedlat)) deallocate(this%dkedlat)
    if (allocated(this%dmfdlon)) deallocate(this%dmfdlon)
    if (allocated(this%dmfdlat)) deallocate(this%dmfdlat)

  end subroutine tend_clear

  subroutine tend_final(this)

    type(tend_type), intent(inout) :: this

    call this%clear()

  end subroutine tend_final

end module tend_mod