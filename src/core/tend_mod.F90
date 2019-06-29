module tend_mod

  use const_mod
  use mesh_mod
  use allocator_mod

  implicit none

  private

  public tend_type
  public tends
  public create_tends

  type tend_type
    type(mesh_type), pointer :: mesh => null()
    real(real_kind), allocatable, dimension(:,:) :: du
    real(real_kind), allocatable, dimension(:,:) :: dv
    real(real_kind), allocatable, dimension(:,:) :: dgd
    ! Individual tendencies
    real(real_kind), allocatable, dimension(:,:) :: qhv
    real(real_kind), allocatable, dimension(:,:) :: qhu
    real(real_kind), allocatable, dimension(:,:) :: dEdlon
    real(real_kind), allocatable, dimension(:,:) :: dEdlat
    real(real_kind), allocatable, dimension(:,:) :: div_mass_flux
    ! Derived variables
    real(real_kind), allocatable, dimension(:,:) :: dpvdlon
    real(real_kind), allocatable, dimension(:,:) :: dpvdlat
  contains
    procedure :: init => tend_init
    procedure :: clear => tend_clear
    final :: tend_final
  end type tend_type

  type(tend_type), allocatable :: tends(:)

contains

  subroutine create_tends()

    integer i

    if (.not. allocated(tends)) then
      allocate(tends(0:2))
      do i = lbound(tends, 1), ubound(tends, 1)
        call tends(i)%init(mesh)
      end do
    end if

  end subroutine create_tends

  subroutine tend_init(this, mesh)

    class(tend_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%du           , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dv           , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dgd          , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%qhu          , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%qhv          , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dEdlon       , half_lon=.true., full_lon=.true.)
    call allocate_array(mesh, this%dEdlat       , full_lon=.true., half_lon=.true.)
    call allocate_array(mesh, this%div_mass_flux, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dpvdlon      , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dpvdlat      , half_lon=.true., full_lat=.true.)

  end subroutine tend_init

  subroutine tend_clear(this)

    class(tend_type), intent(inout) :: this

    if (allocated(this%du           )) deallocate(this%du           )
    if (allocated(this%dv           )) deallocate(this%dv           )
    if (allocated(this%dgd          )) deallocate(this%dgd          )
    if (allocated(this%qhu          )) deallocate(this%qhu          )
    if (allocated(this%qhv          )) deallocate(this%qhv          )
    if (allocated(this%dEdlon       )) deallocate(this%dEdlon       )
    if (allocated(this%dEdlat       )) deallocate(this%dEdlat       )
    if (allocated(this%div_mass_flux)) deallocate(this%div_mass_flux)
    if (allocated(this%dpvdlon      )) deallocate(this%dpvdlon      )
    if (allocated(this%dpvdlat      )) deallocate(this%dpvdlat      )

  end subroutine tend_clear

  subroutine tend_final(this)

    type(tend_type), intent(inout) :: this

    call this%clear()

  end subroutine tend_final

end module tend_mod