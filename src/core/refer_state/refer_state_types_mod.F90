module refer_state_types_mod

  use const_mod
  use allocator_mod
  use mesh_mod
  use static_mod

  implicit none

  private

  public refer_state_type

  type refer_profile_type
    integer num_lev
  end type refer_profile_type

  type refer_state_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:  ) :: phs
    real(r8), allocatable, dimension(:,:,:) :: ph
    real(r8), allocatable, dimension(:,:,:) :: gz
    real(r8), allocatable, dimension(:,:,:) :: t 
    real(r8), allocatable, dimension(:,:,:) :: pt
    real(r8), allocatable, dimension(:,:,:) :: rhod
  contains
    procedure :: init  => refer_state_init
    procedure :: clear => refer_state_clear
    final :: refer_state_final
  end type refer_state_type

contains

  subroutine refer_state_init(this, static)

    class(refer_state_type), intent(inout) :: this
    type(static_type), intent(in) :: static

    call this%clear()

    this%mesh => static%mesh

    call allocate_array(this%mesh, this%phs , full_lon=.true., full_lat=.true.)
    call allocate_array(this%mesh, this%ph  , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%gz  , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%t   , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%pt  , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%rhod, full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine refer_state_init

  subroutine refer_state_clear(this)

    class(refer_state_type), intent(inout) :: this

    if (allocated(this%phs)) deallocate(this%phs)
    if (allocated(this%ph )) deallocate(this%ph )
    if (allocated(this%gz )) deallocate(this%gz )
    if (allocated(this%t  )) deallocate(this%t  )
    if (allocated(this%pt )) deallocate(this%pt )

  end subroutine refer_state_clear

  subroutine refer_state_final(this)

    type(refer_state_type), intent(inout) :: this

    call this%clear()

  end subroutine refer_state_final

end module refer_state_types_mod
