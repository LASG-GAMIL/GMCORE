module state_mod

  use mesh_mod
  use allocator_mod

  implicit none

  private

  public state_type
  public state
  public create_states

  type state_type
    type(mesh_type), pointer :: mesh => null()
    real, allocatable :: u (:,:)
    real, allocatable :: v (:,:)
    real, allocatable :: gd(:,:)
    real, allocatable :: pv(:,:)
  contains
    procedure :: init => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

  type(state_type), allocatable :: state(:)

contains

  subroutine create_states()

    integer i

    if (.not. allocated(state)) then
      allocate(state(0:2))
      do i = lbound(state, 1), ubound(state, 1)
        call state(i)%init(mesh)
      end do
    end if

  end subroutine create_states

  subroutine state_init(this, mesh)

    class(state_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh

    call this%clear()

    this%mesh => mesh

#ifdef STAGGER_V_ON_POLE
    call allocate_array(mesh, this%u , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%v , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%gd, full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%pv, half_lon=.true., full_lat=.true.)
#else
    call allocate_array(mesh, this%u , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%v , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%gd, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv, half_lon=.true., half_lat=.true.)
#endif

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%u )) deallocate(this%u )
    if (allocated(this%v )) deallocate(this%v )
    if (allocated(this%gd)) deallocate(this%gd)
    if (allocated(this%pv)) deallocate(this%pv)

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
