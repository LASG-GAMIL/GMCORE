module static_mod

  use mesh_mod
  use allocator_mod

  implicit none

  type static_type
    type(mesh_type), pointer :: mesh => null()
    real, allocatable :: ghs(:,:)
  contains
    procedure :: init => static_init
    procedure :: clear => static_clear
    final :: static_final
  end type static_type

  type(static_type) static

contains

  subroutine create_static()

    call static%init(mesh)

  end subroutine create_static

  subroutine static_init(this, mesh)

    class(static_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh

    call this%clear()

    this%mesh => mesh

#ifdef STAGGER_V_ON_POLE
    call allocate_array(mesh, this%ghs, full_lon=.true., half_lat=.true.)
#else
    call allocate_array(mesh, this%ghs, full_lon=.true., full_lat=.true.)
#endif

  end subroutine static_init

  subroutine static_clear(this)

    class(static_type), intent(inout) :: this

    if (allocated(this%ghs)) deallocate(this%ghs)

  end subroutine static_clear

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

end module static_mod