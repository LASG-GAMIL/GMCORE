module static_mod

  use const_mod
  use mesh_mod
  use allocator_mod

  implicit none

  private

  public static_type

  type static_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable :: gzs(:,:)
  contains
    procedure :: init => static_init
    procedure :: clear => static_clear
    final :: static_final
  end type static_type

contains

  subroutine static_init(this, mesh)

    class(static_type), intent(inout)         :: this
    type(mesh_type   ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%gzs, full_lon=.true., full_lat=.true.)

  end subroutine static_init

  subroutine static_clear(this)

    class(static_type), intent(inout) :: this

    if (allocated(this%gzs)) deallocate(this%gzs)

  end subroutine static_clear

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

end module static_mod