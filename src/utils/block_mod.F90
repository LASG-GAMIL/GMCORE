module block_mod

  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod

  implicit none

  type block_type
    type(mesh_type), pointer :: mesh => null()
    type(state_type), pointer :: state => null()
    type(static_type), pointer :: static => null()
    type(tend_type), pointer :: tend => null()
  contains
    procedure :: init => block_init
    final :: block_final
  end type block_type

contains

  subroutine block_init(this)

    class(block_type), intent(inout) :: this

  end subroutine block_init

  subroutine block_final(this)

    type(block_type), intent(inout) :: this

  end subroutine block_final

end module block_mod