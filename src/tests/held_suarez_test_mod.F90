module held_suarez_test_mod

  use flogger
  use const_mod
  use block_mod
  use parallel_mod

  implicit none

contains

  subroutine held_suarez_test_set_initial_condition(block)

    type(block_type), intent(inout), target :: block

    integer i, j
    type(mesh_type), pointer :: mesh

    mesh => block%mesh

  end subroutine held_suarez_test_set_initial_condition

end module held_suarez_test_mod
