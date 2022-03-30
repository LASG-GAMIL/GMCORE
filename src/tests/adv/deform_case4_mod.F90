module deform_case4_mod

  use const_mod
  use block_mod
  use transport_mod

  implicit none

  private

contains

  subroutine deform_case4_set_ic(block)

    type(block_type), intent(inout) :: block

  end subroutine deform_case4_set_ic

end module deform_case4_mod