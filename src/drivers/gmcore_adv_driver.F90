program gmcore_adv_driver

  use namelist_mod
  use gmcore_mod
  use deform_case4_mod

  implicit none

  character(256) namelist_path
  integer iblk

  interface
    subroutine set_ic_interface(block)
      import block_type
      type(block_type), intent(inout) :: block
    end subroutine set_ic_interface
  end interface
  procedure(set_ic_interface), pointer :: set_ic

  call get_command_argument(1, namelist_path)
  if (namelist_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call parse_namelist(namelist_path)

  call gmcore_init(namelist_path)

  select case (test_case)
  case ('deform_case4')
    set_ic => deform_case4_set_ic
  case default
    call log_error('Unknown test case ' // trim(test_case) // '!', pid=proc%id)
  end select

  do iblk = 1, size(proc%blocks)
    call set_ic(proc%blocks(iblk))
  end do

  call gmcore_final()

end program gmcore_adv_driver