program gmcore_swm_driver

  use log_mod
  use namelist_mod
  use block_mod
  use process_mod
  use gmcore_mod
  use mountain_zonal_flow_test_mod
  use rossby_haurwitz_wave_test_mod
  use jet_zonal_flow_test_mod
  use steady_geostrophic_flow_test_mod
  use cross_pole_flow_test_mod

  implicit none

  character(256) namelist_file_path
  integer iblk

  interface
    subroutine set_initial_condition_interface(block)
      import block_type
      type(block_type), intent(inout), target :: block
    end subroutine
  end interface
  procedure(set_initial_condition_interface), pointer :: set_initial_condition

  call get_command_argument(1, namelist_file_path)
  if (namelist_file_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call parse_namelist(namelist_file_path)
  call gmcore_init()

  select case (test_case)
  case ('mountain_zonal_flow')
    set_initial_condition => mountain_zonal_flow_test_set_initial_condition
  case ('rossby_haurwitz_wave')
    set_initial_condition => rossby_haurwitz_wave_test_set_initial_condition
  case ('jet_zonal_flow')
    set_initial_condition => jet_zonal_flow_test_set_initial_condition
  case ('steady_geostrophic_flow')
    set_initial_condition => steady_geostrophic_flow_test_set_initial_condition
  case ('cross_pole_flow')
    set_initial_condition => cross_pole_flow_test_set_initial_condition
  case default
    call log_error('Unknown test case ' // trim(test_case) // '!')
  end select

  do iblk = 1, size(process%blocks)
    call set_initial_condition(process%blocks(iblk))
  end do

  call gmcore_run()

  call gmcore_final()

end program gmcore_swm_driver
