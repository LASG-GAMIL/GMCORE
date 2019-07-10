program gmcore_swm_driver

  use log_mod
  use namelist_mod
  use gmcore_mod
  use mountain_zonal_flow_test_mod
  use rossby_haurwitz_wave_test_mod
  use jet_zonal_flow_test_mod

  implicit none

  character(256) namelist_file_path

  call get_command_argument(1, namelist_file_path)
  if (namelist_file_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call parse_namelist(namelist_file_path)
  call gmcore_init()

  select case (test_case)
  case ('mountain_zonal_flow')
    call mountain_zonal_flow_test_set_initial_condition()
  case ('rossby_haurwitz_wave')
    call rossby_haurwitz_wave_test_set_initial_condition()
  case ('jet_zonal_flow')
    call jet_zonal_flow_test_set_initial_condition()
  case default
    call log_error('Unknown test case ' // trim(test_case) // '!')
  end select

  call gmcore_run()

  call gmcore_final()

end program gmcore_swm_driver
