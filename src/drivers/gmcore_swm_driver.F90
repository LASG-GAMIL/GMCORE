program gmcore_swm_driver

  use log_mod
  use namelist_mod
  use gmcore_mod
  use mountain_zonal_flow_test_mod

  implicit none

  call parse_namelist('namelist.gmcore_swm')
  call gmcore_init()

  select case (test_case)
  case ('mountain_zonal_flow')
    call mountain_zonal_flow_test_set_initial_condition()
  case default
    call log_error('Unknown test case ' // trim(test_case) // '!')
  end select

  call gmcore_run()

  call gmcore_final()

end program gmcore_swm_driver
