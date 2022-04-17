program gmcore_swm_driver

  use flogger
  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
  use restart_mod
  use gmcore_mod
  use mountain_zonal_flow_test_mod
  use rossby_haurwitz_wave_test_mod
  use jet_zonal_flow_test_mod
  use steady_geostrophic_flow_test_mod
  use cross_pole_flow_test_mod
  use shallow_water_waves_test_mod
  use vortex_erosion_test_mod

  implicit none

  character(256) namelist_path
  integer iblk

  interface
    subroutine set_ic_interface(block)
      import block_type
      type(block_type), intent(inout), target :: block
    end subroutine
  end interface
  procedure(set_ic_interface), pointer :: set_ic

  call get_command_argument(1, namelist_path)
  if (namelist_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call parse_namelist(namelist_path)

  call const_init(planet)

  call gmcore_init(namelist_path)

  if (restart) then
    call restart_read()
  else
    select case (test_case)
    case ('mountain_zonal_flow')
      set_ic => mountain_zonal_flow_test_set_ic
    case ('rossby_haurwitz_wave')
      set_ic => rossby_haurwitz_wave_test_set_ic
    case ('jet_zonal_flow')
      set_ic => jet_zonal_flow_test_set_ic
    case ('steady_geostrophic_flow')
      set_ic => steady_geostrophic_flow_test_set_ic
    case ('cross_pole_flow')
      set_ic => cross_pole_flow_test_set_ic
    case ('shallow_water_waves')
      set_ic => shallow_water_waves_test_set_ic
    case ('vortex_erosion')
      set_ic => vortex_erosion_test_set_ic
    case default
      call log_error('Unknown test case ' // trim(test_case) // '!', pid=proc%id)
    end select

    do iblk = 1, size(blocks)
      call set_ic(blocks(iblk))
    end do
  end if

  call gmcore_run()

  call gmcore_final()

end program gmcore_swm_driver
