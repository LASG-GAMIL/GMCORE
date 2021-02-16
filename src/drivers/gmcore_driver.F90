program gmcore_driver

  use log_mod
  use namelist_mod
  use block_mod
  use diag_state_mod
  use parallel_mod
  use initial_mod
  use restart_mod
  use gmcore_mod
  use steady_state_test_mod
  use rossby_haurwitz_wave_3d_test_mod
  use mountain_wave_test_mod
  use baroclinic_wave_test_mod
  use held_suarez_test_mod
  use steady_state_pgf_test_mod
  use ksp15_test_mod
  use dcmip31_test_mod

  implicit none

  character(256) namelist_path
  integer iblk

  interface
    subroutine set_ic_interface(block)
      import block_type
      type(block_type), intent(inout), target :: block
    end subroutine set_ic_interface
  end interface
  procedure(set_ic_interface), pointer :: set_ic

  call get_command_argument(1, namelist_path)
  if (namelist_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call parse_namelist(namelist_path)

  if (initial_file == 'N/A' .and. .not. restart) then
    select case (test_case)
    case ('pgf_test')
      call steady_state_pgf_test_set_params()
    case ('ksp15_01', 'ksp15_02')
      call ksp15_test_set_params()
    case ('dcmip31')
      call dcmip31_test_set_params()
      init_diag_state => dcmip31_test_init_diag_state
    end select
  end if

  call gmcore_init(namelist_path)

  if (initial_file /= 'N/A') then
    call initial_read()
  else if (restart) then
    call restart_read()
  else
    select case (test_case)
    case ('steady_state')
      set_ic => steady_state_test_set_ic
    case ('rossby_haurwitz_wave')
      set_ic => rossby_haurwitz_wave_3d_test_set_ic
    case ('mountain_wave')
      set_ic => mountain_wave_test_set_ic
    case ('baroclinic_wave')
      set_ic => baroclinic_wave_test_set_ic
    case ('held_suarez')
      set_ic => held_suarez_test_set_ic
    case ('pgf_test')
      set_ic => steady_state_pgf_test_set_ic
    case ('ksp15_01')
      set_ic => ksp15_01_test_set_ic
    case ('ksp15_02')
      set_ic => ksp15_02_test_set_ic
    case ('dcmip31')
      set_ic => dcmip31_test_set_ic
    case default
      call log_error('Unknown test case ' // trim(test_case) // '!')
    end select

    do iblk = 1, size(proc%blocks)
      call set_ic(proc%blocks(iblk))
    end do
  end if

  call gmcore_run()

  call gmcore_final()

end program gmcore_driver
