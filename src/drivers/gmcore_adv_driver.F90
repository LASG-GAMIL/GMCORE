program gmcore_adv_driver

  use namelist_mod
  use gmcore_mod
  use history_mod
  use process_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use operators_mod, only: calc_qmf
  use time_schemes_mod, only: time_integrator
  use deform_case4_mod

  implicit none

  character(256) namelist_path
  integer iblk

  interface
    subroutine set_ic_interface(block)
      import block_type
      type(block_type), intent(inout) :: block
    end subroutine set_ic_interface
    subroutine set_uv_interface(block, state, time_in_seconds)
      import block_type, state_type
      type(block_type), intent(in   ) :: block
      type(state_type), intent(inout) :: state
      real(8), intent(in) :: time_in_seconds
    end subroutine set_uv_interface
  end interface
  procedure(set_ic_interface), pointer :: set_ic
  procedure(set_uv_interface), pointer :: set_uv

  call get_command_argument(1, namelist_path)
  if (namelist_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call parse_namelist(namelist_path)

  call gmcore_init(namelist_path)

  select case (test_case)
  case ('deform_case4')
    set_ic => deform_case4_set_ic
    set_uv => deform_case4_set_uv
  case default
    call log_error('Unknown test case ' // trim(test_case) // '!', pid=proc%id)
  end select

  do iblk = 1, size(proc%blocks)
    call set_ic(proc%blocks(iblk))
    call set_uv(proc%blocks(iblk), proc%blocks(iblk)%state(old), elapsed_seconds)
    call adv_accum_wind(proc%blocks(iblk), old)
  end do

  call history_setup_h0_adv(proc%blocks)
  call history_write_h0(proc%blocks, old)

  do while (.not. time_is_finished())
    do iblk = 1, size(proc%blocks)
      call set_uv(proc%blocks(iblk), proc%blocks(iblk)%state(new), elapsed_seconds + dt_adv)
      call adv_accum_wind(proc%blocks(iblk), new)
    end do
    do iblk = 1, size(proc%blocks)
      call time_integrator(adv_operator, proc%blocks(iblk), old, new, dt_adv)
    end do
    call time_advance(dt_adv)
    if (is_root_proc() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
    call history_write_h0(proc%blocks, old)
  end do

  call gmcore_final()

contains

  subroutine adv_operator(block, old_state, star_state, new_state, tend1, tend2, dt, pass)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    type(tend_type ), intent(inout) :: tend1
    type(tend_type ), intent(in   ) :: tend2
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    call calc_qmf(block, star_state)

  end subroutine adv_operator

end program gmcore_adv_driver
