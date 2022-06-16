program gmcore_adv_driver

  use flogger
  use namelist_mod
  use gmcore_mod
  use history_mod
  use parallel_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use solid_rotation_test_mod
  use deform_test_mod
  use moving_vortices_test_mod
  use dcmip12_test_mod

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
  case ('solid_rotation')
    call solid_rotation_test_init()
    set_ic => solid_rotation_test_set_ic
    set_uv => solid_rotation_test_set_uv
  case ('deform_case1')
    call deform_test_init(1)
    set_ic => deform_test_set_ic
    set_uv => deform_case1_test_set_uv
  case ('deform_case2')
    call deform_test_init(2)
    set_ic => deform_test_set_ic
    set_uv => deform_case2_test_set_uv
  case ('deform_case3')
    call deform_test_init(3)
    set_ic => deform_test_set_ic
    set_uv => deform_case3_test_set_uv
  case ('deform_case4')
    call deform_test_init(4)
    set_ic => deform_test_set_ic
    set_uv => deform_case4_test_set_uv
  case ('moving_vortices')
    call moving_vortices_test_init()
    set_ic => moving_vortices_test_set_ic
    set_uv => moving_vortices_test_set_uv
  case ('dcmip12')
    call dcmip12_test_init()
    set_ic => dcmip12_test_set_ic
    set_uv => dcmip12_test_set_uv
  case default
    call log_error('Unknown test case ' // trim(test_case) // '!', pid=proc%id)
  end select

  do iblk = 1, size(blocks)
    call set_uv(blocks(iblk), blocks(iblk)%state(old), elapsed_seconds)
    call set_ic(blocks(iblk))
  end do

  call history_setup_h0_adv(blocks)
  call output(old)

  do while (.not. time_is_finished())
    do iblk = 1, size(blocks)
      call set_uv(blocks(iblk), blocks(iblk)%state(new), elapsed_seconds + dt_adv)
      call adv_run(blocks(iblk), new)
    end do
    call diagnose(blocks)
    if (is_root_proc() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
    call time_advance(dt_adv)
    call output(old)
  end do

  call gmcore_final()

contains

  subroutine diagnose(blocks)

    type(block_type), intent(inout) :: blocks(:)

    integer i, j, k, l, iblk
    real(r8) qm

    qm = 0
    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh              , &
                 old  => blocks(iblk)%adv_batches(1)%old, &
                 q    => blocks(iblk)%adv_batches(1)%q  )
      do l = 1, size(blocks(iblk)%adv_batches(1)%tracer_names)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              qm = qm + q(i,j,k,l,old) * mesh%area_cell(j)
            end do
          end do
        end do
        call global_sum(proc%comm, qm)
        call log_add_diag('qm' // to_str(l), qm)
      end do
      end associate
    end do

  end subroutine diagnose

  subroutine output(itime)

    integer, intent(in) :: itime

    real(8), save :: time1 = 0, time2
    integer i, j, k, iblk

    if (time_is_alerted('history_write')) then
      if (time_step == 0) call cpu_time(time1)
      call cpu_time(time2)
      if (time_step /= 0) then
        if (is_root_proc()) call log_notice('Time cost ' // to_str(time2 - time1, 5) // ' seconds.')
        time1 = time2
      end if
      call history_write_h0(blocks, itime)
    end if

  end subroutine output

end program gmcore_adv_driver
