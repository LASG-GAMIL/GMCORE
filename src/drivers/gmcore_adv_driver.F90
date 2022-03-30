program gmcore_adv_driver

  use namelist_mod
  use transport_mod
  use gmcore_mod

  implicit none

  character(256) namelist_path

  call get_command_argument(1, namelist_path)
  if (namelist_path == '') then
    call log_error('You should give a namelist file path!')
  end if

  call parse_namelist(namelist_path)

  call gmcore_init(namelist_path)

  call transport_add_tracer('test1', 'q1', dt_dynamics)
  call transport_add_tracer('test2', 'q2', dt_dynamics)
  call transport_add_tracer('test1', 'q3', dt_dynamics)
  call transport_add_tracer('test3', 'q4', dt_dynamics)
  call transport_add_tracer('test2', 'q5', dt_dynamics)

  call transport_allocate_tracers(proc%blocks(1))

  call gmcore_final()

end program gmcore_adv_driver