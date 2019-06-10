module gmcore_mod

  use parallel_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use history_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod

  implicit none

  private

  public gmcore_init
  public gmcore_run
  public gmcore_final

contains

  subroutine gmcore_init()

    call parallel_init()
    call time_init()
    call create_meshes()
    call create_states()
    call create_static()
    call create_tends()
    call history_init()

  end subroutine gmcore_init

  subroutine gmcore_run()

    call output(state(old))

  end subroutine gmcore_run

  subroutine gmcore_final()

    call parallel_final()

  end subroutine gmcore_final

  subroutine output(state)

    type(state_type), intent(in) :: state

    if (time_is_alerted('history_write')) call history_write(state, static)

  end subroutine output

end module gmcore_mod
