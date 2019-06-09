module gmcore_mod

  use parallel_mod
  use mesh_mod
  use state_mod
  use tend_mod

  implicit none

  private

  public gmcore_init
  public gmcore_run
  public gmcore_final

contains

  subroutine gmcore_init()

    call parallel_init()
    call create_meshes()
    call create_states()
    call create_tends()

  end subroutine gmcore_init

  subroutine gmcore_run()


  end subroutine gmcore_run

  subroutine gmcore_final()

    call parallel_final()

  end subroutine gmcore_final

end module gmcore_mod
