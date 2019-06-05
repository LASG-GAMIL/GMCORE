module gmcore_mod

  use mesh_mod
  use state_mod

  implicit none

  private

  public gmcore_init

contains

  subroutine gmcore_init()

    call create_meshes()
    call create_states()

  end subroutine

end module gmcore_mod
