module nest_mod

  use const_mod
  use namelist_mod
  use block_mod

  implicit none

  private

  public nest_init
  public nest_downward
  public nest_upward
  public nest_max_dom

contains

  subroutine nest_init()

  end subroutine nest_init

  subroutine nest_add_domain()

  end subroutine nest_add_domain

  subroutine nest_downward()

  end subroutine nest_downward

  subroutine nest_coarse_to_fine()

  end subroutine nest_coarse_to_fine

  subroutine nest_upward()

  end subroutine nest_upward

  subroutine nest_fine_to_coarse()

  end subroutine nest_fine_to_coarse

end module nest_mod
