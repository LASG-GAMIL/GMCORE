module moist_mod

  use namelist_mod
  use adv_mod

  implicit none

  private

  public moist_init

contains

  subroutine moist_init()

    call adv_add_tracer('moist', dt_adv, 'qv', 'water vapor')

  end subroutine moist_init

end module moist_mod
