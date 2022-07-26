module moist_mod

  use namelist_mod
  use adv_mod

  implicit none

  private

  public moist_init

contains

  subroutine moist_init()

    call adv_add_tracer('moist', dt_adv, 'qv', 'water vapor' )
    ! call adv_add_tracer('moist', dt_adv, 'ql', 'cloud liquid')
    ! call adv_add_tracer('moist', dt_adv, 'qi', 'cloud ice'   )
    ! call adv_add_tracer('moist', dt_adv, 'qr', 'rain'        )
    ! call adv_add_tracer('moist', dt_adv, 'qs', 'snow'        )
    ! call adv_add_tracer('moist', dt_adv, 'qg', 'graupels'    )
    ! call adv_add_tracer('moist', dt_adv, 'qh', 'hail'        )

  end subroutine moist_init

end module moist_mod
