module damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use zonal_damp_mod
  use meridional_damp_mod
  use div_damp_mod
  use vor_damp_mod

  implicit none

  private

  public damp_init
  public damp_final
  public polar_damp
  public div_damp
  public vor_damp

contains

  subroutine damp_init()

    integer j, jr, k, r

    call zonal_damp_init()
    call meridional_damp_init()
    call div_damp_init()
    call vor_damp_init()

  end subroutine damp_init

  subroutine damp_final()

    call zonal_damp_final()
    call meridional_damp_final()
    call div_damp_final()
    call vor_damp_final()

  end subroutine damp_final

  subroutine polar_damp(block, dt, state)

    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    call zonal_damp_on_lon_edge(block, polar_damp_order, dt, state%u)
    call zonal_damp_on_lat_edge(block, polar_damp_order, dt, state%v)
    if (baroclinic) then
      call zonal_damp_on_cell(block, polar_damp_order, dt, state%pt)
    end if

  end subroutine polar_damp

end module damp_mod