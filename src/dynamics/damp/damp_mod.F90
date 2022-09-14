module damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use filter_mod
  use vor_damp_mod
  use div_damp_mod
  use smag_damp_mod
  use laplace_damp_mod
  use lon_damp_mod
  use lat_damp_mod
  use operators_mod

  implicit none

  private

  public damp_init
  public damp_final
  public damp_run

contains

  subroutine damp_init(blocks)

    type(block_type), intent(in) :: blocks(:)

    call vor_damp_init(blocks)
    call div_damp_init(blocks)
    call smag_damp_init()
    call laplace_damp_init()
    call lon_damp_init()
    call lat_damp_init()

  end subroutine damp_init

  subroutine damp_final()

    call vor_damp_final()
    call div_damp_final()
    call smag_damp_final()
    call laplace_damp_final()
    call lon_damp_final()
    call lat_damp_final()

  end subroutine damp_final

  subroutine damp_run(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    if (use_div_damp) then
      call div_damp_run(block, state)
    end if
    if (use_vor_damp) then
      call vor_damp_run(block, dt, state)
    end if
    if (use_smag_damp) then
      call smag_damp_run(block, dt, tend, state)
    end if

  end subroutine damp_run

end module damp_mod
