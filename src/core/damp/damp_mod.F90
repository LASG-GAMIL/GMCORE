module damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use vor_damp_mod
  use div_damp_mod
  use smag_damp_mod
  use laplace_damp_mod
  use zonal_damp_mod
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
    call laplace_damp_init()
    call zonal_damp_init()

  end subroutine damp_init

  subroutine damp_final()

    call vor_damp_final()
    call div_damp_final()
    call laplace_damp_final()
    call zonal_damp_final()

  end subroutine damp_final

  subroutine damp_run(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    integer cyc

    if (use_div_damp) then
      do cyc = 1, div_damp_cycles
        call div_damp_run(block, dt, state)
      end do
    end if
    if (use_vor_damp) then
      do cyc = 1, div_damp_cycles
        call vor_damp_run(block, dt, state)
      end do
    end if
    if (use_smag_damp) then
      call smag_damp_run(block, dt, tend, state)
    end if
    if (use_zonal_damp) then
      do cyc = 1, zonal_damp_cycles
        call zonal_damp_on_lon_edge(block, zonal_damp_order, dt, state%u_lon)
        call zonal_damp_on_lat_edge(block, zonal_damp_order, dt, state%v_lat)
        if (baroclinic) then
          call zonal_damp_on_cell(block, zonal_damp_order, dt, state%phs)
          call zonal_damp_on_cell(block, zonal_damp_order, dt, state%pt)
        else
          call zonal_damp_on_cell(block, zonal_damp_order, dt, state%gz)
        end if
      end do
    end if

  end subroutine damp_run

end module damp_mod
