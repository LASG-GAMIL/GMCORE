module damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use zonal_damp_mod
  use div_damp_mod
  use vor_damp_mod
  use smag_damp_mod
#ifdef USE_HZD
  use tridiag_hzd_mod
#endif

  implicit none

  private

  public damp_init
  public damp_final
  public damp_run
  public polar_damp_run
  public div_damp_run
  public vor_damp_run
  public smag_damp_run

contains

  subroutine damp_init()

#ifdef USE_HZD
    ! 初始化通信设置
    call zonal_comm_init(proc%comm, proc%id, proc%ngb(east)%id, proc%ngb(west)%id)
    call barrier()
#endif

    call zonal_damp_init()
    call div_damp_init()
    call vor_damp_init()

  end subroutine damp_init

  subroutine damp_final()

    call zonal_damp_final()
    call div_damp_final()
    call vor_damp_final()

  end subroutine damp_final

  subroutine damp_run(dt, new, blocks)

    real(8), intent(in) :: dt
    integer, intent(in) :: new
    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      if (use_div_damp) then
        call div_damp_run(blocks(iblk), dt, blocks(iblk)%state(new))
      end if
      if (use_vor_damp) then
        call vor_damp_run(blocks(iblk), dt, blocks(iblk)%state(new))
      end if
      if (use_polar_damp) then
        call polar_damp_run(blocks(iblk), dt, blocks(iblk)%state(new))
      end if
      if (use_smag_damp) then
        call smag_damp_run(blocks(iblk), dt, blocks(iblk)%tend(new), blocks(iblk)%state(new))
      end if
    end do

  end subroutine damp_run

  subroutine polar_damp_run(block, dt, state)

    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    integer cyc

    do cyc = 1, polar_damp_cycles
      call zonal_damp_on_cell(block, polar_damp_order, dt, state%pt)
      call zonal_damp_on_lon_edge(block, polar_damp_order, dt, state%u)
      call zonal_damp_on_lat_edge(block, polar_damp_order, dt, state%v)
      if (polar_damp_phs) call zonal_damp_on_cell(block, polar_damp_order, dt, state%phs, lat0=polar_damp_phs_lat0)
    end do
    call fill_halo(block, state%pt, full_lon=.true. , full_lat=.true. , full_lev=.true.)
    call fill_halo(block, state%u , full_lon=.false., full_lat=.true. , full_lev=.true.)
    call fill_halo(block, state%v , full_lon=.true. , full_lat=.false., full_lev=.true.)
    if (polar_damp_phs) call fill_halo(block, state%phs, full_lon=.true., full_lat=.true.)
    if (nonhydrostatic) then
      do cyc = 1, polar_damp_cycles
        call zonal_damp_on_lev_edge(block, polar_damp_order, dt, state%w_lev)
      end do
    end if

  end subroutine polar_damp_run

end module damp_mod
