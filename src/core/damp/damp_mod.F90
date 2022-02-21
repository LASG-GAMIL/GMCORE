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
  use filter_mod

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

  end subroutine damp_init

  subroutine damp_final()

    call vor_damp_final()
    call div_damp_final()
    call laplace_damp_final()

  end subroutine damp_final

  subroutine damp_run(dt, new, blocks)

    real(8), intent(in) :: dt
    integer, intent(in) :: new
    type(block_type), intent(inout) :: blocks(:)

    integer iblk, cyc

    do iblk = 1, size(blocks)
      if (use_div_damp) then
        do cyc = 1, div_damp_cycles
          call div_damp_run(blocks(iblk), dt, blocks(iblk)%state(new))
        end do
      end if
      if (use_vor_damp) then
        do cyc = 1, div_damp_cycles
          call vor_damp_run(blocks(iblk), dt, blocks(iblk)%state(new))
        end do
      end if
      if (use_smag_damp) then
        call smag_damp_run(blocks(iblk), dt, blocks(iblk)%tend(new), blocks(iblk)%state(new))
      end if
      if (use_div_damp .or. use_vor_damp) then
        associate (block => blocks(iblk)               , &
                   u     => blocks(iblk)%state(new)%u  , &
                   v     => blocks(iblk)%state(new)%v  , &
                   u_f   => blocks(iblk)%state(new)%u_f, &
                   v_f   => blocks(iblk)%state(new)%v_f)
        call filter_on_lon_edge(block, u, u_f)
        call fill_halo(block, u_f, full_lon=.false., full_lat=.true., full_lev=.true.)
        call filter_on_lat_edge(block, v, v_f)
        call fill_halo(block, v_f, full_lon=.true., full_lat=.false., full_lev=.true.)
        end associate
      end if
    end do

  end subroutine damp_run

end module damp_mod
