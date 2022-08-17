module ref_mod

  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod

  implicit none

  private

  public ref_calc_ps

contains

  subroutine ref_calc_ps()

    real(r8) p0, T0, A
    integer i, j, iblk

    select case (planet)
    case ('earth')
      do iblk = 1, size(blocks)
        associate (block  => blocks(iblk)              , &
                   mesh   => blocks(iblk)%mesh         , &
                   gzs    => blocks(iblk)%static%gzs   , &
                   ref_ps => blocks(iblk)%static%ref_ps)
        p0 = 1.0e5_r8
        T0 = 300.0_r8
        A  = 50.0_r8
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ref_ps(i,j) = p0 * exp(-T0 / A + sqrt((T0 / A)**2 - 2 * gzs(i,j) / A / Rd))
          end do
        end do
        call fill_halo(block, ref_ps, full_lon=.true., full_lat=.true.)
        end associate
      end do
    end select

  end subroutine ref_calc_ps

end module ref_mod
