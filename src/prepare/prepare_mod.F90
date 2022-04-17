module prepare_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use process_mod
  use block_mod
  use topo_mod
  use bkg_mod
  use operators_mod

  implicit none

contains

  subroutine prepare_static(block)

    class(block_type), intent(inout) :: block

    integer i, j

    associate (mesh    => block%mesh          , &
               gzs     => block%static%gzs    , & ! in
               dzsdlon => block%static%dzsdlon, & ! out
               dzsdlat => block%static%dzsdlat)   ! out
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          dzsdlon(i,j) = (gzs(i+1,j) - gzs(i,j)) / g / mesh%de_lon(j)
        end do
      end do
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dzsdlat(i,j) = (gzs(i,j+1) - gzs(i,j)) / g / mesh%de_lat(j)
        end do
      end do
      call fill_halo(block, dzsdlon, full_lon=.false., full_lat=.true.)
      call fill_halo(block, dzsdlat, full_lon=.true., full_lat=.false.)
    end associate

  end subroutine prepare_static

  subroutine prepare_run()

    integer iblk, i

    call topo_read(topo_file)
    do iblk = 1, size(blocks)
      call topo_regrid(blocks(iblk))
    end do
    if (use_topo_smooth) then
      do iblk = 1, size(blocks)
        call topo_smooth(blocks(iblk))
      end do
    end if

    call bkg_read(bkg_type, bkg_file)

    call bkg_regrid_phs()
    call bkg_calc_ph()
    call bkg_regrid_pt()
    call bkg_regrid_u()
    call bkg_regrid_v()

    if (nonhydrostatic) then
      do iblk = 1, size(blocks)
        call diag_ph    (blocks(iblk), blocks(iblk)%state(1))
        call diag_t     (blocks(iblk), blocks(iblk)%state(1))
        call diag_gz_lev(blocks(iblk), blocks(iblk)%state(1))
      end do
    end if

  end subroutine prepare_run

  subroutine prepare_final()

    call topo_final()
    call bkg_final()

  end subroutine prepare_final

end module prepare_mod
