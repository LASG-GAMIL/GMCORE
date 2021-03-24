module prepare_mod

  use const_mod
  use topo_mod
  use bkg_mod
  use process_mod

contains

  subroutine prepare_run(topo_file, bkg_file, bkg_type, &
                         zero_min_lon, zero_max_lon, zero_min_lat, zero_max_lat, &
                         smth_min_lon, smth_max_lon, smth_min_lat, smth_max_lat, &
                         smth_steps)

    character(*), intent(in) :: topo_file
    character(*), intent(in) :: bkg_file
    character(*), intent(in) :: bkg_type
    real(r8), intent(in), optional :: zero_min_lon(:)
    real(r8), intent(in), optional :: zero_max_lon(:)
    real(r8), intent(in), optional :: zero_min_lat(:)
    real(r8), intent(in), optional :: zero_max_lat(:)
    real(r8), intent(in), optional :: smth_min_lon(:)
    real(r8), intent(in), optional :: smth_max_lon(:)
    real(r8), intent(in), optional :: smth_min_lat(:)
    real(r8), intent(in), optional :: smth_max_lat(:)
    integer , intent(in), optional :: smth_steps(:)

    integer iblk

    call topo_read(topo_file)
    do iblk = 1, size(proc%blocks)
      call topo_regrid(proc%blocks(iblk))
    end do
  
    if (present(zero_min_lon) .and. present(zero_max_lon) .and. present(zero_min_lat) .and. present(zero_max_lat)) then
      do iblk = 1, size(proc%blocks)
        do i = 1, size(zero_min_lon)
          if (zero_min_lon(i) /= -1.0e33 .and. zero_max_lon(i) /= -1.0e33 .and. &
              zero_min_lat(i) /= -1.0e33 .and. zero_max_lat(i) /= -1.0e33) then
            call topo_zero(proc%blocks(iblk), zero_min_lon(i), zero_max_lon(i), zero_min_lat(i), zero_max_lat(i))
          end if
        end do
      end do
    end if
    if (present(smth_min_lon) .and. present(smth_max_lon) .and. present(smth_min_lat) .and. present(smth_max_lat) .and. present(smth_steps)) then
      do iblk = 1, size(proc%blocks)
        do i = 1, size(smth_min_lon)
          if (smth_min_lon(i) /= -1.0e33 .and. smth_max_lon(i) /= -1.0e33 .and. &
              smth_min_lat(i) /= -1.0e33 .and. smth_max_lat(i) /= -1.0e33) then
            call topo_smth(proc%blocks(iblk), smth_min_lon(i), smth_max_lon(i), smth_min_lat(i), smth_max_lat(i), smth_steps(i))
          end if
        end do
      end do
    end if
  
    call bkg_read(bkg_type, bkg_file)
  
    call bkg_regrid_phs()
    call bkg_calc_ph()
    call bkg_regrid_pt()
    call bkg_regrid_u()
    call bkg_regrid_v()

  end subroutine prepare_run

  subroutine prepare_final()

    call topo_final()
    call bkg_final()

  end subroutine prepare_final

end module prepare_mod
