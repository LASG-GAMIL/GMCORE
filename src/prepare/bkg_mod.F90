module bkg_mod

  use fiona
  use flogger
  use const_mod
  use block_mod
  use process_mod
  use era5_reader_mod

  implicit none

  private

  public bkg_read
  public bkg_final

contains

  subroutine bkg_read(bkg_type, bkg_file)

    character(*), intent(in) :: bkg_type
    character(*), intent(in) :: bkg_file

    select case (bkg_type)
    case ('era5')
      call era5_reader_run(bkg_file)
    case default
      if (is_root_proc()) call log_error('Unknown bkg_type ' // trim(bkg_type) // '!')
    end select

  end subroutine bkg_read

  subroutine bkg_final()

    call era5_reader_final()

  end subroutine bkg_final

end module bkg_mod
