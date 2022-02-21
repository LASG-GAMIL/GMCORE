module waccm_reader_mod

  use fiona
  use flogger
  use const_mod
  use block_mod
  use process_mod

  implicit none

  integer waccm_nlon
  integer waccm_nlat
  integer waccm_nlev

  real(r8) waccm_p0
  real(r8), allocatable, dimension(:    ) :: waccm_lon
  real(r8), allocatable, dimension(:    ) :: waccm_lat
  real(r8), allocatable, dimension(:    ) :: waccm_hyam
  real(r8), allocatable, dimension(:    ) :: waccm_hybm
  real(r8), allocatable, dimension(:,:,:) :: waccm_u
  real(r8), allocatable, dimension(:,:,:) :: waccm_v
  real(r8), allocatable, dimension(:,:,:) :: waccm_t
  real(r8), allocatable, dimension(:,:,:) :: waccm_p
  real(r8), allocatable, dimension(:,:  ) :: waccm_ps
  real(r8), allocatable, dimension(:,:  ) :: waccm_zs

contains

  subroutine waccm_reader_run(bkg_file)

    character(*), intent(in) :: bkg_file

    integer k
    real(r8), allocatable :: tmp(:,:)

    call waccm_reader_final()

    if (is_root_proc()) call log_notice('Use WACCM ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('waccm', file_path=bkg_file)
    call fiona_get_dim('waccm', 'lon', size=waccm_nlon)
    call fiona_get_dim('waccm', 'lat', size=waccm_nlat)
    call fiona_get_dim('waccm', 'lev', size=waccm_nlev)

    allocate(waccm_lon (waccm_nlon))
    allocate(waccm_lat (waccm_nlat))
    allocate(waccm_hyam(waccm_nlev))
    allocate(waccm_hybm(waccm_nlev))
    allocate(waccm_u   (waccm_nlon,waccm_nlat,waccm_nlev))
    allocate(waccm_v   (waccm_nlon,waccm_nlat,waccm_nlev))
    allocate(waccm_t   (waccm_nlon,waccm_nlat,waccm_nlev))
    allocate(waccm_p   (waccm_nlon,waccm_nlat,waccm_nlev))
    allocate(waccm_ps  (waccm_nlon,waccm_nlat           ))
    allocate(waccm_zs  (waccm_nlon,waccm_nlat           ))

    call fiona_start_input('waccm')
    call fiona_input('waccm', 'lon' , waccm_lon )
    call fiona_input('waccm', 'lat' , waccm_lat )
    call fiona_input('waccm', 'P0'  , waccm_p0  )
    call fiona_input('waccm', 'hyam', waccm_hyam)
    call fiona_input('waccm', 'hybm', waccm_hybm)
    call fiona_input('waccm', 'U'   , waccm_u   )
    call fiona_input('waccm', 'V'   , waccm_v   )
    call fiona_input('waccm', 'T'   , waccm_t   )
    call fiona_input('waccm', 'PS'  , waccm_ps  )
    call fiona_input('waccm', 'PHIS', waccm_zs  )
    call fiona_end_input('waccm')

    waccm_zs = waccm_zs / g

    do k = 1, waccm_nlev
      waccm_p(:,:,k) = waccm_hyam(k) * waccm_p0 + waccm_hybm(k) * waccm_ps(:,:)
    end do

  end subroutine waccm_reader_run

  subroutine waccm_reader_final()

    if (allocated(waccm_lon )) deallocate(waccm_lon )
    if (allocated(waccm_lat )) deallocate(waccm_lat )
    if (allocated(waccm_hyam)) deallocate(waccm_hyam)
    if (allocated(waccm_hybm)) deallocate(waccm_hybm)
    if (allocated(waccm_u   )) deallocate(waccm_u   )
    if (allocated(waccm_v   )) deallocate(waccm_v   )
    if (allocated(waccm_t   )) deallocate(waccm_t   )
    if (allocated(waccm_p   )) deallocate(waccm_p   )
    if (allocated(waccm_ps  )) deallocate(waccm_ps  )
    if (allocated(waccm_zs  )) deallocate(waccm_zs  )

  end subroutine waccm_reader_final

end module waccm_reader_mod
