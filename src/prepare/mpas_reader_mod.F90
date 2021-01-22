module mpas_reader_mod

  use fiona
  use flogger
  use const_mod
  use formula_mod
  use block_mod
  use process_mod

  implicit none

  ! This is not the raw MPAS data.
  integer num_mpas_lon
  integer num_mpas_lat
  integer num_mpas_lev

  real(r8), allocatable, dimension(:    ) :: mpas_lon
  real(r8), allocatable, dimension(:    ) :: mpas_lat
  real(r8), allocatable, dimension(:,:,:) :: mpas_u
  real(r8), allocatable, dimension(:,:,:) :: mpas_v
  real(r8), allocatable, dimension(:,:,:) :: mpas_pt
  real(r8), allocatable, dimension(:,:,:) :: mpas_t
  real(r8), allocatable, dimension(:,:,:) :: mpas_p
  real(r8), allocatable, dimension(:,:  ) :: mpas_ps
  real(r8), allocatable, dimension(:,:  ) :: mpas_zs
  real(r8), allocatable, dimension(:,:  ) :: mpas_zbot

contains

  subroutine mpas_reader_run(bkg_file)

    character(*), intent(in) :: bkg_file

    real(r8), allocatable :: tmp(:,:,:)
    integer k

    call mpas_reader_final()

    if (is_root_proc()) call log_notice('Use MPAS ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('mpas', file_path=bkg_file)
    call fiona_get_dim('mpas', 'lon', size=num_mpas_lon)
    call fiona_get_dim('mpas', 'lat', size=num_mpas_lat)
    call fiona_get_dim('mpas', 'lev', size=num_mpas_lev)

    allocate(tmp      (num_mpas_lon,num_mpas_lat,num_mpas_lev))
    allocate(mpas_lon (num_mpas_lon))
    allocate(mpas_lat (num_mpas_lat))
    allocate(mpas_u   (num_mpas_lon,num_mpas_lat,num_mpas_lev))
    allocate(mpas_v   (num_mpas_lon,num_mpas_lat,num_mpas_lev))
    allocate(mpas_pt  (num_mpas_lon,num_mpas_lat,num_mpas_lev))
    allocate(mpas_t   (num_mpas_lon,num_mpas_lat,num_mpas_lev))
    allocate(mpas_p   (num_mpas_lon,num_mpas_lat,num_mpas_lev))
    allocate(mpas_ps  (num_mpas_lon,num_mpas_lat))
    allocate(mpas_zs  (num_mpas_lon,num_mpas_lat))
    allocate(mpas_zbot(num_mpas_lon,num_mpas_lat))

    call fiona_start_input('mpas')
    call fiona_input('mpas', 'lon', mpas_lon )
    call fiona_input('mpas', 'lat', mpas_lat )
    call fiona_input('mpas', 'u'  , mpas_u   ); tmp = mpas_u (:,:,num_mpas_lev:1:-1); mpas_u  = tmp
    call fiona_input('mpas', 'v'  , mpas_v   ); tmp = mpas_v (:,:,num_mpas_lev:1:-1); mpas_v  = tmp
    call fiona_input('mpas', 'pt' , mpas_pt  ); tmp = mpas_pt(:,:,num_mpas_lev:1:-1); mpas_pt = tmp
    call fiona_input('mpas', 'ph' , mpas_p   ); tmp = mpas_p (:,:,num_mpas_lev:1:-1); mpas_p  = tmp
    call fiona_input('mpas', 'phs', mpas_ps  )
    call fiona_input('mpas', 'ter', mpas_zs  )
    call fiona_input('mpas', 'zs' , mpas_zbot)
    call fiona_end_input('mpas')

    mpas_t = temperature(mpas_pt, mpas_p)

    deallocate(tmp)

  end subroutine mpas_reader_run

  subroutine mpas_reader_final()

    if (allocated(mpas_lon )) deallocate(mpas_lon )
    if (allocated(mpas_lat )) deallocate(mpas_lat )
    if (allocated(mpas_u   )) deallocate(mpas_u   )
    if (allocated(mpas_v   )) deallocate(mpas_v   )
    if (allocated(mpas_pt  )) deallocate(mpas_pt  )
    if (allocated(mpas_p   )) deallocate(mpas_p   )
    if (allocated(mpas_ps  )) deallocate(mpas_ps  )
    if (allocated(mpas_zs  )) deallocate(mpas_zs  )
    if (allocated(mpas_zbot)) deallocate(mpas_zbot)

  end subroutine mpas_reader_final

end module mpas_reader_mod
