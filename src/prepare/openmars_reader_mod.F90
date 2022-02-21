module openmars_reader_mod

  use fiona
  use flogger
  use const_mod
  use block_mod
  use process_mod

  implicit none

  integer nlon_openmars
  integer nlat_openmars
  integer nlev_openmars

  real(r8), allocatable, dimension(:    ) :: openmars_lon
  real(r8), allocatable, dimension(:    ) :: openmars_lat
  real(r8), allocatable, dimension(:    ) :: openmars_lev
  real(r8), allocatable, dimension(:,:,:) :: openmars_u
  real(r8), allocatable, dimension(:,:,:) :: openmars_v
  real(r8), allocatable, dimension(:,:,:) :: openmars_t
  real(r8), allocatable, dimension(:,:,:) :: openmars_p
  real(r8), allocatable, dimension(:,:  ) :: openmars_ps

contains

  subroutine openmars_reader_run(bkg_file)

    character(*), intent(in) :: bkg_file

    integer k
    real(r8), allocatable :: tmp(:,:,:)

    call openmars_reader_final()

    if (is_root_proc()) call log_notice('Use OpenMARS ' // trim(bkg_file) // ' as background.')

    call fiona_open_dataset('openmars', file_path=bkg_file)
    call fiona_get_dim('openmars', 'lon', size=nlon_openmars)
    call fiona_get_dim('openmars', 'lat', size=nlat_openmars)
    call fiona_get_dim('openmars', 'lev', size=nlev_openmars)

    allocate(openmars_lon(nlon_openmars))
    allocate(openmars_lat(nlat_openmars))
    allocate(openmars_lev(nlev_openmars))
    allocate(openmars_u  (nlon_openmars,nlat_openmars,nlev_openmars))
    allocate(openmars_v  (nlon_openmars,nlat_openmars,nlev_openmars))
    allocate(openmars_t  (nlon_openmars,nlat_openmars,nlev_openmars))
    allocate(openmars_p  (nlon_openmars,nlat_openmars,nlev_openmars))
    allocate(openmars_ps (nlon_openmars,nlat_openmars              ))

    call fiona_start_input('openmars')
    call fiona_input('openmars', 'lon' , openmars_lon)
    call fiona_input('openmars', 'lat' , openmars_lat)
    call fiona_input('openmars', 'lev' , openmars_lev)
    call fiona_input('openmars', 'u'   , openmars_u  )
    call fiona_input('openmars', 'v'   , openmars_v  )
    call fiona_input('openmars', 'temp', openmars_t  )
    call fiona_input('openmars', 'ps'  , openmars_ps )
    call fiona_end_input('openmars')

    allocate(tmp(nlon_openmars,nlat_openmars,nlev_openmars))

    ! Reverse latitude order (from South Pole to North Pole).
    tmp(1,:,1) = openmars_lat(nlat_openmars:1:-1); openmars_lat = tmp(1,:,1)
    tmp(:,:,:) = openmars_u (:,nlat_openmars:1:-1,:); openmars_u  = tmp(:,:,:)
    tmp(:,:,:) = openmars_v (:,nlat_openmars:1:-1,:); openmars_v  = tmp(:,:,:)
    tmp(:,:,:) = openmars_t (:,nlat_openmars:1:-1,:); openmars_t  = tmp(:,:,:)
    tmp(:,:,1) = openmars_ps(:,nlat_openmars:1:-1  ); openmars_ps = tmp(:,:,1)

    ! Flip along the Meridian and change longitude to 0-360.
    tmp(1:nlon_openmars/2,1,1) = openmars_lon(nlon_openmars/2+1:nlon_openmars)
    tmp(nlon_openmars/2+1:nlon_openmars,1,1) = 360 + openmars_lon(1:nlon_openmars/2)
    openmars_lon = tmp(:,1,1)
    tmp(1:nlon_openmars/2,:,:) = openmars_u(nlon_openmars/2+1:nlon_openmars,:,:)
    tmp(nlon_openmars/2+1:nlon_openmars,:,:) = openmars_u(1:nlon_openmars/2,:,:)
    openmars_u = tmp
    tmp(1:nlon_openmars/2,:,:) = openmars_v(nlon_openmars/2+1:nlon_openmars,:,:)
    tmp(nlon_openmars/2+1:nlon_openmars,:,:) = openmars_v(1:nlon_openmars/2,:,:)
    openmars_v = tmp
    tmp(1:nlon_openmars/2,:,:) = openmars_t(nlon_openmars/2+1:nlon_openmars,:,:)
    tmp(nlon_openmars/2+1:nlon_openmars,:,:) = openmars_t(1:nlon_openmars/2,:,:)
    openmars_t = tmp
    tmp(1:nlon_openmars/2,:,1) = openmars_ps(nlon_openmars/2+1:nlon_openmars,:)
    tmp(nlon_openmars/2+1:nlon_openmars,:,1) = openmars_ps(1:nlon_openmars/2,:)
    openmars_ps = tmp(:,:,1)

    ! Reverse vertical levels (from top to bottom).
    tmp = openmars_u(:,:,nlev_openmars:1:-1); openmars_u = tmp
    tmp = openmars_v(:,:,nlev_openmars:1:-1); openmars_v = tmp
    tmp = openmars_t(:,:,nlev_openmars:1:-1); openmars_t = tmp
    tmp(1,1,:) = openmars_lev(nlev_openmars:1:-1); openmars_lev = tmp(1,1,:)

    deallocate(tmp)

    ! Calcuate pressure based on sigma formula.
    do k = 1, nlev_openmars
      openmars_p(:,:,k) = openmars_lev(k) * openmars_ps(:,:)
    end do

  end subroutine openmars_reader_run

  subroutine openmars_reader_final()

    if (allocated(openmars_lon)) deallocate(openmars_lon)
    if (allocated(openmars_lat)) deallocate(openmars_lat)
    if (allocated(openmars_lev)) deallocate(openmars_lev)
    if (allocated(openmars_u  )) deallocate(openmars_u  )
    if (allocated(openmars_v  )) deallocate(openmars_v  )
    if (allocated(openmars_t  )) deallocate(openmars_t  )
    if (allocated(openmars_p  )) deallocate(openmars_p  )
    if (allocated(openmars_ps )) deallocate(openmars_ps )

  end subroutine openmars_reader_final

end module openmars_reader_mod
