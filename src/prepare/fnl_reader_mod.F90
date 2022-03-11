#ifdef HAS_ECCODES
module fnl_reader_mod

  use eccodes
  use const_mod
  use process_mod

  implicit none

  integer num_fnl_lon
  integer num_fnl_lat
  integer num_fnl_lev

  real(r8), allocatable, dimension(:    ) :: fnl_lon
  real(r8), allocatable, dimension(:    ) :: fnl_lat
  real(r8), allocatable, dimension(:    ) :: fnl_lev
  real(r8), allocatable, dimension(:,:,:) :: fnl_u
  real(r8), allocatable, dimension(:,:,:) :: fnl_v
  real(r8), allocatable, dimension(:,:,:) :: fnl_t
  real(r8), allocatable, dimension(:,:  ) :: fnl_ps
  real(r8), allocatable, dimension(:,:  ) :: fnl_zs

contains

  subroutine fnl_reader_run(bkg_file)

    character(*), intent(in) :: bkg_file

    integer ifile, igrib, idx, iret, i, j, k, numberOfPoints
    real(r8) start_lon, end_lon, dlon
    real(r8) start_lat, end_lat, dlat
    real(r8) plev(1000)
    real(r8), allocatable :: global_lons_dummy(:)
    real(r8), allocatable :: global_lats_dummy(:)
    real(r8), allocatable :: global_array_dummy(:)
    real(r8), allocatable :: global_array(:,:)
    character(30) shortName, typeOfLevel

    call codes_index_create(idx, bkg_file, 'shortName,typeOfLevel')
    call codes_index_select(idx, 'shortName', 'u')
    call codes_index_select(idx, 'typeOfLevel', 'isobaricInhPa')
    call codes_new_from_index(idx, igrib, iret)
    call codes_get(igrib, 'Ni', num_fnl_lon)
    call codes_get(igrib, 'Nj', num_fnl_lat)
    call codes_get(igrib, 'longitudeOfFirstGridPoint', start_lon); start_lon = start_lon * 1d-6
    call codes_get(igrib, 'longitudeOfLastGridPoint' , end_lon  ); end_lon   = end_lon   * 1d-6
    call codes_get(igrib, 'latitudeOfFirstGridPoint' , start_lat); start_lat = start_lat * 1d-6
    call codes_get(igrib, 'latitudeOfLastGridPoint'  , end_lat  ); end_lat   = end_lat   * 1d-6
    call codes_get(igrib, 'iDirectionIncrement'      , dlon     ); dlon      = dlon      * 1d-6
    call codes_get(igrib, 'jDirectionIncrement'      , dlat     ); dlat      = dlat      * 1d-6
    call codes_get(igrib, 'numberOfPoints', numberOfPoints)
    do while (iret /= CODES_END_OF_INDEX)
      num_fnl_lev = num_fnl_lev + 1
      call codes_get(igrib, 'level', plev(num_fnl_lev))
      call codes_new_from_index(idx, igrib, iret)
    end do
    call codes_release(igrib)

    allocate(fnl_lon(num_fnl_lon))
    allocate(fnl_lat(num_fnl_lat))
    allocate(fnl_lev(num_fnl_lev)); fnl_lev = plev(:num_fnl_lev)
    allocate(fnl_u  (num_fnl_lon,num_fnl_lat,num_fnl_lev))
    allocate(fnl_v  (num_fnl_lon,num_fnl_lat,num_fnl_lev))
    allocate(fnl_t  (num_fnl_lon,num_fnl_lat,num_fnl_lev))
    allocate(fnl_ps (num_fnl_lon,num_fnl_lat))
    allocate(fnl_zs (num_fnl_lon,num_fnl_lat))

    allocate(global_lons_dummy(numberOfPoints))
    allocate(global_lats_dummy(numberOfPoints))
    allocate(global_array_dummy(numberOfPoints))
    allocate(global_array(num_fnl_lon,num_fnl_lat))

    do i = 1, num_fnl_lon
      fnl_lon(i) = start_lon + (i - 1) * dlon
    end do
    ! NOTE: Here we like from South Pole to North Pole.
    do j = 1, num_fnl_lat
      fnl_lat(j) = end_lat + (j - 1) * dlat
    end do

    call codes_index_create(idx, bkg_file, 'shortName,typeOfLevel')
    ! --------------------------------------------------------------------------
    ! u
    call codes_index_select(idx, 'shortName', 'u')
    call codes_index_select(idx, 'typeOfLevel', 'isobaricInhPa')
    call codes_new_from_index(idx, igrib, iret)
    k = 0
    do while (iret /= CODES_END_OF_INDEX)
      call codes_grib_get_data(igrib, global_lats_dummy, global_lons_dummy, global_array_dummy)
      global_array = reshape(global_array_dummy, [num_fnl_lon,num_fnl_lat])
      k = k + 1
      do j = 1, num_fnl_lat
        do i = 1, num_fnl_lon
          fnl_u(i,j,k) = global_array(i,num_fnl_lat-j+1)
        end do
      end do
      call codes_new_from_index(idx, igrib, iret)
    end do
    call codes_release(igrib)
    ! --------------------------------------------------------------------------
    ! v
    call codes_index_select(idx, 'shortName', 'v')
    call codes_index_select(idx, 'typeOfLevel', 'isobaricInhPa')
    call codes_new_from_index(idx, igrib, iret)
    k = 0
    do while (iret /= CODES_END_OF_INDEX)
      call codes_grib_get_data(igrib, global_lats_dummy, global_lons_dummy, global_array_dummy)
      global_array = reshape(global_array_dummy, [num_fnl_lon,num_fnl_lat])
      k = k + 1
      do j = 1, num_fnl_lat
        do i = 1, num_fnl_lon
          fnl_v(i,j,k) = global_array(i,num_fnl_lat-j+1)
        end do
      end do
      call codes_new_from_index(idx, igrib, iret)
    end do
    call codes_release(igrib)
    ! --------------------------------------------------------------------------
    call codes_index_select(idx, 'shortName', 't')
    call codes_index_select(idx, 'typeOfLevel', 'isobaricInhPa')
    call codes_new_from_index(idx, igrib, iret)
    k = 0
    do while (iret /= CODES_END_OF_INDEX)
      call codes_grib_get_data(igrib, global_lats_dummy, global_lons_dummy, global_array_dummy)
      global_array = reshape(global_array_dummy, [num_fnl_lon,num_fnl_lat])
      k = k + 1
      do j = 1, num_fnl_lat
        do i = 1, num_fnl_lon
          fnl_t(i,j,k) = global_array(i,num_fnl_lat-j+1)
        end do
      end do
      call codes_new_from_index(idx, igrib, iret)
    end do
    call codes_release(igrib)
    ! --------------------------------------------------------------------------
    ! sp
    call codes_index_select(idx, 'shortName', 'sp')
    call codes_index_select(idx, 'typeOfLevel', 'surface')
    call codes_new_from_index(idx, igrib, iret)
    call codes_grib_get_data(igrib, global_lats_dummy, global_lons_dummy, global_array_dummy)
    global_array = reshape(global_array_dummy, [num_fnl_lon,num_fnl_lat])
    do j = 1, num_fnl_lat
      do i = 1, num_fnl_lon
        fnl_ps(i,j) = global_array(i,num_fnl_lat-j+1)
      end do
    end do
    call codes_release(igrib)
    ! --------------------------------------------------------------------------
    ! orog
    call codes_index_select(idx, 'shortName', 'orog')
    call codes_index_select(idx, 'typeOfLevel', 'surface')
    call codes_new_from_index(idx, igrib, iret)
    call codes_grib_get_data(igrib, global_lats_dummy, global_lons_dummy, global_array_dummy)
    global_array = reshape(global_array_dummy, [num_fnl_lon,num_fnl_lat])
    do j = 1, num_fnl_lat
      do i = 1, num_fnl_lon
        fnl_zs(i,j) = global_array(i,num_fnl_lat-j+1)
      end do
    end do
    call codes_release(igrib)

    call codes_index_release(idx)

    deallocate(global_lats_dummy, global_lons_dummy, global_array_dummy, global_array)

    ! Change units.
    fnl_lev = fnl_lev * 100.0_r8

  end subroutine fnl_reader_run

  subroutine fnl_reader_final()

    if (allocated(fnl_lon)) deallocate(fnl_lon)
    if (allocated(fnl_lat)) deallocate(fnl_lat)
    if (allocated(fnl_lev)) deallocate(fnl_lev)
    if (allocated(fnl_u  )) deallocate(fnl_u  )
    if (allocated(fnl_v  )) deallocate(fnl_v  )
    if (allocated(fnl_t  )) deallocate(fnl_t  )
    if (allocated(fnl_ps )) deallocate(fnl_ps )
    if (allocated(fnl_zs )) deallocate(fnl_zs )

  end subroutine fnl_reader_final

end module fnl_reader_mod
#endif
