module topo_mod

  use fiona
  use flogger
  use string
  use const_mod
  use block_mod
  use process_mod

  implicit none

  private

  public topo_read
  public topo_regrid
  public topo_final

  integer num_topo_lon
  integer num_topo_lat
  real(r8), allocatable :: topo_lon(:)
  real(r8), allocatable :: topo_lat(:)
  real(r8), allocatable :: topo(:,:)

contains

  subroutine topo_read(topo_file)

    character(*), intent(in) :: topo_file

    call topo_final()

    call log_notice('Use ' // trim(topo_file) // ' as topography.')

    call fiona_open_dataset('topo', file_path=topo_file)
    call fiona_get_dim('topo', 'x', size=num_topo_lon)
    call fiona_get_dim('topo', 'y', size=num_topo_lat)

    call log_notice('num_topo_lon = ' // to_string(num_topo_lon) // ', num_topo_lat = ' // to_string(num_topo_lat))

    allocate(topo_lon(num_topo_lon))
    allocate(topo_lat(num_topo_lat))
    allocate(topo(num_topo_lon,num_topo_lat))

    call fiona_start_input('topo')
    call fiona_input('topo', 'x', topo_lon)
    call fiona_input('topo', 'y', topo_lat)
    call fiona_input('topo', 'z', topo)
    call fiona_end_input('topo')

  end subroutine topo_read

  subroutine topo_regrid(block)

    type(block_type), intent(inout) :: block

    integer ix1(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)
    integer ix2(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)
    integer jy1(block%mesh%full_lat_ibeg:block%mesh%full_lat_iend)
    integer jy2(block%mesh%full_lat_ibeg:block%mesh%full_lat_iend)

    real(8) lon1, lon2
    real(8) lat1, lat2
    real(8) gzs0, landfrac0, sgh0, num_topo_grid, num_topo_grid0

    integer i, j, k, l, i1, i2, j1, j2

    associate (mesh => block%mesh, gzs => block%static%gzs)
      ! To obtain the indices of east-west boundary for boxes of model grids
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon1 = mesh%half_lon_deg(i-1)
        lon2 = mesh%half_lon_deg(i  )

        if (lon1 < topo_lon(1)) then
          lon1 = lon1 + 360.0d0
        else if (lon1 > topo_lon(num_topo_lon)) then
          lon1 = lon1 - 360.0d0
        end if
        if (lon2 < topo_lon(1)) then
          lon2 = lon2 + 360.0d0
        else if (lon2 > topo_lon(num_topo_lon)) then
          lon2 = lon2 - 360.0d0
        end if

        i1 = 0
        i2 = 0
        do k = 1, num_topo_lon
          if (topo_lon(k) >= lon1 .and. i1 == 0) i1 = k
          if (topo_lon(k) >= lon2 .and. i2 == 0) i2 = k - 1
          if (i1 /= 0 .and. i2 /= 0) exit
        end do

        if (i1 == 0 .or. i2 == 0) then
          print *, i, mesh%full_lon_deg(i)
          if (is_root_proc()) call log_error('topo_regrid error!', __FILE__, __LINE__)
        end if

        ix1(i) = i1
        ix2(i) = i2
      end do

      ! To obtain the indices of south-north boundary for boxes of model grids
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        lat1 = mesh%half_lat_deg(j-1)
        lat2 = mesh%half_lat_deg(j  )

        j1 = 0
        j2 = 0
        do k = 1, num_topo_lat
          if (topo_lat(k) >= lat1 .and. j1 == 0) j1 = k
          if (topo_lat(k) >= lat2 .and. j2 == 0) j2 = k - 1
          if (j1 /= 0 .and. j2 /= 0) exit
        end do

        jy1(j) = j1
        jy2(j) = j2
      end do

      if (mesh%has_south_pole()) then
        jy1(mesh%full_lat_ibeg) = 1
        jy2(mesh%full_lat_ibeg) = jy1(mesh%full_lat_ibeg+1) - 1
      end if
      if (mesh%has_north_pole()) then
        jy1(mesh%full_lat_iend) = jy2(mesh%full_lat_iend-1) + 1
        jy2(mesh%full_lat_iend) = num_topo_lat
      end if
 
      ! To calculate the averaged terrain of the model grid box as the model terrain,
      ! the corresponding standard deviation, and the land-sea mask
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        j1 = jy1(j)
        j2 = jy2(j)

        if (ix1(mesh%full_lon_ibeg) > ix2(mesh%full_lon_ibeg)) then
          ! Split into two parts
          i1 = 1
          i2 = ix2(mesh%full_lon_ibeg)
          call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
          gzs(mesh%full_lon_ibeg,j) = gzs0
          num_topo_grid0 = num_topo_grid

          i1 = ix1(mesh%full_lon_ibeg)
          i2 = num_topo_lon
          call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
          gzs(mesh%full_lon_ibeg,j) = gzs(mesh%full_lon_ibeg,j) + gzs0
          num_topo_grid0 = num_topo_grid + num_topo_grid0

          gzs(mesh%full_lon_ibeg,j) = gzs(mesh%full_lon_ibeg,j) / num_topo_grid0
        else
          i1 = ix1(mesh%full_lon_ibeg)
          i2 = ix2(mesh%full_lon_ibeg)
          call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
          gzs(mesh%full_lon_ibeg,j) = gzs0 / num_topo_grid
        end if

        do i = mesh%full_lon_ibeg + 1, mesh%full_lon_iend - 1
          i1 = ix1(i)
          i2 = ix2(i)
          call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
          gzs(i,j) = gzs0 / num_topo_grid
        end do

        if (ix1(mesh%full_lon_iend) > ix2(mesh%full_lon_iend)) then
          i1 = 1
          i2 = ix2(mesh%full_lon_iend)
          call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
          gzs(mesh%full_lon_iend,j) = gzs0
          num_topo_grid0 = num_topo_grid

          i1 = ix1(mesh%full_lon_ibeg)
          i2 = num_topo_lon
          call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
          gzs(mesh%full_lon_iend,j) = gzs(mesh%full_lon_iend,j) + gzs0
          num_topo_grid0 = num_topo_grid + num_topo_grid0

          gzs(mesh%full_lon_iend,j) = gzs(mesh%full_lon_iend,j) / num_topo_grid0
        else
          i1 = ix1(mesh%full_lon_iend)
          i2 = ix2(mesh%full_lon_iend)
          call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
          gzs(mesh%full_lon_iend,j) = gzs0 / num_topo_grid
        end if
      end do

      ! For the poles
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        j1 = jy1(j)
        j2 = jy2(j)
        i1 = 1
        i2 = num_topo_lon
        call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
        gzs0 = gzs0 / num_topo_grid
        sgh0 = sgh0 / num_topo_grid
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gzs(i,j) = gzs0
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        j1 = jy1(j)
        j2 = jy2(j)
        i1 = 1
        i2 = num_topo_lon
        call termask(gzs0, sgh0, landfrac0, num_topo_grid, i1, i2, j1, j2)
        gzs0 = gzs0 / num_topo_grid
        sgh0 = sgh0 / num_topo_grid
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gzs(i,j) = gzs0
        end do
      end if
    end associate

  end subroutine topo_regrid

  subroutine termask(gzs, sgh, landfrac, num_topo_grid, i1, i2, j1, j2)

    integer, intent(in) :: i1, i2, j1, j2
    real(r8), intent(out) :: gzs, sgh, landfrac, num_topo_grid

    real(r8) osmm, oshs, stdgs
    integer k, l
  
    gzs           = 0.0d0
    landfrac      = 0.0d0
    num_topo_grid = 0.0d0
    do k = j1, j2
      do l = i1, i2
        if (topo(l,k) > 0.0d0) then
          gzs      = gzs      + topo(l,k)
          landfrac = landfrac + 1.0d0
        end if
        num_topo_grid = num_topo_grid + 1.0d0
      end do
    end do
    osmm = gzs / num_topo_grid
  
    sgh = 0.0d0
    do k = j1, j2
      do l = i1, i2
        oshs = topo(l,k)
        if (oshs <= 0.0d0) oshs = 0.0
        stdgs = oshs - osmm
        sgh = sgh + stdgs * stdgs
      end do
    end do
  
  end subroutine termask

  subroutine topo_final()

    if (allocated(topo_lon)) deallocate(topo_lon)
    if (allocated(topo_lat)) deallocate(topo_lat)
    if (allocated(topo    )) deallocate(topo    )

  end subroutine topo_final

end module topo_mod
