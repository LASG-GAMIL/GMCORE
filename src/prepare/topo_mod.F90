module topo_mod

  use fiona
  use flogger
  use string
  use const_mod
  use block_mod
  use process_mod
  use parallel_mod

  implicit none

  private

  public topo_read
  public topo_regrid
  public topo_zero
  public topo_smth
  public topo_final

  integer num_topo_lon
  integer num_topo_lat
  real(r8), allocatable, dimension(:  ) :: topo_lon
  real(r8), allocatable, dimension(:  ) :: topo_lat
  real(r8), allocatable, dimension(:,:) :: topo_gzs

contains

  subroutine topo_read(topo_file)

    character(*), intent(in) :: topo_file

    call topo_final()

    call log_notice('Use ' // trim(topo_file) // ' as topography.')

    call fiona_open_dataset('topo', file_path=topo_file)
    call fiona_get_dim('topo', 'x', size=num_topo_lon)
    call fiona_get_dim('topo', 'y', size=num_topo_lat)

    call log_notice('num_topo_lon = ' // to_str(num_topo_lon) // ', num_topo_lat = ' // to_str(num_topo_lat))

    allocate(topo_lon(num_topo_lon))
    allocate(topo_lat(num_topo_lat))
    allocate(topo_gzs(num_topo_lon,num_topo_lat))

    call fiona_start_input('topo')
    call fiona_input('topo', 'x', topo_lon)
    call fiona_input('topo', 'y', topo_lat)
    call fiona_input('topo', 'z', topo_gzs)
    call fiona_end_input('topo')

    topo_gzs = topo_gzs * g

  end subroutine topo_read

  subroutine topo_regrid(block)

    type(block_type), intent(inout) :: block

    integer ix1(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)
    integer ix2(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)
    integer jy1(block%mesh%full_lat_ibeg:block%mesh%full_lat_iend)
    integer jy2(block%mesh%full_lat_ibeg:block%mesh%full_lat_iend)

    real(r8) lon1, lon2
    real(r8) lat1, lat2
    real(r8) gzs0, landfrac0, sgh0, grid_count, grid_count0

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

        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          i1 = ix1(i)
          i2 = ix2(i)
          if (i1 > i2) then
            ! Split into two parts
            i1 = 1
            i2 = ix2(i)
            call count_topo_grids(gzs0, sgh0, landfrac0, grid_count, i1, i2, j1, j2)
            gzs(i,j) = gzs0
            grid_count0 = grid_count
            i1 = ix1(i)
            i2 = num_topo_lon
            call count_topo_grids(gzs0, sgh0, landfrac0, grid_count, i1, i2, j1, j2)
            gzs(i,j) = gzs(i,j) + gzs0
            grid_count0 = grid_count + grid_count0
            gzs(i,j) = gzs(i,j) / grid_count0
          else
            call count_topo_grids(gzs0, sgh0, landfrac0, grid_count, i1, i2, j1, j2)
            gzs(i,j) = gzs0 / grid_count
          end if
        end do
      end do

      ! For the poles
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        j1 = jy1(j)
        j2 = jy2(j)
        i1 = 1
        i2 = num_topo_lon
        call count_topo_grids(gzs0, sgh0, landfrac0, grid_count, i1, i2, j1, j2)
        gzs0 = gzs0 / grid_count
        sgh0 = sgh0 / grid_count
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
        call count_topo_grids(gzs0, sgh0, landfrac0, grid_count, i1, i2, j1, j2)
        gzs0 = gzs0 / grid_count
        sgh0 = sgh0 / grid_count
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gzs(i,j) = gzs0
        end do
      end if
      ! Fill halo outside poles for later smoothing if used.
      if (mesh%has_south_pole()) then
        do j = 1, mesh%lat_halo_width
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            i1 = i + mesh%num_full_lon / 2
            if (i1 > mesh%num_full_lon) i1 = i1 - mesh%num_full_lon
            gzs(i,mesh%full_lat_ibeg-j) = gzs(i1,mesh%full_lat_ibeg+j-1)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        do j = 1, mesh%lat_halo_width
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            i1 = i + mesh%num_full_lon / 2
            if (i1 > mesh%num_full_lon) i1 = i1 - mesh%num_full_lon
            gzs(i,mesh%full_lat_iend+j) = gzs(i1,mesh%full_lat_iend-j+1)
          end do
        end do
      end if
      call fill_halo(block, gzs, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine topo_regrid

  subroutine count_topo_grids(gzs, sgh, landfrac, grid_count, i1, i2, j1, j2)

    integer, intent(in) :: i1, i2, j1, j2
    real(r8), intent(out) :: gzs, sgh, landfrac, grid_count

    real(r8) osmm, oshs, stdgs
    integer k, l
  
    gzs        = 0.0d0
    landfrac   = 0.0d0
    grid_count = 0.0d0
    do k = j1, j2
      do l = i1, i2
        if (topo_gzs(l,k) > 0.0d0) then
          gzs      = gzs      + topo_gzs(l,k)
          landfrac = landfrac + 1.0d0
        end if
        grid_count = grid_count + 1.0d0
      end do
    end do
    osmm = gzs / grid_count
  
    sgh = 0.0d0
    do k = j1, j2
      do l = i1, i2
        oshs = topo_gzs(l,k)
        if (oshs <= 0.0d0) oshs = 0.0
        stdgs = oshs - osmm
        sgh = sgh + stdgs * stdgs
      end do
    end do
  
  end subroutine count_topo_grids

  subroutine topo_final()

    if (allocated(topo_lon)) deallocate(topo_lon)
    if (allocated(topo_lat)) deallocate(topo_lat)
    if (allocated(topo_gzs)) deallocate(topo_gzs)

  end subroutine topo_final

  subroutine topo_zero(block, min_lon, max_lon, min_lat, max_lat)

    type(block_type), intent(inout) :: block
    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat
    real(r8), intent(in) :: max_lat

    integer ibeg, iend, jbeg, jend, i, j
    real(r8) ib, jb, zero_wgt

    associate (mesh => block%mesh, gzs => block%static%gzs)
      do ibeg = mesh%full_lon_ibeg, mesh%full_lon_iend
        if (block%mesh%full_lon_deg(ibeg) >= min_lon) exit
      end do
      do iend = mesh%full_lon_ibeg, mesh%full_lon_iend
        if (block%mesh%full_lon_deg(iend) >= max_lon) exit
      end do
      do jbeg = mesh%full_lat_ibeg, mesh%full_lat_iend
        if (block%mesh%full_lat_deg(jbeg) >= min_lat) exit
      end do
      do jend = mesh%full_lat_ibeg, mesh%full_lat_iend
        if (block%mesh%full_lat_deg(jend) >= max_lat) exit
      end do
      ibeg = max(1, ibeg); iend = min(mesh%full_lon_iend, iend)
      if (ibeg > mesh%full_lon_iend .or. iend > mesh%full_lon_iend .or. &
          jbeg > mesh%full_lat_iend .or. jend > mesh%full_lat_iend) then
        call log_error('Invalid region parameters for topo_zero!', __FILE__, __LINE__)
      end if
      call log_notice('Zero topography in (' // to_str(ibeg) // '-' // to_str(iend) // ',' // &
                                                to_str(jbeg) // '-' // to_str(jend) // ').')
      do j = jbeg, jend
        do i = ibeg, iend
          ib = merge(i - ibeg, iend - i, i < (ibeg + iend) * 0.5)
          jb = merge(j - jbeg, jend - j, j < (jbeg + jend) * 0.5)
          zero_wgt = 1 - max(exp(ib**2 * log(0.1_r8) / ((iend - ibeg) * 0.1_r8)**2), &
                             exp(jb**2 * log(0.1_r8) / ((jend - jbeg) * 0.1_r8)**2))
          gzs(i,j) = (1 - zero_wgt) * gzs(i,j)
        end do
      end do
    end associate

  end subroutine topo_zero

  subroutine topo_smth(block, min_lon, max_lon, min_lat, max_lat, steps)

    type(block_type), intent(inout) :: block
    real(r8), intent(in) :: min_lon
    real(r8), intent(in) :: max_lon
    real(r8), intent(in) :: min_lat
    real(r8), intent(in) :: max_lat
    integer, intent(in) :: steps

    real(r8), allocatable :: smooth_kernel(:,:)
    real(r8), allocatable :: smoothed_gzs(:,:)
    integer ibeg, iend, jbeg, jend, i, j, step
    real(r8) ib, jb, smooth_wgt

    associate (mesh => block%mesh, gzs => block%static%gzs)
      do ibeg = mesh%full_lon_ibeg, mesh%full_lon_iend
        if (block%mesh%full_lon_deg(ibeg) >= min_lon) exit
      end do
      do iend = mesh%full_lon_ibeg, mesh%full_lon_iend
        if (block%mesh%full_lon_deg(iend) >= max_lon) exit
      end do
      do jbeg = mesh%full_lat_ibeg, mesh%full_lat_iend
        if (block%mesh%full_lat_deg(jbeg) >= min_lat) exit
      end do
      do jend = mesh%full_lat_ibeg, mesh%full_lat_iend
        if (block%mesh%full_lat_deg(jend) >= max_lat) exit
      end do
      ibeg = max(1, ibeg); iend = min(mesh%full_lon_iend, iend)
      if (ibeg > mesh%full_lon_iend .or. iend > mesh%full_lon_iend .or. &
          jbeg > mesh%full_lat_iend .or. jend > mesh%full_lat_iend) then
        call log_error('Invalid region parameters for topo_smooth!', __FILE__, __LINE__)
      end if
      call log_notice('Smooth topography in (' // to_str(ibeg) // '-' // to_str(iend) // ',' // &
                                                  to_str(jbeg) // '-' // to_str(jend) // ').')
      allocate(smooth_kernel(-5:5,-5:5))
      do j = -5, 5
        do i = -5, 5
          smooth_kernel(i,j) = exp(-0.5_r8 * (i**2 + j**2))
        end do
      end do
      smooth_kernel = smooth_kernel / sum(smooth_kernel)
      allocate(smoothed_gzs(ibeg:iend,jbeg:jend))
      do step = 1, steps
        do j = jbeg, jend
          do i = ibeg, iend
            if (gzs(i,j) == 0) then
              smoothed_gzs(i,j) = gzs(i,j)
            else
              ib = merge(i - ibeg, iend - i, i < (ibeg + iend) * 0.5)
              jb = merge(j - jbeg, jend - j, j < (jbeg + jend) * 0.5)
              smooth_wgt = 1 - max(exp(ib**2 * log(0.1_r8) / ((iend - ibeg) * 0.1_r8)**2), &
                                   exp(jb**2 * log(0.1_r8) / ((jend - jbeg) * 0.1_r8)**2))
              smoothed_gzs(i,j) = sum(smooth_kernel * gzs(i-5:i+5,j-5:j+5)) * smooth_wgt + (1 - smooth_wgt) * gzs(i,j)
            end if
          end do
        end do
        gzs(ibeg:iend,jbeg:jend) = smoothed_gzs(ibeg:iend,jbeg:jend)
      end do
      deallocate(smooth_kernel, smoothed_gzs)
    end associate

  end subroutine topo_smth

end module topo_mod
