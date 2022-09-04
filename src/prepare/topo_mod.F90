module topo_mod

  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod
  use block_mod
  use process_mod
  use parallel_mod
  use filter_mod
  use laplace_damp_mod

  implicit none

  private

  public topo_read
  public topo_regrid
  public topo_smooth
  public topo_final

  type topo_block_type
    logical :: initialized = .false.
    integer id
    integer lon_ibeg
    integer lon_iend
    integer lat_ibeg
    integer lat_iend
    integer num_lon
    integer num_lat
    real(r8), allocatable :: lon(:)
    real(r8), allocatable :: lat(:)
    real(r8), allocatable :: gzs(:,:)
  contains
    procedure :: init      => topo_block_init
    procedure :: clear     => topo_block_clear
    procedure :: contain   => topo_block_contain
    procedure :: fill_grid => topo_block_fill_grid
    final :: topo_block_final
  end type topo_block_type

  type(topo_block_type) topo, topo1, topo2

contains

  subroutine topo_block_init(this, id, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    class(topo_block_type), intent(inout) :: this
    integer, intent(in) :: id
    integer, intent(in) :: lon_ibeg
    integer, intent(in) :: lon_iend
    integer, intent(in) :: lat_ibeg
    integer, intent(in) :: lat_iend

    call this%clear()

    this%id       = id
    this%lon_ibeg = lon_ibeg
    this%lon_iend = lon_iend
    this%lat_ibeg = lat_ibeg
    this%lat_iend = lat_iend
    this%num_lon  = lon_iend - lon_ibeg + 1
    this%num_lat  = lat_iend - lat_ibeg + 1

    allocate(this%lon(lon_ibeg:lon_iend))
    allocate(this%lat(lat_ibeg:lat_iend))
    allocate(this%gzs(lon_ibeg:lon_iend,lat_ibeg:lat_iend))

    this%initialized = .true.

  end subroutine topo_block_init

  subroutine topo_block_clear(this)

    class(topo_block_type), intent(inout) :: this

    if (allocated(this%lon)) deallocate(this%lon)
    if (allocated(this%lat)) deallocate(this%lat)
    if (allocated(this%gzs)) deallocate(this%gzs)

    this%initialized = .false.

  end subroutine topo_block_clear

  logical function topo_block_contain(this, lon, lat) result(res)

    class(topo_block_type), intent(inout) :: this
    real(r8), intent(in), optional :: lon
    real(r8), intent(in), optional :: lat

    if (.not. this%initialized) then
      res = .false.
      return
    end if

    if (present(lon) .and. present(lat)) then
      res = this%lon(this%lon_ibeg) <= lon .and. lon <= this%lon(this%lon_iend) .and. &
            this%lat(this%lat_ibeg) <= lat .and. lat <= this%lat(this%lat_iend)
    else if (present(lon)) then
      res = this%lon(this%lon_ibeg) <= lon .and. lon <= this%lon(this%lon_iend)
    else if (present(lat)) then
      res = this%lat(this%lat_ibeg) <= lat .and. lat <= this%lat(this%lat_iend)
    end if

  end function topo_block_contain

  subroutine topo_block_fill_grid(this, lon1, lon2, lat1, lat2, gzs, std, lnd, cnt)

    class(topo_block_type), intent(inout) :: this
    real(r8), intent(in) :: lon1
    real(r8), intent(in) :: lon2
    real(r8), intent(in) :: lat1
    real(r8), intent(in) :: lat2
    real(r8), intent(inout) :: gzs
    real(r8), intent(inout) :: std
    real(r8), intent(inout) :: lnd
    integer , intent(inout) :: cnt

    integer is, ie, js, je, i, j

    if (.not. this%initialized) then
      return
    end if

    if (lon1 < this%lon(this%lon_ibeg)) then
      is = this%lon_ibeg
    else if (lon1 >= this%lon(this%lon_ibeg)) then
      do is = this%lon_ibeg, this%lon_iend
        if (lon1 <= this%lon(is)) exit
      end do
    end if
    if (lon2 > this%lon(this%lon_iend)) then
      ie = this%lon_iend
    else if (lon2 <= this%lon(this%lon_iend)) then
      do ie = this%lon_iend, this%lon_ibeg, -1
        if (lon2 >= this%lon(ie)) exit
      end do
    end if
    do js = this%lat_ibeg, this%lat_iend
      if (lat1 <= this%lat(js)) exit
    end do
    do je = this%lat_iend, this%lat_ibeg, -1
      if (lat2 >= this%lat(je)) exit
    end do

    select case (planet)
    case ('earth')
      do j = js, je
        do i = is, ie
          if (this%lon(i) < lon1 .or. this%lon(i) > lon2 .or. this%lat(j) < lat1 .or. this%lat(j) > lat2) then
            stop 999
          end if
          if (this%gzs(i,j) > 0) then
            gzs = gzs + this%gzs(i,j)
            std = std + this%gzs(i,j)**2
          end if
        end do
      end do
      lnd = lnd + count(this%gzs(is:ie,js:je) > 0)
    case ('mars')
      gzs = gzs + sum(this%gzs(is:ie,js:je))
      std = std + sum(this%gzs(is:ie,js:je)**2)
      lnd = lnd + (ie - is + 1) * (je - js + 1)
    end select
    cnt = cnt + (ie - is + 1) * (je - js + 1)

  end subroutine topo_block_fill_grid

  subroutine topo_block_final(this)

    type(topo_block_type), intent(inout) :: this

    call this%clear()

  end subroutine topo_block_final

  subroutine topo_read(topo_file)

    character(*), intent(in) :: topo_file

    integer i, j, nx, ny
    integer is, ie, js, je
    integer is1, ie1, is2, ie2
    real(r8) lon1, lon2, lat1, lat2
    real(r8), allocatable :: x(:), y(:)

    call topo_final()

    if (is_root_proc()) call log_notice('Use ' // trim(topo_file) // ' as topography.')

    lon1 = blocks(1)%mesh%half_lon_deg(blocks(1)%mesh%half_lon_ibeg-1)
    lon2 = blocks(1)%mesh%half_lon_deg(blocks(1)%mesh%half_lon_iend  )
    lat1 = blocks(1)%mesh%half_lat_deg(blocks(1)%mesh%half_lat_ibeg-1)
    lat2 = blocks(1)%mesh%half_lat_deg(blocks(1)%mesh%half_lat_iend+1)
    lat1 = merge(lat1, -90.0_r8, lat1 /= inf)
    lat2 = merge(lat2,  90.0_r8, lat2 /= inf)
    is1 = -999; ie1 = -999; is2 = -999; ie2 = -999

    select case (topo_type)
    case ('etopo1')
      if (planet /= 'earth') call log_error('Topography file ' // trim(topo_file) // ' is used for the Earth!')
      call fiona_open_dataset('topo', file_path=topo_file)
      call fiona_get_dim('topo', 'x', size=nx)
      call fiona_get_dim('topo', 'y', size=ny)

      allocate(x(nx), y(ny))

      call fiona_start_input('topo')
      call fiona_input('topo', 'x', x)
      call fiona_input('topo', 'y', y)

      ! ETOPO1's longitude is from -180 to 180.
      if (lon1 > 180) then
        do i = 1, nx
          if (lon1 <= x(i) + 180) then
            is = i - 1
            exit
          end if
        end do
      else
        do i = 1, nx
          if (lon1 <= x(i)) then
            is = i - 1
            exit
          end if
        end do
      end if
      if (lon2 <= x(nx)) then
        do i = nx, 1, -1
          if (lon2 >= x(i)) then
            ie = i
            exit
          end if
        end do
      else
        ie = nx
        do i = 1, nx
          if (lon2 - 360 <= x(i)) then
            ie2 = i
            exit
          end if
        end do
        is2 = 1
      end if
      do j =  1, ny
        if (lat1 <= y(j)) then
          js = j - 1
          exit
        end if
      end do
      js = max(1, js)
      do j = ny, 1, -1
        if (lat2 >= y(j)) then
          je = j + 1
          exit
        end if
      end do
      je = min(ny, je)
      deallocate(x, y)

      call topo%init(0, is, ie, js, je)
      call fiona_input('topo', 'x', topo%lon(is:ie), start=is, count=ie-is+1)
      call fiona_input('topo', 'y', topo%lat(js:je), start=js, count=je-js+1)
      call fiona_input('topo', 'z', topo%gzs(is:ie,js:je), start=[is,js], count=[ie-is+1,je-js+1])
      if (is2 /= -999 .and. ie2 /= -999) then
        call topo2%init(2, is2, ie2, js, je)
        call fiona_input('topo', 'x', topo2%lon(is2:ie2), start=is2, count=ie2-is2+1)
        call fiona_input('topo', 'y', topo2%lat(js :je ), start=js , count=je -js +1)
        call fiona_input('topo', 'z', topo2%gzs(is2:ie2,js:je), start=[is2,js], count=[ie2-is2+1,je-js+1])
        topo2%lon = topo2%lon + 360
      end if
      call fiona_end_input('topo')
    case ('mola32')
      if (planet /= 'mars') call log_error('Topography file ' // trim(topo_file) // ' is used for the Mars!')
      call fiona_open_dataset('topo', file_path=topo_file)
      call fiona_get_dim('topo', 'longitude', size=nx)
      call fiona_get_dim('topo', 'latitude' , size=ny)

      allocate(x(nx), y(ny))

      call fiona_start_input('topo')
      call fiona_input('topo', 'longitude', x)
      call fiona_input('topo', 'latitude', y)

      if (lon1 >= x(1)) then
        do i = 1, nx
          if (lon1 <= x(i)) then
            is = i - 1
            exit
          end if
        end do
      else
        is = 1
        do i = nx, 1, -1
          if (360 + lon1 >= x(i)) then
            is1 = i
            exit
          end if
        end do
        ie1 = nx
      end if
      if (lon2 <= x(nx)) then
        do i = nx, 1, -1
          if (lon2 >= x(i)) then
            ie = i + 1
            exit
          end if
        end do
      else
        ie = nx
        do i = 1, nx
          if (lon2 - 360 <= x(i)) then
            ie2 = i
            exit
          end if
        end do
        is2 = 1
      end if
      do j = 1, ny
        if (lat1 <= y(j)) then
          js = j - 1
          exit
        end if
      end do
      js = max(1, js)
      do j = ny, 1, -1
        if (lat2 >= y(j)) then
          je = j + 1
          exit
        end if
      end do
      je = min(ny, je)
      deallocate(x, y)

      call topo%init(0, is, ie, js, je)
      call fiona_input('topo', 'longitude', topo%lon(is:ie), start=is, count=ie-is+1)
      call fiona_input('topo', 'latitude' , topo%lat(js:je), start=js, count=je-js+1)
      call fiona_input('topo', 'alt'      , topo%gzs(is:ie,js:je), start=[is,js], count=[ie-is+1,je-js+1])
      if (is1 /= -999 .and. ie1 /= -999) then
        call topo1%init(1, is1, ie1, js, je)
        call fiona_input('topo', 'longitude', topo1%lon(is1:ie1), start=is1, count=ie1-is1+1)
        call fiona_input('topo', 'latitude' , topo1%lat(js :je ), start=js , count=je -js +1)
        call fiona_input('topo', 'alt'      , topo1%gzs(is1:ie1,js:je), start=[is1,js], count=[ie1-is1+1,je-js+1])
        topo1%lon = topo1%lon - 360
      end if
      if (is2 /= -999 .and. ie2 /= -999) then
        call topo2%init(2, is2, ie2, js, je)
        call fiona_input('topo', 'longitude', topo2%lon(is2:ie2), start=is2, count=ie2-is2+1)
        call fiona_input('topo', 'latitude' , topo2%lat(js :je ), start=js , count=je -js +1)
        call fiona_input('topo', 'alt'      , topo2%gzs(is2:ie2,js:je), start=[is2,js], count=[ie2-is2+1,je-js+1])
        topo2%lon = topo2%lon + 360
      end if
      call fiona_end_input('topo')
    case default
      call log_error('Unknown topo_type "' // trim(topo_type) // '"!', pid=proc%id)
    end select

    topo%gzs = topo%gzs * g
    if (topo1%initialized) topo1%gzs = topo1%gzs * g
    if (topo2%initialized) topo2%gzs = topo2%gzs * g

  end subroutine topo_read

  subroutine topo_regrid(block)

    type(block_type), intent(inout) :: block

    real(r8) lon1, lon2, lat1, lat2, pole_gzs, pole_std, pole_lnd
    integer i, j, pole_n
    integer n(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)

    associate (mesh => block%mesh           , &
               lnd  => block%static%landmask, &
               gzs  => block%static%gzs     , &
               std  => block%static%zs_std)
      lnd = 0; gzs = 0; std = 0
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        lat1 = mesh%half_lat_deg(j-1); lat1 = merge(lat1, -90.0_r8, lat1 /= inf)
        lat2 = mesh%half_lat_deg(j  ); lat2 = merge(lat2,  90.0_r8, lat2 /= inf)
        n = 0
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          lon1 = mesh%half_lon_deg(i-1)
          lon2 = mesh%half_lon_deg(i  )
          if (topo%contain(lon=lon1) .or. topo%contain(lon=lon2)) then
            call topo %fill_grid(lon1, lon2, lat1, lat2, gzs(i,j), std(i,j), lnd(i,j), n(i))
          end if
          if (topo1%contain(lon=lon1)) then
            call topo1%fill_grid(lon1, lon2, lat1, lat2, gzs(i,j), std(i,j), lnd(i,j), n(i))
          end if
          if (topo2%contain(lon=lon2)) then
            call topo2%fill_grid(lon1, lon2, lat1, lat2, gzs(i,j), std(i,j), lnd(i,j), n(i))
          end if
          if (.not. mesh%is_pole(j)) then
            gzs(i,j) = gzs(i,j) / n(i)
            std(i,j) = (std(i,j) - 2 * gzs(i,j)**2 * n(i) + gzs(i,j)**2) / n(i) / g
            lnd(i,j) = lnd(i,j) / n(i)
          end if
        end do
        if (mesh%is_pole(j)) then
          call zonal_sum(proc%zonal_circle, gzs(mesh%full_lon_ibeg:mesh%full_lon_iend,j), pole_gzs)
          call zonal_sum(proc%zonal_circle, std(mesh%full_lon_ibeg:mesh%full_lon_iend,j), pole_std)
          call zonal_sum(proc%zonal_circle, lnd(mesh%full_lon_ibeg:mesh%full_lon_iend,j), pole_lnd)
          call zonal_sum(proc%zonal_circle, n, pole_n)
          gzs(mesh%full_lon_ibeg:mesh%full_lon_iend,j) = pole_gzs / pole_n
          std(mesh%full_lon_ibeg:mesh%full_lon_iend,j) = (pole_std - 2 * pole_gzs**2 * pole_n + pole_gzs**2) / pole_n / g
          lnd(mesh%full_lon_ibeg:mesh%full_lon_iend,j) = pole_lnd / pole_n
        end if
      end do
      call fill_halo(block, gzs, full_lon=.true., full_lat=.true.)
      call fill_halo(block, std, full_lon=.true., full_lat=.true.)
      call fill_halo(block, lnd, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine topo_regrid

  subroutine topo_smooth(block)

    type(block_type), intent(inout) :: block

    real(r8) wgt
    integer i, j, cyc

    associate (mesh     => block%mesh        , &
               gzs      => block%static%gzs  , &
               gzs_f    => block%state(1)%phs, & ! Borrow the array.
               landmask => block%static%landmask)
    do cyc = 1, topo_smooth_cycles
      call filter_on_cell(block, gzs, gzs_f)
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (abs(mesh%full_lat_deg(j)) > 60) then
          wgt = sin(pi05 * (1 - (pi05 - abs(mesh%full_lat(j))) / (30 * rad)))
          gzs(:,j) = wgt * gzs_f(:,j) + (1 - wgt) * gzs(:,j)
        end if
      end do
      call fill_halo(block, gzs, full_lon=.true., full_lat=.true.)
      ! do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      !   do i = mesh%full_lon_ibeg, mesh%full_lon_iend
      !     if (landmask(i,j) == 0) gzs(i,j) = 0
      !   end do
      ! end do
    end do
    ! call laplace_damp_on_cell(block, 8, gzs, fill=.true.)
    end associate

  end subroutine topo_smooth

  subroutine topo_final()

    call topo %clear()
    call topo1%clear()
    call topo2%clear()

  end subroutine topo_final

end module topo_mod
