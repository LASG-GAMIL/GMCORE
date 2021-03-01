program gmcore_prepare

  use fiona
  use string
  use topo_mod
  use bkg_mod
  use mesh_mod
  use process_mod
  use block_mod
  use vert_coord_mod
  use initial_mod
  use namelist_mod

  implicit none

  character(256) namelist_file
  character(256) topo_file
  character(256) bkg_file
  character(30) :: bkg_type = 'era5'
  character(30) :: initial_time = '1970-01-01'

  real(r8) :: zero_min_lon(100) = -1.0e33
  real(r8) :: zero_max_lon(100) = -1.0e33
  real(r8) :: zero_min_lat(100) = -1.0e33
  real(r8) :: zero_max_lat(100) = -1.0e33
  real(r8) :: smth_min_lon(100) = -1.0e33
  real(r8) :: smth_max_lon(100) = -1.0e33
  real(r8) :: smth_min_lat(100) = -1.0e33
  real(r8) :: smth_max_lat(100) = -1.0e33
  integer :: smth_steps(100) = 1

  integer iblk, i

  namelist /gmcore_prepare_params/ &
    initial_time                 , &
    topo_file                    , &
    bkg_type                     , &
    bkg_file                     , &
    initial_file                 , &
    num_lon                      , &
    num_lat                      , &
    num_lev                      , &
    coarse_pole_mul              , &
    coarse_pole_decay            , &
    vert_coord_scheme            , &
    vert_coord_template          , &
    zero_min_lon                 , &
    zero_max_lon                 , &
    zero_min_lat                 , &
    zero_max_lat                 , &
    smth_min_lon                 , &
    smth_max_lon                 , &
    smth_min_lat                 , &
    smth_max_lat                 , &
    smth_steps

  call get_command_argument(1, namelist_file)

  open(10, file=namelist_file)
  read(10, nml=gmcore_prepare_params)
  close(10)

  write(*, *) '=================== GMCORE Parameters ==================='
  write(*, *) 'num_lon              = ', to_str(num_lon)
  write(*, *) 'num_lat              = ', to_str(num_lat)
  write(*, *) 'num_lev              = ', to_str(num_lev)
  if (coarse_pole_mul /= 0) then
  write(*, *) 'coarse_pole_mul      = ', to_str(coarse_pole_mul, 2)
  write(*, *) 'coarse_pole_decay    = ', to_str(coarse_pole_decay, 2)
  end if
  write(*, *) 'vert_coord_scheme    = ', trim(vert_coord_scheme)
  write(*, *) 'vert_coord_template  = ', trim(vert_coord_template)
  write(*, *) 'initial_time         = ', trim(initial_time)
  write(*, *) 'namelist_file        = ', trim(namelist_file)
  write(*, *) 'topo_file            = ', trim(topo_file)
  write(*, *) 'bkg_file             = ', trim(bkg_file)
  write(*, *) 'initial_file         = ', trim(initial_file)
  write(*, *) 'bkg_type             = ', trim(bkg_type)
  write(*, *) '========================================================='

  time_scheme = 'N/A'

  call fiona_init(start_time=initial_time, time_units='hours')

  call global_mesh%init_global(num_lon, num_lat, num_lev, lon_halo_width=5, lat_halo_width=5)
  call process_init()
  call vert_coord_init(num_lev, scheme=vert_coord_scheme, template=vert_coord_template)
  call process_create_blocks()

  call topo_read(topo_file)

  do iblk = 1, size(proc%blocks)
    call topo_regrid(proc%blocks(iblk))
    do i = 1, size(zero_min_lon)
      if (zero_min_lon(i) /= -1.0e33 .and. zero_max_lon(i) /= -1.0e33 .and. &
          zero_min_lat(i) /= -1.0e33 .and. zero_max_lat(i) /= -1.0e33) then
        call topo_zero(proc%blocks(iblk), zero_min_lon(i), zero_max_lon(i), zero_min_lat(i), zero_max_lat(i))
      end if
    end do
    do i = 1, size(smth_min_lon)
      if (smth_min_lon(i) /= -1.0e33 .and. smth_max_lon(i) /= -1.0e33 .and. &
          smth_min_lat(i) /= -1.0e33 .and. smth_max_lat(i) /= -1.0e33) then
        call topo_smth(proc%blocks(iblk), smth_min_lon(i), smth_max_lon(i), smth_min_lat(i), smth_max_lat(i), smth_steps(i))
      end if
    end do
  end do

  call bkg_read(bkg_type, bkg_file)

  call bkg_regrid_phs()
  call bkg_calc_ph()
  call bkg_regrid_pt()
  call bkg_regrid_u()
  call bkg_regrid_v()

  call initial_write(initial_file)

  call topo_final()
  call bkg_final()
  call process_final()

end program gmcore_prepare
