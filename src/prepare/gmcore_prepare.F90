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

  integer iblk

  namelist /gmcore_prepare_params/ &
    initial_time                 , &
    topo_file                    , &
    bkg_type                     , &
    bkg_file                     , &
    initial_file                 , &
    num_lon                      , &
    num_lat                      , &
    num_lev                      , &
    coarse_polar_lat0            , &
    coarse_polar_decay           , &
    vert_coord_scheme            , &
    vert_coord_template

  call get_command_argument(1, namelist_file)

  open(10, file=namelist_file)
  read(10, nml=gmcore_prepare_params)
  close(10)

  write(*, *) '=================== GMCORE Parameters ==================='
  write(*, *) 'num_lon              = ', to_str(num_lon)
  write(*, *) 'num_lat              = ', to_str(num_lat)
  write(*, *) 'num_lev              = ', to_str(num_lev)
  write(*, *) 'vert_coord_scheme    = ', trim(vert_coord_scheme)
  write(*, *) 'vert_coord_template  = ', trim(vert_coord_template)
  write(*, *) 'initial_time         = ', trim(initial_time)
  write(*, *) 'namelist_file        = ', trim(namelist_file)
  write(*, *) 'topo_file            = ', trim(topo_file)
  write(*, *) 'bkg_file             = ', trim(bkg_file)
  write(*, *) 'initial_file         = ', trim(initial_file)
  write(*, *) 'bkg_type             = ', trim(bkg_type)
  write(*, *) '========================================================='

  call fiona_init(start_time=initial_time, time_units='hours')

  call global_mesh%init_global(num_lon, num_lat, num_lev)
  call process_init()
  call vert_coord_init(num_lev, scheme=vert_coord_scheme, template=vert_coord_template)
  call process_create_blocks()

  call topo_read(topo_file)

  do iblk = 1, size(proc%blocks)
    call topo_regrid(proc%blocks(iblk))
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
