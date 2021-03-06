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
  use prepare_mod

  implicit none

  character(256) namelist_file
  character(30) :: initial_time = '1970-01-01T00:00:00'

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
    output_group_size            , &
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

  time_scheme = 'N/A'

  call fiona_init()

  call global_mesh%init_global(num_lon, num_lat, num_lev, lon_halo_width=5, lat_halo_width=5)
  call process_init()
  call vert_coord_init(num_lev, scheme=vert_coord_scheme, template=vert_coord_template)
  call process_create_blocks()

  if (is_root_proc()) then
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
    write(*, *) 'output_group_size    = ', to_str(output_group_size)
    write(*, *) 'initial_time         = ', trim(initial_time)
    write(*, *) 'namelist_file        = ', trim(namelist_file)
    write(*, *) 'topo_file            = ', trim(topo_file)
    write(*, *) 'bkg_file             = ', trim(bkg_file)
    write(*, *) 'initial_file         = ', trim(initial_file)
    write(*, *) 'bkg_type             = ', trim(bkg_type)
    write(*, *) '========================================================='
  end if

  call prepare_run(topo_file, bkg_file, bkg_type, &
                   zero_min_lon, zero_max_lon, zero_min_lat, zero_max_lat, &
                   smth_min_lon, smth_max_lon, smth_min_lat, smth_max_lat, &
                   smth_steps)

  call initial_write(initial_file, initial_time)

  call prepare_final()

  call process_final()

end program gmcore_prepare
