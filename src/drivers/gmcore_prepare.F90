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
  use damp_mod
  use prepare_mod

  implicit none

  character(256) namelist_file

  call get_command_argument(1, namelist_file)

  call parse_namelist(namelist_file)

  time_scheme = 'N/A'

  call fiona_init()

  call const_init(planet)
  call global_mesh%init_global(num_lon, num_lat, num_lev, lon_halo_width=2, lat_halo_width=2)
  call process_init()
  call vert_coord_init(num_lev, scheme=vert_coord_scheme, template=vert_coord_template)
  call process_create_blocks()
  call damp_init(proc%blocks)

  if (is_root_proc()) then
    write(*, *) '=================== GMCORE Parameters ==================='
    write(*, *) 'nonhydrostatic       = ', to_str(nonhydrostatic)
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
    write(*, *) 'use_topo_smooth      = ', to_str(use_topo_smooth)
    if (use_topo_smooth) then
    write(*, *) 'topo_smooth_cycles   = ', to_str(topo_smooth_cycles)
    end if
    write(*, *) 'bkg_file             = ', trim(bkg_file)
    write(*, *) 'initial_file         = ', trim(initial_file)
    write(*, *) 'bkg_type             = ', trim(bkg_type)
    write(*, *) '========================================================='
  end if

  call prepare_run()

  call initial_write(initial_file, initial_time)

  call prepare_final()

  call process_final()

end program gmcore_prepare
