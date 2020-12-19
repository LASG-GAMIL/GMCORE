program gmcore_prepare

  use fiona
  use topo_mod
  use bkg_mod
  use mesh_mod
  use process_mod
  use block_mod
  use vert_coord_mod
  use initial_mod

  implicit none

  character(256) namelist_file
  character(256) topo_file
  character(256) bkg_file
  character(256) initial_file
  character(30) bkg_type
  character(30) initial_time
  character(30) vert_coord_template
  integer num_lon
  integer num_lat
  integer num_lev

  integer iblk

  namelist /gmcore_prepare_params/ &
    topo_file                    , &
    bkg_type                     , &
    bkg_file                     , &
    initial_file                 , &
    num_lon                      , &
    num_lat                      , &
    num_lev                      , &
    vert_coord_template

  call get_command_argument(1, namelist_file)

  open(10, file=namelist_file)
  read(10, nml=gmcore_prepare_params)
  close(10)

  if (initial_time == '') then
    initial_time = '1970-01-01'
  end if

  call fiona_init(start_time=initial_time, time_units='hours')

  call global_mesh%init_global(num_lon, num_lat, num_lev)
  call process_init()
  call vert_coord_init(num_lev, template=vert_coord_template)
  call process_create_blocks()

  call topo_read(topo_file)

  do iblk = 1, size(proc%blocks)
    call topo_regrid(proc%blocks(iblk))
  end do

  call bkg_read(bkg_type, bkg_file)

  call initial_write(initial_file)

  call topo_final()
  call bkg_final()
  call process_final()

end program gmcore_prepare
