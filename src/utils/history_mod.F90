module history_mod

  use container
  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod, dt => dt_dyn
  use parallel_mod
  use allocator_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use block_mod
  use diag_state_mod
  use adv_batch_mod

  implicit none

  private

  public history_init
  public history_final
  public history_setup_h0_adv
  public history_write_h0
  public history_write_h1

  character(8), parameter ::    cell_dims_2d(3) = ['lon ', 'lat ',         'time']
  character(8), parameter ::    cell_dims_3d(4) = ['lon ', 'lat ', 'lev ', 'time']
  character(8), parameter ::     vtx_dims_2d(3) = ['ilon', 'ilat',         'time']
  character(8), parameter ::     vtx_dims_3d(4) = ['ilon', 'ilat', 'lev ', 'time']
  character(8), parameter ::     lon_dims_2d(3) = ['ilon', 'lat ',         'time']
  character(8), parameter ::     lon_dims_3d(4) = ['ilon', 'lat ', 'lev ', 'time']
  character(8), parameter ::     lat_dims_2d(3) = ['lon ', 'ilat',         'time']
  character(8), parameter ::     lat_dims_3d(4) = ['lon ', 'ilat', 'lev ', 'time']
  character(8), parameter ::     lev_dims_3d(4) = ['lon ', 'lat ', 'ilev', 'time']
  character(8), parameter :: lon_lev_dims_3d(4) = ['ilon', 'lat ', 'ilev', 'time']
  character(8), parameter :: lat_lev_dims_3d(4) = ['lon ', 'ilat', 'ilev', 'time']

contains

  subroutine history_init()

    character(10) time_value, time_units
    character(4) diag_dims(4)
    real(8) seconds, months

    if (history_interval(1) == 'N/A') call log_error('Parameter history_interval is not set!')
    if (case_name == 'N/A') call log_error('Parameter case_name is not set!')

    time_value = split_string(history_interval(1), ' ', 1)
    time_units = split_string(history_interval(1), ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      if (is_root_proc()) call log_error('Invalid history interval ' // trim(history_interval(1)) // '!')
    end select

    call time_add_alert('history_write', seconds=seconds)

    if (output_h0_new_file == '') then
      call time_add_alert('h0_new_file', seconds=seconds)
      if (is_root_proc()) call log_notice('Output data every ' // trim(history_interval(1)) // '.')
    else if (output_h0_new_file == 'one_file') then
      if (is_root_proc()) call log_notice('Output data in one file.')
    else
      time_value = split_string(output_h0_new_file, ' ', 1)
      time_units = split_string(output_h0_new_file, ' ', 2)
      if (time_units == 'months') then
        read(time_value, *) months
        if (is_root_proc()) call log_notice('Output data every ' // trim(time_value) // ' months.')
        call time_add_alert('h0_new_file', months=months)
      else
        read(time_value, *) seconds
        select case (time_units)
        case ('days')
          if (is_root_proc()) call log_notice('Output data every ' // trim(time_value) // ' days.')
          seconds = seconds * 86400
        case ('hours')
          if (is_root_proc()) call log_notice('Output data every ' // trim(time_value) // ' hours.')
          seconds = seconds * 3600
        case ('minutes')
          if (is_root_proc()) call log_notice('Output data every ' // trim(time_value) // ' minutes.')
          seconds = seconds * 60
        case ('seconds')
          if (is_root_proc()) call log_notice('Output data every ' // trim(time_value) // ' seconds.')
          seconds = seconds
        case default
          if (is_root_proc()) call log_error('Invalid output_h0_new_file ' // trim(output_h0_new_file) // '!')
        end select
        call time_add_alert('h0_new_file', seconds=seconds)
      end if
    end if

    call fiona_init(time_units, start_time_str)

    if (hydrostatic) then
      call history_setup_h0_hydrostatic()
      call history_setup_h1_hydrostatic()
    else if (nonhydrostatic) then
      call history_setup_h0_nonhydrostatic()
      call history_setup_h1_nonhydrostatic()
    else if (.not. advection) then
      call history_setup_h0_swm()
      call history_setup_h1_swm()
    end if

  end subroutine history_init

  subroutine history_final()

    call fiona_final()

  end subroutine history_final

  subroutine history_setup_h0_swm()

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size, split_file=split_h0)
    ! Dimensions
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    ! Variables
    call fiona_add_var('h0', 'tm'   , long_name='total mass'               , units='m'      , dim_names=['time'])
    call fiona_add_var('h0', 'te'   , long_name='total energy'             , units='m4 s-4' , dim_names=['time']    , dtype='r8')
    call fiona_add_var('h0', 'tpe'  , long_name='total potential enstrophy', units='m2 s-5' , dim_names=['time']    , dtype='r8')
    call fiona_add_var('h0', 'u'    , long_name='u wind component'         , units='m s-1'  , dim_names= lon_dims_2d, dtype='r8')
    call fiona_add_var('h0', 'v'    , long_name='v wind component'         , units='m s-1'  , dim_names= lat_dims_2d, dtype='r8')
    call fiona_add_var('h0', 'z'    , long_name='height'                   , units='m'      , dim_names=cell_dims_2d, dtype='r8')
    call fiona_add_var('h0', 'pv'   , long_name='potential vorticity'      , units='s-1 m-1', dim_names= vtx_dims_2d, dtype='r8')
    call fiona_add_var('h0', 'div'  , long_name='divergence'               , units='s-1'    , dim_names=cell_dims_2d, dtype='r8')
    if (test_case == 'stratospheric_vortex_erosion') then
      call fiona_add_var('h0', 'zs', long_name='surface height', units='m', dim_names=cell_dims_2d)
    else
      call fiona_add_var('h0', 'zs', long_name='surface height', units='m', dim_names=['lon', 'lat'])
    end if

  end subroutine history_setup_h0_swm

  subroutine history_setup_h0_adv(blocks)

    type(block_type), intent(in) :: blocks(:)

    integer i, j

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size, split_file=split_h0)
    ! Dimensions
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
    ! Variables
    associate (adv_batches => blocks(1)%adv_batches)
    do i = 1, size(adv_batches)
      do j = 1, size(adv_batches(i)%tracer_names)
        call fiona_add_var('h0', adv_batches(i)%tracer_names(j), &
                           long_name=adv_batches(i)%tracer_long_names(j), &
                           units=adv_batches(i)%tracer_units(j), &
                           dim_names=cell_dims_3d, dtype='r8')
      end do
      call fiona_add_var('h0', 'u'   , long_name='', units='', dim_names= lon_dims_3d)
      call fiona_add_var('h0', 'v'   , long_name='', units='', dim_names= lat_dims_3d)
      call fiona_add_var('h0', 'cflx', long_name='', units='', dim_names= lon_dims_3d)
      call fiona_add_var('h0', 'cfly', long_name='', units='', dim_names= lat_dims_3d)
      call fiona_add_var('h0', 'divx', long_name='', units='', dim_names=cell_dims_3d)
      call fiona_add_var('h0', 'divy', long_name='', units='', dim_names=cell_dims_3d)
      if (global_mesh%num_full_lev > 1) then
        call fiona_add_var('h0', 'z'    , long_name='', units='', dim_names=cell_dims_3d)
        call fiona_add_var('h0', 'z_lev', long_name='', units='', dim_names= lev_dims_3d)
        call fiona_add_var('h0', 'w'    , long_name='', units='', dim_names= lev_dims_3d)
        call fiona_add_var('h0', 'cflz' , long_name='', units='', dim_names= lev_dims_3d)
      end if
    end do
    end associate

  end subroutine history_setup_h0_adv

  subroutine history_setup_h0_hydrostatic()

    integer i, j

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size, split_file=split_h0)
    ! Dimensions
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
    ! Variables
    call fiona_add_var('h0', 'tm'     , long_name='total mass'                  , units='m'     , dim_names=['time'])
    call fiona_add_var('h0', 'te'     , long_name='total energy'                , units='m4 s-4', dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'tpe'    , long_name='total potential enstrophy'   , units='m2 s-5', dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'te_ke'  , long_name='total kinetic energy'        , units=''      , dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'te_ie'  , long_name='total internal energy'       , units=''      , dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'te_pe'  , long_name='total potential energy'      , units=''      , dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'ref_ps' , long_name='reference surface pressure'  , units='Pa'    , dim_names=['lon', 'lat'], dtype='r8')
    call fiona_add_var('h0', 'zs'     , long_name='surface height'              , units='m'     , dim_names=['lon', 'lat'], dtype='r8')
    call fiona_add_var('h0', 'dzsdlon', long_name='zonal zs gradient'           , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'dzsdlat', long_name='meridional zs gradient'      , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'phs'    , long_name='surface hydrostatic pressure', units='Pa'    , dim_names=cell_dims_2d  , dtype='r8')
    call fiona_add_var('h0', 'u'      , long_name='u wind component'            , units='m s-1' , dim_names= lon_dims_3d  , dtype='r8')
    call fiona_add_var('h0', 'v'      , long_name='v wind component'            , units='m s-1' , dim_names= lat_dims_3d  , dtype='r8')
    call fiona_add_var('h0', 'pt'     , long_name='potential temperature'       , units='K'     , dim_names=cell_dims_3d  , dtype='r8')
    call fiona_add_var('h0', 't'      , long_name='temperature'                 , units='K'     , dim_names=cell_dims_3d  , dtype='r8')
    call fiona_add_var('h0', 'z'      , long_name='height'                      , units='m'     , dim_names=cell_dims_3d  , dtype='r8')
    call fiona_add_var('h0', 'ph'     , long_name='hydrostatic pressure'        , units='Pa'    , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'vor'    , long_name='relative vorticity'          , units='s-1'   , dim_names= vtx_dims_3d)
    call fiona_add_var('h0', 'div'    , long_name='divergence'                  , units='s-1'   , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'landmask', long_name='land mask'                  , units=''      , dim_names=['lon', 'lat'])
    associate (adv_batches => blocks(1)%adv_batches)
    do i = 1, size(adv_batches)
      do j = 1, size(adv_batches(i)%tracer_names)
        call fiona_add_var('h0', adv_batches(i)%tracer_names(j), &
                           long_name=adv_batches(i)%tracer_long_names(j), &
                           units=adv_batches(i)%tracer_units(j), &
                           dim_names=cell_dims_3d, dtype='r8')
      end do
    end do
    end associate

  end subroutine history_setup_h0_hydrostatic

  subroutine history_setup_h0_nonhydrostatic()

    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size, split_file=split_h0)
    ! Dimensions
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h0', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
    ! Variables
    call fiona_add_var('h0', 'tm'     , long_name='total mass'                  , units='m'     , dim_names=['time'])
    call fiona_add_var('h0', 'te'     , long_name='total energy'                , units='m4 s-4', dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'tpe'    , long_name='total potential enstrophy'   , units='m2 s-5', dim_names=['time']      , dtype='r8')
    call fiona_add_var('h0', 'zs'     , long_name='surface height'              , units='m'     , dim_names=['lon', 'lat'], dtype='r8')
    call fiona_add_var('h0', 'dzsdlon', long_name='zonal zs gradient'           , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'dzsdlat', long_name='meridional zs gradient'      , units=''      , dim_names=['lon', 'lat'])
    call fiona_add_var('h0', 'phs'    , long_name='surface hydrostatic pressure', units='Pa'    , dim_names=cell_dims_2d, dtype='r8')
    call fiona_add_var('h0', 'u'      , long_name='u wind component'            , units='m s-1' , dim_names= lon_dims_3d, dtype='r8')
    call fiona_add_var('h0', 'v'      , long_name='v wind component'            , units='m s-1' , dim_names= lat_dims_3d, dtype='r8')
    call fiona_add_var('h0', 'pt'     , long_name='potential temperature'       , units='K'     , dim_names=cell_dims_3d, dtype='r8')
    call fiona_add_var('h0', 't'      , long_name='temperature'                 , units='K'     , dim_names=cell_dims_3d, dtype='r8')
    call fiona_add_var('h0', 'z'      , long_name='height'                      , units='m'     , dim_names=cell_dims_3d, dtype='r8')
    call fiona_add_var('h0', 'ph'     , long_name='hydrostatic pressure'        , units='Pa'    , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'vor'    , long_name='relative vorticity'          , units='s-1'   , dim_names= vtx_dims_3d)
    call fiona_add_var('h0', 'div'    , long_name='divergence'                  , units='s-1'   , dim_names=cell_dims_3d)
    call fiona_add_var('h0', 'w'      , long_name='vertical speed'              , units='m s-1' , dim_names= lev_dims_3d)
    call fiona_add_var('h0', 'p'      , long_name='full pressure'               , units='Pa'    , dim_names= lev_dims_3d)
    call fiona_add_var('h0', 'rhod'   , long_name='dry air density'             , units='kg m-3', dim_names=cell_dims_3d)

  end subroutine history_setup_h0_nonhydrostatic

  subroutine history_setup_h1_swm()

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size)
    ! Dimensions
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    ! Variables
    call fiona_add_var('h1', 'dudt'    , long_name='u wind component tendency'                     , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'dvdt'    , long_name='v wind component tendency'                     , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'dgzdt'   , long_name='geopotential tendency'                         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'qhv'     , long_name='nonliear zonal Coriolis force'                 , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'qhu'     , long_name='nonliear meridional Coriolis force'            , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'pgf_lon' , long_name='zonal geopotential energy gradient force'      , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'dkedlon' , long_name='zonal kinetic energy gradient force'           , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'pgf_lat' , long_name='meridional geopotential energy gradient force' , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'dkedlat' , long_name='meridional kinetic energy gradient force'      , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'dmfdlon' , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dmfdlat' , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'mfx_lon' , long_name='normal mass flux on U grid'                    , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'mfy_lon' , long_name='tangent mass flux on U grid'                   , units='', dim_names= lon_dims_2d)
    call fiona_add_var('h1', 'mfy_lat' , long_name='normal mass flux on V grid'                    , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'mfx_lat' , long_name='tangent mass flux on V grid'                   , units='', dim_names= lat_dims_2d)
    call fiona_add_var('h1', 'ke'      , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_2d)

  end subroutine history_setup_h1_swm

  subroutine history_setup_h1_hydrostatic()

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size)
    ! Dimensions
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
    ! Variables
    call fiona_add_var('h1', 'dudt'         , long_name='u wind component tendency'                     , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt'         , long_name='v wind component tendency'                     , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dphsdt'       , long_name='surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dptdt'        , long_name='potential temperature tendency'                , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dptfdlon'     , long_name='zonal potential temperature flux gradient'     , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dptfdlat'     , long_name='meridional potential temperature flux gradient', units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dptfdlev'     , long_name='vertical potential temperature flux gradient'  , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'we_lev'       , long_name='vertical coordinate velocity'                  , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'pgf_lon'      , long_name='zonal pressure gradient force'                 , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'wedudlev'     , long_name='vertical advection of u'                       , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'pgf_lat'      , long_name='meridional pressure gradient force'            , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'wedvdlev'     , long_name='vertical advection of v'                       , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'qhv'          , long_name='nonliear zonal Coriolis force'                 , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'qhu'          , long_name='nonliear meridional Coriolis force'            , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dkedlon'      , long_name='zonal kinetic energy gradient force'           , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dkedlat'      , long_name='meridional kinetic energy gradient force'      , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmfdlon'      , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dmfdlat'      , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'mfx_lon'      , long_name='normal mass flux on U grid'                    , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'mfy_lon'      , long_name='tangent mass flux on U grid'                   , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'mfy_lat'      , long_name='normal mass flux on V grid'                    , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'mfx_lat'      , long_name='tangent mass flux on V grid'                   , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'm'            , long_name='dph on full levels'                            , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'ke'           , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_3d)
    if (use_smag_damp) then
      call fiona_add_var('h1', 'smag_t'     , long_name='Tension for Smagorinsky diffusion'             , units='', dim_names=cell_dims_3d)
      call fiona_add_var('h1', 'smag_s'     , long_name='Shear for Smagorinsky diffusion'               , units='', dim_names= vtx_dims_3d)
      call fiona_add_var('h1', 'kmh_lon'    , long_name='Smagorinsky diffusion coefficient for u'       , units='', dim_names= lon_dims_3d)
      call fiona_add_var('h1', 'kmh_lat'    , long_name='Smagorinsky diffusion coefficient for v'       , units='', dim_names= lat_dims_3d)
    end if

  end subroutine history_setup_h1_hydrostatic

  subroutine history_setup_h1_nonhydrostatic

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size)
    ! Dimensions
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
    ! Variables
    call fiona_add_var('h1', 'dudt'         , long_name='u wind component tendency'                     , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dvdt'         , long_name='v wind component tendency'                     , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dphsdt'       , long_name='surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d)
    call fiona_add_var('h1', 'dptdt'        , long_name='potential temperature tendency'                , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dptfdlon'     , long_name='zonal potential temperature flux gradient'     , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dptfdlat'     , long_name='meridional potential temperature flux gradient', units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dptfdlev'     , long_name='vertical potential temperature flux gradient'  , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'we_lev', long_name='vertical coordinate velocity'                  , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'pgf_lon'      , long_name='zonal pressure gradient force'                 , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'wedudlev'     , long_name='vertical advection of u'                       , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'pgf_lat'      , long_name='meridional pressure gradient force'            , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'wedvdlev'     , long_name='vertical advection of v'                       , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'qhv'          , long_name='nonliear zonal Coriolis force'                 , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'qhu'          , long_name='nonliear meridional Coriolis force'            , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dkedlon'      , long_name='zonal kinetic energy gradient force'           , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'dkedlat'      , long_name='meridional kinetic energy gradient force'      , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'dmfdlon'      , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'dmfdlat'      , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'mfx_lon'      , long_name='normal mass flux on U grid'                    , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'mfy_lon'      , long_name='tangent mass flux on U grid'                   , units='', dim_names= lon_dims_3d)
    call fiona_add_var('h1', 'mfy_lat'      , long_name='normal mass flux on V grid'                    , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'mfx_lat'      , long_name='tangent mass flux on V grid'                   , units='', dim_names= lat_dims_3d)
    call fiona_add_var('h1', 'm'            , long_name='dph on full levels'                            , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'ke'           , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_3d)
    call fiona_add_var('h1', 'adv_gz_lon'   , long_name='advection of gz'                               , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'adv_gz_lat'   , long_name='advection of gz'                               , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'adv_gz_lev'   , long_name='advection of gz'                               , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'adv_w_lon'    , long_name='advection of w'                                , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'adv_w_lat'    , long_name='advection of w'                                , units='', dim_names= lev_dims_3d)
    call fiona_add_var('h1', 'adv_w_lev'    , long_name='advection of w'                                , units='', dim_names= lev_dims_3d)

  end subroutine history_setup_h1_nonhydrostatic

  subroutine history_write_h0_swm(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime 

    integer iblk, is, ie, js, je
    integer start(2), count(2)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh        , &
                 state  => blocks(iblk)%state(itime), &
                 static => blocks(iblk)%static)
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      start = [is,js]
      count = [mesh%num_full_lon,mesh%num_full_lat]
      call fiona_output('h0', 'zs' , static%gzs(is:ie,js:je  ) / g, start=start, count=count)
      call fiona_output('h0', 'z'  ,  state%gz (is:ie,js:je,1) / g, start=start, count=count)
      call fiona_output('h0', 'div',  state%div(is:ie,js:je,1)    , start=start, count=count)
      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      start = [is,js]
      count = [mesh%num_half_lon,mesh%num_full_lat]
      call fiona_output('h0', 'u'  ,  state%u_lon(is:ie,js:je,1)  , start=start, count=count)
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      start = [is,js]
      count = [mesh%num_full_lon,mesh%num_half_lat]
      call fiona_output('h0', 'v'  ,  state%v_lat(is:ie,js:je,1)  , start=start, count=count)
      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      start = [is,js]
      count = [mesh%num_half_lon,mesh%num_half_lat]
      call fiona_output('h0', 'pv' ,  state%pv (is:ie,js:je,1)    , start=start, count=count)
 
      call fiona_output('h0', 'tm' , state %tm)
      call fiona_output('h0', 'te' , state %te)
      call fiona_output('h0', 'tpe', state %tpe)
      end associate
    end do

  end subroutine history_write_h0_swm

  subroutine history_write_h0_adv(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime 

    integer iblk, is, ie, js, je, ks, ke, i, j
    integer start(3), count(3)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    call fiona_output('h0', 'lev' , global_mesh%full_lev(1:global_mesh%num_full_lev))
    call fiona_output('h0', 'ilev', global_mesh%half_lev(1:global_mesh%num_half_lev))

    do iblk = 1, size(blocks)
      associate (mesh        => blocks(iblk)%mesh       , &
                 adv_batches => blocks(iblk)%adv_batches, &
                 state       => blocks(iblk)%state(itime))
      do i = 1, size(adv_batches)
        do j = 1, size(adv_batches(i)%tracer_names)
          is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
          js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
          ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
          start = [is,js,ks]
          count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
          call fiona_output('h0', adv_batches(i)%tracer_names(j), &
                            adv_batches(i)%q(is:ie,js:je,ks:ke,j,adv_batches(i)%old), &
                            start=start, count=count)
        end do
        is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
        js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
        call fiona_output('h0', 'divx', adv_batches(i)%divx(is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output('h0', 'divy', adv_batches(i)%divy(is:ie,js:je,ks:ke), start=start, count=count)
        is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
        js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]
        call fiona_output('h0', 'u'   , adv_batches(i)%u   (is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output('h0', 'cflx', adv_batches(i)%cflx(is:ie,js:je,ks:ke), start=start, count=count)
        is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
        js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
        ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
        start = [is,js,ks]
        count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]
        call fiona_output('h0', 'v'   , adv_batches(i)%v   (is:ie,js:je,ks:ke), start=start, count=count)
        call fiona_output('h0', 'cfly', adv_batches(i)%cfly(is:ie,js:je,ks:ke), start=start, count=count)
        if (global_mesh%num_full_lev > 1) then
          is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
          js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
          ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
          start = [is,js,ks]
          count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
          call fiona_output('h0', 'z'    , state%gz    (is:ie,js:je,ks:ke) / g, start=start, count=count)
          is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
          js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
          ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
          start = [is,js,ks]
          count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]
          call fiona_output('h0', 'z_lev', state%gz_lev(is:ie,js:je,ks:ke) / g, start=start, count=count)
          call fiona_output('h0', 'w'    , adv_batches(i)%we  (is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output('h0', 'cflz' , adv_batches(i)%cflz(is:ie,js:je,ks:ke), start=start, count=count)
        end if
      end do
      end associate
    end do

  end subroutine history_write_h0_adv

  subroutine history_write_h0_hydrostatic(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime 

    integer iblk, is, ie, js, je, ks, ke, i, j
    integer start(3), count(3)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    call fiona_output('h0', 'lev' , global_mesh%full_lev(1:global_mesh%num_full_lev))
    call fiona_output('h0', 'ilev', global_mesh%half_lev(1:global_mesh%num_half_lev))

    do iblk = 1, size(blocks)
      associate (mesh        => blocks(iblk)%mesh        , &
                 state       => blocks(iblk)%state(itime), &
                 static      => blocks(iblk)%static      , &
                 adv_batches => blocks(iblk)%adv_batches )
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
      call fiona_output('h0', 'ref_ps'  , static%ref_ps (is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'zs'      , static%gzs    (is:ie,js:je) / g       , start=start, count=count)
      call fiona_output('h0', 'dzsdlon' , static%dzsdlon(is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'dzsdlat' , static%dzsdlat(is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'landmask', static%landmask(is:ie,js:je)          , start=start, count=count)
      call fiona_output('h0', 'z'       ,  state%gz     (is:ie,js:je,ks:ke) / g , start=start, count=count)
      call fiona_output('h0', 'phs'     ,  state%phs    (is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'ph'      ,  state%ph     (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'pt'      ,  state%pt     (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 't'       ,  state%t      (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'div'     ,  state%div    (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]
      call fiona_output('h0', 'u'       ,  state%u_lon  (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]
      call fiona_output('h0', 'v'       ,  state%v_lat  (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_half_lat,mesh%num_full_lev]
      call fiona_output('h0', 'vor'     ,  state%vor    (is:ie,js:je,ks:ke)     , start=start, count=count)
      
      call fiona_output('h0', 'tm'   , state%tm)
      call fiona_output('h0', 'te'   , state%te)
      call fiona_output('h0', 'tpe'  , state%tpe)
      call fiona_output('h0', 'te_ke', state%te_ke)
      call fiona_output('h0', 'te_ie', state%te_ie)
      call fiona_output('h0', 'te_pe', state%te_pe)

      do i = 1, size(adv_batches)
        do j = 1, size(adv_batches(i)%tracer_names)
          is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
          js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
          ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
          start = [is,js,ks]
          count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
          call fiona_output('h0', adv_batches(i)%tracer_names(j), &
                            adv_batches(i)%q(is:ie,js:je,ks:ke,j,adv_batches(i)%old), &
                            start=start, count=count)
        end do
      end do
      end associate
    end do

  end subroutine history_write_h0_hydrostatic

  subroutine history_write_h0_nonhydrostatic(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime 

    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)

    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    call fiona_output('h0', 'lev' , global_mesh%full_lev(1:global_mesh%num_full_lev))
    call fiona_output('h0', 'ilev', global_mesh%half_lev(1:global_mesh%num_half_lev))

    do iblk = 1, size(blocks)
      associate (mesh   => blocks(iblk)%mesh        , &
                 state  => blocks(iblk)%state(itime), &
                 static => blocks(iblk)%static)
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
      call fiona_output('h0', 'zs'      , static%gzs    (is:ie,js:je) / g       , start=start, count=count)
      call fiona_output('h0', 'dzsdlon' , static%dzsdlon(is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'dzsdlat' , static%dzsdlat(is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'z'       ,  state%gz     (is:ie,js:je,ks:ke) / g , start=start, count=count)
      call fiona_output('h0', 'phs'     ,  state%phs    (is:ie,js:je)           , start=start, count=count)
      call fiona_output('h0', 'ph'      ,  state%ph     (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'pt'      ,  state%pt     (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 't'       ,  state%t      (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'div'     ,  state%div    (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'rhod'    ,  state%rhod   (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]
      call fiona_output('h0', 'u'       ,  state%u_lon  (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]
      call fiona_output('h0', 'v'       ,  state%v_lat  (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_half_lat,mesh%num_full_lev]
      call fiona_output('h0', 'vor'     ,  state%vor    (is:ie,js:je,ks:ke)     , start=start, count=count)
      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]
      call fiona_output('h0', 'w'       , state%w_lev   (is:ie,js:je,ks:ke)     , start=start, count=count)
      call fiona_output('h0', 'p'       , state%p_lev   (is:ie,js:je,ks:ke)     , start=start, count=count)
      
      call fiona_output('h0', 'tm' , state %tm)
      call fiona_output('h0', 'te' , state %te)
      call fiona_output('h0', 'tpe', state %tpe)
      end associate
    end do

  end subroutine history_write_h0_nonhydrostatic

  subroutine history_write_h0(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    real(8) time1, time2

    if (.not. time_has_alert('h0_new_file')) then
      if (restart) then
        restart = .false.
        call log_notice('Skip write state at restart step, and set restart to false.', pid=proc%id)
        return ! Skip the restart step.
      end if
      call fiona_start_output('h0', elapsed_seconds, new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h0', elapsed_seconds, new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h0', elapsed_seconds, new_file=.false.)
    end if

    if (is_root_proc()) then
      call log_notice('Write state.')
      call cpu_time(time1)
    end if

    if (advection) then
      call history_write_h0_adv(blocks, itime)
    else if (hydrostatic) then
      call history_write_h0_hydrostatic(blocks, itime)
    else if (nonhydrostatic) then
      call history_write_h0_nonhydrostatic(blocks, itime)
    else
      call history_write_h0_swm(blocks, itime)
    end if

    call fiona_end_output('h0')

    if (is_root_proc()) then
      call cpu_time(time2)
      call log_notice('Done write state cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine history_write_h0

  subroutine history_write_h1_swm(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    integer is, ie, js, je, ks, ke
    integer start(2), count(2)

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h1', elapsed_seconds, new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h1', elapsed_seconds, new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h1', elapsed_seconds, new_file=.false.)
    end if
    call fiona_output('h1', 'lon'   , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h1', 'lat'   , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h1', 'ilon'  , global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h1', 'ilat'  , global_mesh%half_lat_deg(1:global_mesh%num_half_lat))

    associate (mesh   => blocks(1)%mesh        , &
               state  => blocks(1)%state(itime), &
               tend   => blocks(1)%tend (itime), &
               static => blocks(1)%static)
    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    start = [is,js]
    count = [mesh%num_full_lon,mesh%num_full_lat]
    call fiona_output('h1', 'dmfdlon' ,  tend%dmfdlon (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'dmfdlat' ,  tend%dmfdlat (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'ke'      , state%ke      (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'dgzdt'   ,  tend%dgz     (is:ie,js:je,1), start=start, count=count)

    is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    start = [is,js]
    count = [mesh%num_half_lon,mesh%num_full_lat]
    call fiona_output('h1', 'qhv'     ,  tend%qhv      (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'pgf_lon' ,  tend%pgf_lon  (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'dkedlon' ,  tend%dkedlon  (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'dudt   ' ,  tend%du       (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'mfx_lon' , state%mfx_lon  (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'mfy_lon' , state%mfy_lon  (is:ie,js:je,1), start=start, count=count)

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
    start = [is,js]
    count = [mesh%num_full_lon,mesh%num_half_lat]
    call fiona_output('h1', 'qhu'     ,  tend%qhu      (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'pgf_lat' ,  tend%pgf_lat  (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'dkedlat' ,  tend%dkedlat  (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'dvdt'    ,  tend%dv       (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'mfy_lat' , state%mfy_lat  (is:ie,js:je,1), start=start, count=count)
    call fiona_output('h1', 'mfx_lat' , state%mfx_lat  (is:ie,js:je,1), start=start, count=count)
    end associate

    call fiona_end_output('h1')

  end subroutine history_write_h1_swm

  subroutine history_write_h1_hydrostatic(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    integer is, ie, js, je, ks, ke
    integer start(3), count(3)

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h1', elapsed_seconds, new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h1', elapsed_seconds, new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h1', elapsed_seconds, new_file=.false.)
    end if
    call fiona_output('h1', 'lon'   , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h1', 'lat'   , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h1', 'ilon'  , global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h1', 'ilat'  , global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    call fiona_output('h1', 'lev'   , global_mesh%full_lev(1:global_mesh%num_full_lev))
    call fiona_output('h1', 'ilev'  , global_mesh%half_lev(1:global_mesh%num_half_lev))

    associate (mesh   => blocks(1)%mesh        , &
               state  => blocks(1)%state(itime), &
               tend   => blocks(1)%tend (itime), &
               static => blocks(1)%static)
    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
    call fiona_output('h1', 'dmfdlon' ,  tend%dmfdlon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmfdlat' ,  tend%dmfdlat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'ke'      , state%ke       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'm'       , state%m        (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dphsdt'  ,  tend%dphs     (is:ie,js:je      ), start=start, count=count)
    call fiona_output('h1', 'dptdt'   ,  tend%dpt      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dptfdlon',  tend%dptfdlon (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dptfdlat',  tend%dptfdlat (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dptfdlev',  tend%dptfdlev (is:ie,js:je,ks:ke), start=start, count=count)
    if (use_smag_damp) then
      call fiona_output('h1', 'smag_t', state%smag_t(is:ie,js:je,ks:ke), start=start, count=count)
    end if

    is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
    js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_half_lon,mesh%num_half_lat,mesh%num_full_lev]
    if (use_smag_damp) then
      call fiona_output('h1', 'smag_s', state%smag_s(is:ie,js:je,ks:ke), start=start, count=count)
    end if

    is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]
    call fiona_output('h1', 'qhv'     ,  tend%qhv      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'pgf_lon' ,  tend%pgf_lon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dkedlon' ,  tend%dkedlon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dudt   ' ,  tend%du       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfx_lon', state%mfx_lon   (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfy_lon', state%mfy_lon   (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'wedudlev',  tend%wedudlev (is:ie,js:je,ks:ke), start=start, count=count)
    if (use_smag_damp) then
      call fiona_output('h1', 'kmh_lon' , state%kmh_lon  (is:ie,js:je,ks:ke), start=start, count=count)
    end if

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]
    call fiona_output('h1', 'qhu'     ,  tend%qhu      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'pgf_lat' ,  tend%pgf_lat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dkedlat' ,  tend%dkedlat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dvdt'    ,  tend%dv       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfy_lat', state%mfy_lat   (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfx_lat', state%mfx_lat   (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'wedvdlev',  tend%wedvdlev (is:ie,js:je,ks:ke), start=start, count=count)
    if (use_smag_damp) then
      call fiona_output('h1', 'kmh_lat' , state%kmh_lat  (is:ie,js:je,ks:ke), start=start, count=count)
    end if

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]
    call fiona_output('h1', 'we_lev', state%we_lev(is:ie,js:je,ks:ke), start=start, count=count)
    end associate
    
    call fiona_end_output('h1')

  end subroutine history_write_h1_hydrostatic

  subroutine history_write_h1_nonhydrostatic(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    integer is, ie, js, je, ks, ke
    integer start(3), count(3)

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h1', elapsed_seconds, new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h1', elapsed_seconds, new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h1', elapsed_seconds, new_file=.false.)
    end if
    call fiona_output('h1', 'lon'   , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h1', 'lat'   , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    call fiona_output('h1', 'ilon'  , global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
    call fiona_output('h1', 'ilat'  , global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    call fiona_output('h1', 'lev'   , global_mesh%full_lev(1:global_mesh%num_full_lev))
    call fiona_output('h1', 'ilev'  , global_mesh%half_lev(1:global_mesh%num_half_lev))

    associate (mesh   => blocks(1)%mesh        , &
               state  => blocks(1)%state(itime), &
               tend   => blocks(1)%tend (itime), &
               static => blocks(1)%static)
    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]
    call fiona_output('h1', 'dmfdlon' ,  tend%dmfdlon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmfdlat' ,  tend%dmfdlat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'ke'      , state%ke       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'm'       , state%m        (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dphsdt'  ,  tend%dphs     (is:ie,js:je      ), start=start, count=count)
    call fiona_output('h1', 'dptdt'   ,  tend%dpt      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dptfdlon',  tend%dptfdlon (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dptfdlat',  tend%dptfdlat (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dptfdlev',  tend%dptfdlev (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]
    call fiona_output('h1', 'qhv'     ,  tend%qhv      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'pgf_lon' ,  tend%pgf_lon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dkedlon' ,  tend%dkedlon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dudt   ' ,  tend%du       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfx_lon' , state%mfx_lon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfy_lon' , state%mfy_lon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'wedudlev',  tend%wedudlev (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]
    call fiona_output('h1', 'qhu'     ,  tend%qhu      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'pgf_lat' ,  tend%pgf_lat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dkedlat' ,  tend%dkedlat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dvdt'    ,  tend%dv       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfy_lat' , state%mfy_lat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mfx_lat' , state%mfx_lat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'wedvdlev',  tend%wedvdlev (is:ie,js:je,ks:ke), start=start, count=count)

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]
    call fiona_output('h1', 'we_lev', state%we_lev(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_gz_lon'   ,  tend%adv_gz_lon   (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_gz_lat'   ,  tend%adv_gz_lat   (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_gz_lev'   ,  tend%adv_gz_lev   (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_w_lon'    ,  tend%adv_w_lon    (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_w_lat'    ,  tend%adv_w_lat    (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'adv_w_lev'    ,  tend%adv_w_lev    (is:ie,js:je,ks:ke), start=start, count=count)
    end associate

    call fiona_end_output('h1')

  end subroutine history_write_h1_nonhydrostatic

  subroutine history_write_h1(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    if (hydrostatic) then
      call history_write_h1_hydrostatic(blocks, itime)
    else if (nonhydrostatic) then
      call history_write_h1_nonhydrostatic(blocks, itime)
    else
      call history_write_h1_swm(blocks, itime)
    end if

  end subroutine history_write_h1

end module history_mod
