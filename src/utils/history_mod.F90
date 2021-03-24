module history_mod

  use fiona
  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use allocator_mod
  use time_mod, dt => dt_in_seconds
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use block_mod
  use diag_state_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write_state
  public history_write_debug

contains

  subroutine history_init()

    character(10) time_value, time_units
    character(4) cell_dims(4), cell_dims_2d(3)
    character(4) lon_dims(4), lon_dims_2d(3)
    character(4) lat_dims(4), lat_dims_2d(3)
    character(4) vtx_dims(4), vtx_dims_2d(3)
    character(4) lon_lev_dims(4), lat_lev_dims(4), lev_dims(4)
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

        cell_dims(1) =  'lon';     cell_dims(2) =  'lat';     cell_dims(3) =  'lev';     cell_dims(4) = 'time'
         lon_dims(1) = 'ilon';      lon_dims(2) =  'lat';      lon_dims(3) =  'lev';      lon_dims(4) = 'time'
     lon_lev_dims(1) = 'ilon';  lon_lev_dims(2) =  'lat';  lon_lev_dims(3) = 'ilev';  lon_lev_dims(4) = 'time'
         lat_dims(1) =  'lon';      lat_dims(2) = 'ilat';      lat_dims(3) =  'lev';      lat_dims(4) = 'time'
     lat_lev_dims(1) =  'lon';  lat_lev_dims(2) = 'ilat';  lat_lev_dims(3) = 'ilev';  lat_lev_dims(4) = 'time'
         vtx_dims(1) = 'ilon';      vtx_dims(2) = 'ilat';      vtx_dims(3) =  'lev';      vtx_dims(4) = 'time'
         lev_dims(1) =  'lon';      lev_dims(2) =  'lat';      lev_dims(3) = 'ilev';      lev_dims(4) = 'time'
     cell_dims_2d(1) =  'lon';  cell_dims_2d(2) =  'lat';  cell_dims_2d(3) = 'time'
      lon_dims_2d(1) = 'ilon';   lon_dims_2d(2) =  'lat';   lon_dims_2d(3) = 'time'
      lat_dims_2d(1) =  'lon';   lat_dims_2d(2) = 'ilat';   lat_dims_2d(3) = 'time'
      vtx_dims_2d(1) = 'ilon';   vtx_dims_2d(2) = 'ilat';   vtx_dims_2d(3) = 'time'

    call fiona_init(time_units, start_time_str)
    call fiona_create_dataset('h0', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size, split_file=split_h0)
    call fiona_add_att('h0', 'time_step_size', dt)
    call fiona_add_dim('h0', 'time' , add_var=.true.)
    call fiona_add_dim('h0', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h0', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)

    if (.not. baroclinic .or. .not. diag_state(1)%is_init()) then
      call fiona_add_dim('h0', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
      call fiona_add_dim('h0', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    end if

    ! Common or barotropic variables
    if (baroclinic) then
      call fiona_add_var('h0', 'phs'  , long_name='surface hydrostatic pressure', units='Pa'     , dim_names=cell_dims_2d)
    else
      call fiona_add_var('h0', 'u'    , long_name='u wind component'            , units='m s-1'  , dim_names=lon_dims_2d)
      call fiona_add_var('h0', 'v'    , long_name='v wind component'            , units='m s-1'  , dim_names=lat_dims_2d)
      call fiona_add_var('h0', 'z'    , long_name='height'                      , units='m'      , dim_names=cell_dims_2d)
      call fiona_add_var('h0', 'vor'  , long_name='relative vorticity'          , units='s-1'    , dim_names=vtx_dims_2d)
      call fiona_add_var('h0', 'div'  , long_name='divergence'                  , units='s-1'    , dim_names=cell_dims_2d)
    end if
    call fiona_add_var('h0', 'zs' , long_name='surface height'           , units='m'     , dim_names=cell_dims_2d)
    call fiona_add_var('h0', 'tm' , long_name='total mass'               , units='m'     , dim_names=['time'])
    call fiona_add_var('h0', 'te' , long_name='total energy'             , units='m4 s-4', dim_names=['time'], data_type='r8')
    call fiona_add_var('h0', 'tpe', long_name='total potential enstrophy', units='m2 s-5', dim_names=['time'], data_type='r8')
    call fiona_add_var('h0', 'tpv', long_name='total potential vorticity', units='m2 s-5', dim_names=['time'], data_type='r8')

    if (baroclinic) then
      if (diag_state(1)%is_init()) then
        select case (diag_state(1)%level_type)
        case (height_levels)
          diag_dims(1) = 'lon'; diag_dims(2) = 'lat'; diag_dims(3) = 'zlev'; diag_dims(4) = 'time'
          call fiona_add_dim('h0', 'zlev', size=size(diag_state(1)%levels), add_var=.true., decomp=.false.)
          call fiona_add_var('h0', 'p'   , long_name='pressure'           , units='Pa', dim_names=diag_dims)
        case (pressure_levels)
          diag_dims(1) = 'lon'; diag_dims(2) = 'lat'; diag_dims(3) = 'plev'; diag_dims(4) = 'time'
          call fiona_add_dim('h0', 'plev', size=size(diag_state(1)%levels), add_var=.true., decomp=.false.)
          call fiona_add_var('h0', 'z'   , long_name='geopotential height', units='m', dim_names=diag_dims)
        end select
        call fiona_add_var('h0', 'u'    , long_name='u wind component'      , units='m s-1'  , dim_names=diag_dims)
        call fiona_add_var('h0', 'v'    , long_name='v wind component'      , units='m s-1'  , dim_names=diag_dims)
        call fiona_add_var('h0', 'pt'   , long_name='potential temperature' , units='K'      , dim_names=diag_dims)
        call fiona_add_var('h0', 't'    , long_name='temperature'           , units='K'      , dim_names=diag_dims)
        if (nonhydrostatic) then
          call fiona_add_var('h0', 'w'  , long_name='vertical speed'        , units='m s-1'  , dim_names=diag_dims)
        end if
      else
        call fiona_add_dim('h0', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
        call fiona_add_dim('h0', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
        call fiona_add_var('h0', 'u'    , long_name='u wind component'      , units='m s-1'  , dim_names=lon_dims)
        call fiona_add_var('h0', 'v'    , long_name='v wind component'      , units='m s-1'  , dim_names=lat_dims)
        call fiona_add_var('h0', 'pt'   , long_name='potential temperature' , units='K'      , dim_names=cell_dims)
        call fiona_add_var('h0', 'z'    , long_name='height'                , units='m'      , dim_names=lev_dims)
        call fiona_add_var('h0', 'vor'  , long_name='relative vorticity'    , units='s-1'    , dim_names=vtx_dims)
        call fiona_add_var('h0', 'div'  , long_name='divergence'            , units='s-1'    , dim_names=cell_dims)
        if (nonhydrostatic) then
          call fiona_add_var('h0', 'w', long_name='vertical speed', units='m s-1', dim_names=lev_dims)
          call fiona_add_var('h0', 'p', long_name='full pressure' , units='Pa'   , dim_names=lev_dims)
        end if
      end if
    end if

    call fiona_create_dataset('h1', desc=case_desc, file_prefix=trim(case_name), mpi_comm=proc%comm, group_size=output_group_size)
    call fiona_add_att('h1', 'time_step_size', dt)
    call fiona_add_dim('h1', 'time' , add_var=.true.)
    call fiona_add_dim('h1', 'lon'  , size=global_mesh%num_full_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'lat'  , size=global_mesh%num_full_lat, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilon' , size=global_mesh%num_half_lon, add_var=.true., decomp=.true.)
    call fiona_add_dim('h1', 'ilat' , size=global_mesh%num_half_lat, add_var=.true., decomp=.true.)
    if (baroclinic) then
      call fiona_add_dim('h1', 'lev'  , size=global_mesh%num_full_lev, add_var=.true., decomp=.false.)
      call fiona_add_dim('h1', 'ilev' , size=global_mesh%num_half_lev, add_var=.true., decomp=.false.)
      call fiona_add_var('h1', 'dudt'         , long_name='u wind component tendency'                     , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'dvdt'         , long_name='v wind component tendency'                     , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'dphsdt'       , long_name='surface hydrostatic pressure tendency'         , units='', dim_names=cell_dims_2d)
      call fiona_add_var('h1', 'dptdt'        , long_name='potential temperature tendency'                , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'dptfdlon'     , long_name='zonal potential temperature flux gradient'     , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'dptfdlat'     , long_name='meridional potential temperature flux gradient', units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'dptfdlev'     , long_name='vertical potential temperature flux gradient'  , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'wedphdlev_lev', long_name='vertical coordinate velocity'                  , units='', dim_names=lev_dims)
      call fiona_add_var('h1', 'pgf_lon'      , long_name='zonal pressure gradient force'                 , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'wedudlev'     , long_name='vertical advection of u'                       , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'pgf_lat'      , long_name='meridional pressure gradient force'            , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'wedvdlev'     , long_name='vertical advection of v'                       , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'qhv'          , long_name='nonliear zonal Coriolis force'                 , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'qhu'          , long_name='nonliear meridional Coriolis force'            , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'dkedlon'      , long_name='zonal kinetic energy gradient force'           , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'dkedlat'      , long_name='meridional kinetic energy gradient force'      , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'dmfdlon'      , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'dmfdlat'      , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'mf_lon_n'     , long_name='normal mass flux on U grid'                    , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'mf_lon_t'     , long_name='tangent mass flux on U grid'                   , units='', dim_names=lon_dims)
      call fiona_add_var('h1', 'mf_lat_n'     , long_name='normal mass flux on V grid'                    , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'mf_lat_t'     , long_name='tangent mass flux on V grid'                   , units='', dim_names=lat_dims)
      call fiona_add_var('h1', 'm'            , long_name='dph on full levels'                            , units='', dim_names=cell_dims)
      call fiona_add_var('h1', 'ke'           , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims)
      if (nonhydrostatic) then
        call fiona_add_var('h1', 'adv_gz_lon', long_name='advection of gz', units='', dim_names=lev_dims)
        call fiona_add_var('h1', 'adv_gz_lat', long_name='advection of gz', units='', dim_names=lev_dims)
        call fiona_add_var('h1', 'adv_gz_lev', long_name='advection of gz', units='', dim_names=lev_dims)
        call fiona_add_var('h1', 'adv_w_lon' , long_name='advection of w' , units='', dim_names=lev_dims)
        call fiona_add_var('h1', 'adv_w_lat' , long_name='advection of w' , units='', dim_names=lev_dims)
        call fiona_add_var('h1', 'adv_w_lev' , long_name='advection of w' , units='', dim_names=lev_dims)
        call fiona_add_var('h1', 'dzsdlon', long_name='zonal zs gradient', units='', dim_names=cell_dims_2d)
        call fiona_add_var('h1', 'dzsdlat', long_name='meridional zs gradient', units='', dim_names=cell_dims_2d)
      end if
    else
      call fiona_add_var('h1', 'dudt'         , long_name='u wind component tendency'                     , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'dvdt'         , long_name='v wind component tendency'                     , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'dgzdt'        , long_name='geopotential tendency'                         , units='', dim_names=cell_dims_2d)
      call fiona_add_var('h1', 'qhv'          , long_name='nonliear zonal Coriolis force'                 , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'qhu'          , long_name='nonliear meridional Coriolis force'            , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'pgf_lon'      , long_name='zonal geopotential energy gradient force'      , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'dkedlon'      , long_name='zonal kinetic energy gradient force'           , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'pgf_lat'      , long_name='meridional geopotential energy gradient force' , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'dkedlat'      , long_name='meridional kinetic energy gradient force'      , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'dmfdlon'      , long_name='zonal mass flux divergence'                    , units='', dim_names=cell_dims_2d)
      call fiona_add_var('h1', 'dmfdlat'      , long_name='meridional mass flux divergence'               , units='', dim_names=cell_dims_2d)
      call fiona_add_var('h1', 'mf_lon_n'     , long_name='normal mass flux on U grid'                    , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'mf_lon_t'     , long_name='tangent mass flux on U grid'                   , units='', dim_names=lon_dims_2d)
      call fiona_add_var('h1', 'mf_lat_n'     , long_name='normal mass flux on V grid'                    , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'mf_lat_t'     , long_name='tangent mass flux on V grid'                   , units='', dim_names=lat_dims_2d)
      call fiona_add_var('h1', 'ke'           , long_name='kinetic energy on cell grid'                   , units='', dim_names=cell_dims_2d)
    end if

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

  end subroutine history_init

  subroutine history_final()

  end subroutine history_final

  subroutine history_write_state(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    integer iblk, is, ie, js, je, ks, ke
    integer start(3), count(3)
    real(8) time1, time2

    if (is_root_proc()) then
      call log_notice('Write state.')
      call cpu_time(time1)
    end if

    if (.not. time_has_alert('h0_new_file')) then
      call fiona_start_output('h0', elapsed_seconds, new_file=time_step==0)
    else if (time_is_alerted('h0_new_file')) then
      call fiona_start_output('h0', elapsed_seconds, new_file=.true., tag=curr_time%format('%Y-%m-%d_%H_%M'))
    else
      call fiona_start_output('h0', elapsed_seconds, new_file=.false.)
    end if
    call fiona_output('h0', 'lon' , global_mesh%full_lon_deg(1:global_mesh%num_full_lon))
    call fiona_output('h0', 'lat' , global_mesh%full_lat_deg(1:global_mesh%num_full_lat))
    if (baroclinic) then
      if (diag_state(1)%is_init()) then
        select case (diag_state(1)%level_type)
        case (height_levels)
          call fiona_output('h0', 'zlev', diag_state(1)%levels)
        case (pressure_levels)
          call fiona_output('h0', 'plev', diag_state(1)%levels)
        end select
      else
        call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
        call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
        call fiona_output('h0', 'lev' , global_mesh%full_lev)
        call fiona_output('h0', 'ilev', global_mesh%half_lev)
      end if
    else
      call fiona_output('h0', 'ilon', global_mesh%half_lon_deg(1:global_mesh%num_half_lon))
      call fiona_output('h0', 'ilat', global_mesh%half_lat_deg(1:global_mesh%num_half_lat))
    end if

    do iblk = 1, size(blocks)
      mesh => blocks(iblk)%mesh
      state => blocks(iblk)%state(itime)
      static => blocks(iblk)%static

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

      call fiona_output('h0', 'zs' , static%gzs(is:ie,js:je) / g, start=start, count=count)

      if (baroclinic) then
        call fiona_output('h0', 'phs', state%phs(is:ie,js:je), start=start, count=count)
        if (diag_state(1)%is_init()) then
          count(3) = size(diag_state(iblk)%levels)
          call fiona_output('h0', 'u'  , diag_state(iblk)%u (is:ie,js:je,:), start=start, count=count)
          call fiona_output('h0', 'v'  , diag_state(iblk)%v (is:ie,js:je,:), start=start, count=count)
          call fiona_output('h0', 'pt' , diag_state(iblk)%pt(is:ie,js:je,:), start=start, count=count)
          call fiona_output('h0', 't'  , diag_state(iblk)%t (is:ie,js:je,:), start=start, count=count)
          if (nonhydrostatic) then
            call fiona_output('h0', 'w', diag_state(iblk)%w (is:ie,js:je,:), start=start, count=count)
          end if
          select case (diag_state(1)%level_type)
          case (height_levels)
            call fiona_output('h0', 'p', diag_state(iblk)%p (is:ie,js:je,:), start=start, count=count)
          case (pressure_levels)
            call fiona_output('h0', 'z', diag_state(iblk)%z (is:ie,js:je,:), start=start, count=count)
          end select
        else
          call fiona_output('h0', 'pt' , state%pt (is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output('h0', 'div', state%div(is:ie,js:je,ks:ke), start=start, count=count)
        end if
      else
        call fiona_output('h0', 'z'  , state%gz (is:ie,js:je,ks:ke) / g, start=start, count=count)
        call fiona_output('h0', 'div', state%div(is:ie,js:je,ks:ke), start=start, count=count)
      end if

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

      if (.not. diag_state(1)%is_init()) then
        call fiona_output('h0', 'u', state%u(is:ie,js:je,ks:ke), start=start, count=count)
      end if

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

      if (.not. diag_state(1)%is_init()) then
        call fiona_output('h0', 'v', state%v(is:ie,js:je,ks:ke), start=start, count=count)
      end if

      is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
      js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
      ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
      start = [is,js,ks]
      count = [mesh%num_half_lon,mesh%num_half_lat,mesh%num_full_lev]

      if (.not. diag_state(1)%is_init()) then
        call fiona_output('h0', 'vor', state%vor(is:ie,js:je,ks:ke), start=start, count=count)
      end if

      is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
      js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
      ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
      start = [is,js,ks]
      count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]

      if (nonhydrostatic .and. .not. diag_state(1)%is_init()) then
        call fiona_output('h0', 'z', state%gz_lev(is:ie,js:je,ks:ke) / g, start=start, count=count)
        if (nonhydrostatic) then
          call fiona_output('h0', 'w', state%w_lev(is:ie,js:je,ks:ke), start=start, count=count)
          call fiona_output('h0', 'p', state%p_lev(is:ie,js:je,ks:ke), start=start, count=count)
        end if
      end if

      call fiona_output('h0', 'tm' , state %tm)
      call fiona_output('h0', 'te' , state %te)
      call fiona_output('h0', 'tpe', state %tpe)
      call fiona_output('h0', 'tpv', state %tav)
    end do
    call fiona_end_output('h0')

    if (is_root_proc()) then
      call cpu_time(time2)
      call log_notice('Done write state cost ' // to_str(time2 - time1, 5) // ' seconds.')
    end if

  end subroutine history_write_state

  subroutine history_write_debug(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(tend_type), pointer :: tend

    integer is, ie, js, je, ks, ke
    integer start(3), count(3)

    mesh => blocks(1)%mesh
    state => blocks(1)%state(itime)
    tend => blocks(1)%tend(itime)

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
    if (baroclinic) then
      call fiona_output('h1', 'lev' , global_mesh%full_lev)
      call fiona_output('h1', 'ilev', global_mesh%half_lev)
    end if

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_full_lev]

    call fiona_output('h1', 'dmfdlon' , tend%dmfdlon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dmfdlat' , tend%dmfdlat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'ke'      , state%ke      (is:ie,js:je,ks:ke), start=start, count=count)
    if (baroclinic) then
      call fiona_output('h1', 'm'       , state%m      (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'dphsdt'  , tend%dphs    (is:ie,js:je      ), start=start, count=count)
      call fiona_output('h1', 'dptdt'   , tend%dpt     (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'dptfdlon', tend%dptfdlon(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'dptfdlat', tend%dptfdlat(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'dptfdlev', tend%dptfdlev(is:ie,js:je,ks:ke), start=start, count=count)
    else
      call fiona_output('h1', 'dgzdt'   , tend%dgz     (is:ie,js:je,ks:ke), start=start, count=count)
    end if

    is = mesh%half_lon_ibeg; ie = mesh%half_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_half_lon,mesh%num_full_lat,mesh%num_full_lev]

    call fiona_output('h1', 'qhv'     , tend%qhv      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'pgf_lon' , tend%pgf_lon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dkedlon' , tend%dkedlon  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dudt   ' , tend%du       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mf_lon_n', state%mf_lon_n(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mf_lon_t', state%mf_lon_t(is:ie,js:je,ks:ke), start=start, count=count)

    if (baroclinic) then
      call fiona_output('h1', 'wedudlev', tend%wedudlev(is:ie,js:je,ks:ke), start=start, count=count)
    end if

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%half_lat_ibeg; je = mesh%half_lat_iend
    ks = mesh%full_lev_ibeg; ke = mesh%full_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_half_lat,mesh%num_full_lev]

    call fiona_output('h1', 'qhu'     , tend%qhu      (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'pgf_lat' , tend%pgf_lat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dkedlat' , tend%dkedlat  (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'dvdt'    , tend%dv       (is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mf_lat_n', state%mf_lat_n(is:ie,js:je,ks:ke), start=start, count=count)
    call fiona_output('h1', 'mf_lat_t', state%mf_lat_t(is:ie,js:je,ks:ke), start=start, count=count)

    if (baroclinic) then
      call fiona_output('h1', 'wedvdlev', tend%wedvdlev(is:ie,js:je,ks:ke), start=start, count=count)
    end if

    is = mesh%full_lon_ibeg; ie = mesh%full_lon_iend
    js = mesh%full_lat_ibeg; je = mesh%full_lat_iend
    ks = mesh%half_lev_ibeg; ke = mesh%half_lev_iend
    start = [is,js,ks]
    count = [mesh%num_full_lon,mesh%num_full_lat,mesh%num_half_lev]

    if (baroclinic) then
      call fiona_output('h1', 'wedphdlev_lev', state%wedphdlev_lev(is:ie,js:je,ks:ke), start=start, count=count)
    end if
    if (nonhydrostatic) then
      call fiona_output('h1', 'adv_gz_lon', tend%adv_gz_lon(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'adv_gz_lat', tend%adv_gz_lat(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'adv_gz_lev', tend%adv_gz_lev(is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'adv_w_lon' , tend%adv_w_lon (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'adv_w_lat' , tend%adv_w_lat (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'adv_w_lev' , tend%adv_w_lev (is:ie,js:je,ks:ke), start=start, count=count)
      call fiona_output('h1', 'dzsdlon', blocks(1)%static%dzsdlon(is:ie,js:je), start=start, count=count)
      call fiona_output('h1', 'dzsdlat', blocks(1)%static%dzsdlat(is:ie,js:je), start=start, count=count)
    end if

    call fiona_end_output('h1')

  end subroutine history_write_debug

end module history_mod
