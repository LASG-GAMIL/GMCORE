module bkg_mod

  use fiona
  use flogger
  use namelist_mod
  use const_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use process_mod
  use parallel_mod
  use era5_reader_mod
#ifdef HAS_ECCODES
  use fnl_reader_mod
#endif
  use mpas_reader_mod
  use waccm_reader_mod
  use openmars_reader_mod
  use latlon_interp_mod
  use vert_interp_mod
  use interp_mod

  implicit none

  private

  public bkg_read
  public bkg_final
  public bkg_regrid_phs
  public bkg_calc_ph
  public bkg_regrid_pt
  public bkg_regrid_u
  public bkg_regrid_v

contains

  subroutine bkg_read(bkg_type_, bkg_file)

    character(*), intent(in) :: bkg_type_
    character(*), intent(in) :: bkg_file

    bkg_type = bkg_type_

    select case (bkg_type)
    case ('era5')
      call era5_reader_run(bkg_file)
#ifdef HAS_ECCODES
    case ('fnl')
      call fnl_reader_run(bkg_file)
#endif
    case ('mpas')
      call mpas_reader_run(bkg_file)
    case ('waccm')
      call waccm_reader_run(bkg_file)
    case ('openmars')
      call openmars_reader_run(bkg_file)
    case default
      if (is_root_proc()) call log_error('Unknown bkg_type ' // trim(bkg_type) // '!')
    end select

  end subroutine bkg_read

  subroutine bkg_final()

    call era5_reader_final()
#ifdef HAS_ECCODES
    call fnl_reader_final()
#endif
    call mpas_reader_final()
    call waccm_reader_final()
    call openmars_reader_final()

  end subroutine bkg_final

  subroutine bkg_regrid_phs()

    real(r8), allocatable, dimension(:,:) :: p0, t0, z0, t0_p
    real(r8) lapse_kappa
    integer iblk, i, j

    logical do_hydrostatic_correct, do_drymass_correct

    if (is_root_proc()) call log_notice('Regrid mean sea level pressure and calculate surface pressure based on pressure-height formula.')

    lapse_kappa = lapse_rate * Rd_o_g
    do_drymass_correct = .false.
    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%mesh        , &
                 gzs   => blocks(iblk)%static%gzs  , &
                 phs   => blocks(iblk)%state(1)%phs)
        allocate(p0  (mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
        allocate(t0  (mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
        allocate(z0  (mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
        allocate(t0_p(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))

        select case (bkg_type)
        case ('era5')
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_mslp, mesh, p0)
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_t(:,:,num_era5_lev), mesh, t0)
          t0_p = era5_lev(num_era5_lev)
          z0 = 0.0_r8
          do_hydrostatic_correct = .true.
#ifdef HAS_ECCODES
        case ('fnl')
          call latlon_interp_bilinear_cell(fnl_lon, fnl_lat, fnl_ps, mesh, phs)
          do_hydrostatic_correct = .false.
#endif
        case ('mpas')
          call latlon_interp_bilinear_cell(mpas_lon, mpas_lat, mpas_ps, mesh, phs)
          do_hydrostatic_correct = .false.
        case ('waccm')
          call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_ps, mesh, p0)
          call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_zs, mesh, z0)
          call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_p(:,:,waccm_nlev), mesh, t0_p)
          call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_t(:,:,waccm_nlev), mesh, t0  )
          do_hydrostatic_correct = .true.
        case ('openmars')
          call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_ps, mesh, phs)
          do_hydrostatic_correct = .false.
        end select
        ! According to pressure-height formula based on hydrostatic assumption.
        if (do_hydrostatic_correct) then
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              t0(i,j) = t0(i,j) * (p0(i,j) / t0_p(i,j))**lapse_kappa
              phs(i,j) = p0(i,j) * (1.0_r8 - lapse_rate * (gzs(i,j) / g - z0(i,j)) / t0(i,j))**(1.0_r8 / lapse_kappa)
            end do
          end do
        end if
        ! Evaluate dry mass surface pressure
        if (do_drymass_correct) then
          call calc_dry_air_ps(era5_lon, era5_lat, era5_lev, era5_ps, era5_q, era5_psd)
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_psd, mesh, phs)
        end if

        call fill_halo(block, phs, full_lon=.true., full_lat=.true.)
        deallocate(p0, t0, z0, t0_p)
      end associate
    end do

  end subroutine bkg_regrid_phs

  subroutine calc_dry_air_ps(lon, lat, lev, ps, q, psd)

    real(r8), intent(in) :: lon(num_era5_lon), &
                            lat(num_era5_lat), &
                            lev(num_era5_lev)
    real(r8), intent(in) :: ps(num_era5_lon,num_era5_lat), &
                             q(num_era5_lon,num_era5_lat,num_era5_lev)
    real(r8), intent(out) :: psd(num_era5_lon,num_era5_lat)

    integer i, j, k, nlev_eff
    real(r8) pfull(num_era5_lev  ), hpfull(num_era5_lev  ), dpfull(num_era5_lev)
    real(r8) pface(num_era5_lev+1), hpface(num_era5_lev+1)

    do j = 1, num_era5_lat
      do i = 1, num_era5_lon
        pfull = lev(:)
        ! Check effective full level number.
        nlev_eff = 99999
        do k = 1, num_era5_lev
          if (pfull(k) >= ps(i,j)) then
            nlev_eff = k
            exit
          end if 
        end do
        if (nlev_eff == 99999) nlev_eff = num_era5_lev 

        pface(2:nlev_eff) = (pfull(1:nlev_eff-1) + pfull(2:nlev_eff)) * 0.5_r8
        pface(1) = 0.0_r8
        pface(nlev_eff+1) = ps(i,j)
 
        dpfull(1:nlev_eff) = pface(2:nlev_eff+1) - pface(1:nlev_eff)
        hpface(1) = pface(1)

        do k = 2, nlev_eff+1
          hpface(k) = hpface(k-1) + dpfull(k-1) * (1.0_r8 - q(i,j,k-1))
        end do 
        hpfull(1:nlev_eff) = (hpface(1:nlev_eff) + hpface(2:nlev_eff+1)) * 0.5_r8
        psd(i,j) = hpface(nlev_eff+1)
      end do
    end do

  end subroutine calc_dry_air_ps

  subroutine bkg_calc_ph()

    integer iblk, i, j, k

    if (is_root_proc()) call log_notice('Calculate pressure on each grid.')

    do iblk = 1, size(blocks)
      associate (mesh => blocks(iblk)%mesh        , &
                 phs  => blocks(iblk)%state(1)%phs, &
                 ph   => blocks(iblk)%state(1)%ph)
        ! Calculate pressure on GMCORE grids.
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              ph(i,j,k) = vert_coord_calc_ph(k, phs(i,j))
            end do
          end do
        end do
      end associate
    end do

  end subroutine bkg_calc_ph

  subroutine bkg_regrid_pt()

    real(r8), allocatable, dimension(:,:,:) :: t1, pt1, p1
    integer iblk, i, j, k

    if (is_root_proc()) call log_notice('Regrid temperature and calculate potential temperature.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)             , &
                 mesh  => blocks(iblk)%mesh        , &
                 phs   => blocks(iblk)%state(1)%phs, &
                 ph    => blocks(iblk)%state(1)%ph , &
                 t     => blocks(iblk)%state(1)%t  , &
                 pt    => blocks(iblk)%state(1)%pt)
        select case (bkg_type)
        case ('era5')
          allocate(t1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_era5_lev))
          do k = 1, num_era5_lev
            call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_t(:,:,k), mesh, t1(:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_log_linear(era5_lev, t1(i,j,:), ph(i,j,1:mesh%num_full_lev), t(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
              pt(i,j,:) = potential_temperature(t(i,j,:), ph(i,j,:))
            end do
          end do
          deallocate(t1)
#ifdef HAS_ECCODES
        case ('fnl')
          allocate(t1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_fnl_lev))
          do k = 1, num_fnl_lev
            call latlon_interp_bilinear_cell(fnl_lon, fnl_lat, fnl_t(:,:,k), mesh, t1(:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_log_linear(fnl_lev, t1(i,j,:), ph(i,j,1:mesh%num_full_lev), t(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
              pt(i,j,:) = potential_temperature(t(i,j,:), ph(i,j,:))
            end do
          end do
          deallocate(t1)
#endif
        case ('mpas')
          allocate(pt1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_mpas_lev))
          allocate(p1 (mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_mpas_lev))
          do k = 1, num_mpas_lev
            call latlon_interp_bilinear_cell(mpas_lon, mpas_lat, mpas_pt(:,:,k), mesh, pt1(:,:,k))
            call latlon_interp_bilinear_cell(mpas_lon, mpas_lat, mpas_p (:,:,k), mesh, p1 (:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_log_linear(p1(i,j,:), pt1(i,j,:), ph(i,j,1:mesh%num_full_lev), pt(i,j,1:mesh%num_full_lev), allow_extrap=.false.)
            end do
          end do
          deallocate(pt1, p1)
        case ('waccm')
          allocate(t1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,waccm_nlev))
          allocate(p1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,waccm_nlev))
          do k = 1, waccm_nlev
            call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_t(:,:,k), mesh, t1(:,:,k))
            call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), ph(i,j,1:mesh%num_full_lev), t(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
              pt(i,j,:) = potential_temperature(t(i,j,:), ph(i,j,:))
            end do
          end do
          deallocate(t1, p1)
        case ('openmars')
          allocate(t1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,nlev_openmars))
          allocate(p1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,nlev_openmars))
          do k = 1, nlev_openmars
            call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_t(:,:,k), mesh, t1(:,:,k))
            call latlon_interp_bilinear_cell(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_log_linear(p1(i,j,:), t1(i,j,:), ph(i,j,1:mesh%num_full_lev), t(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
              pt(i,j,:) = potential_temperature(t(i,j,:), ph(i,j,:))
            end do
          end do
          deallocate(t1, p1)
        end select
        call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      end associate
    end do

  end subroutine bkg_regrid_pt

  subroutine bkg_regrid_u()

    real(r8), allocatable, dimension(:,:,:) :: u1, p1
    integer iblk, i, j, k

    if (is_root_proc()) call log_notice('Regrid u wind component.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)            , &
                 mesh  => blocks(iblk)%mesh       , &
                 ph    => blocks(iblk)%state(1)%ph, &
                 u     => blocks(iblk)%state(1)%u_lon)
        select case (bkg_type)
        case ('era5')
          allocate(u1(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_era5_lev))
          do k = 1, num_era5_lev
            call latlon_interp_bilinear_lon_edge(era5_lon, era5_lat, era5_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              call vert_interp_linear(era5_lev, u1(i,j,:), ph(i,j,1:mesh%num_full_lev), u(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
            end do
          end do
          deallocate(u1)
#ifdef HAS_ECCODES
        case ('fnl')
          allocate(u1(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_fnl_lev))
          do k = 1, num_fnl_lev
            call latlon_interp_bilinear_lon_edge(fnl_lon, fnl_lat, fnl_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              call vert_interp_linear(fnl_lev, u1(i,j,:), ph(i,j,1:mesh%num_full_lev), u(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
            end do
          end do
          deallocate(u1)
#endif
        case ('mpas')
          allocate(u1(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_mpas_lev))
          allocate(p1(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_mpas_lev))
          do k = 1, num_mpas_lev
            call latlon_interp_bilinear_lon_edge(mpas_lon, mpas_lat, mpas_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
            call latlon_interp_bilinear_lon_edge(mpas_lon, mpas_lat, mpas_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              call vert_interp_linear(p1(i,j,:), u1(i,j,:), ph(i,j,1:mesh%num_full_lev), u(i,j,1:mesh%num_full_lev), allow_extrap=.false.)
            end do
          end do
          deallocate(u1, p1)
        case ('waccm')
          allocate(u1(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,waccm_nlev))
          allocate(p1(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,waccm_nlev))
          do k = 1, waccm_nlev
            call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
            call latlon_interp_bilinear_cell(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                call vert_interp_linear(p1(i,j,:), u1(i,j,:), ph(i,j,1:mesh%num_full_lev), u(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
            end do
          end do
          deallocate(u1, p1)
        case ('openmars')
          allocate(u1(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,nlev_openmars))
          allocate(p1(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,nlev_openmars))
          do k = 1, nlev_openmars
            call latlon_interp_bilinear_lon_edge(openmars_lon, openmars_lat, openmars_u(:,:,k), mesh, u1(:,:,k), zero_pole=.true.)
            call latlon_interp_bilinear_lon_edge(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              call vert_interp_linear(p1(i,j,:), u1(i,j,:), ph(i,j,1:mesh%num_full_lev), u(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
            end do
          end do
          deallocate(u1, p1)
        end select
        call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
      end associate
    end do

  end subroutine bkg_regrid_u

  subroutine bkg_regrid_v()

    real(r8), allocatable, dimension(:,:,:) :: v1, p1
    integer iblk, i, j, k

    if (is_root_proc()) call log_notice('Regrid v wind component.')

    do iblk = 1, size(blocks)
      associate (block => blocks(iblk)            , &
                 mesh  => blocks(iblk)%mesh       , &
                 ph    => blocks(iblk)%state(1)%ph, &
                 v     => blocks(iblk)%state(1)%v_lat)
        select case (bkg_type)
        case ('era5')
          allocate(v1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,num_era5_lev))
          do k = 1, num_era5_lev
            call latlon_interp_bilinear_lat_edge(era5_lon, era5_lat, era5_v(:,:,k), mesh, v1(:,:,k))
          end do
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_linear(era5_lev, v1(i,j,:), ph(i,j,1:mesh%num_full_lev), v(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
            end do
          end do
          deallocate(v1)
#ifdef HAS_ECCODES
        case ('fnl')
          allocate(v1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,num_fnl_lev))
          do k = 1, num_fnl_lev
            call latlon_interp_bilinear_lat_edge(fnl_lon, fnl_lat, fnl_v(:,:,k), mesh, v1(:,:,k))
          end do
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_linear(fnl_lev, v1(i,j,:), ph(i,j,1:mesh%num_full_lev), v(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
            end do
          end do
          deallocate(v1)
#endif
        case ('mpas')
          allocate(v1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,num_mpas_lev))
          allocate(p1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,num_mpas_lev))
          do k = 1, num_mpas_lev
            call latlon_interp_bilinear_lat_edge(mpas_lon, mpas_lat, mpas_v(:,:,k), mesh, v1(:,:,k))
            call latlon_interp_bilinear_lat_edge(mpas_lon, mpas_lat, mpas_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_linear(p1(i,j,:), v1(i,j,:), ph(i,j,1:mesh%num_full_lev), v(i,j,1:mesh%num_full_lev), allow_extrap=.false.)
            end do
          end do
          deallocate(v1, p1)
        case ('waccm')
          allocate(v1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,waccm_nlev))
          allocate(p1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,waccm_nlev))
          do k = 1, waccm_nlev
            call latlon_interp_bilinear_lat_edge(waccm_lon, waccm_lat, waccm_v(:,:,k), mesh, v1(:,:,k))
            call latlon_interp_bilinear_lat_edge(waccm_lon, waccm_lat, waccm_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_linear(p1(i,j,:), v1(i,j,:), ph(i,j,1:mesh%num_full_lev), v(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
            end do
          end do
          deallocate(v1, p1)
        case ('openmars')
          allocate(v1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,nlev_openmars))
          allocate(p1(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,nlev_openmars))
          do k = 1, nlev_openmars
            call latlon_interp_bilinear_lat_edge(openmars_lon, openmars_lat, openmars_v(:,:,k), mesh, v1(:,:,k))
            call latlon_interp_bilinear_lat_edge(openmars_lon, openmars_lat, openmars_p(:,:,k), mesh, p1(:,:,k))
          end do
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_linear(p1(i,j,:), v1(i,j,:), ph(i,j,1:mesh%num_full_lev), v(i,j,1:mesh%num_full_lev), allow_extrap=.true.)
            end do
          end do
          deallocate(v1, p1)
        end select
        call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
      end associate
    end do

  end subroutine bkg_regrid_v

end module bkg_mod
