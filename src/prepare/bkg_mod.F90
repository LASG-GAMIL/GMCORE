module bkg_mod

  use fiona
  use flogger
  use const_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use process_mod
  use parallel_mod
  use era5_reader_mod
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

  character(30) :: bkg_type = 'N/A'

contains

  subroutine bkg_read(bkg_type_, bkg_file)

    character(*), intent(in) :: bkg_type_
    character(*), intent(in) :: bkg_file

    bkg_type = bkg_type_

    select case (bkg_type)
    case ('era5')
      call era5_reader_run(bkg_file)
    case default
      if (is_root_proc()) call log_error('Unknown bkg_type ' // trim(bkg_type) // '!')
    end select

  end subroutine bkg_read

  subroutine bkg_final()

    call era5_reader_final()

  end subroutine bkg_final

  subroutine bkg_regrid_phs()

    real(r8), allocatable, dimension(:,:) :: p0, t0, z0
    real(r8) t0_p
    real(r8), parameter :: lapse = 0.006_r8 ! K m-1
    real(r8), parameter :: kappa = Rd / g
    real(r8), parameter :: lapse_kappa = lapse * kappa
    integer iblk, i, j

    call log_notice('Regrid mean sea level pressure and calculate surface pressure based on pressure-height formula.')

    do iblk = 1, size(proc%blocks)
      associate (mesh => proc%blocks(iblk)%mesh, &
                 gzs => proc%blocks(iblk)%static%gzs, &
                 phs => proc%blocks(iblk)%state(1)%phs)
        allocate(p0(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
        allocate(t0(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))
        allocate(z0(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub))

        select case (bkg_type)
        case ('era5')
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_mslp, mesh, p0)
          call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_t(:,:,num_era5_lev), mesh, t0)
          t0_p = era5_lev(num_era5_lev)
          z0 = 0.0_r8
        end select
        ! According to pressure-height formula based on hydrostatic assumption.
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            t0(i,j) = t0(i,j) * (p0(i,j) / t0_p)**lapse_kappa
            phs(i,j) = p0(i,j) * (1.0_r8 - lapse * (gzs(i,j) / g - z0(i,j)) / t0(i,j))**(1.0_r8 / lapse_kappa)
          end do
        end do

        deallocate(p0, t0, z0)
      end associate
    end do

  end subroutine bkg_regrid_phs

  subroutine bkg_calc_ph()

    integer iblk, i, j, k

    call log_notice('Calculate pressure on each grid.')

    do iblk = 1, size(proc%blocks)
      associate (mesh => proc%blocks(iblk)%mesh        , &
                 phs  => proc%blocks(iblk)%state(1)%phs, &
                 ph   => proc%blocks(iblk)%state(1)%ph)
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

    real(r8), allocatable :: tmp(:,:,:)
    integer iblk, i, j, k

    call log_notice('Regrid temperature and calculate potential temperature.')

    do iblk = 1, size(proc%blocks)
      associate (block => proc%blocks(iblk)             , &
                 mesh  => proc%blocks(iblk)%mesh        , &
                 phs   => proc%blocks(iblk)%state(1)%phs, &
                 ph    => proc%blocks(iblk)%state(1)%ph , &
                 t     => proc%blocks(iblk)%state(1)%t  , &
                 pt    => proc%blocks(iblk)%state(1)%pt)
        select case (bkg_type)
        case ('era5')
          allocate(tmp(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_era5_lev))
          do k = 1, num_era5_lev
            call latlon_interp_bilinear_cell(era5_lon, era5_lat, era5_t(:,:,k), mesh, tmp(:,:,k))
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_log_linear(era5_lev, tmp(i,j,:), ph(i,j,:), t(i,j,:), allow_extrap=.true.)
              pt(i,j,:) = potential_temperature(t(i,j,:), ph(i,j,:))
            end do
          end do
        end select
        deallocate(tmp)
        call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      end associate
    end do

  end subroutine bkg_regrid_pt

  subroutine bkg_regrid_u()

    real(r8), allocatable :: tmp(:,:,:)
    integer iblk, i, j, k

    call log_notice('Regrid u wind component.')

    do iblk = 1, size(proc%blocks)
      associate (block => proc%blocks(iblk)             , &
                 mesh  => proc%blocks(iblk)%mesh        , &
                 ph    => proc%blocks(iblk)%state(1)%ph , &
                 u     => proc%blocks(iblk)%state(1)%u)
        select case (bkg_type)
        case ('era5')
          allocate(tmp(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,num_era5_lev))
          do k = 1, num_era5_lev
            call latlon_interp_bilinear_lon_edge(era5_lon, era5_lat, era5_u(:,:,k), mesh, tmp(:,:,k), zero_pole=.true.)
          end do
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              call vert_interp_linear(era5_lev, tmp(i,j,:), ph(i,j,:), u(i,j,:), allow_extrap=.true.)
            end do
          end do
          deallocate(tmp)
          call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
        end select
      end associate
    end do

  end subroutine bkg_regrid_u

  subroutine bkg_regrid_v()

    real(r8), allocatable :: tmp(:,:,:)
    integer iblk, i, j, k

    call log_notice('Regrid v wind component.')

    do iblk = 1, size(proc%blocks)
      associate (block => proc%blocks(iblk)             , &
                 mesh  => proc%blocks(iblk)%mesh        , &
                 ph    => proc%blocks(iblk)%state(1)%ph , &
                 v     => proc%blocks(iblk)%state(1)%v)
        select case (bkg_type)
        case ('era5')
          allocate(tmp(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,num_era5_lev))
          do k = 1, num_era5_lev
            call latlon_interp_bilinear_lat_edge(era5_lon, era5_lat, era5_v(:,:,k), mesh, tmp(:,:,k))
          end do
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              call vert_interp_linear(era5_lev, tmp(i,j,:), ph(i,j,:), v(i,j,:), allow_extrap=.true.)
            end do
          end do
          deallocate(tmp)
          call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
        end select
      end associate
    end do

  end subroutine bkg_regrid_v

end module bkg_mod
