module dcmip31_test_mod

  use flogger
  use namelist_mod
  use const_mod
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod
  use operators_mod
  use diag_state_mod

  implicit none

  private

  public dcmip31_test_set_params
  public dcmip31_test_set_ic
  public dcmip31_test_init_diag_state

  real(r8), parameter :: X    = 125.0_r8
  real(r8), parameter :: ptop = 273.919e2_r8 ! Pa
  real(r8), parameter :: u0   = 20.0_r8      ! m s-1
  real(r8), parameter :: teq  = 300.0_r8     ! K
  real(r8), parameter :: peq  = 1.0e5_r8     ! Pa
  real(r8), parameter :: d    = 5000.0_r8    ! m
  real(r8), parameter :: d2   = d**2
  real(r8), parameter :: lonc = 2 * pi / 3
  real(r8), parameter :: latc = 0.0_r8
  real(r8), parameter :: dpt  = 1.0_r8       ! K
  real(r8), parameter :: lz   = 20000.0_r8   ! m
  real(r8), parameter :: N    = 0.01         ! s-1
  real(r8), parameter :: N2   = N**2
  real(r8), parameter :: t0   = g**2 / N**2 / cp
  real(r8), parameter :: p0   = 1.0e5_r8     ! Pa

contains

  subroutine dcmip31_test_set_params()

    omega = 0.0_r8
    radius = radius / X

  end subroutine dcmip31_test_set_params

  subroutine dcmip31_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, k
    real(r8) ts, cos_2lat, r

    associate (mesh   => block%mesh           , &
               u      => block%state(1)%u     , &
               v      => block%state(1)%v     , &
               w      => block%state(1)%w     , &
               phs    => block%state(1)%phs   , &
               ph_lev => block%state(1)%ph_lev, &
               ph     => block%state(1)%ph    , &
               pt     => block%state(1)%pt    , &
               t      => block%state(1)%t     , &
               gz     => block%state(1)%gz    , &
               gzs    => block%static%gzs)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            u(i,j,k) = u0 * mesh%full_cos_lat(j)
          end do
        end do
      end do
      call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

      v = 0.0_r8
      gzs = 0.0_r8

      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        cos_2lat = cos(2 * mesh%full_lat(j))
        ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          phs(i,j) = peq * exp(u0 / (4 * t0 * Rd) * (u0 + 2 * omega * radius) * (cos_2lat - 1)) * &
                     (ts / teq)**(1 / Rd_o_cp)
        end do
      end do
      call fill_halo(block, phs, full_lon=.true., full_lat=.true.)

      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, phs(i,j))
          end do
        end do
      end do
      call fill_halo(block, ph_lev, full_lon=.true., full_lat=.true., full_lev=.false.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ph(i,j,k) = 0.5d0 * (ph_lev(i,j,k) + ph_lev(i,j,k+1))
          end do
        end do
      end do
      call fill_halo(block, ph, full_lon=.true., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          cos_2lat = cos(2 * mesh%full_lat(j))
          ts = t0 + (teq - t0) * exp(-u0 * N2 / (4 * g**2) * (u0 + 2 * omega * radius) * (cos_2lat - 1))
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pt(i,j,k) = ts * (p0 / phs(i,j))**Rd_o_cp / (         &
              ts / t0 * ((ph(i,j,k) / phs(i,j))**Rd_o_cp - 1) + 1 &
            )
            t(i,j,k) = temperature(pt(i,j,k), ph(i,j,k))
          end do
        end do
      end do
      call fill_halo(block, t, full_lon=.true., full_lat=.true., full_lev=.true.)

      if (nonhydrostatic) then
        w = 0.0_r8
        call diag_gz_lev(block, block%state(1))
        call interp_gz(block, block%state(1))
      end if

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ! Perturbation
            r = radius * acos(sin(latc) * mesh%full_sin_lat(j) + cos(latc) * mesh%full_cos_lat(j) * cos(mesh%full_lon(i) - lonc))
            pt(i,j,k) = pt(i,j,k) + dpt * d2 / (d2 + r**2) * sin(pi2 * gz(i,j,k) / g / lz)
          end do
        end do
      end do
      call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate
  
  end subroutine dcmip31_test_set_ic

  subroutine dcmip31_test_init_diag_state(blocks)

    type(block_type), intent(in) :: blocks(:)

    real(r8) z(11)
    integer iblk

    z = [0.0_r8, 1.0e3_r8, 2.0e3_r8, 3.0e3_r8, 4.0e3_r8, 5.0e3_r8, 6.0e3_r8, 7.0e3_r8, 8.0e3_r8, 9.0e3_r8, 10.0e3_r8]

    do iblk = 1, size(blocks)
      call diag_state(iblk)%init_height_levels(blocks(iblk), z, instance)
    end do

  end subroutine dcmip31_test_init_diag_state
  
end module dcmip31_test_mod
