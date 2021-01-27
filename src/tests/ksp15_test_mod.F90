module ksp15_test_mod

  use namelist_mod
  use const_mod
  use parallel_mod
  use block_mod
  use vert_coord_mod
  use formula_mod

  implicit none

  private

  public ksp15_01_test_set_ic

  real(r8), parameter :: teq  = 300.0_r8     ! K
  real(r8), parameter :: peq  = 1.0e5_r8     ! Pa
  real(r8), parameter :: ueq  = 20.0_r8      ! m s-1
  real(r8), parameter :: X    = 166.7_r8
  real(r8), parameter :: h0   = 250.0_r8     ! m
  real(r8), parameter :: d0   = 5000.0_r8    ! m
  real(r8), parameter :: xi0  = 4000.0_r8    ! m
  real(r8), parameter :: lonc = pi

contains

  subroutine ksp15_01_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    real(r8) dlon, r0
    integer i, j, k

    omega = 0.0_r8

    associate (mesh   => block%mesh           , &
               u      => block%state(1)%u     , &
               v      => block%state(1)%v     , &
               t      => block%state(1)%t     , &
               pt     => block%state(1)%pt    , &
               gzs    => block%static  %gzs   , &
               phs    => block%state(1)%phs   , &
               ph     => block%state(1)%ph    , &
               ph_lev => block%state(1)%ph_lev, &
               gz_lev => block%state(1)%gz_lev)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            u(i,j,k) = ueq * mesh%full_cos_lat(j)
          end do
        end do
      end do
      call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dlon = abs(mesh%full_lon(i) - lonc)
          dlon = min(pi2 - dlon, dlon)
          r0 = radius / X * dlon
          gzs(i,j) = g * h0 * exp(-r0**2 / d0**2) * cos(pi * r0 / xi0)**2 * mesh%full_cos_lat(j)
        end do
      end do
      call fill_halo(block, gzs, full_lon=.true., full_lat=.true.)

      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          phs(i,j) = peq * exp(-0.5_r8 * ueq**2 / Rd / teq * mesh%full_sin_lat(j)**2 - gzs(i,j) / Rd / teq)
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
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            t(i,j,k) = teq
            pt(i,j,k) = potential_temperature(t(i,j,k), ph(i,j,k))
          end do
        end do
      end do
      call fill_halo(block, t , full_lon=.true., full_lat=.true., full_lev=.true.)
      call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)

      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            gz_lev(i,j,k) = Rd * teq * log(peq / ph_lev(i,j,k)) - 0.5_r8 * ueq**2 * mesh%full_sin_lat(j)**2
          end do
        end do
      end do
      call fill_halo(block, gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
    end associate

  end subroutine ksp15_01_test_set_ic

end module ksp15_test_mod
