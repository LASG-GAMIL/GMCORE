module ffsl_mod

  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
  use adv_batch_mod

  implicit none

  private

  public ffsl_init
  public ffsl_calc_mass_flux
  public ffsl_calc_tracer_flux

  interface
    subroutine flux_interface(block, batch, u, v, mx, my, mfx, mfy)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: mx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: my (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    end subroutine flux_interface
    pure real(r8) function slope_interface(f)
      import r8
      real(r8), intent(in) :: f(-1:1)
    end function slope_interface
  end interface

  procedure(flux_interface ), pointer :: flux  => null()
  procedure(slope_interface), pointer :: slope => null()

contains

  subroutine ffsl_init()

    select case (ffsl_flux_type)
    case ('van_leer')
      flux => flux_van_leer
    case ('ppm')
      flux => flux_ppm
    end select

    select case (limiter_type)
    case ('none')
      slope => slope_simple
    case ('mono')
      slope => slope_mono
    case ('pd')
      slope => slope_pd
    end select

  end subroutine ffsl_init

  subroutine ffsl_calc_mass_flux(block, batch, m, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: m  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k
    real(r8) work(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend,block%mesh%num_full_lev)
    real(r8) pole(block%mesh%num_full_lev)

    associate (mesh => block%mesh, &
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               mx   => batch%qx  , & ! work array
               my   => batch%qy)     ! work array
    ! Run inner advective operators.
    call flux(block, batch, u, v, m, m, mfx, mfy)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ! Subtract divergence terms from flux to form advective operators.
          mx(i,j,k) = m(i,j,k) - 0.5_r8 * (mfx(i,j,k) - mfx(i-1,j,k) - divx(i,j,k) * m(i,j,k))
          my(i,j,k) = m(i,j,k) - 0.5_r8 * (mfy(i,j,k) - mfy(i,j-1,k) - divy(i,j,k) * m(i,j,k))
        end do
      end do
    end do
    ! Handle the Pole boundary conditions.
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = mfy(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          mx(i,j,k) = m(i,j,k)
          my(i,j,k) = m(i,j,k) - 0.5_r8 * (pole(k) - divy(i,j,k) * m(i,j,k))
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = mfy(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          mx(i,j,k) = m(i,j,k)
          my(i,j,k) = m(i,j,k) + 0.5_r8 * (pole(k) - divy(i,j,k) * m(i,j,k))
        end do
      end do
    end if
    call fill_halo(block, mx, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, my, full_lon=.true., full_lat=.true., full_lev=.true.)
    ! Run outer flux form operators.
    call flux(block, batch, u, v, my, mx, mfx, mfy)
    end associate

  end subroutine ffsl_calc_mass_flux

  subroutine ffsl_calc_tracer_flux(block, batch, q, qmfx, qmfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: q    (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qmfx (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qmfy (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k
    real(r8) work(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend,block%mesh%num_full_lev)
    real(r8) pole(block%mesh%num_full_lev)

    associate (mesh => block%mesh, &
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               mfx  => batch%mfx , & ! in
               mfy  => batch%mfy , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               qx   => batch%qx  , & ! work array
               qy   => batch%qy)     ! work array
    ! Run inner advective operators.
    call flux(block, batch, u, v, q, q, qmfx, qmfy)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ! Subtract divergence terms from flux to form advective operators.
          qx(i,j,k) = q(i,j,k) - 0.5_r8 * (qmfx(i,j,k) - qmfx(i-1,j,k) - divx(i,j,k) * q(i,j,k))
          qy(i,j,k) = q(i,j,k) - 0.5_r8 * (qmfy(i,j,k) - qmfy(i,j-1,k) - divy(i,j,k) * q(i,j,k))
        end do
      end do
    end do
    ! Handle the Pole boundary conditions.
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = qmfy(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          qx(i,j,k) = q(i,j,k)
          qy(i,j,k) = q(i,j,k) - 0.5_r8 * (pole(k) - divy(i,j,k) * q(i,j,k))
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = qmfy(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          qx(i,j,k) = q(i,j,k)
          qy(i,j,k) = q(i,j,k) + 0.5_r8 * (pole(k) - divy(i,j,k) * q(i,j,k))
        end do
      end do
    end if
    call fill_halo(block, qx, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, qy, full_lon=.true., full_lat=.true., full_lev=.true.)
    ! Run outer flux form operators.
    call flux(block, batch, mfx, mfy, qy, qx, qmfx, qmfy)
    end associate

  end subroutine ffsl_calc_tracer_flux

  subroutine flux_van_leer(block, batch, u, v, mx, my, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: mx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: my (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k, iu, ju, ci, i1, i2
    real(r8) cf, sm, dm

    associate (mesh => block %mesh, &
               cflx => batch %cflx, & ! in
               cfly => batch %cfly)   ! in
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          i1 = i - merge(ci, 0, ci > 0) + 1
          i2 = i - merge(0, ci, ci > 0)
          sm = sum(mx(i1:i2,j,k))
          iu = merge(i - ci, i - ci + 1, cf > 0)
          dm = slope(mx(iu-1:iu+1,j,k))
          mfx(i,j,k) = cf * (mx(iu,j,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf)) + sm
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          cf = cfly(i,j,k)
          ju = merge(j, j + 1, cf > 0)
          dm = slope(my(i,ju-1:ju+1,k))
          mfy(i,j,k) = cf * (my(i,ju,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf))
        end do
      end do
    end do
    call fill_halo(block, mfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    end associate

  end subroutine flux_van_leer

  subroutine flux_ppm(block, batch, u, v, mx, my, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: mx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: my (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k, iu, ju, ci, i1, i2
    real(r8) cf, s1, s2, ds1, ds2, ds3, sm

    associate (mesh => block %mesh, &
               cflx => batch %cflx, & ! in
               cfly => batch %cfly, & ! in
               mxl  => batch %qxl , & ! work array
               myl  => batch %qyl , & ! work array
               dmx  => batch %dqx , & ! work array
               dmy  => batch %dqy , & ! work array
               mx6  => batch %qx6 , & ! work array
               my6  => batch %qy6 )   ! work array
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          call ppm(mx(i-2:i+2,j,k), mxl(i,j,k), dmx(i,j,k), mx6(i,j,k))
          call ppm(my(i,j-2:j+2,k), myl(i,j,k), dmy(i,j,k), my6(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, mxl, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, dmx, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mx6, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, myl, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, dmy, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, my6, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          i1 = i - merge(ci, 0, ci > 0) + 1
          i2 = i - merge(0, ci, ci > 0)
          sm = sum(mx(i1:i2,j,k))
          iu = merge(i - ci, i - ci + 1, cf > 0)
          s1 = merge(1 - abs(cf), 0.0_r8, cf >= 0)
          s2 = merge(1.0_r8, abs(cf), cf >= 0)
          ds1 = s2    - s1
          ds2 = s2**2 - s1**2
          ds3 = s2**3 - s1**3
          mfx(i,j,k) = sign(mxl(iu,j,k) * ds1 + 0.5_r8 * dmx(iu,j,k) * ds2 + mx6(iu,j,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8), cf) + sm
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          cf = cfly(i,j,k)
          ju = merge(j, j + 1, cf > 0)
          s1 = merge(1 - abs(cf), 0.0_r8, cf >= 0)
          s2 = merge(1.0_r8, abs(cf), cf >= 0)
          ds1 = s2    - s1
          ds2 = s2**2 - s1**2
          ds3 = s2**3 - s1**3
          mfy(i,j,k) = sign(myl(i,ju,k) * ds1 + 0.5_r8 * dmy(i,ju,k) * ds2 + my6(i,ju,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8), cf)
        end do
      end do
    end do
    call fill_halo(block, mfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    end associate

  end subroutine flux_ppm

  subroutine ppm(f, fl, df, f6)

    real(r8), intent(in ) :: f(-2:2)
    real(r8), intent(out) :: fl
    real(r8), intent(out) :: df
    real(r8), intent(out) :: f6

    real(r8) dfl, dfr, fr

    ! Calculate values at left and right cell interfaces.
    dfl = slope(f(-2:0))
    df  = slope(f(-1:1))
    dfr = slope(f( 0:2))
    ! Why (B2) in Lin (2004) divide (dfl - df) and (df - dfr) by 3?
    fl = 0.5_r8 * (f(-1) + f(0)) + (dfl - df) / 6.0_r8
    fr = 0.5_r8 * (f( 1) + f(0)) + (df - dfr) / 6.0_r8
    ! Why (B3) and (B4) in Lin (2004) multiply df by 2?
    fl = f(0) - sign(min(abs(df), abs(fl - f(0))), df)
    fr = f(0) + sign(min(abs(df), abs(fr - f(0))), df)
    f6 = 6 * f(0) - 3 * (fl + fr)
    df = fr - fl

  end subroutine ppm

  pure real(r8) function slope_simple(f) result(res)

    real(r8), intent(in) :: f(-1:1) 

    res = (f(1) - f(-1)) * 0.5_r8

  end function slope_simple

  pure real(r8) function slope_mono(f) result(res)

    real(r8), intent(in) :: f(-1:1)

    real(r8) df, df_min, df_max

    df = (f(1) - f(-1)) * 0.5_r8 ! Initial guess
    df_min = 2 * (f(0) - minval(f))
    df_max = 2 * (maxval(f) - f(0))
    res = sign(min(abs(df), df_min, df_max), df)

  end function slope_mono

  pure real(r8) function slope_pd(f) result(res)

    real(r8), intent(in) :: f(-1:1)

    real(r8) df

    df = (f(1) - f(-1)) * 0.5_r8 ! Initial guess
    res = sign(min(abs(df), 2 * f(0)), df)

  end function slope_pd

end module ffsl_mod
