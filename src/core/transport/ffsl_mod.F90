module ffsl_mod

  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
  use tracer_mod

  implicit none

  private

  public ffsl_init

  interface
    subroutine ffsl_interface(block, tracer, dt, mfx, mfy, u, v, qx, qy, qfx, qfy)
      import block_type, tracer_type, r8
      type(block_type ), intent(in   ) :: block
      type(tracer_type), intent(inout) :: tracer
      real(r8), intent(in ) :: dt
      real(r8), intent(in ) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, & ! Mass flux along x-axis
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, & ! Mass flux along y-axis
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: qx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, & ! Tracer mixing ratio transported along x-axis
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: qy (block%mesh%full_lon_lb:block%mesh%full_lon_ub, & ! Tracer mixing ratio transported along y-axis
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: qfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, & ! Tracer mass flux along x-axis
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: qfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, & ! Tracer mass flux along y-axis
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    end subroutine ffsl_interface
    pure real(r8) function slope_interface(f)
      import r8
      real(r8), intent(in) :: f(-1:1)
    end function slope_interface
  end interface

  procedure(ffsl_interface), pointer :: ffsl => null()
  procedure(slope_interface), pointer :: slope => null()

contains

  subroutine ffsl_init()

    select case (ffsl_flux_type)
    case ('van_leer')
      ffsl => ffsl_van_leer
    case ('ppm')
      ffsl => ffsl_ppm
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

  subroutine ffsl_run(block, tracer, dt, mfx, mfy, u, v)

    type(block_type ), intent(in   ) :: block
    type(tracer_type), intent(inout) :: tracer
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k

    associate (mesh => block%mesh, &
               q    => tracer%q  , &
               qx   => tracer%qx , &
               qy   => tracer%qy , &
               qfx  => tracer%qfx, &
               qfy  => tracer%qfy)
    ! --------------------------------------------------------------------------
    ! Run inner operators.
    call ffsl(block, tracer, dt, mfx, mfy, u, v, q , q , qfx, qfy)
    ! --------------------------------------------------------------------------
    ! Run outer operators.
    call ffsl(block, tracer, dt, mfx, mfy, u, v, qy, qx, qfx, qfy)
    end associate

  end subroutine ffsl_run

  subroutine ffsl_van_leer(block, tracer, dt, mfx, mfy, u, v, qx, qy, qfx, qfy)

    type(block_type ), intent(in   ) :: block
    type(tracer_type), intent(inout) :: tracer
    real(r8), intent(in ) :: dt
    real(r8), intent(in ) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: qx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: qy (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    real(r8) c, dq
    integer i, j, k, iu, ju

    associate (mesh => block%mesh)
    ! Along x-axis
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          iu = merge(i - 1, i, u(i,j,k) > 0)
          c  = u(i,j,k) * dt / mesh%de_lon(j)
          dq = slope(qx(iu-1:iu+1,j,k))
          qfx(i,j,k) = mfx(i,j,k) * (qx(iu,j,k) + dq * 0.5_r8 * (sign(1.0_r8, c) - c))
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ju = merge(j - 1, j, v(i,j,k) > 0)
          c  = v(i,j,k) * dt / mesh%de_lat(ju)
          dq = slope(qy(i,ju-1:ju+1,k))
          qfy(i,j,k) = mfy(i,j,k) * (qy(i,ju,k) + dq * 0.5_r8 * (sign(1.0_r8, c) - c))
        end do
      end do
    end do
    call fill_halo(block, qfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, qfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    end associate

  end subroutine ffsl_van_leer

  subroutine ffsl_ppm(block, tracer, dt, mfx, mfy, u, v, qx, qy, qfx, qfy)

    type(block_type ), intent(in   ) :: block
    type(tracer_type), intent(inout) :: tracer
    real(r8), intent(in ) :: dt
    real(r8), intent(in ) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: qx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: qy (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    real(r8) c, s1, s2, ds1, ds2, ds3
    integer i, j, k, iu, ju

    associate (mesh => block%mesh, &
               qxl  => tracer%qxl, &
               qyl  => tracer%qyl, &
               dqx  => tracer%dqx, &
               dqy  => tracer%dqy, &
               qx6  => tracer%qx6, &
               qy6  => tracer%qy6)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          call ppm(qx(i-2:i+2,j,k), qxl(i,j), dqx(i,j), qx6(i,j))
          call ppm(qy(i,j-2:j+2,k), qyl(i,j), dqy(i,j), qy6(i,j))
        end do
      end do
      call fill_halo(block, qxl, full_lon=.true., full_lat=.true., south_halo=.false., north_halo=.false.)
      call fill_halo(block, dqx, full_lon=.true., full_lat=.true., south_halo=.false., north_halo=.false.)
      call fill_halo(block, qx6, full_lon=.true., full_lat=.true., south_halo=.false., north_halo=.false.)
      call fill_halo(block, qyl, full_lon=.true., full_lat=.true.,  west_halo=.false.,  east_halo=.false.)
      call fill_halo(block, dqy, full_lon=.true., full_lat=.true.,  west_halo=.false.,  east_halo=.false.)
      call fill_halo(block, qy6, full_lon=.true., full_lat=.true.,  west_halo=.false.,  east_halo=.false.)
      ! Along x-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          c = u(i,j,k) * dt / mesh%de_lon(j)
          iu = merge(i - 1, i, c > 0)
          s1 = merge(1 - abs(c), 0.0_r8, c >= 0)
          s2 = merge(1.0_r8, abs(c), c >= 0)
          ds1 = s2    - s1
          ds2 = s2**2 - s1**2
          ds3 = s2**3 - s1**3
          qfx(i,j,k) = sign(qxl(iu,j) * ds1 + 0.5_r8 * dqx(iu,j) * ds2 + qx6(iu,j) * (0.5_r8 * ds2 - ds3 / 3.0_r8), c)
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          c = v(i,j,k) * dt / mesh%de_lat(j)
          ju = merge(j - 1, j, c > 0)
          s1 = merge(1 - abs(c), 0.0_r8, c >= 0)
          s2 = merge(1.0_r8, abs(c), c >= 0)
          ds1 = s2    - s1
          ds2 = s2**2 - s1**2
          ds3 = s2**3 - s1**3
          qfy(i,j,k) = sign(qyl(i,ju) * ds1 + 0.5_r8 * dqy(i,ju) * ds2 + qy6(i,ju) * (0.5_r8 * ds2 - ds3 / 3.0_r8), c)
        end do
      end do
    end do
    call fill_halo(block, qfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, qfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    end associate

  end subroutine ffsl_ppm

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
