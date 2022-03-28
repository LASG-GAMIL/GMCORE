module ffsl_mod

  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
  use tracer_mod
  use transport_batch_mod

  implicit none

  private

  public ffsl_init

  interface
    subroutine flux_interface(block, batch, tracer)
      import block_type, transport_batch_type, tracer_type
      type(block_type          ), intent(in   ) :: block
      type(transport_batch_type), intent(in   ) :: batch
      type(tracer_type         ), intent(inout) :: tracer
    end subroutine flux_interface
    pure real(r8) function slope_interface(f)
      import r8
      real(r8), intent(in) :: f(-1:1)
    end function slope_interface
  end interface

  procedure(flux_interface ), pointer :: mass_flux   => null()
  procedure(flux_interface ), pointer :: tracer_flux => null()
  procedure(slope_interface), pointer :: slope       => null()

contains

  subroutine ffsl_init()

    select case (ffsl_flux_type)
    case ('van_leer')
      mass_flux   => mass_flux_van_leer
      tracer_flux => tracer_flux_van_leer
    case ('ppm')
      mass_flux   => mass_flux_ppm
      tracer_flux => tracer_flux_ppm
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

  subroutine ffsl_calc_mass_flux(block, state, batch)

    type(block_type          ), intent(in) :: block
    type(state_type          ), intent(in) :: state
    type(transport_batch_type), intent(in) :: batch

    associate (mfx => batch%mfx, mfy => batch%mfy)
    end associate

  end subroutine ffsl_calc_mass_flux

  subroutine mass_flux_van_leer(block, batch, tracer)

    type(block_type          ), intent(in   ) :: block
    type(transport_batch_type), intent(in   ) :: batch
    type(tracer_type         ), intent(inout) :: tracer

    integer i, j, k, iu, ju, ci, i1, i2
    real(r8) cf, dm

    associate (mesh => block%mesh, &
               cx   => batch%cx  , & ! in
               cy   => batch%cy  , & ! in
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               mx   => tracer%qx , & ! in
               my   => tracer%qy , & ! in
               mfx  => tracer%mfx, & ! out
               mfy  => tracer%mfy)   ! out
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ci = int(cx(i,j,k))
          cf = cx(i,j,k) - ci
          i1 = i - merge(0, ci, ci > 0) + 1
          i2 = i + merge(ci, 0, ci > 0)
          mfx(i,j,k) = u(i,j,k) * sum(mx(i1:i2,j,k))
          iu = merge(i - ci, i + 1 + ci, cf > 0)
          dm = slope(mx(iu-1:iu+1,j,k))
          mfx(i,j,k) = mfx(i,j,k) + u(i,j,k) * (mx(iu,j,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf))
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          cf = cy(i,j,k)
          ju = merge(j, j + 1, cf > 0)
          dm = slope(my(i,ju-1:ju+1,k))
          mfx(i,j,k) = v(i,j,k) * (my(i,ju,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf))
        end do
      end do
    end do
    call fill_halo(block, mfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    end associate

  end subroutine mass_flux_van_leer

  subroutine tracer_flux_van_leer(block, batch, tracer)

    type(block_type          ), intent(in   ) :: block
    type(transport_batch_type), intent(in   ) :: batch
    type(tracer_type         ), intent(inout) :: tracer

    integer i, j, k, iu, ju, ci, i1, i2
    real(r8) cf, dq

    associate (mesh => block%mesh, &
               cx   => batch%cx  , & ! in
               cy   => batch%cy  , & ! in
               mfx  => batch%mfx , & ! in
               mfy  => batch%mfy , & ! in
               qx   => tracer%qx , & ! in
               qy   => tracer%qy , & ! in
               qmfx => tracer%mfx, & ! out
               qmfy => tracer%mfy)   ! out
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ci = int(cx(i,j,k))
          cf = cx(i,j,k) - ci
          i1 = i - merge(0, ci, ci > 0) + 1
          i2 = i + merge(ci, 0, ci > 0)
          qmfx(i,j,k) = sum(qx(i1:i2,j,k))
          iu = merge(i - ci, i + 1 + ci, cf > 0)
          dq = slope(qx(iu-1:iu+1,j,k))
          qmfx(i,j,k) = mfx(i,j,k) * (qx(iu,j,k) + dq * 0.5_r8 * (sign(1.0_r8, cf) - cf))
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          cf = cy(i,j,k)
          ju = merge(j, j + 1, cf > 0)
          dq = slope(qy(i,ju-1:ju+1,k))
          qmfy(i,j,k) = mfy(i,j,k) * (qy(i,ju,k) + dq * 0.5_r8 * (sign(1.0_r8, cf) - cf))
        end do
      end do
    end do
    call fill_halo(block, qmfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, qmfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    end associate

  end subroutine tracer_flux_van_leer

  subroutine mass_flux_ppm(block, batch, tracer)

    type(block_type          ), intent(in   ) :: block
    type(transport_batch_type), intent(in   ) :: batch
    type(tracer_type         ), intent(inout) :: tracer

    integer i, j, k, iu, ju, ci, i1, i2
    real(r8) cf, s1, s2, ds1, ds2, ds3

    associate (mesh => block%mesh, &
               cx   => batch%cx  , & ! in
               cy   => batch%cy  , & ! in
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               mx   => tracer%qx , & ! in
               my   => tracer%qy , & ! in
               mxl  => tracer%qxl, & ! out
               myl  => tracer%qyl, & ! out
               dmx  => tracer%dqx, & ! out
               dmy  => tracer%dqy, & ! out
               mx6  => tracer%qx6, & ! out
               my6  => tracer%qy6, & ! out
               mfx  => tracer%mfx, & ! out
               mfy  => tracer%mfy)   ! out
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
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
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend

        end do
      end do
    end do
    end associate
  end subroutine mass_flux_ppm

  subroutine tracer_flux_ppm(block, batch, tracer)

    type(block_type          ), intent(in   ) :: block
    type(transport_batch_type), intent(in   ) :: batch
    type(tracer_type         ), intent(inout) :: tracer

    integer i, j, k, iu, ju, ci, i1, i2
    real(r8) cf, s1, s2, ds1, ds2, ds3

    associate (mesh => block%mesh, &
               cx   => batch%cx  , & ! in
               cy   => batch%cy  , & ! in
               mfx  => tracer%mfx, & ! in
               mfy  => tracer%mfy, & ! in
               qx   => tracer%qx , & ! in
               qy   => tracer%qy , & ! in
               qxl  => tracer%qxl, & ! out
               qyl  => tracer%qyl, & ! out
               dqx  => tracer%dqx, & ! out
               dqy  => tracer%dqy, & ! out
               qx6  => tracer%qx6, & ! out
               qy6  => tracer%qy6, & ! out
               qmfx => tracer%mfx, & ! out
               qmfy => tracer%mfy)   ! out
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          call ppm(qx(i-2:i+2,j,k), qxl(i,j,k), dqx(i,j,k), qx6(i,j,k))
          call ppm(qy(i,j-2:j+2,k), qyl(i,j,k), dqy(i,j,k), qy6(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, qxl, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, dqx, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, qx6, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, qyl, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, dqy, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, qy6, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        end do
      end do
    end do
    call fill_halo(block, mfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    end associate

  end subroutine tracer_flux_ppm

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
