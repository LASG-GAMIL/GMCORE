module moving_vortices_test_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use parallel_mod
  use history_mod
  use block_mod
  use adv_mod

  implicit none

  private

  public moving_vortices_test_init
  public moving_vortices_test_set_ic
  public moving_vortices_test_set_uv

  real(r8), parameter :: period = 12 * 86400    ! Rotation period
  real(r8), parameter :: rho0   = 3             ! 
  real(r8), parameter :: gamma  = 5             ! Field stiffness parameter
  real(r8), parameter :: alpha  = pi05          ! Angle between Rotation axis and the equator
  real(r8), parameter :: lonp0  = pi            ! Rotation north pole longitude
  real(r8), parameter :: latp0  = pi05 - alpha  ! Rotation north pole latitude
  real(r8), parameter :: lonv0  = pi * 1.5_r8   ! Initial vortex longitude
  real(r8), parameter :: latv0  = 0             ! Initial vortex latitude
  real(r8) u0                                   ! Velocity amplitude
  real(r8) lonvr0                               ! Initial vortex longitude in rotated coordinate
  real(r8) lonvr                                ! Vortex longitude in rotated coordinate
  real(r8) latvr                                ! Vortex latitude  in rotated coordinate
  real(r8) lonv                                 ! Vortex longitude
  real(r8) latv                                 ! Vortex latitude

contains

  subroutine moving_vortices_test_init()

    u0 = pi2 * radius / (12.0_r8 * 86400.0_r8)
    
    call rotation_transform(lonp0, latp0, lonv0, latv0, lonvr, latvr)
    lonv = lonv0
    latv = latv0
    lonvr0 = lonvr

  end subroutine moving_vortices_test_init

  subroutine moving_vortices_test_set_ic(block)

    type(block_type), intent(inout) :: block

    integer i, j
    real(r8) lon, lat, lonr, latr

    call adv_add_tracer('moving_vortices', dt, 'q0', 'background tracer')
    call adv_add_tracer('moving_vortices', dt, 'q1', 'vortex tracer')

    call adv_allocate_tracers(block)

    associate (mesh => block%mesh              , &
               old  => block%adv_batches(1)%old, &
               q    => block%adv_batches(1)%q  )
    ! Background tracer
    q(:,:,:,1,old) = 1
    ! Vortex tracer 
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        call rotation_transform(lonv0, latv0, lon, lat, lonr, latr)
        q(i,j,1,2,old) = 1 - tanh(rho(latr) / gamma * sin(lonr))
      end do
    end do
    call fill_halo(block, q(:,:,:,2,old), full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine moving_vortices_test_set_ic

  subroutine moving_vortices_test_set_uv(block, state, time_in_seconds)

    type(block_type), intent(in   ) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: time_in_seconds

    integer i, j
    real(r8) lon, lat, dlon, latr

    lonvr = lonvr0 + u0 / radius * time_in_seconds
    if (lonvr > pi2) lonvr = lonvr - pi2
    call inverse_rotation_transform(lonp0, latp0, lonv, latv, lonvr, latvr)

    associate (mesh => block%mesh   , &
               m    => state%m      , &
               u    => state%u_lon  , &
               v    => state%v_lat  , &
               mfx  => state%mfx_lon, &
               mfy  => state%mfy_lat)
    m = 1
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      lat = mesh%full_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        lon = mesh%half_lon(i)
        dlon = lon - lonv
        call rotation_transform(lonv, latv, lon, lat, lat_r=latr)
        u(i,j,1) = u0 * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha)) + &
                   a_omg(latr) * (sin(latv) * cos(lat) - cos(latv) * cos(dlon) * sin(lat))
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    mfx = u
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      lat = mesh%half_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        dlon = lon - lonv
        call rotation_transform(lonv, latv, lon, lat, lat_r=latr)
        v(i,j,1) = -u0 * sin(lon) * sin(alpha) + a_omg(latr) * cos(latv) * sin(dlon)
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    mfy = v
    end associate

  end subroutine moving_vortices_test_set_uv

  real(r8) function rho(lat) result(res)

    real(r8), intent(in) :: lat

    res = rho0 * cos(lat)

  end function rho

  real(r8) function a_omg(latr) result(res)

    real(r8), intent(in) :: latr

    real(r8), parameter :: c = 1.5_r8 * sqrt(3.0_r8)
    real(r8) r

    r = rho(latr)
    if (abs(r) < eps) then
      res = 0
    else
      res = u0 * c / cosh(r)**2 * tanh(r) / r
    end if

  end function a_omg

end module moving_vortices_test_mod
