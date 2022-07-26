module solid_rotation_test_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use parallel_mod
  use block_mod
  use adv_mod

  implicit none

  public solid_rotation_test_init
  public solid_rotation_test_set_ic
  public solid_rotation_test_set_uv

  real(r8), parameter :: period = 12 * 86400
  real(r8), parameter :: h0     = 1000 ! m
  real(r8), parameter :: lon0   = 3 * pi / 2.0_r8
  real(r8), parameter :: lat0   = 0
  real(r8), parameter :: alpha  = 90.0_r8 * rad
  real(r8) u0

contains

  subroutine solid_rotation_test_init()

    u0 = pi2 * radius / period

  end subroutine solid_rotation_test_init

  subroutine solid_rotation_test_set_ic(block)

    type(block_type), intent(inout) :: block

    integer i, j
    real(r8) lon, lat, r, R0

    call adv_add_tracer('solid_rotation', dt, 'q0', 'background tracer')
    call adv_add_tracer('solid_rotation', dt, 'q1', 'cosine bell tracer')

    call adv_allocate_tracers(block)

    R0 = radius / 3.0_r8

    associate (mesh => block%mesh              , &
               old  => block%adv_batches(1)%old, &
               q    => block%adv_batches(1)%q  )
    ! Background
    q(:,:,:,1,old) = 1
    ! Cosine bell
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        r = calc_distance(lon0, lat0, lon, lat)
        if (r < R0) then
          q(i,j,1,2,old) = h0 / 2.0_r8 * (1 + cos(pi * r / R0))
        else
          q(i,j,1,2,old) = 0
        end if
      end do
    end do
    call fill_halo(block, q(:,:,:,2,old), full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine solid_rotation_test_set_ic

  subroutine solid_rotation_test_set_uv(block, state, time_in_seconds)

    type(block_type), intent(in   ) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: time_in_seconds

    integer i, j
    real(r8) lon, lat

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
        u(i,j,1) = u0 * (cos(lat) * cos(alpha) + sin(lat) * cos(lon) * sin(alpha))
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    mfx = u
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      lat = mesh%half_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        v(i,j,1) = -u0 * sin(lon) * sin(alpha)
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    mfy = v
    end associate

  end subroutine solid_rotation_test_set_uv

end module solid_rotation_test_mod
