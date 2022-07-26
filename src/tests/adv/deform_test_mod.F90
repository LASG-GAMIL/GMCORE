module deform_test_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use parallel_mod
  use history_mod
  use block_mod
  use adv_mod

  implicit none

  private

  public deform_test_init
  public deform_test_set_ic
  public deform_case1_test_set_uv
  public deform_case2_test_set_uv
  public deform_case3_test_set_uv
  public deform_case4_test_set_uv

  real(r8), parameter :: period = 12 * 86400
  real(r8) lon1, lat1, lon2, lat2

contains

  subroutine deform_test_init(case_n)

    integer, intent(in) :: case_n

    select case (case_n)
    case (1)
      lon1 = pi
      lat1 = pi / 3.0_r8
      lon2 = pi
      lat2 = -pi / 3.0_r8
    case (3)
      lon1 = pi * 3.0_r8 / 4.0_r8
      lat1 = 0
      lon2 = pi * 5.0_r8 / 4.0_r8
      lat2 = 0
    case (2, 4)
      lon1 = pi * 5.0_r8 / 6.0_r8
      lat1 = 0
      lon2 = pi * 7.0_r8 / 6.0_r8
      lat2 = 0
    end select

  end subroutine deform_test_init

  subroutine deform_test_set_ic(block)

    type(block_type), intent(inout) :: block

    integer i, j
    real(r8) lon, lat, r, r1, r2, qmax, qmin, c
    real(r8) x(3), x1(3), x2(3)

    call adv_add_tracer('deform_test', dt, 'q0', 'background tracer'       )
    call adv_add_tracer('deform_test', dt, 'q1', 'cosine hills tracer'     )
    call adv_add_tracer('deform_test', dt, 'q2', 'slotted cylinders tracer')
    call adv_add_tracer('deform_test', dt, 'q3', 'gaussian hills tracer'   )

    call adv_allocate_tracers(block)

    call cartesian_transform(lon1, lat1, x1(1), x1(2), x1(3)); x1 = x1 / radius
    call cartesian_transform(lon2, lat2, x2(1), x2(2), x2(3)); x2 = x2 / radius

    associate (mesh => block%mesh              , &
               old  => block%adv_batches(1)%old, &
               q    => block%adv_batches(1)%q  )
    ! Background
    q(:,:,:,1,old) = 1
    ! Cosine hills
    qmax = 1.0_r8; qmin = 0.1_r8; c = 0.9_r8; r = radius * 0.5_r8
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        r1 = calc_distance(lon1, lat1, lon, lat)
        r2 = calc_distance(lon2, lat2, lon, lat)
        if (r1 < r) then
          q(i,j,1,2,old) = qmin + c * qmax * 0.5_r8 * (1 + cos(pi * r1 / r))
        else if (r2 < r) then
          q(i,j,1,2,old) = qmin + c * qmax * 0.5_r8 * (1 + cos(pi * r2 / r))
        else
          q(i,j,1,2,old) = qmin
        end if
      end do
    end do
    call fill_halo(block, q(:,:,:,2,old), full_lon=.true., full_lat=.true., full_lev=.true.)
    ! Slotted cylinders
    qmax = 1.0_r8; qmin = 0.1_r8; r = radius * 0.5_r8
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        r1 = calc_distance(lon1, lat1, lon, lat)
        r2 = calc_distance(lon2, lat2, lon, lat)
        if ((r1 <= r .and. abs(lon - lon1) >= r / radius / 6.0_r8) .or. &
            (r2 <= r .and. abs(lon - lon2) >= r / radius / 6.0_r8)) then
          q(i,j,1,3,old) = qmax
        else if (r1 <= r .and. abs(lon - lon1) < r / radius / 6.0_r8 .and. lat - lat1 < -5.0_r8 / 12.0_r8 * (r / radius)) then
          q(i,j,1,3,old) = qmax
        else if (r2 <= r .and. abs(lon - lon2) < r / radius / 6.0_r8 .and. lat - lat2 >  5.0_r8 / 12.0_r8 * (r / radius)) then
          q(i,j,1,3,old) = qmax
        else
          q(i,j,1,3,old) = qmin
        end if
      end do
    end do
    call fill_halo(block, q(:,:,:,3,old), full_lon=.true., full_lat=.true., full_lev=.true.)
    ! Gaussian hills
    qmax = 0.95_r8; c = 5.0_r8
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        call cartesian_transform(lon, lat, x(1), x(2), x(3))
        x = x / radius
        q(i,j,1,4,old) = qmax * (exp(-c * dot_product(x - x1, x - x1)) + exp(-c * dot_product(x - x2, x - x2)))
      end do
    end do
    call fill_halo(block, q(:,:,:,4,old), full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine deform_test_set_ic

  subroutine deform_case1_test_set_uv(block, state, time_in_seconds)

    type(block_type), intent(in   ) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: time_in_seconds

    integer i, j
    real(r8) lon, lat, k, cos_t

    associate (mesh => block%mesh   , &
               u    => state%u_lon  , &
               v    => state%v_lat  , &
               m    => state%m      , &
               mfx  => state%mfx_lon, &
               mfy  => state%mfy_lat)
    m = 1
    k = 10.0_r8 * radius / period
    cos_t = cos(pi * time_in_seconds / period)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      lat = mesh%full_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        lon = mesh%half_lon(i)
        u(i,j,1) = k * sin(lon / 2.0_r8)**2 * sin(2 * lat) * cos_t
        mfx(i,j,1) = u(i,j,1)
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      lat = mesh%half_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        v(i,j,1) = k / 2.0_r8 * sin(lon) * cos(lat) * cos_t
        mfy(i,j,1) = v(i,j,1)
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine deform_case1_test_set_uv

  subroutine deform_case2_test_set_uv(block, state, time_in_seconds)

    type(block_type), intent(in   ) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: time_in_seconds

    integer i, j
    real(r8) lon, lat, k, cos_t

    associate (mesh => block%mesh   , &
               m    => state%m      , &
               u    => state%u_lon  , &
               v    => state%v_lat  , &
               mfx  => state%mfx_lon, &
               mfy  => state%mfy_lat)
    m = 1
    k = 10.0_r8 * radius / period
    cos_t = cos(pi * time_in_seconds / period)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      lat = mesh%full_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        lon = mesh%half_lon(i)
        u(i,j,1) = k * sin(lon)**2 * sin(2 * lat) * cos_t
        mfx(i,j,1) = u(i,j,1)
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      lat = mesh%half_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        v(i,j,1) = k * sin(2 * lon) * cos(lat) * cos_t
        mfy(i,j,1) = v(i,j,1)
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine deform_case2_test_set_uv

  subroutine deform_case3_test_set_uv(block, state, time_in_seconds)

    type(block_type), intent(in   ) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: time_in_seconds

    integer i, j
    real(r8) lon, lat, k, cos_t

    associate (mesh => block%mesh   , &
               m    => state%m      , &
               u    => state%u_lon  , &
               v    => state%v_lat  , &
               mfx  => state%mfx_lon, &
               mfy  => state%mfy_lat)
    m = 1
    k = 5.0_r8 * radius / period
    cos_t = cos(pi * time_in_seconds / period)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      lat = mesh%full_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        lon = mesh%half_lon(i)
        u(i,j,1) = -k * sin(lon / 2.0_r8)**2 * sin(2 * lat) * cos(lat)**2 * cos_t
        mfx(i,j,1) = u(i,j,1)
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      lat = mesh%half_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        v(i,j,1) = k / 2.0_r8 * sin(lon) * cos(lat)**3 * cos_t
        mfy(i,j,1) = v(i,j,1)
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine deform_case3_test_set_uv

  subroutine deform_case4_test_set_uv(block, state, time_in_seconds)

    type(block_type), intent(in   ) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: time_in_seconds

    integer i, j
    real(r8) lon, lat, k, c1, c2, cos_t

    associate (mesh => block%mesh   , &
               m    => state%m      , &
               u    => state%u_lon  , &
               v    => state%v_lat  , &
               mfx  => state%mfx_lon, &
               mfy  => state%mfy_lat)
    m = 1
    k = 10.0_r8 * radius / period
    c1 = pi2 * time_in_seconds / period
    c2 = pi2 * radius / period
    cos_t = cos(pi * time_in_seconds / period)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      lat = mesh%full_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        lon = mesh%half_lon(i) - c1
        u(i,j,1) = k * sin(lon)**2 * sin(2 * lat) * cos_t + c2 * cos(lat)
        mfx(i,j,1) = u(i,j,1)
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      lat = mesh%half_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i) - c1
        v(i,j,1) = k * sin(2 * lon) * cos(lat) * cos_t
        mfy(i,j,1) = v(i,j,1)
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine deform_case4_test_set_uv

end module deform_test_mod
