module deform_case4_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use parallel_mod
  use history_mod
  use block_mod
  use adv_mod

  implicit none

  private

  public deform_case4_set_ic
  public deform_case4_set_uv

  real(r8), parameter :: period = 12 * 86400
  real(r8), parameter :: lon1   = pi * 5.0_r8 / 6.0_r8
  real(r8), parameter :: lat1   = 0
  real(r8), parameter :: lon2   = pi * 7.0_r8 / 6.0_r8
  real(r8), parameter :: lat2   = 0

contains

  subroutine deform_case4_set_ic(block)

    type(block_type), intent(inout) :: block

    integer i, j
    real(r8) lon, lat, r, r1, r2, qmax, qmin, c
    real(r8) x(3), x1(3), x2(3)

    call adv_add_tracer('deform_case4', dt, 'q1', 'cosine hills tracer'     )
    call adv_add_tracer('deform_case4', dt, 'q2', 'slotted cylinders tracer')
    call adv_add_tracer('deform_case4', dt, 'q3', 'gaussian hills tracer'   )

    call adv_allocate_tracers(block)

    call cartesian_transform(lon1, lat1, x1(1), x1(2), x1(3)); x1 = x1 / radius
    call cartesian_transform(lon2, lat2, x2(1), x2(2), x2(3)); x2 = x2 / radius

    associate (mesh => block%mesh, state => block%state(1))
    ! Cosine hills
    qmax = 1.0_r8; qmin = 0.1_r8; c = 0.9_r8; r = radius * 0.5_r8
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        r1 = calc_distance(lon1, lat1, lon, lat)
        r2 = calc_distance(lon2, lat2, lon, lat)
        if (r1 < r) then
          state%q(i,j,1,1) = qmin + c * qmax * 0.5_r8 * (1 + cos(pi * r1 / r))
        else if (r2 < r) then
          state%q(i,j,1,1) = qmin + c * qmax * 0.5_r8 * (1 + cos(pi * r2 / r))
        else
          state%q(i,j,1,1) = qmin
        end if
      end do
    end do
    call fill_halo(block, state%q(:,:,:,1), full_lon=.true., full_lat=.true., full_lev=.true.)
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
          state%q(i,j,1,2) = qmax
        else if (r1 <= r .and. abs(lon - lon1) < r / radius / 6.0_r8 .and. lat - lat1 < -5.0_r8 / 12.0_r8 * (r / radius)) then
          state%q(i,j,1,2) = qmax
        else if (r2 <= r .and. abs(lon - lon2) < r / radius / 6.0_r8 .and. lat - lat2 >  5.0_r8 / 12.0_r8 * (r / radius)) then
          state%q(i,j,1,2) = qmax
        else
          state%q(i,j,1,2) = qmin
        end if
      end do
    end do
    call fill_halo(block, state%q(:,:,:,2), full_lon=.true., full_lat=.true., full_lev=.true.)
    ! Gaussian hills
    qmax = 0.95_r8; c = 5.0_r8
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        call cartesian_transform(lon, lat, x(1), x(2), x(3))
        x = x / radius
        state%q(i,j,1,3) = qmax * (exp(-c * dot_product(x - x1, x - x1)) + exp(-c * dot_product(x - x2, x - x2)))
      end do
    end do
    call fill_halo(block, state%q(:,:,:,3), full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine deform_case4_set_ic

  subroutine deform_case4_set_uv(block, state, time_in_seconds)

    type(block_type), intent(in   ) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: time_in_seconds

    integer i, j
    real(r8) lon, lat, c0, c1, c2, cos_t

    associate (mesh => block%mesh, u => state%u, v => state%v)
    c0 = 10.0_r8 * radius / period
    c1 = pi2 * time_in_seconds / period
    c2 = pi2 * radius / period
    cos_t = cos(pi * time_in_seconds / period)
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      lat = mesh%full_lat(j)
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        lon = mesh%half_lon(i) - c1
        u(i,j,1) = c0 * sin(lon)**2  * sin(lat**2) * cos_t + c2 * cos(lat)
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      lat = mesh%half_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i) - c1
        v(i,j,1) = c0 * sin(lon * 2) * cos(lat) * cos_t
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine deform_case4_set_uv

end module deform_case4_mod
