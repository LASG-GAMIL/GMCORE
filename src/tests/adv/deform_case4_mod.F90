module deform_case4_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use sphere_geometry_mod
  use parallel_mod
  use block_mod
  use adv_mod

  implicit none

  private

  public deform_case4_set_ic

  real(r8), parameter :: period = 5
  real(r8), parameter :: lon1   = pi * 5.0_r8 / 6.0_r8
  real(r8), parameter :: lat1   = 0
  real(r8), parameter :: lon2   = pi * 7.0_r8 / 6.0_r8
  real(r8), parameter :: lat2   = 0

contains

  subroutine deform_case4_set_ic(block)

    type(block_type), intent(inout) :: block

    integer i, j
    real(r8) lon, lat, r, r1, r2, qmax, qmin, c
    type(tracer_type), pointer :: tracer

    call adv_add_tracer('deform_case4', 'q1', dt)
    call adv_add_tracer('deform_case4', 'q2', dt)

    call adv_allocate_tracers(block)

    associate (mesh => block%mesh)
    ! Cosine hills
    tracer => block%state(1)%adv_batches(1)%get_tracer('q1')
    qmax = 1; qmin = 0.1; c = 0.9; r = radius * 0.5
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        r1 = calc_distance(lon1, lat1, lon, lat)
        r2 = calc_distance(lon2, lat2, lon, lat)
        if (r1 < r) then
          tracer%q(i,j,1) = qmin + c * qmax * 0.5 * (1 + cos(pi * r1 / r))
        else if (r2 < r) then
          tracer%q(i,j,1) = qmin + c * qmax * 0.5 * (1 + cos(pi * r2 / r))
        else
          tracer%q(i,j,1) = qmin
        end if
      end do
    end do
    call fill_halo(block, tracer%q, full_lon=.true., full_lat=.true., full_lev=.true.)
    ! Slotted cylinders
    tracer => block%state(1)%adv_batches(1)%get_tracer('q2')
    qmax = 1; qmin = 0.1; r = 0.5
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      lat = mesh%full_lat(j)
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        lon = mesh%full_lon(i)
        r1 = calc_distance(lon1, lat1, lon, lat)
        r2 = calc_distance(lon2, lat2, lon, lat)
        if ((r1 <= r .and. abs(lon - lon1) >= r / 6.0_r8) .or. &
            (r2 <= r .and. abs(lon - lon2) >= r / 6.0_r8)) then
          tracer%q(i,j,1) = qmax
        else if (r1 <= r .and. abs(lon - lon1) < r / 6.0_r8 .and. lat - lat1 < -5.0_r8 / 12.0_r8 * r) then
          tracer%q(i,j,1) = qmax
        else if (r2 <= r .and. abs(lon - lon2) < r / 6.0_r8 .and. lat - lat2 >  5.0_r8 / 12.0_r8 * r) then
          tracer%q(i,j,1) = qmax
        else
          tracer%q(i,j,1) = qmin
        end if
      end do
    end do
    call fill_halo(block, tracer%q, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine deform_case4_set_ic

end module deform_case4_mod