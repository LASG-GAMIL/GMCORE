module jet_zonal_flow_test_mod

  use flogger
  use string
  use const_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public jet_zonal_flow_test_set_ic

  real(r8), parameter :: u_max = 80.0_r8
  real(r8), parameter :: lat0 = pi / 7.0_r8
  real(r8), parameter :: lat1 = pi / 2.0_r8 - lat0
  real(r8), parameter :: en = exp(-4.0_r8 / (lat1 - lat0)**2_r8)
  real(r8)            :: gh0
  real(r8)            :: ghd
  real(r8), parameter :: lat2 = pi / 4.0_r8
  real(r8), parameter :: alpha = 1.0_r8 / 3.0_r8
  real(r8), parameter :: beta = 1.0_r8 / 15.0_r8

contains

  subroutine jet_zonal_flow_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, neval, ierr
    real(r8) abserr

    associate (mesh   => block%mesh          , &
               u      => block%state(1)%u_lon, &
               v      => block%state(1)%v_lat, &
               gz     => block%state(1)%gz   , &
               gzs    => block%static%gzs)
    gh0 = g * 1.0d4
    ghd = g * 120.0d0

    gzs = 0

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        u(i,j,1) = u_function(mesh%full_lat(j))
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true.)

    v = 0

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      i = mesh%half_lon_ibeg
      if (j == mesh%full_lat_ibeg) then
        gz(i,j,1) = gh0
      else
        call qags(gh_integrand, -0.5*pi, mesh%full_lat(j), 1.0d-12, 1.0d-3, gz(i,j,1), abserr, neval, ierr)
        if (ierr /= 0) then
          call log_error('Failed to calculate integration at (' // to_str(i) // ',' // to_str(j) // ')!', __FILE__, __LINE__)
        end if
        gz(i,j,1) = gh0 - gz(i,j,1)
      end if
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        gz(i,j,1) = gz(mesh%half_lon_ibeg,j,1)
        ! Add perturbation.
        gz(i,j,1) = gz(i,j,1) + ghd * &
          cos(mesh%full_lat(j)) * &
          exp(-(merge(mesh%full_lon(i) - 2*pi, mesh%full_lon(i), mesh%full_lon(i) > pi)  / alpha)**2) * &
          exp(-((lat2 - mesh%full_lat(j)) / beta)**2)
      end do
    end do
    call fill_halo(block, gz, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine jet_zonal_flow_test_set_ic

  real(r8) function gh_integrand(lat) result(res)

    real(r8), intent(in) :: lat

    real(r8) u, f

    u = u_function(lat)
    f = 2 * omega * sin(lat)
    res = radius * u * (f + tan(lat) / radius * u)

  end function gh_integrand

  real(r8) function u_function(lat) result(res)

    real(r8), intent(in) :: lat

    if (lat <= lat0 .or. lat >= lat1) then
      res = 0.0
    else
      res = u_max / en * exp(1 / (lat - lat0) / (lat - lat1))
    end if

  end function u_function

end module jet_zonal_flow_test_mod
