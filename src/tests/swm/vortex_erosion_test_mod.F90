module vortex_erosion_test_mod
  
  ! J. R. Bates & Yong Li (1997) Simulation of Stratospheric Vortex Erosion 
  ! Using Three Different Global Shallow Water Numerical Models.

  use flogger
  use string
  use const_mod
  use time_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public vortex_erosion_test_set_ic
  public vortex_erosion_test_apply_forcing

contains
  
  subroutine vortex_erosion_test_set_ic(block)

    type(block_type), intent(inout), target :: block

    integer i, j, neval, ierr
    real(r8) abserr
    real(r8) gh0
    
    associate(mesh => block%mesh          , &
              u    => block%state(1)%u_lon, &
              v    => block%state(1)%v_lat, &
              gz   => block%state(1)%gz   , &
              gzs  => block%static%gzs)
    gh0 = g * 6.0e3_r8

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      u(:,j,1) = u_function(mesh%full_lat(j))
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true.)
    
    v = 0.0_r8 

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      i = mesh%half_lon_ibeg
      if (j == mesh%full_lat_ibeg) then
        gz(i,j,1) = gh0
      else
        call qags(gh_integrand, -0.5*pi, mesh%full_lat(j), 1.0e-12, 1.0e-3, gz(i,j,1), abserr, neval, ierr)
        if (ierr /= 0) then
          call log_error('Failed to calculate integration at (' // to_str(i) // ',' // to_str(j) // '!')
        end if 
        gz(i,j,1) = gh0 - gz(i,j,1)
      end if
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        gz(i,j,1) = gz(mesh%half_lon_ibeg,j,1)
      end do 
    end do 

    call fill_halo(block, gz, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine vortex_erosion_test_set_ic

  subroutine vortex_erosion_test_apply_forcing(block, static)

    type(block_type), intent(in) :: block
    type(static_type), intent(inout) :: static
    integer i, j, k
    real(r8) hs, elapsed_days, at, b_lat, y
    
    hs = 720.0_r8 ! m 
    elapsed_days = elapsed_seconds / 86400._r8
    elapsed_days = mod(elapsed_days, 20.0_r8)

    associate (mesh => block%mesh, &
               gzs  => static%gzs)
    if (elapsed_days < 4.0_r8) then
      at = 0.5_r8 * (1.0_r8 - cos(pi * elapsed_days / 4.0_r8))
    else if (elapsed_days < 16.0_r8) then
      at = 1.0_r8
    else if (elapsed_days < 20.0_r8) then
      at = 0.5_r8 * (1.0_r8 + cos(pi * (elapsed_days - 16._r8) / 4.0_r8))
    end if

    b_lat = 0.0_r8
    gzs = 0.0_r8
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      if (mesh%full_lat(j) > 0.0_r8) then
        y = (tan(pi * 0.25_r8) / tan(mesh%full_lat(j)))**2
        b_lat = y * exp(1.0_r8 - y)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gzs(i,j) = hs * at * b_lat * mesh%full_sin_lon(i) * g
        end do
      end if
    end do
    call fill_halo(block, gzs, full_lon=.true., full_lat=.true.)
    end associate

  end subroutine vortex_erosion_test_apply_forcing

  real(r8) function gh_integrand(lat) result(res)

    real(r8), intent(in) :: lat

    real(r8) u, f

    u = u_function(lat)
    f = 2 * omega * sin(lat)
    res = radius * u * (f + tan(lat) / radius * u)

  end function gh_integrand

  real(r8) function u_function(lat) result(res)

    real(r8), intent(in) :: lat
    real(r8) lat_deg

    lat_deg = lat * deg

    if (lat_deg <= 0.0_r8) then
      res = -lat_deg / 9.0_r8 - 10.0_r8
    else if (lat_deg < 60.0_r8) then
      res = lat_deg - 10.0_r8
    else
      res = -5.0_r8 / 3.0_r8 * lat_deg + 150
    end if 

  end function u_function

end module vortex_erosion_test_mod
