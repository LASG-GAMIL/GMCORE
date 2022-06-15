module dcmip12_test_mod

  use const_mod
  use namelist_mod, dt => dt_adv
  use parallel_mod
  use history_mod
  use block_mod
  use vert_coord_mod
  use operators_mod
  use interp_mod
  use adv_mod

  implicit none

  private

  public dcmip12_test_init
  public dcmip12_test_set_ic
  public dcmip12_test_set_uv

  real(r8), parameter :: T0   = 300
  real(r8), parameter :: p0   = 1.0e5_r8
  real(r8), parameter :: K0   = 5
  real(r8), parameter :: u0   = 40
  real(r8), parameter :: w0   = 0.15
  real(r8), parameter :: z1   = 2000
  real(r8), parameter :: z2   = 5000
  real(r8), parameter :: z0   = 0.5_r8 * (z1 + z2)
  real(r8), parameter :: ztop = 12000
  real(r8), parameter :: tau  = 86400

  real(r8) rho0

contains

  subroutine dcmip12_test_init()

    ptop = 25494.4
    rho0 = p0 / (Rd * T0)

  end subroutine dcmip12_test_init

  subroutine dcmip12_test_set_ic(block)

    type(block_type), intent(inout) :: block

    real(r8) z
    integer i, j, k

    call adv_add_tracer('dcmip12', dt, 'q0', 'background tracer')
    call adv_add_tracer('dcmip12', dt, 'q1', 'test tracer'      )

    call adv_allocate_tracers(block)

    associate (mesh => block%mesh              , &
               old  => block%adv_batches(1)%old, &
               gz   => block%state(1)%gz       , &
               q    => block%adv_batches(1)%q  )
    ! Background
    q(:,:,:,1,old) = 1
    ! Test
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          z = gz(i,j,k) / g
          if (z1 < z .and. z < z2) then
            q(i,j,k,2,old) = 0.5_r8 * (1 + cos(pi2 * (z - z0) / (z2 - z1)))
          else
            q(i,j,k,2,old) = 0
          end if
        end do
      end do
    end do
    call fill_halo(block, q(:,:,:,2,old), full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine dcmip12_test_set_ic

  subroutine dcmip12_test_set_uv(block, state, time_in_seconds)

    type(block_type), intent(in   ) :: block
    type(state_type), intent(inout) :: state
    real(8), intent(in) :: time_in_seconds

    integer i, j, k
    real(r8) lon, lat, rho, cos_t, dphdlev

    cos_t = cos(pi * time_in_seconds / tau)

    associate (mesh     => block%mesh    , &
               phs      => state%phs     , &
               ph_lev   => state%ph_lev  , &
               ph       => state%ph      , &
               m_lon    => state%m_lon   , &
               m_lat    => state%m_lat   , &
               m_lev    => state%m_lev   , &
               t        => state%t       , &
               gz_lev   => state%gz_lev  , &
               gz       => state%gz      , &
               u        => state%u_lon   , &
               v        => state%v_lat   , &
               mf_lon_n => state%mf_lon_n, &
               mf_lat_n => state%mf_lat_n, &
               we       => state%we_lev)
    phs = p0
    t   = T0
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ph_lev(i,j,k) = vert_coord_calc_ph_lev(k, phs(i,j))
        end do
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ph(i,j,k) = 0.5d0 * (ph_lev(i,j,k) + ph_lev(i,j,k+1))
        end do
      end do
    end do
    call calc_m(block, state)
    call calc_gz_lev(block, state)
    call interp_lev_edge_to_cell(mesh, gz_lev, gz)
    m_lon = 1; m_lat = 1; m_lev = 1
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        lat = mesh%full_lat(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          lon = mesh%half_lon(i)
          u(i,j,k) = u0 * cos(lat)
          mf_lon_n(i,j,k) = u(i,j,k) * m_lon(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        lat = mesh%half_lat(j)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          lon = mesh%full_lon(i)
          rho = ph(i,j,k) / (Rd * T0)
          v(i,j,k) = -radius * w0 * pi * rho0 / (K0 * ztop * rho) * &
            cos(lat) * sin(K0 * lat) * cos(pi * gz(i,j,k) / (g * ztop)) * cos_t
          mf_lat_n(i,j,k) = v(i,j,k) * m_lat(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        lat = mesh%full_lat(j)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          lon = mesh%full_lon(i)
          rho = ph_lev(i,j,k) / (Rd * T0)
          dphdlev = m_lev(i,j,k) / mesh%half_dlev(k)
          we(i,j,k) = -dphdlev * g * w0 * rho0 / K0 * (                   &
            -2 * sin(K0 * lat) * sin(lat) + K0 * cos(lat) * cos(K0 * lat) &
          ) * sin(pi * gz_lev(i,j,k) / (g * ztop)) * cos_t / p0
        end do
      end do
    end do
    end associate

  end subroutine dcmip12_test_set_uv

end module dcmip12_test_mod
