module nh_mod

  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
  use process_mod
  use interp_mod
  use reduce_mod
  use math_mod
  use debug_mod
  use rayleigh_damp_mod

  implicit none

  private

  public nh_solve

contains

  subroutine nh_solve(block, tend, last_state, old_state, new_state, dt)

    type(block_type), intent(inout) :: block
    type(tend_type ), intent(inout) :: tend
    type(state_type), intent(in) :: last_state
    type(state_type), intent(inout) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    logical, save :: first_call = .true.

    call interp_mf            (block, old_state)
    call interp_gz            (block, old_state)
    call apply_bc_w_lev       (block, old_state)
    call interp_w             (block, old_state)
    call interp_wedphdlev     (block, old_state)
    if (first_call) then
      call diag_m_lev         (block, old_state)
      call diag_rhod          (block, old_state)
      call diag_p             (block, old_state)
      call reduce_run         (block, old_state, dt, nh_pass)
      first_call = .false.
    end if
    call diag_m_lev           (block, new_state)

    call calc_adv_gz          (block, old_state, tend)
    call calc_adv_w           (block, old_state, tend)
    call implicit_w_solver    (block, tend, last_state, old_state, new_state, dt)

    call diag_rhod            (block, new_state)
    call diag_p               (block, new_state)
    call reduce_run           (block, new_state, dt, nh_pass)

  end subroutine nh_solve

  subroutine diag_m_lev(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j, k

    associate (mesh => block%mesh)
      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%m_lev(i,j,k) = state%ph(i,j,k) - state%ph(i,j,k-1)
          end do
        end do
      end do
      ! Top boundary
      k = mesh%half_lev_ibeg
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%m_lev(i,j,k) = state%ph(i,j,k) - state%ph_lev(i,j,k)
        end do
      end do
      ! Bottom boundary
      k = mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%m_lev(i,j,k) = state%ph_lev(i,j,k) - state%ph(i,j,k-1)
        end do
      end do
    end associate

  end subroutine diag_m_lev

  subroutine interp_wedphdlev(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lev_edge_to_cell(block%mesh, state%wedphdlev_lev, state%wedphdlev)

  end subroutine interp_wedphdlev

  subroutine interp_mf(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lon_edge_to_lev_lon_edge(block%mesh, state%mf_lon_n, state%mf_lev_lon_n, handle_top_bottom=.true.)
    call interp_lat_edge_to_lev_lat_edge(block%mesh, state%mf_lat_n, state%mf_lev_lat_n, handle_top_bottom=.true.)
    call fill_halo(block, state%mf_lev_lon_n, full_lon=.false., full_lat=.true. , full_lev=.false.)
    call fill_halo(block, state%mf_lev_lat_n, full_lon=.true. , full_lat=.false., full_lev=.false.)

  end subroutine interp_mf

  subroutine interp_gz(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lev_edge_to_cell        (block%mesh, state%gz_lev, state%gz        )
    call interp_lev_edge_to_lev_lon_edge(block%mesh, state%gz_lev, state%gz_lev_lon, u=state%mf_lev_lon_n)
    call interp_lev_edge_to_lev_lat_edge(block%mesh, state%gz_lev, state%gz_lev_lat, v=state%mf_lev_lat_n)
    call fill_halo(block, state%gz_lev_lon, full_lon=.false., full_lat=.true. , full_lev=.false., south_halo=.false., north_halo=.false.)
    call fill_halo(block, state%gz_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false., west_halo=.false., east_halo=.false.)

  end subroutine interp_gz

  subroutine interp_w(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lev_edge_to_cell        (block%mesh, state%w_lev, state%w        )
    call interp_lev_edge_to_lev_lon_edge(block%mesh, state%w_lev, state%w_lev_lon, u=state%mf_lev_lon_n)
    call interp_lev_edge_to_lev_lat_edge(block%mesh, state%w_lev, state%w_lev_lat, v=state%mf_lev_lat_n)
    call fill_halo(block, state%w_lev_lon, full_lon=.false., full_lat=.true. , full_lev=.false., south_halo=.false., north_halo=.false.)
#ifdef V_POLE
    call fill_halo(block, state%w_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false., west_halo=.false., east_halo=.false., south_halo=.false.)
#else
    call fill_halo(block, state%w_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false., west_halo=.false., east_halo=.false., north_halo=.false.)
#endif

  end subroutine interp_w

  subroutine calc_adv_gz(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j, k, move
    real(r8) a, b, dwedphdlevgz, dwedphdlev
    real(r8) work(state%mesh%full_lon_ibeg:state%mesh%full_lon_iend,state%mesh%num_half_lev)
    real(r8) pole(state%mesh%num_half_lev)

    associate (mesh          => block%mesh         , &
               reduced_mesh  => block%reduced_mesh , &
               reduced_tend  => block%reduced_tend , &
               reduced_state => block%reduced_state, &
               mf_lev_lon_n  => state%mf_lev_lon_n , &
               mf_lev_lat_n  => state%mf_lev_lat_n , &
               gz_lev_lon    => state%gz_lev_lon   , &
               gz_lev_lat    => state%gz_lev_lat   , &
               gz_lev        => state%gz_lev       , &
               gz            => state%gz           , &
               m_lev         => state%m_lev        , &
               wedphdlev     => state%wedphdlev    , &
               wedphdlev_lev => state%wedphdlev_lev, &
               adv_gz_lon    => tend%adv_gz_lon    , &
               adv_gz_lat    => tend%adv_gz_lat    , &
               adv_gz_lev    => tend%adv_gz_lev)
      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (reduced_mesh(j)%reduce_factor > 1) then
            adv_gz_lon(:,j,k) = 0.0_r8
            do move = 1, reduced_mesh(j)%reduce_factor
              do i = reduced_mesh(j)%full_lon_ibeg, reduced_mesh(j)%full_lon_iend
                reduced_tend(j)%adv_gz_lon(i,k) = (               &
                  reduced_state(j)%mf_lev_lon_n(k,i  ,0,move) * ( &
                    reduced_state(j)%gz_lev_lon(k,i  ,0,move) -   &
                    reduced_state(j)%gz_lev    (k,i  ,0,move)     &
                  ) -                                             &
                  reduced_state(j)%mf_lev_lon_n(k,i-1,0,move) * ( &
                    reduced_state(j)%gz_lev_lon(k,i-1,0,move) -   &
                    reduced_state(j)%gz_lev    (k,i  ,0,move)     &
                  )                                               &
                ) / reduced_mesh(j)%de_lon(0) / reduced_state(j)%m_lev(k,i,0,move)
              end do
              call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%adv_gz_lon(:,k), mesh, adv_gz_lon(:,j,k))
            end do
            call overlay_inner_halo(block, adv_gz_lon(:,j,k), west_halo=.true.)
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              adv_gz_lon(i,j,k) = (                                             &
                mf_lev_lon_n(i  ,j,k) * (gz_lev_lon(i  ,j,k) - gz_lev(i,j,k)) - &
                mf_lev_lon_n(i-1,j,k) * (gz_lev_lon(i-1,j,k) - gz_lev(i,j,k))   &
              ) / mesh%de_lon(j) / m_lev(i,j,k)
            end do
          end if
        end do
      end do

      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            adv_gz_lat(i,j,k) = (                                             &
              mf_lev_lat_n(i,j+1,k) * (gz_lev_lat(i,j+1,k) - gz_lev(i,j,k)) - &
              mf_lev_lat_n(i,j  ,k) * (gz_lev_lat(i,j  ,k) - gz_lev(i,j,k))   &
            ) / mesh%de_lat(j) / m_lev(i,j,k)
#else
            adv_gz_lat(i,j,k) = (                                             &
              mf_lev_lat_n(i,j  ,k) * (gz_lev_lat(i,j  ,k) - gz_lev(i,j,k)) - &
              mf_lev_lat_n(i,j-1,k) * (gz_lev_lat(i,j-1,k) - gz_lev(i,j,k))   &
            ) / mesh%de_lat(j) / m_lev(i,j,k)
#endif
          end do
        end do
      end do
#ifndef V_POLE
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            work(i,k) = mf_lev_lat_n(i,j,k) * (gz_lev_lat(i,j,k) - gz_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_gz_lat(i,j,k) = pole(k) / m_lev(i,j,k) / global_mesh%num_full_lon
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            work(i,k) = - mf_lev_lat_n(i,j-1,k) * (gz_lev_lat(i,j-1,k) - gz_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_gz_lat(i,j,k) = pole(k) / m_lev(i,j,k) / global_mesh%num_full_lon
          end do
        end do
      end if
#endif

      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
        a = mesh%half_dlev_upper(k) / mesh%half_dlev_lower(k)
        b = mesh%half_dlev_lower(k) / mesh%half_dlev_upper(k)
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dwedphdlevgz = a * (wedphdlev(i,j,k  ) * gz(i,j,k  ) - wedphdlev_lev(i,j,k) * gz_lev(i,j,k)) - &
                           b * (wedphdlev(i,j,k-1) * gz(i,j,k-1) - wedphdlev_lev(i,j,k) * gz_lev(i,j,k))
            dwedphdlev = a * (wedphdlev(i,j,k  ) - wedphdlev_lev(i,j,k)) - &
                         b * (wedphdlev(i,j,k-1) - wedphdlev_lev(i,j,k))
            adv_gz_lev(i,j,k) = (dwedphdlevgz - gz_lev(i,j,k) * dwedphdlev) / m_lev(i,j,k)
          end do
        end do
      end do
      k = mesh%half_lev_ibeg
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          adv_gz_lev(i,j,k) = wedphdlev(i,j,k) * (gz(i,j,k) - gz_lev(i,j,k)) / m_lev(i,j,k)
        end do
      end do
      ! Bottom gz is static topography, so no tendency for it (i.e., gzs).
      k = mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          adv_gz_lev(i,j,k) = wedphdlev(i,j,k-1) * (gz_lev(i,j,k) - gz(i,j,k-1)) / m_lev(i,j,k)
        end do
      end do
#ifndef V_POLE
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_gz_lev(i,j,k) = adv_gz_lev(i,j,k) / global_mesh%num_full_lon
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_gz_lev(i,j,k) = adv_gz_lev(i,j,k) / global_mesh%num_full_lon
          end do
        end do
      end if
#endif
    end associate

  end subroutine calc_adv_gz

  subroutine calc_adv_w(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j, k, move
    real(r8) a, b, dwedphdlevw, dwedphdlev
    real(r8) work(state%mesh%full_lon_ibeg:state%mesh%full_lon_iend,state%mesh%num_half_lev)
    real(r8) pole(state%mesh%num_half_lev)

    associate (mesh          => block%mesh         , &
               reduced_mesh  => block%reduced_mesh , &
               reduced_tend  => block%reduced_tend , &
               reduced_state => block%reduced_state, &
               mf_lev_lon_n  => state%mf_lev_lon_n , &
               mf_lev_lat_n  => state%mf_lev_lat_n , &
               w_lev_lon     => state%w_lev_lon    , &
               w_lev_lat     => state%w_lev_lat    , &
               w_lev         => state%w_lev        , &
               w             => state%w            , &
               m_lev         => state%m_lev        , &
               wedphdlev     => state%wedphdlev    , &
               wedphdlev_lev => state%wedphdlev_lev, &
               adv_w_lon     => tend%adv_w_lon     , &
               adv_w_lat     => tend%adv_w_lat     , &
               adv_w_lev     => tend%adv_w_lev)
      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (reduced_mesh(j)%reduce_factor > 1) then
            adv_w_lon(:,j,k) = 0.0_r8
            do move = 1, reduced_mesh(j)%reduce_factor
              do i = reduced_mesh(j)%full_lon_ibeg, reduced_mesh(j)%full_lon_iend
                reduced_tend(j)%adv_w_lon(i,k) = (                &
                  reduced_state(j)%mf_lev_lon_n(k,i  ,0,move) * ( &
                    reduced_state(j)%w_lev_lon (k,i  ,0,move) -   &
                    reduced_state(j)%w_lev     (k,i  ,0,move)     &
                  ) -                                             &
                  reduced_state(j)%mf_lev_lon_n(k,i-1,0,move) * ( &
                    reduced_state(j)%w_lev_lon (k,i-1,0,move) -   &
                    reduced_state(j)%w_lev     (k,i  ,0,move)     &
                  )                                               &
                ) / reduced_mesh(j)%de_lon(0) / reduced_state(j)%m_lev(k,i,0,move)
              end do
              call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%adv_w_lon(:,k), mesh, adv_w_lon(:,j,k))
            end do
            call overlay_inner_halo(block, adv_w_lon(:,j,k), west_halo=.true.)
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              adv_w_lon(i,j,k) = (                                            &
                mf_lev_lon_n(i  ,j,k) * (w_lev_lon(i  ,j,k) - w_lev(i,j,k)) - &
                mf_lev_lon_n(i-1,j,k) * (w_lev_lon(i-1,j,k) - w_lev(i,j,k))   &
              ) / mesh%de_lon(j) / m_lev(i,j,k)
            end do
          end if
        end do
      end do

      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            adv_w_lat(i,j,k) = (                                            &
              mf_lev_lat_n(i,j+1,k) * (w_lev_lat(i,j+1,k) - w_lev(i,j,k)) - &
              mf_lev_lat_n(i,j  ,k) * (w_lev_lat(i,j  ,k) - w_lev(i,j,k))   &
            ) / mesh%de_lat(j) / m_lev(i,j,k)
#else
            adv_w_lat(i,j,k) = (                                            &
              mf_lev_lat_n(i,j  ,k) * (w_lev_lat(i,j  ,k) - w_lev(i,j,k)) - &
              mf_lev_lat_n(i,j-1,k) * (w_lev_lat(i,j-1,k) - w_lev(i,j,k))   &
            ) / mesh%de_lat(j) / m_lev(i,j,k)
#endif
          end do
        end do
      end do
#ifndef V_POLE
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            work(i,k) = mf_lev_lat_n(i,j,k) * (w_lev_lat(i,j,k) - w_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_w_lat(i,j,k) = pole(k) / m_lev(i,j,k) / global_mesh%num_full_lon
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            work(i,k) = - mf_lev_lat_n(i,j-1,k) * (w_lev_lat(i,j-1,k) - w_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_w_lat(i,j,k) = pole(k) / m_lev(i,j,k) / global_mesh%num_full_lon
          end do
        end do
      end if
#endif

      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
        a = mesh%half_dlev_upper(k) / mesh%half_dlev_lower(k)
        b = mesh%half_dlev_lower(k) / mesh%half_dlev_upper(k)
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dwedphdlevw = a * (wedphdlev(i,j,k  ) * w(i,j,k  ) - wedphdlev_lev(i,j,k) * w_lev(i,j,k)) - &
                          b * (wedphdlev(i,j,k-1) * w(i,j,k-1) - wedphdlev_lev(i,j,k) * w_lev(i,j,k))
            dwedphdlev = a * (wedphdlev(i,j,k  ) - wedphdlev_lev(i,j,k)) - &
                         b * (wedphdlev(i,j,k-1) - wedphdlev_lev(i,j,k))
            adv_w_lev(i,j,k) = (dwedphdlevw - state%w_lev(i,j,k) * dwedphdlev) / m_lev(i,j,k)
          end do
        end do
      end do
      ! Top w is fixed to be zero.
      k = mesh%half_lev_ibeg
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          adv_w_lev(i,j,k) = wedphdlev(i,j,k) * (w(i,j,k) - w_lev(i,j,k)) / m_lev(i,j,k)
        end do
      end do
      ! Bottom w is from boundary condition, so no tendency for it.
      k = mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          adv_w_lev(i,j,k) = wedphdlev(i,j,k-1) * (w_lev(i,j,k) - w(i,j,k-1)) / m_lev(i,j,k)
        end do
      end do
#ifndef V_POLE
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_w_lev(i,j,k) = adv_w_lev(i,j,k) / global_mesh%num_full_lon
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_w_lev(i,j,k) = adv_w_lev(i,j,k) / global_mesh%num_full_lon
          end do
        end do
      end if
#endif
    end associate

  end subroutine calc_adv_w

  subroutine diag_rhod(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    integer i, j, k, kk

    ! Diagnose dry air density from hydrostatic equation.
    associate (mesh     => block%mesh    , & ! in
               gz_lev   => state%gz_lev  , & ! in
               m        => state%m       , & ! in
               rhod     => state%rhod    , & ! out
               rhod_lon => state%rhod_lon, & ! out
               rhod_lat => state%rhod_lat)   ! out
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            rhod(i,j,k) = - m(i,j,k) / (gz_lev(i,j,k+1) - gz_lev(i,j,k))
            if (rhod(i,j,k) <= 0) then
              print *, i, j, k, rhod(i,j,k), m(i,j,k), gz_lev(i,j,k:k+1)
              do kk = 27, 1, -1
                print *, kk, gz_lev(i,j,kk)
              end do
              stop 'negative rhod!'
            end if
          end do
        end do
      end do
      call fill_halo(block, rhod, full_lon=.true. , full_lat=.true., full_lev=.true.)
      call interp_cell_to_lon_edge(mesh, rhod, rhod_lon)
      call interp_cell_to_lat_edge(mesh, rhod, rhod_lat)
      call fill_halo(block, rhod_lon, full_lon=.false., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false., north_halo=.false.)
    end associate

  end subroutine diag_rhod

  subroutine diag_p(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    real(r8), parameter :: p0 = 1.0e5_r8
    integer i, j, k

    associate (mesh => block%mesh          , & ! in
               rhod => state%rhod          , & ! in
               pt => state%pt              , & ! in
               p => state%p                , & ! out
               p_lev => state%p_lev        , & ! out
               p_lev_lon => state%p_lev_lon, & ! out
               p_lev_lat => state%p_lev_lat)   ! out
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            p(i,j,k) = p0 * (Rd * pt(i,j,k) * rhod(i,j,k) / p0)**cp_o_cv
            if (debug_is_inf(p(i,j,k))) then
              print *, i, j, k, pt(i,j,k), rhod(i,j,k)
              stop 'NaN p!'
            end if
          end do
        end do
      end do

      call interp_cell_to_lev_edge(mesh, p, p_lev, handle_top_bottom=.true.)
      call fill_halo(block, p_lev, full_lon=.true. , full_lat=.true., full_lev=.false.)
      call interp_lev_edge_to_lev_lon_edge(mesh, p_lev, p_lev_lon)
      call interp_lev_edge_to_lev_lat_edge(mesh, p_lev, p_lev_lat)
      call fill_halo(block, p_lev_lon, full_lon=.false., full_lat=.true., full_lev=.false., west_halo=.false., south_halo=.false., north_halo=.false.)
    end associate

  end subroutine diag_p

  subroutine diag_linearized_p(block, old_state, new_state)

    type(block_type), intent(in) :: block
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    integer i, j, k

    associate (mesh       => block%mesh      , &
               old_p      => old_state%p     , &
               new_p      => new_state%p     , &
               old_gz_lev => old_state%gz_lev, &
               new_gz_lev => new_state%gz_lev, &
               old_pt     => old_state%pt    , &
               new_pt     => new_state%pt)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_p(i,j,k) = old_p(i,j,k) * (1.0_r8 + cp_o_cv * ( &
              new_pt(i,j,k) / old_pt(i,j,k) -                   &
              (new_gz_lev(i,j,k+1) - new_gz_lev(i,j,k)) /       &
              (old_gz_lev(i,j,k+1) - old_gz_lev(i,j,k))         &
            ))
          end do
        end do
      end do

      call interp_cell_to_lev_edge(mesh, new_p, new_state%p_lev, handle_top_bottom=.true.)
      call fill_halo(block, new_state%p_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      call interp_lev_edge_to_lev_lon_edge(mesh, new_state%p_lev, new_state%p_lev_lon)
      call interp_lev_edge_to_lev_lat_edge(mesh, new_state%p_lev, new_state%p_lev_lat)
    end associate

  end subroutine diag_linearized_p

  subroutine apply_bc_w_lev(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    real(r8) x1, x2, a, b, us, vs
    integer i, j, k

    associate (mesh    => block%mesh          , &
               u       => state%u             , &
               v       => state%v             , &
               dzsdlon => block%static%dzsdlon, &
               dzsdlat => block%static%dzsdlat)
      k = mesh%half_lev_iend
      x1 = mesh%half_lev(k) - mesh%full_lev(k-1)
      x2 = mesh%half_lev(k) - mesh%full_lev(k-2)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          us = a * (u(i-1,j,k-1) + u(i,j,k-1)) + b * (u(i-1,j,k-2) + u(i,j,k-2))
#ifdef V_POLE
          vs = a * (v(i,j+1,k-1) + v(i,j,k-1)) + b * (v(i,j+1,k-2) + v(i,j,k-2))
#else
          vs = a * (v(i,j-1,k-1) + v(i,j,k-1)) + b * (v(i,j-1,k-2) + v(i,j,k-2))
#endif
          state%w_lev(i,j,k) = us * dzsdlon(i,j) + vs * dzsdlat(i,j)
        end do
      end do
    end associate

  end subroutine apply_bc_w_lev

  subroutine implicit_w_solver(block, tend, last_state, old_state, new_state, dt)

    type(block_type), intent(in) :: block
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: last_state
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    real(r8), parameter :: implicit_w_beta = 0.75_r8
    real(r8) w1 (block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) gz1(block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) dgz(block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8) dp1, gdtbeta, gdt1mbeta, gdtbeta2gam
    real(r8) a(block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) b(block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) c(block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) d(block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    integer i, j, k

    call apply_bc_w_lev(block, new_state)

    !
    ! ϕ¹ = ϕⁿ - Δt adv_ϕ* + g Δt (1 - β) w*
    !
    !                                                   
    ! w¹ = wⁿ - Δt adv_w* - g Δt + g Δt (1 - β) (∂p/∂π)*
    !                                                   
    ! Linearized state of ideal gas
    !
    ! 𝜹pⁿ⁺¹ ≈ 𝜹pⁿ + 𝜹(𝜸 pⁿ (𝜹𝜋 θ)ⁿ⁺¹ / (𝜹𝜋 θ)ⁿ) - 𝜹(𝜸 pⁿ 𝜹ϕ¹ / 𝜹ϕⁿ) - 𝜹(𝜸 pⁿ g Δt β 𝜹wⁿ⁺¹ / 𝜹φⁿ)
    !         -----------------------------------------------------
    !                                dp1
    !
    associate (mesh        => block%mesh       , &
               beta        => implicit_w_beta  , &
               adv_gz_lon  => tend%adv_gz_lon  , & ! FIXME: After test success, merge advection tends togethor.
               adv_gz_lat  => tend%adv_gz_lat  , & !
               adv_gz_lev  => tend%adv_gz_lev  , & !
               adv_w_lon   => tend%adv_w_lon   , & !
               adv_w_lat   => tend%adv_w_lat   , & !
               adv_w_lev   => tend%adv_w_lev   , & !
               last_p      => last_state%p     , &
               old_p       => old_state%p      , &
               old_p_lev   => old_state%p_lev  , &
               last_w_lev  => last_state%w_lev , &
               old_w_lev   => old_state%w_lev  , &
               new_w_lev   => new_state%w_lev  , &
               old_m_lev   => old_state%m_lev  , &
               new_m_lev   => new_state%m_lev  , &
               last_gz_lev => last_state%gz_lev, &
               old_gz_lev  => old_state%gz_lev , &
               new_gz_lev  => new_state%gz_lev , &
               last_m      => last_state%m     , &
               new_m       => new_state%m      , &
               last_pt     => last_state%pt    , &
               new_pt      => new_state%pt)
      ! last: n, old: *, new: n + 1
      gdtbeta     = g * dt * beta
      gdt1mbeta   = g * dt * (1 - beta)
      gdtbeta2gam = (g * dt * beta)**2 * cp_o_cv
      ! FIXME: Two Poles may skip the duplicate calculation?
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            dgz(k) = last_gz_lev(i,j,k+1) - last_gz_lev(i,j,k)
          end do
          gz1 = last_gz_lev(i,j,:) - dt * (adv_gz_lon(i,j,:) + adv_gz_lat(i,j,:) + adv_gz_lev(i,j,:)) + gdt1mbeta * old_w_lev(i,j,:)
          w1  = last_w_lev (i,j,:) - dt * (adv_w_lon (i,j,:) + adv_w_lat (i,j,:) + adv_w_lev (i,j,:)) - g * dt
          do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
            w1(k) = w1(k) + gdt1mbeta * (old_p(i,j,k) - old_p(i,j,k-1)) / old_m_lev(i,j,k)
          end do
          ! Top boundary
          k = mesh%half_lev_ibeg
          w1(k) = w1(k) + gdt1mbeta * (old_p(i,j,k) - old_p_lev(i,j,k)) / old_m_lev(i,j,k)
          ! Bottom boundary
          k = mesh%half_lev_iend
          w1(k) = w1(k) + gdt1mbeta * (old_p_lev(i,j,k) - old_p(i,j,k-1)) / old_m_lev(i,j,k)
          ! Use linearized state of ideal gas to calculate the first part of ∂pⁿ⁺¹ (i.e. dp1).
          do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
            dp1 = (last_p(i,j,k) - last_p(i,j,k-1)) + cp_o_cv * ((                                      &
              last_p(i,j,k  ) * new_m(i,j,k  ) * new_pt(i,j,k  ) / last_m(i,j,k  ) / last_pt(i,j,k  ) - &
              last_p(i,j,k-1) * new_m(i,j,k-1) * new_pt(i,j,k-1) / last_m(i,j,k-1) / last_pt(i,j,k-1)   &
            ) - (                                                                                       &
              last_p(i,j,k  ) * (gz1(k+1) - gz1(k  )) / dgz(k  ) -                                      &
              last_p(i,j,k-1) * (gz1(k  ) - gz1(k-1)) / dgz(k-1)                                        &
            ))
            w1(k) = w1(k) + gdtbeta * dp1 / new_m_lev(i,j,k)
          end do
          ! Set coefficients for implicit solver.
          a(1) = 0.0_r8
          b(1) = 1.0_r8
          c(1) = 0.0_r8
          d(1) = 0.0_r8 ! Top w is set to zero.
          do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
            a(k) = gdtbeta2gam * last_p(i,j,k-1) / dgz(k-1)
            b(k) = new_m_lev(i,j,k) - gdtbeta2gam * (last_p(i,j,k) / dgz(k) + last_p(i,j,k-1) / dgz(k-1))
            c(k) = gdtbeta2gam * last_p(i,j,k  ) / dgz(k  )
            d(k) = new_m_lev(i,j,k) * w1(k)
          end do
          a(mesh%num_half_lev) = 0.0_r8
          b(mesh%num_half_lev) = 1.0_r8
          c(mesh%num_half_lev) = 0.0_r8
          d(mesh%num_half_lev) = new_w_lev(i,j,mesh%num_half_lev)
          call tridiag_thomas(a, b, c, d, new_w_lev(i,j,:))

          call rayleigh_damp_w(dt, old_gz_lev(i,j,:), new_w_lev(i,j,:))

          ! Update gz after w is solved.
          do k = mesh%half_lev_ibeg, mesh%half_lev_iend - 1
            new_gz_lev(i,j,k) = gz1(k) + gdtbeta * new_w_lev(i,j,k)
          end do
        end do
      end do
      call fill_halo(block, new_w_lev , full_lon=.true., full_lat=.true., full_lev=.false.)
      call fill_halo(block, new_gz_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
    end associate

  end subroutine implicit_w_solver

  subroutine rayleigh_damp_w(dt, gz, w)

    real(8) , intent(in   ) :: dt
    real(r8), intent(in   ) :: gz(:)
    real(r8), intent(inout) :: w (:)

    real(r8), parameter :: rayleigh_damp_w_coef = 0.2_r8
    real(r8), parameter :: gzd = 10.0e3_r8
    real(r8) c
    integer k

    do k = 2, size(w) - 1
      if (gz(k) >= gz(1) - gzd) then
        c = rayleigh_damp_w_coef * sin(pi05 * (1 - (gz(1) - gz(k)) / gzd))**2
        w(k) = w(k) / (1 + c * dt)
      end if
    end do

  end subroutine rayleigh_damp_w

end module nh_mod
