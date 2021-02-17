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

  public nh_prepare
  public nh_solve
  public interp_gz

contains

  subroutine nh_prepare(blocks)

    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(proc%blocks)
      associate (mesh      => blocks(iblk)%mesh              , &
                 ph        => blocks(iblk)%state(1)%ph       , &
                 ph_lev    => blocks(iblk)%state(1)%ph_lev   , &
                 p         => blocks(iblk)%state(1)%p        , &
                 p_lev     => blocks(iblk)%state(1)%p_lev    , &
                 p_lev_lon => blocks(iblk)%state(1)%p_lev_lon, &
                 p_lev_lat => blocks(iblk)%state(1)%p_lev_lat)
        p = ph; p_lev = ph_lev
        call interp_lev_edge_to_lev_lon_edge(mesh, p_lev, p_lev_lon)
        call interp_lev_edge_to_lev_lat_edge(mesh, p_lev, p_lev_lat)
        call fill_halo(blocks(iblk), p_lev_lon, full_lon=.false., full_lat=.true., full_lev=.false., west_halo=.false., south_halo=.false., north_halo=.false.)
      end associate
      call diag_m_lev(blocks(iblk), blocks(iblk)%state(1))
    end do

  end subroutine nh_prepare

  subroutine nh_solve(block, tend, old_state, star_state, new_state, dt)

    type(block_type), intent(inout) :: block
    type(tend_type ), intent(inout) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    call interp_mf            (block, star_state)
    call interp_gz            (block, star_state)
    call interp_w             (block, star_state)
    call interp_wedphdlev     (block, star_state)
    call reduce_run           (block, star_state, dt, nh_pass_1)
    call diag_m_lev           (block, new_state)

    call calc_adv_gz          (block, star_state, tend)
    call calc_adv_w           (block, star_state, tend)
    call implicit_w_solver    (block, tend, old_state, star_state, new_state, dt)

    call diag_rhod            (block, new_state)
    call diag_linearized_p    (block, old_state, new_state)
    call interp_p             (block, new_state)
    call reduce_run           (block, new_state, dt, nh_pass_2)

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

    associate (mesh         => block%mesh        , &
               u            => state%u           , & ! in
               v            => state%v           , & ! in
               m_lev        => state%m_lev       , & ! in
               u_lev_lon    => state%u_lev_lon   , & ! out
               v_lev_lat    => state%v_lev_lat   , & ! out
               mf_lev_lon_n => state%mf_lev_lon_n, & ! out
               mf_lev_lat_n => state%mf_lev_lat_n)   ! out
      call interp_lon_edge_to_lev_lon_edge(mesh, u, u_lev_lon, handle_top_bottom=.true.)
      call interp_lat_edge_to_lev_lat_edge(mesh, v, v_lev_lat, handle_top_bottom=.true.)
      call interp_lev_edge_to_lev_lon_edge(mesh, m_lev, mf_lev_lon_n)
      call interp_lev_edge_to_lev_lat_edge(mesh, m_lev, mf_lev_lat_n)
      mf_lev_lon_n = mf_lev_lon_n * u_lev_lon
      mf_lev_lat_n = mf_lev_lat_n * v_lev_lat
      call fill_halo(block, u_lev_lon   , full_lon=.false., full_lat=.true. , full_lev=.false.)
      call fill_halo(block, v_lev_lat   , full_lon=.true. , full_lat=.false., full_lev=.false.)
      call fill_halo(block, mf_lev_lon_n, full_lon=.false., full_lat=.true. , full_lev=.false.)
      call fill_halo(block, mf_lev_lat_n, full_lon=.true. , full_lat=.false., full_lev=.false.)
    end associate

  end subroutine interp_mf

  subroutine interp_gz(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    associate (mesh         => block%mesh        , &
               mf_lev_lon_n => state%mf_lev_lon_n, & ! in
               mf_lev_lat_n => state%mf_lev_lat_n, & ! in
               gz_lev       => state%gz_lev      , & ! in
               gz           => state%gz          , & ! out
               gz_lev_lon   => state%gz_lev_lon  , & ! out
               gz_lev_lat   => state%gz_lev_lat)     ! out
      call interp_lev_edge_to_cell        (mesh, gz_lev, gz        )
      call interp_lev_edge_to_lev_lon_edge(mesh, gz_lev, gz_lev_lon, u=mf_lev_lon_n)
      call interp_lev_edge_to_lev_lat_edge(mesh, gz_lev, gz_lev_lat, v=mf_lev_lat_n)
      call fill_halo(block, gz_lev_lon, full_lon=.false., full_lat=.true. , full_lev=.false., south_halo=.false., north_halo=.false.)
      call fill_halo(block, gz_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false., west_halo=.false., east_halo=.false.)
    end associate

  end subroutine interp_gz

  subroutine interp_w(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    associate (mesh         => block%mesh        , &
               mf_lev_lon_n => state%mf_lev_lon_n, & ! in
               mf_lev_lat_n => state%mf_lev_lat_n, & ! in
               w_lev        => state%w_lev       , & ! in
               w            => state%w           , & ! out
               w_lev_lon    => state%w_lev_lon   , & ! out
               w_lev_lat    => state%w_lev_lat)      ! out
      call interp_lev_edge_to_cell        (mesh, w_lev, w        )
      call interp_lev_edge_to_lev_lon_edge(mesh, w_lev, w_lev_lon, u=mf_lev_lon_n)
      call interp_lev_edge_to_lev_lat_edge(mesh, w_lev, w_lev_lat, v=mf_lev_lat_n)
      call fill_halo(block, w_lev_lon, full_lon=.false., full_lat=.true. , full_lev=.false., south_halo=.false., north_halo=.false.)
#ifdef V_POLE
      call fill_halo(block, w_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false., west_halo=.false., east_halo=.false., south_halo=.false.)
#else
      call fill_halo(block, w_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false., west_halo=.false., east_halo=.false., north_halo=.false.)
#endif
    end associate

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
               mf_lev_lon_n  => state%mf_lev_lon_n , & ! in
               mf_lev_lat_n  => state%mf_lev_lat_n , & ! in
               gz_lev_lon    => state%gz_lev_lon   , & ! in
               gz_lev_lat    => state%gz_lev_lat   , & ! in
               gz_lev        => state%gz_lev       , & ! in
               gz            => state%gz           , & ! in
               m_lev         => state%m_lev        , & ! in
               wedphdlev     => state%wedphdlev    , & ! in
               wedphdlev_lev => state%wedphdlev_lev, & ! in
               adv_gz_lon    => tend%adv_gz_lon    , & ! out
               adv_gz_lat    => tend%adv_gz_lat    , & ! out
               adv_gz_lev    => tend%adv_gz_lev)       ! out
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
            adv_gz_lat(i,j,k) = pole(k) / m_lev(i,j,k)
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
            adv_gz_lat(i,j,k) = pole(k) / m_lev(i,j,k)
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
      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
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

      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
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
        do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            work(i,k) = mf_lev_lat_n(i,j,k) * (w_lev_lat(i,j,k) - w_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_w_lat(i,j,k) = pole(k) / m_lev(i,j,k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            work(i,k) = - mf_lev_lat_n(i,j-1,k) * (w_lev_lat(i,j-1,k) - w_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            adv_w_lat(i,j,k) = pole(k) / m_lev(i,j,k)
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
      ! Bottom w is from boundary condition, so no tendency for it.
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
              do kk = mesh%half_lev_iend, mesh%half_lev_ibeg, -1
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

    associate (mesh => block%mesh, & ! in
               rhod => state%rhod, & ! in
               pt   => state%pt  , & ! in
               p    => state%p)      ! out
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
    end associate

  end subroutine diag_p

  subroutine diag_linearized_p(block, old_state, new_state)

    type(block_type), intent(in) :: block
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    integer i, j, k

    associate (mesh       => block%mesh      , &
               old_m      => old_state%m     , & ! in
               new_m      => new_state%m     , & ! in
               old_p      => old_state%p     , & ! in
               new_p      => new_state%p     , & ! out
               old_gz_lev => old_state%gz_lev, & ! in
               new_gz_lev => new_state%gz_lev, & ! in
               old_pt     => old_state%pt    , & ! in
               new_pt     => new_state%pt)       ! in
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_p(i,j,k) = old_p(i,j,k) * (1.0_r8 + cp_o_cv * ( &
              (new_m(i,j,k) * new_pt(i,j,k)) /                  &
              (old_m(i,j,k) * old_pt(i,j,k)) -                  &
              (new_gz_lev(i,j,k+1) - new_gz_lev(i,j,k)) /       &
              (old_gz_lev(i,j,k+1) - old_gz_lev(i,j,k))         &
            ))
          end do
        end do
      end do
    end associate

  end subroutine diag_linearized_p

  subroutine interp_p(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    associate (mesh       => block%mesh      , &
               p          => state%p         , & ! in
               p_lev      => state%p_lev     , & ! out
               p_lev_lon  => state%p_lev_lon , & ! out
               p_lev_lat  => state%p_lev_lat)    ! out
      call interp_cell_to_lev_edge(mesh, p, p_lev, handle_top_bottom=.true.)
      call fill_halo(block, p_lev, full_lon=.true., full_lat=.true., full_lev=.false.)
      call interp_lev_edge_to_lev_lon_edge(mesh, p_lev, p_lev_lon)
      call interp_lev_edge_to_lev_lat_edge(mesh, p_lev, p_lev_lat)
      call fill_halo(block, p_lev_lon, full_lon=.false., full_lat=.true., full_lev=.false., west_halo=.false., south_halo=.false., north_halo=.false.)
    end associate

  end subroutine interp_p

  subroutine apply_bc_w_lev(block, star_state, new_state)

    type(block_type), intent(in) :: block
    type(state_type), intent(in) :: star_state
    type(state_type), intent(inout) :: new_state

    real(r8) us_dzsdlon, vs_dzsdlat
    integer i, j, k

    associate (mesh      => block%mesh          , &
               u_lev_lon => star_state%u_lev_lon, & ! in
               v_lev_lat => star_state%v_lev_lat, & ! in
               dzsdlon   => block%static%dzsdlon, & ! in
               dzsdlat   => block%static%dzsdlat, & ! in
               w_lev     => new_state%w_lev)        ! out
      k = mesh%half_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          us_dzsdlon = (u_lev_lon(i-1,j,k) * dzsdlon(i-1,j) + &
                        u_lev_lon(i  ,j,k) * dzsdlon(i  ,j)) * 0.5_r8
#ifdef V_POLE
          vs_dzsdlat = (v_lev_lat(i,j+1,k) * dzsdlat(i,j+1) + &
                        v_lev_lat(i,j  ,k) * dzsdlat(i,j  )) * 0.5_r8
#else
          vs_dzsdlat = (v_lev_lat(i,j-1,k) * dzsdlat(i,j-1) + &
                        v_lev_lat(i,j  ,k) * dzsdlat(i,j  )) * 0.5_r8
#endif
          w_lev(i,j,k) = us_dzsdlon + vs_dzsdlat
        end do
      end do
    end associate

  end subroutine apply_bc_w_lev

  subroutine implicit_w_solver(block, tend, old_state, star_state, new_state, dt)

    type(block_type), intent(in) :: block
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(in) :: star_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    real(r8) w1 (block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) gz1(block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) dgz(block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8) dp1, gdtbeta, gdt1mbeta, gdtbeta2gam
    real(r8) a(global_mesh%num_half_lev)
    real(r8) b(global_mesh%num_half_lev)
    real(r8) c(global_mesh%num_half_lev)
    real(r8) d(global_mesh%num_half_lev)
    integer i, j, k

    call apply_bc_w_lev(block, star_state, new_state)

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
               beta        => implicit_w_wgt   , &
               adv_gz_lon  => tend%adv_gz_lon  , & ! FIXME: After test success, merge advection tends togethor.
               adv_gz_lat  => tend%adv_gz_lat  , & !
               adv_gz_lev  => tend%adv_gz_lev  , & !
               adv_w_lon   => tend%adv_w_lon   , & !
               adv_w_lat   => tend%adv_w_lat   , & !
               adv_w_lev   => tend%adv_w_lev   , & !
               old_p       => old_state%p      , &
               star_p      => star_state%p     , &
               star_p_lev  => star_state%p_lev , &
               old_w_lev   => old_state%w_lev  , &
               star_w_lev  => star_state%w_lev , &
               new_w_lev   => new_state%w_lev  , &
               star_m_lev  => star_state%m_lev , &
               new_m_lev   => new_state%m_lev  , &
               old_gz_lev  => old_state%gz_lev , &
               star_gz_lev => star_state%gz_lev, &
               new_gz_lev  => new_state%gz_lev , &
               old_m       => old_state%m      , &
               new_m       => new_state%m      , &
               old_pt      => old_state%pt     , &
               new_pt      => new_state%pt)
      ! last: n, old: *, new: n + 1
      gdtbeta     = g * dt * beta
      gdt1mbeta   = g * dt * (1 - beta)
      gdtbeta2gam = (g * dt * beta)**2 * cp_o_cv
      ! FIXME: Two Poles may skip the duplicate calculation?
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            dgz(k) = old_gz_lev(i,j,k+1) - old_gz_lev(i,j,k)
          end do
          gz1 = old_gz_lev(i,j,:) - dt * (adv_gz_lon(i,j,:) + adv_gz_lat(i,j,:) + adv_gz_lev(i,j,:)) + gdt1mbeta * star_w_lev(i,j,:)
          w1  = old_w_lev (i,j,:) - dt * (adv_w_lon (i,j,:) + adv_w_lat (i,j,:) + adv_w_lev (i,j,:)) - g * dt
          do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
            w1(k) = w1(k) + gdt1mbeta * (star_p(i,j,k) - star_p(i,j,k-1)) / star_m_lev(i,j,k)
          end do
          ! Top boundary
          k = mesh%half_lev_ibeg
          w1(k) = w1(k) + gdt1mbeta * (star_p(i,j,k) - star_p_lev(i,j,k)) / star_m_lev(i,j,k)
          ! Bottom boundary
          k = mesh%half_lev_iend
          w1(k) = w1(k) + gdt1mbeta * (star_p_lev(i,j,k) - star_p(i,j,k-1)) / star_m_lev(i,j,k)
          ! Use linearized state of ideal gas to calculate the first part of ∂pⁿ⁺¹ (i.e. dp1).
          do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
            dp1 = (old_p(i,j,k) - old_p(i,j,k-1)) + cp_o_cv * ((                                     &
              old_p(i,j,k  ) * new_m(i,j,k  ) * new_pt(i,j,k  ) / old_m(i,j,k  ) / old_pt(i,j,k  ) - &
              old_p(i,j,k-1) * new_m(i,j,k-1) * new_pt(i,j,k-1) / old_m(i,j,k-1) / old_pt(i,j,k-1)   &
            ) - (                                                                                    &
              old_p(i,j,k  ) * (gz1(k+1) - gz1(k  )) / dgz(k  ) -                                    &
              old_p(i,j,k-1) * (gz1(k  ) - gz1(k-1)) / dgz(k-1)                                      &
            ))
            w1(k) = w1(k) + gdtbeta * dp1 / new_m_lev(i,j,k)
          end do
          ! Set coefficients for implicit solver.
          a(1) = 0.0_r8
          b(1) = 1.0_r8
          c(1) = 0.0_r8
          d(1) = 0.0_r8 ! Top w is set to zero.
          do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
            a(k) = gdtbeta2gam * old_p(i,j,k-1) / dgz(k-1)
            b(k) = new_m_lev(i,j,k) - gdtbeta2gam * (old_p(i,j,k) / dgz(k) + old_p(i,j,k-1) / dgz(k-1))
            c(k) = gdtbeta2gam * old_p(i,j,k  ) / dgz(k  )
            d(k) = new_m_lev(i,j,k) * w1(k)
          end do
          a(mesh%num_half_lev) = 0.0_r8
          b(mesh%num_half_lev) = 1.0_r8
          c(mesh%num_half_lev) = 0.0_r8
          d(mesh%num_half_lev) = new_w_lev(i,j,mesh%num_half_lev)
          call tridiag_thomas(a, b, c, d, new_w_lev(i,j,:))

          call rayleigh_damp_w(dt, star_gz_lev(i,j,:), new_w_lev(i,j,:))

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

    real(r8), parameter :: gzd = 10.0e3_r8 * g
    real(r8) c
    integer k

    do k = 2, size(w) - 1
      if (gz(k) > gz(1) - gzd) then
        c = rayleigh_damp_w_coef * sin(pi05 * (1 - (gz(1) - gz(k)) / gzd))**2
        w(k) = w(k) / (1 + c * dt)
      else
        return
      end if
    end do

  end subroutine rayleigh_damp_w

end module nh_mod
