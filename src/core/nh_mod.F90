module nh_mod

  use const_mod
  use block_mod
  use parallel_mod
  use process_mod
  use interp_mod
  use reduce_mod
  use tridiag_mod

  implicit none

  private

  public calc_m_lev
  public calc_mf_lev_lon_n_mf_lev_lat_n
  public calc_gz_lev_lon_gz_lev_lat
  public calc_w_w_lev_lon_w_lev_lat
  public calc_adv_gz
  public calc_adv_w
  public calc_p

contains

  subroutine calc_m_lev(block, state)

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
      k = mesh%half_lev_ibeg
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%m_lev(i,j,k) = state%ph(i,j,k) - state%ph_lev(i,j,k)
        end do
      end do
      k = mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%m_lev(i,j,k) = state%ph_lev(i,j,k) - state%ph(i,j,k-1)
        end do
      end do
    end associate

  end subroutine calc_m_lev

  subroutine calc_wedphdlev(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lev_edge_to_cell(block%mesh, state%wedphdlev_lev, state%wedphdlev)

  end subroutine calc_wedphdlev

  subroutine calc_mf_lev_lon_n_mf_lev_lat_n(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lon_edge_to_lev_lon_edge(block%mesh, state%mf_lon_n, state%mf_lev_lon_n, handle_top_bottom=.true.)
    call interp_lat_edge_to_lev_lat_edge(block%mesh, state%mf_lat_n, state%mf_lev_lat_n, handle_top_bottom=.true.)
    call fill_halo(block, state%mf_lev_lon_n, full_lon=.false., full_lat=.true. , full_lev=.false.)
    call fill_halo(block, state%mf_lev_lat_n, full_lon=.true. , full_lat=.false., full_lev=.false.)

  end subroutine calc_mf_lev_lon_n_mf_lev_lat_n

  subroutine calc_gz_lev_lon_gz_lev_lat(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lev_edge_to_cell        (block%mesh, state%gz_lev, state%gz        )
    call interp_lev_edge_to_lev_lon_edge(block%mesh, state%gz_lev, state%gz_lev_lon)
    call interp_lev_edge_to_lev_lat_edge(block%mesh, state%gz_lev, state%gz_lev_lat)
    call fill_halo(block, state%gz        , full_lon=.true. , full_lat=.true. , full_lev=.true. )
    call fill_halo(block, state%gz_lev_lon, full_lon=.false., full_lat=.true. , full_lev=.false.)
    call fill_halo(block, state%gz_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false.)

  end subroutine calc_gz_lev_lon_gz_lev_lat

  subroutine calc_w_w_lev_lon_w_lev_lat(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lev_edge_to_cell        (block%mesh, state%w_lev, state%w        )
    call interp_lev_edge_to_lev_lon_edge(block%mesh, state%w_lev, state%w_lev_lon)
    call interp_lev_edge_to_lev_lat_edge(block%mesh, state%w_lev, state%w_lev_lat)
    call fill_halo(block, state%w        , full_lon=.true. , full_lat=.true. , full_lev=.true.)
    call fill_halo(block, state%w_lev_lon, full_lon=.false., full_lat=.true. , full_lev=.false.)
    call fill_halo(block, state%w_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false.)

  end subroutine calc_w_w_lev_lon_w_lev_lat

  subroutine calc_adv_gz(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    integer i, j, k, move
    real(r8) a, b, dwedphdlevgz, dwedphdlev
    real(r8) pole(state%mesh%num_full_lev)

    associate (mesh          => block%mesh,         &
               reduced_mesh  => block%reduced_mesh, &
               reduced_tend  => block%reduced_tend, &
               reduced_state => block%reduced_state)
      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (reduced_mesh(j)%reduce_factor > 1) then
            tend%adv_gz(:,j,k) = 0.0_r8
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
              call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%adv_gz_lon(:,k), mesh, tend%adv_gz(:,j,k))
            end do
            call overlay_inner_halo(block, tend%adv_gz(:,j,k), west_halo=.true.)
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%adv_gz(i,j,k) = (                                                              &
                state%mf_lev_lon_n(i  ,j,k) * (state%gz_lev_lon(i  ,j,k) - state%gz_lev(i,j,k)) - &
                state%mf_lev_lon_n(i-1,j,k) * (state%gz_lev_lon(i-1,j,k) - state%gz_lev(i,j,k))   &
              ) / mesh%de_lon(j) / state%m_lev(i,j,k)
            end do
          end if
        end do
      end do

      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            tend%adv_gz(i,j,k) = tend%adv_gz(i,j,k) + (                                         &
              state%mf_lev_lat_n(i,j+1,k) * (state%gz_lev_lat(i,j+1,k) - state%gz_lev(i,j,k)) - &
              state%mf_lev_lat_n(i,j  ,k) * (state%gz_lev_lat(i,j  ,k) - state%gz_lev(i,j,k))   &
            ) / mesh%de_lat(j) / state%m_lev(i,j,k)
#else
            tend%adv_gz(i,j,k) = tend%adv_gz(i,j,k) + (                                         &
              state%mf_lev_lat_n(i,j  ,k) * (state%gz_lev_lat(i,j  ,k) - state%gz_lev(i,j,k)) - &
              state%mf_lev_lat_n(i,j-1,k) * (state%gz_lev_lat(i,j-1,k) - state%gz_lev(i,j,k))   &
            ) / mesh%de_lat(j) / state%m_lev(i,j,k)
#endif
          end do
        end do
      end do
#ifndef V_POLE
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pole(k) = pole(k) + state%mf_lev_lat_n(i,j,k) * (state%gz_lev_lat(i,j,k) - state%gz_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%adv_gz(i,j,k) = pole(k) / state%m_lev(i,j,k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pole(k) = pole(k) - state%mf_lev_lat_n(i,j-1,k) * (state%gz_lev_lat(i,j-1,k) - state%gz_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%adv_gz(i,j,k) = pole(k) / state%m_lev(i,j,k)
          end do
        end do
      end if
#endif

      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
        a = mesh%half_dlev_upper(k) / mesh%half_dlev_lower(k)
        b = mesh%half_dlev_lower(k) / mesh%half_dlev_upper(k)
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dwedphdlevgz = a * state%wedphdlev(i,j,k+1) * state%gz(i,j,k+1) - &
                           b * state%wedphdlev(i,j,k  ) * state%gz(i,j,k+1) - &
                           (a - b) * state%wedphdlev_lev(i,j,k) * state%gz_lev(i,j,k)
            dwedphdlev = a * state%wedphdlev(i,j,k+1) - &
                         b * state%wedphdlev(i,j,k  ) - &
                         (a - b) * state%wedphdlev_lev(i,j,k)
            tend%adv_gz(i,j,k) = tend%adv_gz(i,j,k) + (dwedphdlevgz - state%gz_lev(i,j,k) * dwedphdlev) / state%m_lev(i,j,k)
          end do
        end do
      end do
      k = mesh%half_lev_ibeg
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%adv_gz(i,j,k) = tend%adv_gz(i,j,k) + state%wedphdlev(i,j,k) * (state%gz(i,j,k) - state%gz_lev(i,j,k)) / state%m_lev(i,j,k)
        end do
      end do
      k = mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%adv_gz(i,j,k) = tend%adv_gz(i,j,k) + state%wedphdlev(i,j,k) * (state%gz_lev(i,j,k) - state%gz(i,j,k)) / state%m_lev(i,j,k)
        end do
      end do
    end associate

  end subroutine calc_adv_gz

  subroutine calc_adv_w(block, state, tend, dt)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real(8), intent(in) :: dt

    integer i, j, k, move
    real(r8) a, b, dwedphdlevw, dwedphdlev
    real(r8) pole(state%mesh%num_full_lev)

    associate (mesh          => block%mesh,         &
               reduced_mesh  => block%reduced_mesh, &
               reduced_tend  => block%reduced_tend, &
               reduced_state => block%reduced_state)
      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (reduced_mesh(j)%reduce_factor > 1) then
            tend%adv_w(:,j,k) = 0.0_r8
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
              call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%adv_w_lon(:,k), mesh, tend%adv_w(:,j,k))
            end do
            call overlay_inner_halo(block, tend%adv_w(:,j,k), west_halo=.true.)
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%adv_w(i,j,k) = tend%adv_w(i,j,k) + (                                         &
                state%mf_lev_lon_n(i  ,j,k) * (state%w_lev_lon(i  ,j,k) - state%w_lev(i,j,k)) - &
                state%mf_lev_lon_n(i-1,j,k) * (state%w_lev_lon(i-1,j,k) - state%w_lev(i,j,k))   &
              ) / mesh%de_lon(j) / state%m_lev(i,j,k)
            end do
          end if
        end do
      end do

      do k = mesh%half_lev_ibeg, mesh%half_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            tend%adv_w(i,j,k) = tend%adv_w(i,j,k) + (                                         &
              state%mf_lev_lat_n(i,j+1,k) * (state%w_lev_lat(i,j+1,k) - state%w_lev(i,j,k)) - &
              state%mf_lev_lat_n(i,j  ,k) * (state%w_lev_lat(i,j  ,k) - state%w_lev(i,j,k))   &
            ) / mesh%de_lat(j) / state%m_lev(i,j,k)
#else
            tend%adv_w(i,j,k) = tend%adv_w(i,j,k) + (                                         &
              state%mf_lev_lat_n(i,j  ,k) * (state%w_lev_lat(i,j  ,k) - state%w_lev(i,j,k)) - &
              state%mf_lev_lat_n(i,j-1,k) * (state%w_lev_lat(i,j-1,k) - state%w_lev(i,j,k))   &
            ) / mesh%de_lat(j) / state%m_lev(i,j,k)
#endif
          end do
        end do
      end do
#ifndef V_POLE
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pole(k) = pole(k) + state%mf_lev_lat_n(i,j,k) * (state%w_lev_lat(i,j,k) - state%w_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%adv_w(i,j,k) = pole(k) / state%m_lev(i,j,k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        pole = 0.0_r8
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pole(k) = pole(k) - state%mf_lev_lat_n(i,j-1,k) * (state%w_lev_lat(i,j-1,k) - state%w_lev(i,j,k))
          end do
        end do
        call zonal_sum(proc%zonal_comm, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tend%adv_w(i,j,k) = pole(k) / state%m_lev(i,j,k)
          end do
        end do
      end if
#endif

      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
        a = mesh%half_dlev_upper(k) / mesh%half_dlev_lower(k)
        b = mesh%half_dlev_lower(k) / mesh%half_dlev_upper(k)
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dwedphdlevw = a * state%wedphdlev(i,j,k+1) * state%w(i,j,k+1) - &
                          b * state%wedphdlev(i,j,k  ) * state%w(i,j,k+1) - &
                          (a - b) * state%wedphdlev_lev(i,j,k) * state%w_lev(i,j,k)
            dwedphdlev = a * state%wedphdlev(i,j,k+1) - &
                         b * state%wedphdlev(i,j,k  ) - &
                         (a - b) * state%wedphdlev_lev(i,j,k)
            tend%adv_w(i,j,k) = tend%adv_w(i,j,k) + (dwedphdlevw - state%w_lev(i,j,k) * dwedphdlev) / state%m_lev(i,j,k)
          end do
        end do
      end do
      k = mesh%half_lev_ibeg
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%adv_w(i,j,k) = tend%adv_w(i,j,k) + state%wedphdlev(i,j,k) * (state%w(i,j,k) - state%w_lev(i,j,k)) / state%m_lev(i,j,k)
        end do
      end do
      k = mesh%half_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%adv_w(i,j,k) = tend%adv_w(i,j,k) + state%wedphdlev(i,j,k) * (state%w_lev(i,j,k) - state%w(i,j,k)) / state%m_lev(i,j,k)
        end do
      end do
    end associate

  end subroutine calc_adv_w

  subroutine calc_p(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

  end subroutine calc_p

  subroutine implicit_w_solver(block, tend, old_state, new_state, dt)

    type(block_type), intent(in) :: block
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    real(r8), parameter :: implicit_w_beta = 0.75_r8
    real(r8) w1 (block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) gz1(block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) dp1(block%mesh%half_lev_lb:block%mesh%half_lev_ub)
    real(r8) p1, p2

    integer i, j, k

    !
    ! ϕ¹ = ϕⁿ - Δt adv_ϕⁿ + g Δt (1 - β) wⁿ⁺¹
    !
    !                                                   
    ! w¹ = wⁿ - Δt adv_wⁿ - g Δt + g Δt (1 - β) (∂p/∂π)ⁿ⁺¹
    !                                                   
    ! Linearized state of ideal gas
    !
    ! 𝜹pⁿ⁺¹ ≈ 𝜹pⁿ + 𝜹(𝜸 pⁿ (𝜹𝜋 θ)ⁿ⁺¹ / (𝜹𝜋 θ)ⁿ) - 𝜹(𝜸 pⁿ 𝜹ϕ¹ / 𝜹ϕⁿ) - 𝜹(𝜸 pⁿ g Δt β 𝜹wⁿ⁺¹ / 𝜹φⁿ)
    !         -----------------------------------------------------
    !                                dp1
    !
    associate (mesh => block%mesh, beta => implicit_w_beta)
      ! FIXME: Two Poles may skip the duplicate calculation?
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          gz1 = old_state%gz_lev(i,j,:) - dt * tend%adv_gz(i,j,:) + g * dt * (1 - beta) * old_state%w_lev(i,j,:)
          w1  = old_state%w_lev (i,j,:) - dt * tend%adv_w (i,j,:) - g * dt
          do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
            w1(k) = w1(k) + g * dt * (1 - beta) * (old_state%p(i,j,k) - old_state%p(i,j,k-1)) / old_state%m_lev(i,j,k)
          end do
          ! Top boundary
          k = mesh%half_lev_ibeg
          w1(k) = w1(k) + g * dt * (1 - beta) * (old_state%p(i,j,k) - old_state%p_lev(i,j,k)) / old_state%m_lev(i,j,k)
          ! Bottom boundary
          k = mesh%half_lev_iend
          w1(k) = w1(k) + g * dt * (1 - beta) * (old_state%p_lev(i,j,k) - old_state%p(i,j,k-1)) / old_state%m_lev(i,j,k)
          ! Use linearized state of ideal gas to calculate the first part of ∂pⁿ⁺¹ (i.e. dp1).
          do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
            p1 = old_state%p_lev(i,j,k-1)
            p2 = old_state%p_lev(i,j,k  )
            dp1(k) = (old_state%p(i,j,k) - old_state%p(i,j,k-1)) + cp_o_cv * ((                                  &
              p2 * new_state%m(i,j,k  ) * new_state%pt(i,j,k  ) / old_state%m(i,j,k  ) / old_state%pt(i,j,k  ) - &
              p1 * new_state%m(i,j,k-1) * new_state%pt(i,j,k-1) / old_state%m(i,j,k-1) / old_state%pt(i,j,k-1)   &
            ) - (                                                                                                &
              p2 * (gz1(k+1) - gz1(k  )) / (old_state%gz_lev(i,j,k+1) - old_state%gz_lev(i,j,k  )) -             &
              p1 * (gz1(k  ) - gz1(k-1)) / (old_state%gz_lev(i,j,k  ) - old_state%gz_lev(i,j,k-1))               &
            ))
          end do
        end do
      end do
    end associate

  end subroutine implicit_w_solver

end module nh_mod
