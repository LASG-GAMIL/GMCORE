module gmcore_mod

  use flogger
  use const_mod
  use namelist_mod
  use parallel_mod
  use time_mod, dt => dt_in_seconds, old => old_time_idx, new => new_time_idx
  use history_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use operators_mod
  use reduce_mod
  use debug_mod
  use damp_mod

  implicit none

  private

  public gmcore_init
  public gmcore_run
  public gmcore_final

  integer, parameter :: csp2      = 1
  integer, parameter :: all_pass  = 0
  integer, parameter :: slow_pass = 1
  integer, parameter :: fast_pass = 2

  interface
    subroutine integrator_interface(dt, static, tends, states, old, new, pass)
      import r8, static_type, tend_type, state_type
      real(r8)         , intent(in   ) :: dt
      type(static_type), intent(in   ) :: static
      type(tend_type  ), intent(inout) :: tends (:)
      type(state_type ), intent(inout) :: states(:)
      integer          , intent(in   ) :: old
      integer          , intent(in   ) :: new
      integer          , intent(in   ) :: pass
    end subroutine integrator_interface

    subroutine splitter_interface(dt, static, tends, states)
      import r8, static_type, tend_type, state_type
      real(r8)         , intent(in   ) :: dt
      type(static_type), intent(in   ) :: static
      type(tend_type  ), intent(inout) :: tends (:)
      type(state_type ), intent(inout) :: states(:)
    end subroutine splitter_interface
  end interface

  procedure(integrator_interface), pointer :: integrator
  procedure(splitter_interface), pointer :: splitter

  real(r8) :: damp_t0 = -1

contains

  subroutine gmcore_init()

    call log_init()
    call parallel_init()
    call time_init()
    call mesh_init_root()
    call state_init_root()
    call static_init_root()
    call tend_init_root()
    call history_init()
    call reduce_init()

    select case (time_scheme)
    case ('pc2')
      integrator => predict_correct
    case ('rk3')
      integrator => runge_kutta_3rd
    case ('rk4')
      integrator => runge_kutta_4th
    case default
      integrator => predict_correct
      call log_notice('Use pc2 integrator.')
    end select

    select case (split_scheme)
    case ('csp2')
      splitter => csp2_splitting
    case default
      splitter => no_splitting
      call log_notice('No fast-slow split.')
    end select

    call time_add_alert('print', hours=1.0_r8)

  end subroutine gmcore_init

  subroutine gmcore_run()

    call calc_m_lon_m_lat(states(old))
    call calc_m_vtx(states(old))
    call operators_prepare(states(old))
    call diagnose(states(old))
    call output(states(old), tends(old))
    call log_print_diag(curr_time%isoformat())

    do while (.not. time_is_finished())
      call time_integrate(dt, static, tends, states)
      if (time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      call time_advance()
      call diagnose(states(old))
      if (states(old)%total_pe > states(new)%total_pe) damp_t0 = elapsed_seconds
      call output(states(old), tends(old))
    end do

  end subroutine gmcore_run

  subroutine gmcore_final()

    call parallel_final()
    call reduce_final()

  end subroutine gmcore_final

  subroutine output(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type ), intent(in) :: tend

    if (time_is_alerted('history_write')) then
      call history_write_state(static, state)
      call history_write_debug(static, state, tend)
    end if

  end subroutine output

  subroutine diagnose(state)

    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j
    ! real(r8) cfl

    mesh => state%mesh

    state%total_m = 0.0_r8
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%total_m = state%total_m + state%gd(i,j) * mesh%cell_area(j)
      end do
    end do

    state%total_ke = 0.0_r8
    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%total_ke = state%total_ke + state%mf_lon_n(i,j) * state%u(i,j) * mesh%lon_edge_area(j) * 2
      end do
    end do
    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%total_ke = state%total_ke + state%mf_lat_n(i,j) * state%v(i,j) * mesh%lat_edge_area(j) * 2
      end do
    end do
    state%total_e = state%total_ke
    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        state%total_e = state%total_e + (state%gd(i,j)**2 / g + state%gd(i,j) * static%ghs(i,j) / g) * mesh%cell_area(j)
      end do
    end do

    state%total_av = 0.0_r8
    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%total_av = state%total_av + state%m_vtx(i,j) * state%pv(i,j) * mesh%vertex_area(j)
      end do
    end do

    state%total_pe = 0.0_r8
    do j = mesh%half_lat_start_idx, mesh%half_lat_end_idx
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        state%total_pe = state%total_pe + state%m_vtx(i,j) * state%pv(i,j)**2 * 0.5_r8 * mesh%vertex_area(j)
      end do
    end do

    ! do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
    !   cfl = 0.0_r8
    !   if (reduced_mesh(j)%reduce_factor > 0) then
    !     do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
    !       cfl = max(dt / fast_cycles * sqrt(state%gd(i,j) + static%ghs(i,j)) / reduced_mesh(j)%de_lon(0), cfl)
    !     end do
    !   else
    !     do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
    !       cfl = max(dt / fast_cycles * sqrt(state%gd(i,j) + static%ghs(i,j)) / mesh%de_lon(j), cfl)
    !     end do
    !   end if
    !   if (cfl > 0.5) then
    !     print *, j, cfl
    !   end if
    ! end do

    call log_add_diag('total_m' , state%total_m )
    call log_add_diag('total_e' , state%total_e )
    call log_add_diag('total_pe', state%total_pe)

  end subroutine diagnose

  subroutine space_operators(static, state, tend, dt, pass)

    type(static_type), intent(in   ) :: static
    type(state_type) , intent(inout) :: state
    type(tend_type)  , intent(inout) :: tend
    real(r8)         , intent(in   ) :: dt
    integer          , intent(in   ) :: pass

    type(mesh_type), pointer :: mesh
    integer i, j

    call operators_prepare(state)
    call reduce_run(state, dt)

    mesh => state%mesh

    select case (pass)
    case (all_pass)
      call calc_qhu_qhv(state, tend, dt)
      call calc_dkedlon_dkedlat(state, tend, dt)
      call calc_dpedlon_dpedlat(static, state, tend, dt)
      call calc_dmfdlon_dmfdlat(state, tend, dt)

      do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
        do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
          tend%du(i,j) =   tend%qhv(i,j) - tend%dpedlon(i,j) - tend%dkedlon(i,j)
        end do
      end do

      do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%dv(i,j) = - tend%qhu(i,j) - tend%dpedlat(i,j) - tend%dkedlat(i,j)
        end do
      end do

      do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%dgd(i,j) = - (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * g
        end do
      end do
    case (slow_pass)
      call calc_qhu_qhv(state, tend, dt)

      do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
        do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
          tend%du(i,j) =   tend%qhv(i,j)
        end do
      end do

      do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%dv(i,j) = - tend%qhu(i,j)
        end do
      end do

      tend%dgd = 0.0_r8
    case (fast_pass)
!$omp sections
!$omp section
      call calc_dkedlon_dkedlat(state, tend, dt)
      call calc_dpedlon_dpedlat(static, state, tend, dt)
      do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
        do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
          tend%du(i,j) = - tend%dpedlon(i,j) - tend%dkedlon(i,j)
        end do
      end do
      do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%dv(i,j) = - tend%dpedlat(i,j) - tend%dkedlat(i,j)
        end do
      end do
!$omp section
      call calc_dmfdlon_dmfdlat(state, tend, dt)
      do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
        do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
          tend%dgd(i,j) = - (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * g
        end do
      end do
!$omp end sections
    end select

    ! call debug_check_space_operators(static, state, tend)

  end subroutine space_operators

  subroutine time_integrate(dt, static, tends, states)

    real(r8)         , intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (:)
    type(state_type ), intent(inout) :: states(:)

    call splitter(dt, static, tends, states)

  end subroutine time_integrate

  subroutine csp2_splitting(dt, static, tends, states)

    real(r8)         , intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (:)
    type(state_type ), intent(inout) :: states(:)

    real(r8) fast_dt
    integer subcycle, t1, t2

    fast_dt = dt / fast_cycles
    t1 = 3
    t2 = old

    call integrator(0.5_r8 * dt, static, tends, states, old, t1, slow_pass)
    do subcycle = 1, fast_cycles
      call integrator(fast_dt, static, tends, states, t1, t2, fast_pass)
      call time_swap_indices(t1, t2)
    end do 
    call integrator(0.5_r8 * dt, static, tends, states, t1, new, slow_pass)

  end subroutine csp2_splitting

  subroutine no_splitting(dt, static, tends, states)

    real(r8)         , intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (:)
    type(state_type ), intent(inout) :: states(:)

    call integrator(dt, static, tends, states, old, new, all_pass)

  end subroutine no_splitting

  subroutine predict_correct(dt, static, tends, states, old, new, pass)

    real(r8)         , intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (:)
    type(state_type ), intent(inout) :: states(:)
    integer          , intent(in   ) :: old
    integer          , intent(in   ) :: new
    integer          , intent(in   ) :: pass

    ! Do first predict step.
    call space_operators(static, states(old), tends(old), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, tends(old), states(old), states(new))

    ! Do second predict step.
    call space_operators(static, states(new), tends(old), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, tends(old), states(old), states(new))

    ! Do correct stepe
    call space_operators(static, states(new), tends(new),          dt, pass)
    call update_state(         dt, tends(new), states(old), states(new))

  end subroutine predict_correct

  subroutine runge_kutta_3rd(dt, static, tends, states, old, new, pass)

    real(r8)         , intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (:) ! -1:2
    type(state_type ), intent(inout) :: states(:) ! -1:2
    integer          , intent(in   ) :: old
    integer          , intent(in   ) :: new
    integer          , intent(in   ) :: pass

    integer s1, s2, s3

    s1 = 3
    s2 = 4
    s3 = new

    call space_operators(static, states(old), tends(s1), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, tends(s1), states(old), states(s1))

    call space_operators(static, states(s1) , tends(s2), 2.0_r8 * dt, pass)
    call update_state(        -dt, tends(s1), states(old), states(s2))
    call update_state(2.0_r8 * dt, tends(s2), states(s2) , states(s2))

    call space_operators(static, states(s2) , tends(s3),          dt, pass)
    tends(old)%du  = (tends(s1)%du  + 4.0_r8 * tends(s2)%du  + tends(s3)%du ) / 6.0_r8
    tends(old)%dv  = (tends(s1)%dv  + 4.0_r8 * tends(s2)%dv  + tends(s3)%dv ) / 6.0_r8
    tends(old)%dgd = (tends(s1)%dgd + 4.0_r8 * tends(s2)%dgd + tends(s3)%dgd) / 6.0_r8
    call update_state(         dt, tends(old), states(old), states(new))

  end subroutine runge_kutta_3rd

  subroutine runge_kutta_4th(dt, static, tends, states, old, new, pass)

    real(r8)         , intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (:) ! -1:2
    type(state_type ), intent(inout) :: states(:) ! -1:2
    integer          , intent(in   ) :: old
    integer          , intent(in   ) :: new
    integer          , intent(in   ) :: pass

    integer s1, s2, s3, s4

    s1 = 3
    s2 = 4
    s3 = 5
    s4 = new

    call space_operators(static, states(old), tends(s1), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, tends(s1), states(old), states(s1))

    call space_operators(static, states(s1) , tends(s2), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, tends(s2), states(old), states(s2))

    call space_operators(static, states(s2) , tends(s3),          dt, pass)
    call update_state(         dt, tends(s3), states(old), states(s3))

    call space_operators(static, states(s3) , tends(s4),          dt, pass)
    tends(old)%du  = (tends(s1)%du  + 2.0_r8 * tends(s2)%du  + 2.0_r8 * tends(s3)%du  + tends(s4)%du ) / 6.0_r8
    tends(old)%dv  = (tends(s1)%dv  + 2.0_r8 * tends(s2)%dv  + 2.0_r8 * tends(s3)%dv  + tends(s4)%dv ) / 6.0_r8
    tends(old)%dgd = (tends(s1)%dgd + 2.0_r8 * tends(s2)%dgd + 2.0_r8 * tends(s3)%dgd + tends(s4)%dgd) / 6.0_r8
    call update_state(         dt, tends(old), states(old), states(new))

  end subroutine runge_kutta_4th

  subroutine update_state(dt, tend, old_state, new_state)

    real(r8)        , intent(in   ) :: dt
    type(tend_type ), intent(in   ) :: tend
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: new_state

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => old_state%mesh

    do j = mesh%full_lat_start_idx, mesh%full_lat_end_idx
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
      end do
    end do

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      do i = mesh%half_lon_start_idx, mesh%half_lon_end_idx
        new_state%u(i,j) = old_state%u(i,j) + dt * tend%du(i,j)
      end do
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
      do i = mesh%full_lon_start_idx, mesh%full_lon_end_idx
        new_state%v(i,j) = old_state%v(i,j) + dt * tend%dv(i,j)
      end do
    end do

    call damp_state(new_state)

    call parallel_fill_halo(mesh, new_state%gd(:,:))
    call parallel_fill_halo(mesh, new_state%u (:,:))
    call parallel_fill_halo(mesh, new_state%v (:,:))

    ! Do not forget to synchronize the mass on edge and vertex for diagnosing!
    call calc_m_lon_m_lat(new_state)
    call calc_mf_lon_n_mf_lat_n(new_state)
    call calc_mf_lon_t_mf_lat_t(new_state)
    call calc_m_vtx(new_state)

    if (pv_scheme == 4) call diagnose(new_state)

  end subroutine update_state

  subroutine damp_state(state)

    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer j
    real(r8) wgt

    mesh => state%mesh

    if (adaptive_damp) then
      wgt = exp(- (elapsed_seconds - damp_t0) / 1800.0_r8)
      call log_add_diag('damp_wgt', wgt)
    else
      wgt = 1.0_r8
    end if

    do j = mesh%full_lat_start_idx_no_pole, mesh%full_lat_end_idx_no_pole
      if (reduced_mesh(j-1)%reduce_factor > 0 .or. reduced_mesh(j)%reduce_factor > 0 .or. reduced_mesh(j+1)%reduce_factor > 0) then
        call damp_run(damp_order, dt, mesh%de_lon(j), wgt, mesh%half_lon_lb, mesh%half_lon_ub, mesh%num_half_lon, state%u(:,j))
      end if
    end do

    do j = mesh%half_lat_start_idx_no_pole, mesh%half_lat_end_idx_no_pole
#ifdef STAGGER_V_ON_POLE
      if (reduced_mesh(j)%reduce_factor > 0 .or. reduced_mesh(j-1)%reduce_factor > 0) then
        call damp_run(damp_order, dt, mesh%le_lat(j), wgt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%num_full_lon, state%v(:,j))
      end if
#else
      if (reduced_mesh(j)%reduce_factor > 0 .or. reduced_mesh(j+1)%reduce_factor > 0) then
        call damp_run(damp_order, dt, mesh%le_lat(j), wgt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%num_full_lon, state%v(:,j))
      end if
#endif
    end do
#ifndef STAGGER_V_ON_POLE
    ! NOTE: These special treatment could be remove in future.
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_start_idx_no_pole
      state%v(:,j) = 0.5 * state%v(:,j+1) + 0.5 * state%v(:,j)
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_end_idx_no_pole
      state%v(:,j) = 0.5 * state%v(:,j-1) + 0.5 * state%v(:,j)
    end if
#endif

  end subroutine damp_state

end module gmcore_mod
