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
  use debug_mod

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
      import real_kind, static_type, tend_type, state_type
      real(real_kind  ), intent(in   ) :: dt
      type(static_type), intent(in   ) :: static
      type(tend_type  ), intent(inout) :: tends (0:2)
      type(state_type ), intent(inout) :: states(0:2)
      integer          , intent(in   ) :: old
      integer          , intent(in   ) :: new
      integer          , intent(in   ) :: pass
    end subroutine integrator_interface

    subroutine splitter_interface(dt, static, tends, states)
      import real_kind, static_type, tend_type, state_type
      real(real_kind  ), intent(in   ) :: dt
      type(static_type), intent(in   ) :: static
      type(tend_type  ), intent(inout) :: tends (0:2)
      type(state_type ), intent(inout) :: states(0:2)
    end subroutine splitter_interface
  end interface

  procedure(integrator_interface), pointer :: integrator
  procedure(splitter_interface), pointer :: splitter

contains

  subroutine gmcore_init()

    call log_init()
    call parallel_init()
    call time_init()
    call create_meshes()
    call create_states()
    call create_static()
    call create_tends()
    call history_init()

    select case (time_scheme)
    case ('predict_correct')
      integrator => predict_correct
    case default
      integrator => predict_correct
      call log_notice('Use predict_correct integrator.')
    end select

    select case (split_scheme)
    case ('csp2')
      splitter => csp2_splitting
    case default
      splitter => no_splitting
      call log_notice('No fast-slow split.')
    end select

  end subroutine gmcore_init

  subroutine gmcore_run()

    call operators_prepare(states(old))
    call diagnose(states(old))
    call output  (states(old), tends(old))
    call log_print_diag(curr_time%isoformat())

    do while (.not. time_is_finished())
      call time_integrate(dt, static, tends, states)
      call time_advance()
      call diagnose(states(old))
      call output  (states(old), tends(old))
      call log_print_diag(curr_time%isoformat())
    end do

  end subroutine gmcore_run

  subroutine gmcore_final()

    call parallel_final()

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

    integer i, j

    state%total_mass = 0.0d0
    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%total_mass = state%total_mass + state%hd(i,j) * state%mesh%cell_area(j)
      end do
    end do

    state%total_energy = 0.0d0
    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%total_energy = state%total_energy + g * state%hd(i,j) * (state%hd(i,j) / 2.0d0 + static%hs(i,j)) * state%mesh%cell_area(j)
      end do
    end do
    do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%total_energy = state%total_energy + state%mass_lon(i,j) * state%u(i,j)**2 / 2.0d0 * state%mesh%lon_edge_area(j)
      end do
    end do
    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
        state%total_energy = state%total_energy + state%mass_lat(i,j) * state%v(i,j)**2 / 2.0d0 * state%mesh%lat_edge_area(j)
      end do
    end do

    state%total_absolute_vorticity = 0.0d0
    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%total_absolute_vorticity = state%total_absolute_vorticity + &
                                         state%mass_vertex(i,j) * state%pv(i,j) * &
                                         state%mesh%vertex_area(j)
      end do
    end do

    state%total_potential_enstrophy = 0.0d0
    do j = state%mesh%half_lat_start_idx, state%mesh%half_lat_end_idx
      do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
        state%total_potential_enstrophy = state%total_potential_enstrophy + &
                                          state%mass_vertex(i,j) * state%pv(i,j)**2 * &
                                          state%mesh%vertex_area(j)
      end do
    end do

    call log_add_diag('total_mass', state%total_mass)
    call log_add_diag('total_energy', state%total_energy)
    call log_add_diag('total_absolute_vorticity', state%total_absolute_vorticity)
    call log_add_diag('total_potential_enstrophy', state%total_potential_enstrophy)

  end subroutine diagnose

  subroutine space_operators(static, state, tend, pass)

    type(static_type), intent(in   ) :: static
    type(state_type) , intent(inout) :: state
    type(tend_type)  , intent(inout) :: tend
    integer          , intent(in   ) :: pass

    integer i, j

    call operators_prepare(state)

    select case (pass)
    case (all_pass)
      call nonlinear_coriolis_operator(state, tend)
      call energy_gradient_operator(static, state, tend)
      call mass_flux_divergence_operator(state, tend)

      do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
        do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
          tend%du(i,j) =   tend%qhv(i,j) - tend%dEdlon(i,j)
        end do
      end do

      do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
        do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
          tend%dv(i,j) = - tend%qhu(i,j) - tend%dEdlat(i,j)
        end do
      end do

      do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
        do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
          tend%dhd(i,j) = - tend%div_mass_flux(i,j)
        end do
      end do
    case (slow_pass)
      call nonlinear_coriolis_operator(state, tend)

      do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
        do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
          tend%du(i,j) =   tend%qhv(i,j)
        end do
      end do

      do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
        do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
          tend%dv(i,j) = - tend%qhu(i,j)
        end do
      end do

      tend%dEdlon = 0.0d0
      tend%dEdlat = 0.0d0
      tend%div_mass_flux = 0.0d0
      tend%dhd = 0.0d0
    case (fast_pass)
      call energy_gradient_operator(static, state, tend)
      call mass_flux_divergence_operator(state, tend)

      do j = state%mesh%full_lat_start_idx_no_pole, state%mesh%full_lat_end_idx_no_pole
        do i = state%mesh%half_lon_start_idx, state%mesh%half_lon_end_idx
          tend%du(i,j) = - tend%dEdlon(i,j)
        end do
      end do

      do j = state%mesh%half_lat_start_idx_no_pole, state%mesh%half_lat_end_idx_no_pole
        do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
          tend%dv(i,j) = - tend%dEdlat(i,j)
        end do
      end do

      do j = state%mesh%full_lat_start_idx, state%mesh%full_lat_end_idx
        do i = state%mesh%full_lon_start_idx, state%mesh%full_lon_end_idx
          tend%dhd(i,j) = - tend%div_mass_flux(i,j)
        end do
      end do

      tend%qhv = 0.0d0
      tend%qhu = 0.0d0
    end select

    call debug_check_space_operators(static, state, tend)

  end subroutine space_operators

  subroutine time_integrate(dt, static, tends, states)

    real(real_kind  ), intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (0:2)
    type(state_type ), intent(inout) :: states(0:2)

    call splitter(dt, static, tends, states)

  end subroutine time_integrate

  subroutine csp2_splitting(dt, static, tends, states)

    real(real_kind  ), intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (0:2)
    type(state_type ), intent(inout) :: states(0:2)

    real(real_kind) fast_dt
    integer subcycle, t1, t2

    fast_dt = dt / fast_cycles
    t1 = 0 
    t2 = old

    call integrator(0.5 * dt, static, tends, states, old, t1, slow_pass)
    do subcycle = 1, fast_cycles
      call integrator(fast_dt, static, tends, states, t1, t2, fast_pass)
      call time_swap_indices(t1, t2)
    end do 
    call integrator(0.5 * dt, static, tends, states, t1, new, slow_pass)

  end subroutine csp2_splitting

  subroutine no_splitting(dt, static, tends, states)

    real(real_kind  ), intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (0:2)
    type(state_type ), intent(inout) :: states(0:2)

    call integrator(dt, static, tends, states, old, new, all_pass)

  end subroutine no_splitting

  subroutine predict_correct(dt, static, tends, states, old, new, pass)

    real(real_kind  ), intent(in   ) :: dt
    type(static_type), intent(in   ) :: static
    type(tend_type  ), intent(inout) :: tends (0:2)
    type(state_type ), intent(inout) :: states(0:2)
    integer          , intent(in   ) :: old
    integer          , intent(in   ) :: new
    integer          , intent(in   ) :: pass

    ! Do first predict step.
    call space_operators(static, states(old), tends(old), pass)
    call update_state(0.5d0 * dt, tends(old), states(old), states(new))

    ! Do second predict step.
    call space_operators(static, states(new), tends(old), pass)
    call update_state(0.5d0 * dt, tends(old), states(old), states(new))

    ! Do correct stepe
    call space_operators(static, states(new), tends(new), pass)
    call update_state(        dt, tends(new), states(old), states(new))

  end subroutine predict_correct

  subroutine update_state(dt, tend, old_state, new_state)

    real(real_kind ), intent(in   ) :: dt
    type(tend_type ), intent(in   ) :: tend
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: new_state

    integer i, j

    do j = new_state%mesh%full_lat_start_idx, new_state%mesh%full_lat_end_idx
      do i = new_state%mesh%full_lon_start_idx, new_state%mesh%full_lon_end_idx
        new_state%hd(i,j) = old_state%hd(i,j) + dt * tend%dhd(i,j)
      end do
    end do

    do j = new_state%mesh%full_lat_start_idx_no_pole, new_state%mesh%full_lat_end_idx_no_pole
      do i = new_state%mesh%half_lon_start_idx, new_state%mesh%half_lon_end_idx
        new_state%u(i,j) = old_state%u(i,j) + dt * tend%du(i,j)
      end do
    end do

    do j = new_state%mesh%half_lat_start_idx_no_pole, new_state%mesh%half_lat_end_idx_no_pole
      do i = new_state%mesh%full_lon_start_idx, new_state%mesh%full_lon_end_idx
        new_state%v(i,j) = old_state%v(i,j) + dt * tend%dv(i,j)
      end do
    end do

    call parallel_fill_halo(new_state%mesh, new_state%hd(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%mesh, new_state%u (:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%mesh, new_state%v (:,:), all_halo=.true.)

  end subroutine update_state

end module gmcore_mod
