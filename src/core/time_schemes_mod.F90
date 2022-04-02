module time_schemes_mod

  use flogger
  use const_mod
  use namelist_mod
  use tend_mod
  use block_mod
  use operators_mod
  use parallel_mod
  use nh_mod
  use filter_mod

  implicit none

  private

  public time_scheme_init
  public time_integrator
  public update_state
  public space_operators_interface

  interface
    subroutine space_operators_interface(block, old_state, star_state, new_state, tend1, tend2, dt, pass)
      import block_type, state_type, tend_type
      type(block_type), intent(inout) :: block
      type(state_type), intent(in   ) :: old_state
      type(state_type), intent(inout) :: star_state
      type(state_type), intent(inout) :: new_state
      type(tend_type ), intent(inout) :: tend1
      type(tend_type ), intent(in   ) :: tend2
      real(8), intent(in) :: dt
      integer, intent(in) :: pass
    end subroutine space_operators_interface

    subroutine step_interface(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)
      import space_operators_interface, block_type, state_type, tend_type
      procedure(space_operators_interface) space_operators
      type(block_type), intent(inout) :: block
      type(state_type), intent(in   ) :: old_state
      type(state_type), intent(inout) :: star_state
      type(state_type), intent(inout) :: new_state
      type(tend_type ), intent(inout) :: tend1
      type(tend_type ), intent(inout) :: tend2
      real(8), intent(in) :: dt
    end subroutine step_interface

    subroutine time_integrator_interface(space_operators, block, old, new, dt)
      import block_type, tend_type, state_type, space_operators_interface
      procedure(space_operators_interface) space_operators
      type(block_type), intent(inout) :: block
      integer, intent(in) :: old
      integer, intent(in) :: new
      real(8), intent(in) :: dt
    end subroutine time_integrator_interface
  end interface

  procedure(step_interface), pointer :: step
  procedure(time_integrator_interface), pointer :: time_integrator

contains

  subroutine time_scheme_init()

    select case (time_scheme)
    case ('euler')
      time_integrator => euler
    case ('pc2')
      time_integrator => predict_correct
    case ('wrfrk3')
      time_integrator => wrf_runge_kutta_3rd
    case default
      time_integrator => predict_correct
    end select

    if (advection) then
      step => step_adv
    else
      step => step_forward_backward
    end if

  end subroutine time_scheme_init

  subroutine step_all(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    type(tend_type ), intent(inout) :: tend1
    type(tend_type ), intent(inout) :: tend2
    real(8), intent(in) :: dt

    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, all_pass)
    call update_state(block, tend1, old_state, new_state, dt)

  end subroutine step_all

  subroutine step_forward_backward(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    type(tend_type ), intent(inout) :: tend1
    type(tend_type ), intent(inout) :: tend2
    real(8), intent(in) :: dt

    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, forward_pass)
    call update_state(block, tend1, old_state, new_state, dt)
    call space_operators(block, old_state, star_state, new_state, tend2, tend1, dt, backward_pass)
    call update_state(block, tend2, old_state, new_state, dt)

  end subroutine step_forward_backward

  subroutine step_adv(space_operators, block, old_state, star_state, new_state, tend1, tend2, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    type(tend_type ), intent(inout) :: tend1
    type(tend_type ), intent(inout) :: tend2
    real(8), intent(in) :: dt

    integer i, j, k, l, m, n
    real(r8) work(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend,block%mesh%num_full_lev)
    real(r8) pole(block%mesh%num_full_lev)

    call space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, all_pass)

    associate (mesh => block%mesh)
    do m = 1, size(block%adv_batches)
      do n = 1, size(block%adv_batches(m)%tracer_names)
        l = block%adv_batches(m)%tracer_idx(n)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%q(i,j,k,l) = old_state%q(i,j,k,l) - (                 &
                star_state%qmf_lon(i,j,k,l) - star_state%qmf_lon(i-1,j,k,l) + &
                star_state%qmf_lat(i,j,k,l) - star_state%qmf_lat(i,j-1,k,l)   &
              )
            end do
          end do
        end do
        if (mesh%has_south_pole()) then
          j = mesh%full_lat_ibeg
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              work(i,k) = star_state%qmf_lat(i,j,k,l)
            end do
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%q(i,j,k,l) = old_state%q(i,j,k,l) - pole(k)
            end do
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_lat_iend
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              work(i,k) = star_state%qmf_lat(i,j-1,k,l)
            end do
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%q(i,j,k,l) = old_state%q(i,j,k,l) + pole(k)
            end do
          end do
        end if
        call fill_halo(block, new_state%q(:,:,:,l), full_lon=.true., full_lat=.true., full_lev=.true.)
      end do
    end do
    end associate

  end subroutine step_adv

  subroutine update_state(block, tend, old_state, new_state, dt)

    type(block_type), intent(inout) :: block
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    integer i, j, k

    associate (mesh => block%mesh)
    if (baroclinic) then
      if (tend%update_phs) then
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%phs(i,j) = old_state%phs(i,j) + dt * tend%dphs(i,j)
          end do
        end do
        call fill_halo(block, new_state%phs, full_lon=.true., full_lat=.true.)
        call filter_on_cell(block, new_state%phs, new_state%phs_f)
        call fill_halo(block, new_state%phs_f, full_lon=.true., full_lat=.true.)

        call diag_ph(block, new_state)
        call diag_m (block, new_state)
        if (nonhydrostatic) then
          call diag_m_lev(block, new_state)
        end if
      else if (tend%copy_phs) then
        new_state%phs    = old_state%phs
        new_state%phs_f  = old_state%phs_f
        new_state%ph_lev = old_state%ph_lev
        new_state%ph     = old_state%ph
        new_state%m      = old_state%m
      end if

      if (tend%update_pt) then
        if (.not. tend%update_phs .and. .not. tend%copy_phs .and. is_root_proc()) call log_error('Mass is not updated or copied!')
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%pt(i,j,k) = (old_state%pt(i,j,k) * old_state%m(i,j,k) + dt * tend%dpt(i,j,k)) / new_state%m(i,j,k)
            end do
          end do
        end do
        call fill_halo(block, new_state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
        ! NOTE: When we use FFSL to transport pt, it may be not needed to filter pt.
        call filter_on_cell(block, new_state%pt, new_state%pt_f)
        call fill_halo(block, new_state%pt_f, full_lon=.true., full_lat=.true., full_lev=.true.)
      else if (tend%copy_pt) then
        new_state%pt   = old_state%pt
        new_state%pt_f = old_state%pt_f
      end if
    else
      if (tend%update_gz) then
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%gz(i,j,k) = old_state%gz(i,j,k) + dt * tend%dgz(i,j,k)
            end do
          end do
        end do
        call fill_halo(block, new_state%gz, full_lon=.true., full_lat=.true.)
        call filter_on_cell(block, new_state%gz, new_state%gz_f)
        call fill_halo(block, new_state%gz_f, full_lon=.true., full_lat=.true.)
      else if (tend%copy_gz) then
        new_state%gz   = old_state%gz
        new_state%gz_f = old_state%gz_f
      end if
    end if

    if (tend%update_u .and. tend%update_v) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            new_state%u(i,j,k) = old_state%u(i,j,k) + dt * tend%du(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%v(i,j,k) = old_state%v(i,j,k) + dt * tend%dv(i,j,k)
          end do
        end do
      end do
      call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)
      call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        new_state%v(:,j,:) = 0.4_r8 * new_state%v(:,j,:) + 0.6_r8 * new_state%v(:,j+1,:)
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        new_state%v(:,j,:) = 0.4_r8 * new_state%v(:,j,:) + 0.6_r8 * new_state%v(:,j-1,:)
      end if
      call filter_on_lon_edge(block, new_state%u, new_state%u_f)
      call filter_on_lat_edge(block, new_state%v, new_state%v_f)
      call fill_halo(block, new_state%u_f, full_lon=.false., full_lat=.true., full_lev=.true.)
      call fill_halo(block, new_state%v_f, full_lon=.true., full_lat=.false., full_lev=.true.)
    end if

    if (time_is_alerted('filter_reset')) then
      if (tend%update_phs) new_state%phs = new_state%phs_f
      if (tend%update_pt ) new_state%pt  = new_state%pt_f
      if (tend%update_u  ) new_state%u   = new_state%u_f
      if (tend%update_v  ) new_state%v   = new_state%v_f
    end if
    end associate

  end subroutine update_state

  subroutine predict_correct(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt

    associate (state => block%state, tend => block%tend)
    call step(space_operators, block, state(old), state(old), state(new), tend(old), tend(new), dt / 2.0_r8)
    call step(space_operators, block, state(old), state(new), state(3  ), tend(old), tend(new), dt / 2.0_r8)
    call step(space_operators, block, state(old), state(3  ), state(new), tend(old), tend(new), dt         )
    end associate

  end subroutine predict_correct

  subroutine wrf_runge_kutta_3rd(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt

    associate (state => block%state, tend => block%tend)
    call step(space_operators, block, state(old), state(old), state(new), tend(old), tend(new), dt / 3.0_r8)
    call step(space_operators, block, state(old), state(new), state(3  ), tend(old), tend(new), dt / 2.0_r8)
    call step(space_operators, block, state(old), state(3  ), state(new), tend(old), tend(new), dt         )
    end associate

  end subroutine wrf_runge_kutta_3rd

  subroutine euler(space_operators, block, old, new, dt)

    procedure(space_operators_interface) space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt

    associate (state => block%state, tend => block%tend)
    call step(space_operators, block, state(old), state(old), state(new), tend(old), tend(new), dt)
    end associate

  end subroutine euler

end module time_schemes_mod
