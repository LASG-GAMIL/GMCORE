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
  public time_scheme_final
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

    call time_scheme_final()

    select case (time_scheme)
    case ('euler')
      time_integrator => euler
    case ('pc2')
      time_integrator => pc2
    case ('wrfrk3')
      time_integrator => wrfrk3
    case default
      time_integrator => pc2
    end select

    step => step_forward_backward

  end subroutine time_scheme_init

  subroutine time_scheme_final()

  end subroutine time_scheme_final

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

  subroutine update_state(block, tend, old_state, new_state, dt)

    type(block_type), intent(inout) :: block
    type(tend_type ), intent(inout) :: tend
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    integer i, j, k
    real(r8) wgt
    real(r8) tmp1(global_mesh%num_full_lev)

    associate (mesh => block%mesh)
    if (baroclinic) then
      if (tend%update_phs) then
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%phs(i,j) = old_state%phs(i,j) + dt * tend%dphs(i,j)
          end do
        end do
        call fill_halo(block, new_state%phs, full_lon=.true., full_lat=.true., south_halo=.false., north_halo=.false.)
        call filter_on_cell(block%big_filter, new_state%phs, new_state%phs_f)
        call fill_halo(block, new_state%phs, full_lon=.true., full_lat=.true.)

        call calc_ph(block, new_state)
        call calc_m (block, new_state)
      else if (tend%copy_phs) then
        new_state%phs    = old_state%phs
        new_state%ph_lev = old_state%ph_lev
        new_state%ph     = old_state%ph
        new_state%m      = old_state%m
      end if

      if (tend%update_pt) then
        if (.not. tend%update_phs .and. .not. tend%copy_phs .and. is_root_proc()) call log_error('Mass is not updated or copied!')
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%pt(i,j,k) = old_state%pt(i,j,k) * old_state%m(i,j,k) + dt * tend%dpt(i,j,k)
            end do
          end do
        end do
        ! ----------------------------------------------------------------------
        call fill_halo(block, new_state%pt, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
        call filter_on_cell(block%small_filter1, new_state%pt)
        ! ----------------------------------------------------------------------
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%pt(i,j,k) = new_state%pt(i,j,k) / new_state%m(i,j,k)
            end do
          end do
        end do
        call fill_halo(block, new_state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)
      else if (tend%copy_pt) then
        new_state%pt = old_state%pt
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
        call fill_halo(block, new_state%gz, full_lon=.true., full_lat=.true., south_halo=.false., north_halo=.false.)
        call filter_on_cell(block%big_filter, new_state%gz, new_state%gz_f)
        call fill_halo(block, new_state%gz_f, full_lon=.true., full_lat=.true.)
        ! ----------------------------------------------------------------------
        call filter_on_cell(block%small_filter1, new_state%gz)
        ! ----------------------------------------------------------------------
        call fill_halo(block, new_state%gz, full_lon=.true., full_lat=.true.)
      else if (tend%copy_gz) then
        new_state%gz = old_state%gz
      end if
    end if

    if (tend%update_u .and. tend%update_v) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            new_state%u_lon(i,j,k) = old_state%u_lon(i,j,k) + dt * tend%du(i,j,k)
          end do
        end do
      end do
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%v_lat(i,j,k) = old_state%v_lat(i,j,k) + dt * tend%dv(i,j,k)
          end do
        end do
      end do
      ! ----------------------------------------------------------------------
      if (baroclinic) then
        call fill_halo(block, new_state%u_lon, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
        call filter_on_lon_edge(block%small_filter2, new_state%u_lon)
        call fill_halo(block, new_state%u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
        call fill_halo(block, new_state%v_lat, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
        call filter_on_lat_edge(block%small_filter2, new_state%v_lat)
        call fill_halo(block, new_state%v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
      end if
      ! ----------------------------------------------------------------------
      call fill_halo(block, new_state%u_lon, full_lon=.false., full_lat=.true., full_lev=.true.)
      call fill_halo(block, new_state%v_lat, full_lon=.true., full_lat=.false., full_lev=.true.)
      wgt = 0.6_r8
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        new_state%v_lat(:,j,:) = (1 - wgt) * new_state%v_lat(:,j,:) + wgt * new_state%v_lat(:,j+1,:)
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        new_state%v_lat(:,j,:) = (1 - wgt) * new_state%v_lat(:,j,:) + wgt * new_state%v_lat(:,j-1,:)
      end if
      call fill_halo(block, new_state%v_lat, full_lon=.true., full_lat=.false., full_lev=.true., west_halo=.false., east_halo=.false.)
      call filter_on_lon_edge(block%big_filter, new_state%u_lon, new_state%u_f)
      call filter_on_lat_edge(block%big_filter, new_state%v_lat, new_state%v_f)
      call fill_halo(block, new_state%u_f, full_lon=.false., full_lat=.true., full_lev=.true.)
      call fill_halo(block, new_state%v_f, full_lon=.true., full_lat=.false., full_lev=.true.)
    end if
    end associate

  end subroutine update_state

  subroutine pc2(space_operators, block, old, new, dt)

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

  end subroutine pc2

  subroutine wrfrk3(space_operators, block, old, new, dt)

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

  end subroutine wrfrk3

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
