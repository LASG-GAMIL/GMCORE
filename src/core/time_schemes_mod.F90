module time_schemes_mod

  use flogger
  use const_mod
  use namelist_mod
  use tend_mod
  use block_mod
  use operators_mod
  use parallel_mod

  implicit none

  private

  public time_scheme_init
  public time_integrator
  public predict_correct
  public runge_kutta_3rd
  public runge_kutta_4th
  public euler
  public ssp_runge_kutta_3rd
  public update_state

  interface
    subroutine space_operators_interface(block, old_state, star_state, new_state, tend, dt, pass)
      import block_type, state_type, tend_type
      type(block_type), intent(inout) :: block
      type(state_type), intent(in) :: old_state
      type(state_type), intent(inout) :: star_state
      type(state_type), intent(inout) :: new_state
      type(tend_type), intent(inout) :: tend
      real(8), intent(in) :: dt
      integer, intent(in) :: pass
    end subroutine space_operators_interface

    subroutine time_integrator_interface(space_operators, block, old, new, dt, pass)
      import block_type, tend_type, state_type, space_operators_interface
      procedure(space_operators_interface), intent(in), pointer :: space_operators
      type(block_type), intent(inout) :: block
      integer, intent(in) :: old
      integer, intent(in) :: new
      real(8), intent(in) :: dt
      integer, intent(in) :: pass
    end subroutine time_integrator_interface
  end interface

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
    case ('rk3')
      time_integrator => runge_kutta_3rd
    case ('rk4')
      time_integrator => runge_kutta_4th
    case ('ssprk3')
      time_integrator => ssp_runge_kutta_3rd
    case default
      time_integrator => predict_correct
    end select

  end subroutine time_scheme_init

  subroutine update_state(block, tend, old_state, new_state, dt, pass)

    type(block_type), intent(inout) :: block
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => old_state%mesh

    if (baroclinic) then
      if (tend%update_phs) then
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%phs(i,j) = old_state%phs(i,j) + dt * tend%dphs(i,j)
          end do
        end do
        call fill_halo(block, new_state%phs, full_lon=.true., full_lat=.true.)

        call diag_ph(block, new_state)
        call diag_m (block, new_state)
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
              new_state%pt(i,j,k) = (old_state%pt(i,j,k) * old_state%m(i,j,k) + dt * tend%dpt(i,j,k)) / new_state%m(i,j,k)
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
        call fill_halo(block, new_state%gz, full_lon=.true., full_lat=.true.)
      else if (tend%copy_gz) then
        new_state%gz = old_state%gz
      end if
    end if

    if (tend%update_u) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            new_state%u(i,j,k) = old_state%u(i,j,k) + dt * tend%du(i,j,k)
          end do
        end do
      end do
      call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)
    end if

    if (tend%update_v) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%v(i,j,k) = old_state%v(i,j,k) + dt * tend%dv(i,j,k)
          end do
        end do
      end do
      if (limit_pole_v) then
        if (mesh%has_south_pole()) then
          j = mesh%half_lat_ibeg
          new_state%v(:,j,:) = (1 - limit_pole_v_wgt) * new_state%v(:,j,:) + limit_pole_v_wgt * new_state%v(:,j+1,:)
        end if
        if (mesh%has_north_pole()) then
          j = mesh%half_lat_iend
          new_state%v(:,j,:) = (1 - limit_pole_v_wgt) * new_state%v(:,j,:) + limit_pole_v_wgt * new_state%v(:,j-1,:)
        end if
      end if
      call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end if

    call operators_prepare(block, new_state, dt, pass)

  end subroutine update_state

  subroutine predict_correct(space_operators, block, old, new, dt, pass)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    call space_operators(block, block%state(old), block%state(old), block%state(new), block%tend(old), 0.5_r8 * dt, pass)
    call update_state(block, block%tend(old), block%state(old), block%state(new), 0.5_r8 * dt, pass)

    call space_operators(block, block%state(old), block%state(new), block%state(3), block%tend(new), 0.5_r8 * dt, pass)
    call update_state(block, block%tend(new), block%state(old), block%state(3), 0.5_r8 * dt, pass)

    call space_operators(block, block%state(old), block%state(3), block%state(new), block%tend(3), dt, pass)
    call update_state(block, block%tend(3), block%state(old), block%state(new), dt, pass)

  end subroutine predict_correct

  subroutine wrf_runge_kutta_3rd(space_operators, block, old, new, dt, pass)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    call space_operators(block, block%state(old), block%state(old), block%state(new), block%tend(old), dt / 3.0_r8, pass)
    call update_state(block, block%tend(old), block%state(old), block%state(new), dt / 3.0_r8, pass)

    call space_operators(block, block%state(old), block%state(new), block%state(3), block%tend(new), 0.5_r8 * dt, pass)
    call update_state(block, block%tend(new), block%state(old), block%state(3), 0.5_r8 * dt, pass)

    call space_operators(block, block%state(old), block%state(3), block%state(new), block%tend(3), dt, pass)
    call update_state(block, block%tend(3), block%state(old), block%state(new), dt, pass)

  end subroutine wrf_runge_kutta_3rd

  subroutine runge_kutta_3rd(space_operators, block, old, new, dt, pass)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    integer s1, s2, s3

    s1 = 3
    s2 = 4
    s3 = new

    call space_operators(block, block%state(old), block%state(old), block%state(s1), block%tend(s1), 0.5_r8 * dt, pass)
    call update_state(block, block%tend(s1), block%state(old), block%state(s1), 0.5_r8 * dt, pass)

    call space_operators(block, block%state(old), block%state(s1), block%state(s2), block%tend(s2), 2.0_r8 * dt, pass)
    call update_state(block, block%tend(s1), block%state(old), block%state(s2), -dt, pass)
    call update_state(block, block%tend(s2), block%state(s2) , block%state(s2), 2.0_r8 * dt, pass)

    call space_operators(block, block%state(old), block%state(s2), block%state(new), block%tend(s3), dt, pass)
    block%tend(old) = (block%tend(s1) + 4.0_r8 * block%tend(s2) + block%tend(s3)) / 6.0_r8
    call update_state(block, block%tend(old), block%state(old), block%state(new), dt, pass)

  end subroutine runge_kutta_3rd

  subroutine runge_kutta_4th(space_operators, block, old, new, dt, pass)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    integer s1, s2, s3, s4

    s1 = 3
    s2 = 4
    s3 = 5
    s4 = new

    call space_operators(block, block%state(old), block%state(old), block%state(s1), block%tend(s1), 0.5_r8 * dt, pass)
    call update_state(block, block%tend(s1), block%state(old), block%state(s1), 0.5_r8 * dt, pass)

    call space_operators(block, block%state(old), block%state(s1), block%state(s2), block%tend(s2), 0.5_r8 * dt, pass)
    call update_state(block, block%tend(s2), block%state(old), block%state(s2), 0.5_r8 * dt, pass)

    call space_operators(block, block%state(old), block%state(s2), block%state(s3), block%tend(s3), dt, pass)
    call update_state(block, block%tend(s3), block%state(old), block%state(s3), dt, pass)

    call space_operators(block, block%state(old), block%state(s3), block%state(new), block%tend(s4), dt, pass)
    block%tend(old) = (block%tend(s1) + 2.0_r8 * block%tend(s2) + 2.0_r8 * block%tend(s3) + block%tend(s4)) / 6.0_r8
    call update_state(block, block%tend(old), block%state(old), block%state(new), dt, pass)

  end subroutine runge_kutta_4th

  subroutine euler(space_operators, block, old, new, dt, pass)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    call space_operators(block, block%state(old), block%state(old), block%state(new), block%tend(old), dt, pass)
    call update_state(block, block%tend(old), block%state(old), block%state(new), dt, pass)

  end subroutine euler

  subroutine ssp_runge_kutta_3rd(space_operators, block, old, new, dt, pass)

    procedure(space_operators_interface), intent(in), pointer :: space_operators
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    integer s1, s2, s3

    s1 = 3
    s2 = 4
    s3 = new

    call space_operators(block, block%state(old), block%state(old), block%state(s1), block%tend(s1), dt, pass)
    call update_state(block, block%tend(s1), block%state(old), block%state(s1), dt, pass)

    call space_operators(block, block%state(old), block%state(s1), block%state(s2), block%tend(s2), dt, pass)
    block%tend(s3) = block%tend(s1) + block%tend(s2)
    call update_state(block, block%tend(s3), block%state(old), block%state(s2), 0.25_r8 * dt, pass)

    call space_operators(block, block%state(old), block%state(s2), block%state(new), block%tend(s3), 0.5_r8 * dt, pass)
    block%tend(old) = (block%tend(s1) + block%tend(s2) + 4.0_r8 * block%tend(s3)) / 6.0_r8
    call update_state(block, block%tend(old), block%state(old), block%state(new), dt, pass)

  end subroutine ssp_runge_kutta_3rd

end module time_schemes_mod
