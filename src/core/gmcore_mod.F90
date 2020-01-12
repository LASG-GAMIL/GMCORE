module gmcore_mod

  use flogger
  use const_mod
  use namelist_mod
  use process_mod
  use parallel_mod
  use time_mod, dt => dt_in_seconds, old => old_time_idx, new => new_time_idx
  use history_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use block_mod
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
    subroutine integrator_interface(dt, block, old, new, pass)
      import r8, block_type, tend_type, state_type
      real(r8), intent(in) :: dt
      type(block_type), intent(inout) :: block
      integer, intent(in) :: old
      integer, intent(in) :: new
      integer, intent(in) :: pass
    end subroutine integrator_interface

    subroutine splitter_interface(dt, block)
      import r8, block_type
      real(r8), intent(in) :: dt
      type(block_type), intent(inout) :: block
    end subroutine splitter_interface
  end interface

  procedure(integrator_interface), pointer :: integrator
  procedure(splitter_interface), pointer :: splitter

  real(r8) :: damp_t0 = -1

contains

  subroutine gmcore_init()

    call log_init()
    call global_mesh%init(num_lon, num_lat)
    call debug_check_areas()
    call process_init()
    call time_init()
    call history_init()
    call reduce_init(proc%blocks)

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

    call operators_prepare(proc%blocks, old)
    call diagnose(proc%blocks, old)
    call output(proc%blocks, old)
    call log_print_diag(curr_time%isoformat())

    do while (.not. time_is_finished())
      call time_integrate(dt, proc%blocks)
      if (time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      call time_advance()
      call operators_prepare(proc%blocks, old)
      call diagnose(proc%blocks, old)
      call output(proc%blocks, old)
    end do

  end subroutine gmcore_run

  subroutine gmcore_final()

    call process_final()

  end subroutine gmcore_final

  subroutine output(blocks, itime)

    type(block_type), intent(in) :: blocks(:)
    integer, intent(in) :: itime

    if (time_is_alerted('history_write')) then
      call history_write_state(blocks, itime)
      call history_write_debug(blocks, itime)
    end if

  end subroutine output

  subroutine diagnose(blocks, itime)

    type(block_type), intent(inout), target :: blocks(:)
    integer, intent(in) :: itime

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    integer i, j, iblk
    real(r8) tm, te, tav, tpe

    do iblk = 1, size(blocks)
      mesh => blocks(iblk)%mesh
      state => blocks(iblk)%state(itime)
      static => blocks(iblk)%static

      tm = 0.0_r8
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tm = tm + state%gd(i,j) * mesh%cell_area(j)
        end do
      end do

      te = 0.0_r8
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          te = te + state%mf_lon_n(i,j) * 0.5_r8 * state%u(i,j) * mesh%lon_edge_area(j) * 2
        end do
      end do
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          te = te + state%mf_lat_n(i,j) * 0.5_r8 * state%v(i,j) * mesh%lat_edge_area(j) * 2
        end do
      end do
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          te = te + (state%gd(i,j)**2 / g * 0.5_r8 + state%gd(i,j) * static%ghs(i,j) / g) * mesh%cell_area(j)
        end do
      end do

      tav = 0.0_r8
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tav = tav + state%m_vtx(i,j) * state%pv(i,j) * mesh%vertex_area(j)
        end do
      end do

      tpe = 0.0_r8
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tpe = tpe + state%m_vtx(i,j) * state%pv(i,j)**2 * 0.5_r8 * mesh%vertex_area(j)
        end do
      end do
    end do

    call log_add_diag('total_m' , tm )
    call log_add_diag('total_e' , te )
    call log_add_diag('total_pe', tpe)

  end subroutine diagnose

  subroutine space_operators(block, state, tend, dt, pass)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass

    type(mesh_type), pointer :: mesh
    integer i, j

    call state%update_all()
    call reduce_run(block, state, dt)

    mesh => state%mesh

    select case (pass)
    case (all_pass)
      call calc_qhu_qhv(block, state, tend, dt)
      call calc_dkedlon_dkedlat(block, state, tend, dt)
      call calc_dpedlon_dpedlat(block, state, tend, dt)
      call calc_dmfdlon_dmfdlat(block, state, tend, dt)

      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%du(i,j) =   tend%qhv(i,j) - tend%dpedlon(i,j) - tend%dkedlon(i,j)
        end do
      end do

      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dv(i,j) = - tend%qhu(i,j) - tend%dpedlat(i,j) - tend%dkedlat(i,j)
        end do
      end do

      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dgd(i,j) = - (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * g
        end do
      end do
    case (slow_pass)
      call calc_qhu_qhv(block, state, tend, dt)

      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%du(i,j) =   tend%qhv(i,j)
        end do
      end do

      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dv(i,j) = - tend%qhu(i,j)
        end do
      end do

      tend%dgd = 0.0_r8
    case (fast_pass)
!$omp sections
!$omp section
      call calc_dkedlon_dkedlat(block, state, tend, dt)
      call calc_dpedlon_dpedlat(block, state, tend, dt)
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          tend%du(i,j) = - tend%dpedlon(i,j) - tend%dkedlon(i,j)
        end do
      end do
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dv(i,j) = - tend%dpedlat(i,j) - tend%dkedlat(i,j)
        end do
      end do
!$omp section
      call calc_dmfdlon_dmfdlat(block, state, tend, dt)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          tend%dgd(i,j) = - (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * g
        end do
      end do
!$omp end sections
    end select

    ! call debug_check_space_operators(static, state, tend)

  end subroutine space_operators

  subroutine time_integrate(dt, blocks)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      call splitter(dt, blocks(iblk))
    end do

  end subroutine time_integrate

  subroutine csp2_splitting(dt, block)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block

    real(r8) fast_dt
    integer subcycle, t1, t2

    fast_dt = dt / fast_cycles
    t1 = 3
    t2 = old

    call integrator(0.5_r8 * dt, block, old, t1, slow_pass)
    do subcycle = 1, fast_cycles
      call integrator(fast_dt, block, t1, t2, fast_pass)
      call time_swap_indices(t1, t2)
    end do
    call integrator(0.5_r8 * dt, block, t1, new, slow_pass)

  end subroutine csp2_splitting

  subroutine no_splitting(dt, block)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block

    call integrator(dt, block, old, new, all_pass)

  end subroutine no_splitting

  subroutine predict_correct(dt, block, old, new, pass)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass

    ! Do first predict step.
    call space_operators(block, block%state(old), block%tend(old), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, block, block%tend(old), block%state(old), block%state(new))

    ! Do second predict step.
    call space_operators(block, block%state(new), block%tend(old), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, block, block%tend(old), block%state(old), block%state(new))

    ! Do correct stepe
    call space_operators(block, block%state(new), block%tend(new),          dt, pass)
    call update_state(         dt, block, block%tend(new), block%state(old), block%state(new))

  end subroutine predict_correct

  subroutine runge_kutta_3rd(dt, block, old, new, pass)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass

    integer s1, s2, s3

    s1 = 3
    s2 = 4
    s3 = new

    call space_operators(block, block%state(old), block%tend(s1), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, block, block%tend(s1), block%state(old), block%state(s1))

    call space_operators(block, block%state(s1) , block%tend(s2), 2.0_r8 * dt, pass)
    call update_state(        -dt, block, block%tend(s1), block%state(old), block%state(s2))
    call update_state(2.0_r8 * dt, block, block%tend(s2), block%state(s2) , block%state(s2))

    call space_operators(block, block%state(s2) , block%tend(s3),          dt, pass)
    block%tend(old)%du  = (block%tend(s1)%du  + 4.0_r8 * block%tend(s2)%du  + block%tend(s3)%du ) / 6.0_r8
    block%tend(old)%dv  = (block%tend(s1)%dv  + 4.0_r8 * block%tend(s2)%dv  + block%tend(s3)%dv ) / 6.0_r8
    block%tend(old)%dgd = (block%tend(s1)%dgd + 4.0_r8 * block%tend(s2)%dgd + block%tend(s3)%dgd) / 6.0_r8
    call update_state(         dt, block, block%tend(old), block%state(old), block%state(new))

  end subroutine runge_kutta_3rd

  subroutine runge_kutta_4th(dt, block, old, new, pass)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass

    integer s1, s2, s3, s4

    s1 = 3
    s2 = 4
    s3 = 5
    s4 = new

    call space_operators(block, block%state(old), block%tend(s1), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, block, block%tend(s1), block%state(old), block%state(s1))

    call space_operators(block, block%state(s1) , block%tend(s2), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, block, block%tend(s2), block%state(old), block%state(s2))

    call space_operators(block, block%state(s2) , block%tend(s3),          dt, pass)
    call update_state(         dt, block, block%tend(s3), block%state(old), block%state(s3))

    call space_operators(block, block%state(s3) , block%tend(s4),          dt, pass)
    block%tend(old)%du  = (block%tend(s1)%du  + 2.0_r8 * block%tend(s2)%du  + 2.0_r8 * block%tend(s3)%du  + block%tend(s4)%du ) / 6.0_r8
    block%tend(old)%dv  = (block%tend(s1)%dv  + 2.0_r8 * block%tend(s2)%dv  + 2.0_r8 * block%tend(s3)%dv  + block%tend(s4)%dv ) / 6.0_r8
    block%tend(old)%dgd = (block%tend(s1)%dgd + 2.0_r8 * block%tend(s2)%dgd + 2.0_r8 * block%tend(s3)%dgd + block%tend(s4)%dgd) / 6.0_r8
    call update_state(         dt, block, block%tend(old), block%state(old), block%state(new))

  end subroutine runge_kutta_4th

  subroutine update_state(dt, block, tend, old_state, new_state)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    type(mesh_type), pointer :: mesh
    integer i, j

    mesh => old_state%mesh

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        new_state%gd(i,j) = old_state%gd(i,j) + dt * tend%dgd(i,j)
      end do
    end do
    call fill_halo(mesh, new_state%gd)

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        new_state%u(i,j) = old_state%u(i,j) + dt * tend%du(i,j)
      end do
    end do
    call fill_halo(mesh, new_state%u)

    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        new_state%v(i,j) = old_state%v(i,j) + dt * tend%dv(i,j)
      end do
    end do
    call fill_halo(mesh, new_state%v)

    call damp_state(block, new_state)

  end subroutine update_state

  subroutine damp_state(block, state)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer j, damp_order
    real(r8) wgt

    mesh => state%mesh

    if (adaptive_damp) then
      wgt = exp(- (elapsed_seconds - damp_t0) / 1800.0_r8)
      call log_add_diag('damp_wgt', wgt)
    else
      wgt = 1.0_r8
    end if

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      damp_order = block%reduced_mesh(j)%damp_order
      if (damp_order > 0) then
        call damp_run(damp_order, dt, mesh%de_lon(j), wgt, mesh%half_lon_lb, mesh%half_lon_ub, mesh%num_half_lon, state%u(:,j))
      end if
    end do

    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
#ifdef V_POLE
      damp_order = max(block%reduced_mesh(j)%damp_order, block%reduced_mesh(j-1)%damp_order)
#else
      damp_order = max(block%reduced_mesh(j)%damp_order, block%reduced_mesh(j+1)%damp_order)
#endif
      if (damp_order > 0) then
        call damp_run(damp_order, dt, mesh%le_lat(j), wgt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%num_full_lon, state%v(:,j))
      end if
    end do

  end subroutine damp_state

end module gmcore_mod
