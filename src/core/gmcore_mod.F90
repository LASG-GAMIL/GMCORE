module gmcore_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use process_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use history_mod
  use restart_mod
  use block_mod
  use vert_coord_mod
  use operators_mod
  use interp_mod
  use reduce_mod
  use debug_mod
  use damp_mod

  implicit none

  private

  public gmcore_init
  public gmcore_run
  public gmcore_final

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

contains

  subroutine gmcore_init(namelist_path)

    character(*), intent(in) :: namelist_path

    character(10) time_value, time_units
    real(r8) seconds

    call log_init()
    call global_mesh%init_global(num_lon, num_lat, num_lev, lon_halo_width=max(1, maxval(reduce_factors) - 1), lat_halo_width=2)
    call debug_check_areas()
    call process_init()
    call vert_coord_init(num_lev, namelist_path)
    call process_create_blocks()
    call time_init()
    call history_init()
    call restart_init()
    call reduce_init(proc%blocks)
    call damp_init()

    select case (time_scheme)
    case ('debug')
      integrator => euler_debug
    case ('pc2')
      integrator => predict_correct
    case ('pc2+fb')
      integrator => predict_correct_with_forward_backward
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
      if (is_root_proc()) call log_notice('No fast-slow split.')
    end select

    time_value = split_string(print_interval, ' ', 1)
    time_units = split_string(print_interval, ' ', 2)
    read(time_value, *) seconds
    select case (time_units)
    case ('days')
      seconds = seconds * 86400
    case ('hours')
      seconds = seconds * 3600
    case ('minutes')
      seconds = seconds * 60
    case ('seconds')
      seconds = seconds
    case default
      call log_error('Invalid print interval ' // trim(print_interval) // '!')
    end select

    call time_add_alert('print', seconds=seconds)

  end subroutine gmcore_init

  subroutine gmcore_run()

    call operators_prepare(proc%blocks, old, dt_in_seconds)
    call diagnose(proc%blocks, old)
    call output(proc%blocks, old)

    do while (.not. time_is_finished())
      call time_integrate(dt_in_seconds, proc%blocks)
      if (proc%id == 0 .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      call time_advance(dt_in_seconds)
      call operators_prepare(proc%blocks, old, dt_in_seconds)
      call diagnose(proc%blocks, old)
      call output(proc%blocks, old)
    end do

  end subroutine gmcore_run

  subroutine gmcore_final()

    call damp_final()
    call history_final()
    call process_final()

  end subroutine gmcore_final

  subroutine output(blocks, itime)

    type(block_type), intent(in), target :: blocks(:)
    integer, intent(in) :: itime

    type(state_type), pointer :: state
    real(r8), save :: time1 = 0, time2
    integer i, j, k, iblk

    if (time_is_alerted('history_write')) then
      if (time_step == 0) call cpu_time(time1)
      call cpu_time(time2)
      if (time_step /= 0) then
        if (is_root_proc()) call log_notice('Time cost ' // to_string(time2 - time1, 5) // ' seconds.')
        time1 = time2
      end if
      ! Interpolate onto isobaric layers.
      do iblk = 1, size(blocks)
        state => blocks(iblk)%state(itime)
        call interp_lon_edge_to_isobaric_level(state%mesh, state%ph, state%u, 85000.0_r8, state%u850, logp=.true.)
        call interp_lon_edge_to_isobaric_level(state%mesh, state%ph, state%u, 70000.0_r8, state%u700, logp=.true.)
        call interp_lat_edge_to_isobaric_level(state%mesh, state%ph, state%v, 85000.0_r8, state%v850, logp=.true.)
        call interp_lat_edge_to_isobaric_level(state%mesh, state%ph, state%v, 70000.0_r8, state%v700, logp=.true.)
        call interp_cell_to_isobaric_level(state%mesh, state%ph, state%t, 85000.0_r8, state%t850, logp=.true.)
        call interp_cell_to_isobaric_level(state%mesh, state%ph, state%t, 70000.0_r8, state%t700, logp=.true.)
      end do
      call history_write_state(blocks, itime)
      call history_write_debug(blocks, itime)
    end if
    if (time_is_alerted('restart_write')) then
      call restart_write(blocks, itime)
    end if

  end subroutine output

  subroutine diagnose(blocks, itime)

    type(block_type), intent(inout), target :: blocks(:)
    integer, intent(in) :: itime

    type(mesh_type), pointer :: mesh
    type(state_type), pointer :: state
    type(static_type), pointer :: static
    integer i, j, k, iblk
    real(r8) tm, te, tav, tpe

    tm = 0.0_r8
    te = 0.0_r8
    tav = 0.0_r8
    tpe = 0.0_r8
    do iblk = 1, size(blocks)
      mesh => blocks(iblk)%mesh
      state => blocks(iblk)%state(itime)
      static => blocks(iblk)%static

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tm = tm + state%m(i,j,k) * mesh%area_cell(j)
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            te = te + state%mf_lon_n(i,j,k) * 0.5_r8 * state%u(i,j,k) * mesh%area_lon(j) * 2
          end do
        end do
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            te = te + state%mf_lat_n(i,j,k) * 0.5_r8 * state%v(i,j,k) * mesh%area_lat(j) * 2
          end do
        end do
      end do
      if (baroclinic .and. hydrostatic) then
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              te = te + state%m(i,j,k) * cp * state%t(i,j,k) * mesh%area_cell(j)
            end do
          end do
        end do
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            te = te + static%gzs(i,j) * state%phs(i,j) * mesh%area_cell(j)
          end do
        end do
      else if (.not. baroclinic) then
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            te = te + (state%m(i,j,1)**2 * g * 0.5_r8 + state%m(i,j,1) * static%gzs(i,j)) * mesh%area_cell(j)
          end do
        end do
      end if

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            tav = tav + state%m_vtx(i,j,k) * state%pv(i,j,k) * mesh%area_vtx(j)
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            tpe = tpe + state%m_lon(i,j,k) * state%pv_lon(i,j,k)**2 * 0.5_r8 * mesh%area_lon(j)
          end do
        end do
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tpe = tpe + state%m_lat(i,j,k) * state%pv_lat(i,j,k)**2 * 0.5_r8 * mesh%area_lat(j)
          end do
        end do
      end do
    end do
    call global_sum(proc%comm, tm)
    call global_sum(proc%comm, te)
    call global_sum(proc%comm, tav)
    call global_sum(proc%comm, tpe)

    do iblk = 1, size(blocks)
      blocks(iblk)%state(itime)%tm  = tm
      blocks(iblk)%state(itime)%te  = te
      blocks(iblk)%state(itime)%tpe = tpe
    end do

    call log_add_diag('tm' , tm )
    call log_add_diag('te' , te )
    call log_add_diag('tpe', tpe)

  end subroutine diagnose

  subroutine space_operators(block, state, tend, dt, pass)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    real(r8), intent(in) :: dt
    integer, intent(in) :: pass

    type(mesh_type), pointer :: mesh
    integer i, j, k

    call wait_halo(state%async(async_gz))
    call wait_halo(state%async(async_u))
    call wait_halo(state%async(async_v))

    call operators_prepare(block, state, dt, pass)
    call reduce_run(block, state, dt, pass)

    mesh => state%mesh

    call tend%reset_flags()

    select case (pass)
    case (all_pass)
      if (baroclinic .and. hydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, state, tend, dt)
        call calc_dphs             (block, state, tend, dt)
        call calc_wedphdlev        (block, state, tend, dt)
        call calc_wedudlev_wedvdlev(block, state, tend, dt)
        call calc_dptfdlon_dptfdlat(block, state, tend, dt)
        call calc_dptfdlev         (block, state, tend, dt)
        call calc_qhu_qhv          (block, state, tend, dt)
        call calc_dkedlon_dkedlat  (block, state, tend, dt)
        call calc_dpedlon_dpedlat  (block, state, tend, dt)
        call calc_dpdlon_dpdlat    (block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) =   tend%qhv(i,j,k) - tend%dpedlon(i,j,k) - tend%dkedlon(i,j,k) - tend%dpdlon(i,j,k) - tend%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%qhu(i,j,k) - tend%dpedlat(i,j,k) - tend%dkedlat(i,j,k) - tend%dpdlat(i,j,k) - tend%wedvdlev(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dpt(i,j,k) = - tend%dptfdlon(i,j,k) - tend%dptfdlat(i,j,k) - tend%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend%updated_du   = .true.
        tend%updated_dv   = .true.
        tend%updated_dphs = .true.
        tend%updated_dpt  = .true.
      else
        call calc_dmfdlon_dmfdlat(block, state, tend, dt)
        call calc_qhu_qhv        (block, state, tend, dt)
        call calc_dkedlon_dkedlat(block, state, tend, dt)
        call calc_dpedlon_dpedlat(block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) =   tend%qhv(i,j,k) - tend%dpedlon(i,j,k) - tend%dkedlon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%qhu(i,j,k) - tend%dpedlat(i,j,k) - tend%dkedlat(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dgz(i,j,k) = - (tend%dmfdlon(i,j,k) + tend%dmfdlat(i,j,k)) * g
            end do
          end do
        end do

        tend%updated_du  = .true.
        tend%updated_dv  = .true.
        tend%updated_dgz = .true.
      end if
    case (all_pass + forward_pass)
      if (baroclinic .and. hydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, state, tend, dt)
        call calc_dphs             (block, state, tend, dt)
        call calc_wedphdlev        (block, state, tend, dt)
        call calc_dptfdlon_dptfdlat(block, state, tend, dt)
        call calc_dptfdlev         (block, state, tend, dt)
        call calc_wedudlev_wedvdlev(block, state, tend, dt)
        call calc_qhu_qhv          (block, state, tend, dt)
        call calc_dkedlon_dkedlat  (block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) =   tend%qhv(i,j,k) - tend%dkedlon(i,j,k) - tend%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%qhu(i,j,k) - tend%dkedlat(i,j,k) - tend%wedvdlev(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dpt(i,j,k) = - tend%dptfdlon(i,j,k) - tend%dptfdlat(i,j,k) - tend%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend%updated_dphs = .true.
        tend%updated_dpt  = .true.
      else
        if (is_root_proc()) call log_error('Unfinished branch!', __FILE__, __LINE__)
      end if
    case (all_pass + backward_pass)
      if (baroclinic .and. hydrostatic) then
        call calc_dpedlon_dpedlat  (block, state, tend, dt)
        call calc_dpdlon_dpdlat    (block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) = tend%du(i,j,k) - tend%dpedlon(i,j,k) - tend%dpdlon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = tend%dv(i,j,k) - tend%dpedlat(i,j,k) - tend%dpdlat(i,j,k)
            end do
          end do
        end do

        tend%updated_du = .true.
        tend%updated_dv = .true.
      else
        if (is_root_proc()) call log_error('Unfinished branch!', __FILE__, __LINE__)
      end if
    case (slow_pass)
      if (baroclinic .and. hydrostatic) then
        call calc_qhu_qhv          (block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) =   tend%qhv(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%qhu(i,j,k)
            end do
          end do
        end do

        tend%updated_du   = .true.
        tend%updated_dv   = .true.
        tend%copy_pt      = .true.
        tend%copy_phs     = .true.
      else
        call calc_qhu_qhv          (block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) =   tend%qhv(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%qhu(i,j,k)
            end do
          end do
        end do

        tend%updated_du = .true.
        tend%updated_dv = .true.
        tend%copy_gz    = .true.
      end if
    case (fast_pass)
      if (baroclinic .and. hydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, state, tend, dt)
        call calc_dphs             (block, state, tend, dt)
        call calc_wedphdlev        (block, state, tend, dt)
        call calc_wedudlev_wedvdlev(block, state, tend, dt)
        call calc_dptfdlev         (block, state, tend, dt)
        call calc_dptfdlon_dptfdlat(block, state, tend, dt)
        call calc_dpedlon_dpedlat  (block, state, tend, dt)
        call calc_dkedlon_dkedlat  (block, state, tend, dt)
        call calc_dpdlon_dpdlat    (block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) = - tend%wedudlev(i,j,k) - tend%dkedlon(i,j,k) - tend%dpedlon(i,j,k) - tend%dpdlon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%wedvdlev(i,j,k) - tend%dkedlat(i,j,k) - tend%dpedlat(i,j,k) - tend%dpdlat(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dpt(i,j,k) = - tend%dptfdlon(i,j,k) - tend%dptfdlat(i,j,k) - tend%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend%updated_du   = .true.
        tend%updated_dv   = .true.
        tend%updated_dpt  = .true.
        tend%updated_dphs = .true.
      else
        call calc_dkedlon_dkedlat(block, state, tend, dt)
        call calc_dpedlon_dpedlat(block, state, tend, dt)
        call calc_dmfdlon_dmfdlat(block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) = - tend%dpedlon(i,j,k) - tend%dkedlon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%dpedlat(i,j,k) - tend%dkedlat(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dgz(i,j,k) = - (tend%dmfdlon(i,j,k) + tend%dmfdlat(i,j,k)) * g
            end do
          end do
        end do

        tend%updated_du  = .true.
        tend%updated_dv  = .true.
        tend%updated_dgz = .true.
      end if
    end select

    ! call debug_check_space_operators(block, state, tend)

  end subroutine space_operators

  subroutine time_integrate(dt, blocks)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      call splitter(dt, blocks(iblk))

      if (use_polar_damp) then
        call latlon_damp_lon(blocks(iblk), polar_damp_order, dt, blocks(iblk)%state(new)%u)
        call latlon_damp_lat(blocks(iblk), polar_damp_order, dt, blocks(iblk)%state(new)%v)
      end if
      if (use_div_damp) then
        call div_damp(blocks(iblk), blocks(iblk)%state(old), blocks(iblk)%state(new), dt)
      end if
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

  subroutine euler_debug(dt, block, old, new, pass)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass

    call space_operators(block, block%state(old), block%tend(old), dt, pass)
    call update_state(dt, block, block%tend(old), block%state(old), block%state(new), pass)

  end subroutine euler_debug

  subroutine predict_correct(dt, block, old, new, pass)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass

    ! Do first predict step.
    call space_operators(block, block%state(old), block%tend(old), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, block, block%tend(old), block%state(old), block%state(new), pass)

    ! Do second predict step.
    call space_operators(block, block%state(new), block%tend(old), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, block, block%tend(old), block%state(old), block%state(new), pass)

    ! Do correct step.
    call space_operators(block, block%state(new), block%tend(new),          dt, pass)
    call update_state(         dt, block, block%tend(new), block%state(old), block%state(new), pass)

  end subroutine predict_correct

  subroutine predict_correct_with_forward_backward(dt, block, old, new, pass)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass

    call space_operators(block, block%state(old), block%tend(old), 0.5_r8 * dt, pass + forward_pass)
    call update_state(0.5_r8 * dt, block, block%tend(old), block%state(old), block%state(new), pass) ! Update phs, pt, du, dv
    call space_operators(block, block%state(new), block%tend(old), 0.5_r8 * dt, pass + backward_pass)
    call update_state(0.5_r8 * dt, block, block%tend(old), block%state(old), block%state(new), pass) ! Update u, v

    call space_operators(block, block%state(new), block%tend(old), 0.5_r8 * dt, pass + forward_pass)
    call update_state(0.5_r8 * dt, block, block%tend(old), block%state(old), block%state(new), pass) ! Update phs, pt, du, dv
    call space_operators(block, block%state(new), block%tend(old), 0.5_r8 * dt, pass + backward_pass)
    call update_state(0.5_r8 * dt, block, block%tend(old), block%state(old), block%state(new), pass) ! Update u, v

    call space_operators(block, block%state(new), block%tend(new),          dt, pass + forward_pass)
    call update_state(         dt, block, block%tend(new), block%state(old), block%state(new), pass) ! Update phs, pt, du, dv
    call space_operators(block, block%state(new), block%tend(new),          dt, pass + backward_pass)
    call update_state(         dt, block, block%tend(new), block%state(old), block%state(new), pass) ! Update u, v

  end subroutine predict_correct_with_forward_backward

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
    call update_state(0.5_r8 * dt, block, block%tend(s1), block%state(old), block%state(s1), pass)

    call space_operators(block, block%state(s1) , block%tend(s2), 2.0_r8 * dt, pass)
    call update_state(        -dt, block, block%tend(s1), block%state(old), block%state(s2), pass)
    call update_state(2.0_r8 * dt, block, block%tend(s2), block%state(s2) , block%state(s2), pass)

    call space_operators(block, block%state(s2) , block%tend(s3),          dt, pass)
    block%tend(old)%du = (block%tend(s1)%du + 4.0_r8 * block%tend(s2)%du + block%tend(s3)%du) / 6.0_r8
    block%tend(old)%dv = (block%tend(s1)%dv + 4.0_r8 * block%tend(s2)%dv + block%tend(s3)%dv) / 6.0_r8
    if (baroclinic) then
      block%tend(old)%dpt  = (block%tend(s1)%dpt  + 4.0_r8 * block%tend(s2)%dpt  + block%tend(s3)%dpt ) / 6.0_r8
      block%tend(old)%dphs = (block%tend(s1)%dphs + 4.0_r8 * block%tend(s2)%dphs + block%tend(s3)%dphs) / 6.0_r8
    else
      block%tend(old)%dgz = (block%tend(s1)%dgz + 4.0_r8 * block%tend(s2)%dgz + block%tend(s3)%dgz) / 6.0_r8
    end if
    call update_state(         dt, block, block%tend(old), block%state(old), block%state(new), pass)

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
    call update_state(0.5_r8 * dt, block, block%tend(s1), block%state(old), block%state(s1), pass)

    call space_operators(block, block%state(s1) , block%tend(s2), 0.5_r8 * dt, pass)
    call update_state(0.5_r8 * dt, block, block%tend(s2), block%state(old), block%state(s2), pass)

    call space_operators(block, block%state(s2) , block%tend(s3),          dt, pass)
    call update_state(         dt, block, block%tend(s3), block%state(old), block%state(s3), pass)

    call space_operators(block, block%state(s3) , block%tend(s4),          dt, pass)
    block%tend(old)%du = (block%tend(s1)%du + 2.0_r8 * block%tend(s2)%du + 2.0_r8 * block%tend(s3)%du + block%tend(s4)%du) / 6.0_r8
    block%tend(old)%dv = (block%tend(s1)%dv + 2.0_r8 * block%tend(s2)%dv + 2.0_r8 * block%tend(s3)%dv + block%tend(s4)%dv) / 6.0_r8
    if (baroclinic) then
      block%tend(old)%dpt  = (block%tend(s1)%dpt  + 2.0_r8 * block%tend(s2)%dpt  + 2.0_r8 * block%tend(s3)%dpt  + block%tend(s4)%dpt ) / 6.0_r8
      block%tend(old)%dphs = (block%tend(s1)%dphs + 2.0_r8 * block%tend(s2)%dphs + 2.0_r8 * block%tend(s3)%dphs + block%tend(s4)%dphs) / 6.0_r8
    else
      block%tend(old)%dgz = (block%tend(s1)%dgz + 2.0_r8 * block%tend(s2)%dgz + 2.0_r8 * block%tend(s3)%dgz + block%tend(s4)%dgz) / 6.0_r8
    end if
    call update_state(         dt, block, block%tend(old), block%state(old), block%state(new), pass)

  end subroutine runge_kutta_4th

  subroutine update_state(dt, block, tend, old_state, new_state, pass)

    real(r8), intent(in) :: dt
    type(block_type), intent(inout) :: block
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state
    integer, intent(in) :: pass

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => old_state%mesh

    if (baroclinic) then
      if (tend%updated_dphs) then
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%phs(i,j) = old_state%phs(i,j) + dt * tend%dphs(i,j)
          end do
        end do
        call fill_halo(block, new_state%phs, full_lon=.true., full_lat=.true.)

        call calc_ph_lev_ph(block, new_state)
        call calc_m        (block, new_state)
      else if (tend%copy_phs) then
        new_state%phs    = old_state%phs
        new_state%ph_lev = old_state%ph_lev
        new_state%ph     = old_state%ph
        new_state%m      = old_state%m
      end if

      if (tend%updated_dpt) then
        if (.not. tend%updated_dphs .and. .not. tend%copy_phs .and. is_root_proc()) call log_error('Mass is not updated or copied!')
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
      if (tend%updated_dgz) then
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

    if (tend%updated_du) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            new_state%u(i,j,k) = old_state%u(i,j,k) + dt * tend%du(i,j,k)
          end do
        end do
      end do
      call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)
    end if

    if (tend%updated_dv) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            new_state%v(i,j,k) = old_state%v(i,j,k) + dt * tend%dv(i,j,k)
          end do
        end do
      end do
      call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end if

  end subroutine update_state

  subroutine damp_state(block, state)

    type(block_type), intent(inout) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k, damp_order

    mesh => state%mesh

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        damp_order = block%reduced_mesh(j)%damp_order
        if (damp_order > 0) then
          if (damp_order == 1) cycle ! User can choose not to damp except for cases when potential enstrophy increases.
          call zonal_damp(block, damp_order, dt_in_seconds, mesh%de_lon(j),      &
                          mesh%half_lon_lb, mesh%half_lon_ub, mesh%num_half_lon, &
                          state%u(:,j,k))
        end if
      end do

      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
#ifdef V_POLE
        damp_order = max(block%reduced_mesh(j)%damp_order, block%reduced_mesh(j-1)%damp_order)
#else
        damp_order = max(block%reduced_mesh(j)%damp_order, block%reduced_mesh(j+1)%damp_order)
#endif
        if (damp_order > 0) then
          if (damp_order == 1) cycle ! User can choose not to damp except for cases when potential enstrophy increases.
          call zonal_damp(block, damp_order, dt_in_seconds, mesh%le_lat(j),      &
                          mesh%full_lon_lb, mesh%full_lon_ub, mesh%num_full_lon, &
                          state%v(:,j,k))
        end if
      end do
    end do

  end subroutine damp_state

end module gmcore_mod
