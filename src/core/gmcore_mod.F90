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
  use tend_mod
  use block_mod
  use vert_coord_mod
  use time_schemes_mod
  use operators_mod
  use interp_mod
  use reduce_mod
  use debug_mod
  use pgf_mod
  use damp_mod
  use test_forcing_mod

  implicit none

  private

  public gmcore_init
  public gmcore_run
  public gmcore_final

  interface
    subroutine splitter_interface(dt, block)
      import block_type
      real(8), intent(in) :: dt
      type(block_type), intent(inout) :: block
    end subroutine splitter_interface
  end interface

  procedure(splitter_interface), pointer :: splitter

contains

  subroutine gmcore_init(namelist_path)

    character(*), intent(in) :: namelist_path

    character(10) time_value, time_units
    real(8) seconds

    call log_init()
    call global_mesh%init_global(num_lon, num_lat, num_lev, lon_halo_width=max(2, maxval(reduce_factors) - 1), lat_halo_width=2)
    !call debug_check_areas()
    call process_init()
    call vert_coord_init(num_lev, namelist_path)
    call process_create_blocks()
    call time_init()
    call history_init()
    call restart_init()
    call reduce_init(proc%blocks)
    call time_scheme_init()
    call pgf_init()
    call damp_init()

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

    if (is_root_proc()) call print_namelist()

  end subroutine gmcore_init

  subroutine gmcore_run()

    call operators_prepare(proc%blocks, old, dt_in_seconds)
    call diagnose(proc%blocks, old)
    call output(old)

    do while (.not. time_is_finished())
      call time_integrate(dt_in_seconds, proc%blocks)
      if (is_root_proc() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      call time_advance(dt_in_seconds)
      call operators_prepare(proc%blocks, old, dt_in_seconds)
      call diagnose(proc%blocks, old)
      call output(old)
    end do

  end subroutine gmcore_run

  subroutine gmcore_final()

    call damp_final()
    call history_final()
    call process_final()

  end subroutine gmcore_final

  subroutine output(itime)

    integer, intent(in) :: itime

    real(8), save :: time1 = 0, time2
    integer i, j, k, iblk

    if (time_is_alerted('history_write')) then
      if (time_step == 0) call cpu_time(time1)
      call cpu_time(time2)
      if (time_step /= 0) then
        if (is_root_proc()) call log_notice('Time cost ' // to_string(time2 - time1, 5) // ' seconds.')
        time1 = time2
      end if
      if (baroclinic) then
        ! Interpolate onto isobaric layers.
        do iblk = 1, size(proc%blocks)
          associate (state => proc%blocks(iblk)%state(itime))
            call interp_lon_edge_to_isobaric_level(state%mesh, state%ph, state%u, 85000.0_r8, state%u850, logp=.true.)
            call interp_lon_edge_to_isobaric_level(state%mesh, state%ph, state%u, 70000.0_r8, state%u700, logp=.true.)
            call interp_lat_edge_to_isobaric_level(state%mesh, state%ph, state%v, 85000.0_r8, state%v850, logp=.true.)
            call interp_lat_edge_to_isobaric_level(state%mesh, state%ph, state%v, 70000.0_r8, state%v700, logp=.true.)
            call interp_cell_to_isobaric_level(state%mesh, state%ph, state%t, 85000.0_r8, state%t850, logp=.true.)
            call interp_cell_to_isobaric_level(state%mesh, state%ph, state%t, 70000.0_r8, state%t700, logp=.true.)
          end associate
        end do
      end if
      call history_write_state(proc%blocks, itime)
      call history_write_debug(proc%blocks, itime)
    end if
    if (time_is_alerted('restart_write')) then
      call restart_write(itime)
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
            tav = tav + state%pv(i,j,k) * mesh%area_vtx(j)
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
      blocks(iblk)%state(itime)%tav = tav
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
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    type(mesh_type), pointer :: mesh
    integer i, j, k

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
        call pgf_run               (block, state, tend)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) =   tend%qhv(i,j,k) - tend%pgf_lon(i,j,k) - tend%dkedlon(i,j,k) - tend%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%qhu(i,j,k) - tend%pgf_lat(i,j,k) - tend%dkedlat(i,j,k) - tend%wedvdlev(i,j,k)
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
        call pgf_run             (block, state, tend    )

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) =   tend%qhv(i,j,k) - tend%pgf_lon(i,j,k) - tend%dkedlon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%qhu(i,j,k) - tend%pgf_lat(i,j,k) - tend%dkedlat(i,j,k)
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
        call calc_dkedlon_dkedlat  (block, state, tend, dt)
        call pgf_run               (block, state, tend    )

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) = - tend%wedudlev(i,j,k) - tend%dkedlon(i,j,k) - tend%pgf_lon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%wedvdlev(i,j,k) - tend%dkedlat(i,j,k) - tend%pgf_lat(i,j,k)
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
        call pgf_run             (block, state, tend    )
        call calc_dmfdlon_dmfdlat(block, state, tend, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%du(i,j,k) = - tend%pgf_lon(i,j,k) - tend%dkedlon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dv(i,j,k) = - tend%pgf_lat(i,j,k) - tend%dkedlat(i,j,k)
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

    real(8), intent(in) :: dt
    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      call splitter(dt, blocks(iblk))

      if (use_div_damp) then
        call div_damp(blocks(iblk), dt, blocks(iblk)%state(new))
      end if
      if (use_vor_damp) then
        call vor_damp(blocks(iblk), dt, blocks(iblk)%state(new))
      end if
      if (use_polar_damp) then
        call polar_damp(blocks(iblk), dt, blocks(iblk)%state(new))
      end if
      call test_forcing_run(blocks(iblk), dt, blocks(iblk)%state(new))
      if (use_div_damp .or. use_vor_damp .or. use_vor_damp) then
        call operators_prepare(blocks(iblk), blocks(iblk)%state(new), dt, all_pass)
      end if
    end do

  end subroutine time_integrate

  subroutine csp2_splitting(dt, block)

    real(8), intent(in) :: dt
    type(block_type), intent(inout) :: block

    real(r8) fast_dt
    integer subcycle, t1, t2

    fast_dt = dt / fast_cycles
    t1 = 3
    t2 = old

    call time_integrator(space_operators, 0.5_r8 * dt, block, old, t1, slow_pass)
    do subcycle = 1, fast_cycles
      call time_integrator(space_operators, fast_dt, block, t1, t2, fast_pass)
      call time_swap_indices(t1, t2)
    end do
    call time_integrator(space_operators, 0.5_r8 * dt, block, t1, new, slow_pass)

  end subroutine csp2_splitting

  subroutine no_splitting(dt, block)

    real(8), intent(in) :: dt
    type(block_type), intent(inout) :: block

    call time_integrator(space_operators, dt, block, old, new, all_pass)

  end subroutine no_splitting

end module gmcore_mod
