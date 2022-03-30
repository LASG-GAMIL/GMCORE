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
  use debug_mod
  use adv_mod
  use pgf_mod
  use damp_mod
  use diag_state_mod
  use test_forcing_mod
  use filter_mod

  implicit none

  private

  public gmcore_init
  public gmcore_run
  public gmcore_final

  public block_type
  public proc

  procedure(space_operators_interface), pointer :: operators

contains

  subroutine gmcore_init(namelist_path, comm)

    character(*), intent(in) :: namelist_path
    integer, intent(in), optional :: comm

    character(10) time_value, time_units
    real(8) seconds

    call log_init()
    call global_mesh%init_global(num_lon, num_lat, num_lev, &
                                 lon_halo_width=2, &
                                 lat_halo_width=2)
    call process_init(comm)
    call vert_coord_init(num_lev, namelist_path)
    call process_create_blocks()
    call time_init(dt_dyn)
    call diag_state_init(proc%blocks)
    call history_init()
    call restart_init()
    call time_scheme_init()
    call adv_init()
    call pgf_init()
    call interp_init()
    call damp_init(proc%blocks)

    operators => space_operators

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
    if (filter_reset_interval > 0) call time_add_alert('filter_reset', seconds=filter_reset_interval)

    if (is_root_proc()) call print_namelist()

  end subroutine gmcore_init

  subroutine gmcore_run()

    integer j, iblk, itime

    do iblk = 1, size(proc%blocks)
      associate (block => proc%blocks(iblk)     , &
                 mesh  => proc%blocks(iblk)%mesh, &
                 state => proc%blocks(iblk)%state(old))
      if (baroclinic) then 
        call prepare_static(block)
        ! Ensure bottom gz_lev is the same as gzs.
        do itime = lbound(block%state, 1), ubound(block%state, 1)
          block%state(itime)%gz_lev(:,:,global_mesh%half_lev_iend) = block%static%gzs
        end do
      end if
        associate (state => block%state(old))
        call filter_on_lon_edge(block, state%u, state%u_f)
        call fill_halo(block, state%u_f, full_lon=.false., full_lat=.true., full_lev=.true.)
        call filter_on_lat_edge(block, state%v, state%v_f)
        call fill_halo(block, state%v_f, full_lon=.true., full_lat=.false., full_lev=.true.)
        if (baroclinic) then
          call filter_on_cell(block, state%phs, state%phs_f)
          call fill_halo(block, state%phs_f, full_lon=.true., full_lat=.true.)
          call filter_on_cell(block, state%pt, state%pt_f)
          call fill_halo(block, state%pt_f, full_lon=.true., full_lat=.true., full_lev=.true.)
        else
          state%gz_f = state%gz
        end if
        end associate
      end associate
    end do

    call operators_prepare(proc%blocks, old, dt_dyn)
    if (nonhydrostatic) call nh_prepare(proc%blocks)
    call diagnose(proc%blocks, old)
    call output(old)

    do while (.not. time_is_finished())
      call time_integrate(dt_dyn, proc%blocks)
      if (is_root_proc() .and. time_is_alerted('print')) call log_print_diag(curr_time%isoformat())
      call time_advance(dt_dyn)
      call operators_prepare(proc%blocks, old, dt_dyn)
      call diagnose(proc%blocks, old)
      call output(old)
    end do

  end subroutine gmcore_run

  subroutine gmcore_final()

    call log_final()
    call time_final()
    call interp_final()
    call adv_final()
    call damp_final()
    call diag_state_final()
    call history_final()
    call process_final()

  end subroutine gmcore_final

  subroutine prepare_static(block)

    class(block_type), intent(inout) :: block

    integer i, j

    associate (mesh    => block%mesh          , &
               gzs     => block%static%gzs    , & ! in
               dzsdlon => block%static%dzsdlon, & ! out
               dzsdlat => block%static%dzsdlat)   ! out
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          dzsdlon(i,j) = (gzs(i+1,j) - gzs(i,j)) / g / mesh%de_lon(j)
        end do
      end do
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          dzsdlat(i,j) = (gzs(i,j+1) - gzs(i,j)) / g / mesh%de_lat(j)
        end do
      end do
      call fill_halo(block, dzsdlon, full_lon=.false., full_lat=.true.)
      call fill_halo(block, dzsdlat, full_lon=.true., full_lat=.false.)
    end associate

  end subroutine prepare_static

  subroutine output(itime)

    integer, intent(in) :: itime

    real(8), save :: time1 = 0, time2
    integer i, j, k, iblk

    if (time_is_alerted('history_write')) then
      if (time_step == 0) call cpu_time(time1)
      call cpu_time(time2)
      if (time_step /= 0) then
        if (is_root_proc()) call log_notice('Time cost ' // to_str(time2 - time1, 5) // ' seconds.')
        time1 = time2
      end if
      if (output_h0) call history_write_state(proc%blocks, itime)
      if (output_h1) call history_write_debug(proc%blocks, itime)
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
    real(r8) tm, te, tav, tpe, tpt, max_w

    tm = 0.0_r8
    te = 0.0_r8
    tav = 0.0_r8
    tpe = 0.0_r8
    tpt = 0.0_r8
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
      if (hydrostatic) then
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
      else if (nonhydrostatic) then
      else
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

      if (baroclinic) then
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tpt = tpt + state%m(i,j,k) * state%pt(i,j,k) * mesh%area_cell(j)
            end do
          end do
        end do
      end if
    end do
    call global_sum(proc%comm, tm)
    call global_sum(proc%comm, te)
    call global_sum(proc%comm, tav)
    call global_sum(proc%comm, tpe)
    if (baroclinic) call global_sum(proc%comm, tpt)

    do iblk = 1, size(blocks)
      blocks(iblk)%state(itime)%tm  = tm
      blocks(iblk)%state(itime)%te  = te
      blocks(iblk)%state(itime)%tav = tav
      blocks(iblk)%state(itime)%tpe = tpe
      if (diag_state(iblk)%is_init()) call diag_state(iblk)%run(proc%blocks(iblk)%state(itime))
    end do

    call log_add_diag('tm' , tm )
    if (baroclinic) call log_add_diag('tpt', tpt)
    call log_add_diag('te' , te )
    call log_add_diag('tpe', tpe)

    if (nonhydrostatic) then
      max_w = 0
      do iblk = 1, size(blocks)
        max_w = max(max_w, maxval(abs(blocks(iblk)%state(itime)%w)))
      end do
      call global_max(proc%comm, max_w)
      call log_add_diag('w', max_w)
    end if

  end subroutine diagnose

  subroutine space_operators(block, old_state, star_state, new_state, tend1, tend2, dt, pass)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in   ) :: old_state
    type(state_type), intent(inout) :: star_state
    type(state_type), intent(inout) :: new_state
    type(tend_type ), intent(inout) :: tend1
    type(tend_type ), intent(in   ) :: tend2
    real(8), intent(in) :: dt
    integer, intent(in) :: pass

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => star_state%mesh

    call tend1%reset_flags()

    select case (pass)
    case (all_pass)
      call operators_prepare(block, star_state, dt, pass)
      if (hydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, star_state, tend1, dt)
        call calc_dphs             (block, star_state, tend1, dt)
        call calc_wedphdlev_lev    (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call calc_dptfdlon_dptfdlat(block, star_state, tend1, dt)
        call calc_dptfdlev         (block, star_state, tend1, dt)
        call calc_qhu_qhv          (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat  (block, star_state, tend1, dt)
        call pgf_run               (block, star_state, tend1)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%pgf_lon(i,j,k) - tend1%dkedlon(i,j,k) - tend1%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%pgf_lat(i,j,k) - tend1%dkedlat(i,j,k) - tend1%wedvdlev(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dpt(i,j,k) = - tend1%dptfdlon(i,j,k) - tend1%dptfdlat(i,j,k) - tend1%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
        tend1%update_phs = .true.
        tend1%update_pt  = .true.
      else if (nonhydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, star_state, tend1, dt)
        call calc_dphs             (block, star_state, tend1, dt)
        call calc_wedphdlev_lev    (block, star_state, tend1, dt)
        call calc_dptfdlon_dptfdlat(block, star_state, tend1, dt)
        call calc_dptfdlev         (block, star_state, tend1, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dpt(i,j,k) = - tend1%dptfdlon(i,j,k) - tend1%dptfdlat(i,j,k) - tend1%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend1%update_phs = .true.
        tend1%update_pt  = .true.
        call update_state(block, tend1, old_state, new_state, dt)

        call nh_solve(block, tend1, old_state, star_state, new_state, dt)

        call calc_qhu_qhv          (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat  (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call pgf_run               (block, new_state, tend1)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%pgf_lon(i,j,k) - tend1%dkedlon(i,j,k) - tend1%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%pgf_lat(i,j,k) - tend1%dkedlat(i,j,k) - tend1%wedvdlev(i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
      else
        call calc_dmfdlon_dmfdlat(block, star_state, tend1, dt)
        call calc_qhu_qhv        (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat(block, star_state, tend1, dt)
        call pgf_run             (block, star_state, tend1    )

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%pgf_lon(i,j,k) - tend1%dkedlon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%pgf_lat(i,j,k) - tend1%dkedlat(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dgz(i,j,k) = - (tend1%dmfdlon(i,j,k) + tend1%dmfdlat(i,j,k)) * g
            end do
          end do
        end do

        tend1%update_u  = .true.
        tend1%update_v  = .true.
        tend1%update_gz = .true.
      end if
    case (forward_pass)
      call operators_prepare(block, star_state, dt, pass)
      if (hydrostatic) then
        call calc_dmfdlon_dmfdlat  (block, star_state, tend1, dt)
        call calc_dphs             (block, star_state, tend1, dt)
        call calc_wedphdlev_lev    (block, star_state, tend1, dt)
        call calc_wedudlev_wedvdlev(block, star_state, tend1, dt)
        call calc_dptfdlon_dptfdlat(block, star_state, tend1, dt)
        call calc_dptfdlev         (block, star_state, tend1, dt)
        call calc_qhu_qhv          (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat  (block, star_state, tend1, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%dkedlon(i,j,k) - tend1%wedudlev(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%dkedlat(i,j,k) - tend1%wedvdlev(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dpt(i,j,k) = - tend1%dptfdlon(i,j,k) - tend1%dptfdlat(i,j,k) - tend1%dptfdlev(i,j,k)
            end do
          end do
        end do

        tend1%update_phs = .true.
        tend1%update_pt  = .true.
      else if (nonhydrostatic) then

      else
        call calc_dmfdlon_dmfdlat(block, star_state, tend1, dt)
        call calc_qhu_qhv        (block, star_state, tend1, dt)
        call calc_dkedlon_dkedlat(block, star_state, tend1, dt)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(i,j,k) =   tend1%qhv(i,j,k) - tend1%dkedlon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(i,j,k) = - tend1%qhu(i,j,k) - tend1%dkedlat(i,j,k)
            end do
          end do

          do j = mesh%full_lat_ibeg, mesh%full_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dgz(i,j,k) = - (tend1%dmfdlon(i,j,k) + tend1%dmfdlat(i,j,k)) * g
            end do
          end do
        end do

        tend1%update_gz = .true.
      end if
    case (backward_pass)
      call operators_prepare(block, new_state, dt, pass)
      if (hydrostatic) then
        call pgf_run(block, new_state, tend1)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(i,j,k) = tend2%du(i,j,k) - tend1%pgf_lon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(i,j,k) = tend2%dv(i,j,k) - tend1%pgf_lat(i,j,k)
            end do
          end do
        end do

        tend1%update_u   = .true.
        tend1%update_v   = .true.
      else if (nonhydrostatic) then

      else
        call pgf_run(block, new_state, tend1)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend1%du(i,j,k) = tend2%du(i,j,k) - tend1%pgf_lon(i,j,k)
            end do
          end do

          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend1%dv(i,j,k) = tend2%dv(i,j,k) - tend1%pgf_lat(i,j,k)
            end do
          end do
        end do

        tend1%update_u  = .true.
        tend1%update_v  = .true.
      end if
    end select

  end subroutine space_operators

  subroutine time_integrate(dt, blocks)

    real(8), intent(in) :: dt
    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      call time_integrator(operators, blocks(iblk), old, new, dt)
      call test_forcing_run(blocks(iblk), dt, blocks(iblk)%static, blocks(iblk)%state(new))
    end do
    call damp_run(dt, new, blocks)

  end subroutine time_integrate

end module gmcore_mod
