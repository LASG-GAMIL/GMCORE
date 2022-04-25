module adv_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
  use adv_batch_mod
  use ffsl_mod
  use upwind_mod
  use weno_mod
  use tvd_mod

  implicit none

  private

  public adv_init
  public adv_final
  public adv_add_tracer
  public adv_allocate_tracers
  public adv_accum_wind
  public adv_calc_mass_hflx_cell
  public adv_calc_mass_hflx_vtx
  public adv_calc_tracer_hflx_cell
  public adv_calc_tracer_hflx_vtx
  public adv_calc_tracer_hval_vtx

  public upwind1
  public upwind3
  public weno3
  public weno5
  public tvd

  interface
    subroutine calc_hflx_cell_interface(block, batch, m, mfx, mfy, dt)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: m  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(8), intent(in), optional :: dt
    end subroutine calc_hflx_cell_interface
    subroutine calc_hflx_vtx_interface(block, batch, m, mfx, mfy, dt)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: m  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(8), intent(in), optional :: dt
    end subroutine calc_hflx_vtx_interface
    subroutine calc_hval_vtx_interface(block, batch, m, mvx, mvy, dt)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: m  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mvx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mvy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(8), intent(in), optional :: dt
    end subroutine calc_hval_vtx_interface
    subroutine calc_vflx_cell_interface(block, batch, m, mfz, dt)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: m  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfz(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%half_lev_lb:block%mesh%half_lev_ub)
      real(8), intent(in), optional :: dt
    end subroutine calc_vflx_cell_interface
  end interface

  procedure(calc_hflx_cell_interface), pointer :: adv_calc_mass_hflx_cell   => null()
  procedure(calc_hflx_vtx_interface ), pointer :: adv_calc_mass_hflx_vtx    => null()
  procedure(calc_hflx_cell_interface), pointer :: adv_calc_tracer_hflx_cell => null()
  procedure(calc_hflx_vtx_interface ), pointer :: adv_calc_tracer_hflx_vtx  => null()
  procedure(calc_hval_vtx_interface ), pointer :: adv_calc_tracer_hval_vtx  => null()

  integer ntracer
  character(30), allocatable :: batch_names(:)
  character(30), allocatable :: tracer_names(:)
  character(30), allocatable :: tracer_long_names(:)
  character(30), allocatable :: tracer_units(:)
  real(8), allocatable :: tracer_dt(:)

contains

  subroutine adv_init()

    call adv_final()

    select case (adv_scheme)
    case ('ffsl')
      call ffsl_init()
      adv_calc_mass_hflx_cell   => ffsl_calc_mass_hflx_cell
      adv_calc_mass_hflx_vtx    => ffsl_calc_mass_hflx_vtx
      adv_calc_tracer_hflx_cell => ffsl_calc_tracer_hflx_cell
      adv_calc_tracer_hflx_vtx  => ffsl_calc_tracer_hflx_vtx
      adv_calc_tracer_hval_vtx  => ffsl_calc_tracer_hval_vtx
    end select

    ntracer = 0
    allocate(batch_names      (1000))
    allocate(tracer_names     (1000))
    allocate(tracer_long_names(1000))
    allocate(tracer_units     (1000))
    allocate(tracer_dt        (1000))

  end subroutine adv_init

  subroutine adv_add_tracer(batch_name, dt, name, long_name, units)

    character(*), intent(in) :: batch_name
    real(r8), intent(in) :: dt
    character(*), intent(in) :: name
    character(*), intent(in), optional :: long_name
    character(*), intent(in), optional :: units

    ntracer = ntracer + 1
    if (ntracer > size(batch_names)) then
      call log_error('Insufficient character array!', __FILE__, __LINE__, pid=proc%id)
    end if
    batch_names (ntracer) = batch_name
    tracer_names(ntracer) = name
    tracer_dt   (ntracer) = dt
    if (present(long_name)) tracer_long_names(ntracer) = long_name
    if (present(units    )) tracer_units     (ntracer) = units

  end subroutine adv_add_tracer

  subroutine adv_allocate_tracers(block)

    type(block_type), intent(inout) :: block

    integer nbatch, nbatch_tracer, i, j, k, is, ie, js, je, ks, ke
    logical found
    character(30) unique_batch_names(100)
    real(r8) unique_tracer_dt(100)

    nbatch = 0
    do i = 1, ntracer
      found = .false.
      do j = 1, i - 1
        if (i /= j .and. batch_names(i) == batch_names(j)) then
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        nbatch = nbatch + 1
        unique_batch_names(nbatch) = batch_names(i)
        unique_tracer_dt  (nbatch) = tracer_dt  (i)
      end if
    end do
    if (is_root_proc()) then
      call log_notice('There are ' // to_str(nbatch) // ' advection batches.')
      do i = 1, nbatch
        write(*, *) '- ', trim(unique_batch_names(i)), real(unique_tracer_dt(i))
      end do
    end if

    ! Initialize advection batches in block objects and allocate tracer arrays in state objects.
    if (advection) then
      allocate(block%adv_batches(nbatch))
      do i = 1, nbatch
        call block%adv_batches(i)%init(block%mesh, 'cell', unique_batch_names(i), unique_tracer_dt(i))
      end do
    else
      call block%adv_batch_mass%init(block%mesh, 'cell', 'mass', dt_dyn)
      call block%adv_batch_pv  %init(block%mesh, 'vtx' , 'pv'  , dt_dyn)
      if (nbatch > 0) then
        allocate(block%adv_batches(nbatch))
        do i = 1, nbatch
          call block%adv_batches(i)%init(block%mesh, 'cell', unique_batch_names(i), unique_tracer_dt(i))
        end do
      end if
    end if
    if (ntracer > 0) then
      associate (mesh => block%mesh)
      do i = 1, size(block%state)
        allocate(block%state(i)%q      (mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub,ntracer))
        allocate(block%state(i)%qmf_lon(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub,ntracer))
        allocate(block%state(i)%qmf_lat(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub,ntracer))
        allocate(block%state(i)%qmf_lev(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub,ntracer))
      end do
      end associate
    end if

    ! Record tracer information in advection batches.
    do i = 1, nbatch
      nbatch_tracer = 0
      do j = 1, ntracer
        if (batch_names(j) == block%adv_batches(i)%alert_key) then
          nbatch_tracer = nbatch_tracer + 1
        end if
      end do
      call block%adv_batches(i)%allocate_tracers(nbatch_tracer)
      k = 0
      do j = 1, ntracer
        if (batch_names(j) == block%adv_batches(i)%alert_key) then
          k = k + 1
          block%adv_batches(i)%tracer_idx       (k) =                   j
          block%adv_batches(i)%tracer_names     (k) = tracer_names     (j)
          block%adv_batches(i)%tracer_long_names(k) = tracer_long_names(j)
          block%adv_batches(i)%tracer_units     (k) = tracer_units     (j)
        end if
      end do
    end do

  end subroutine adv_allocate_tracers

  subroutine adv_accum_wind(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    real(r8) work(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend,block%mesh%num_full_lev)
    real(r8) pole(block%mesh%num_full_lev)
    integer i, j, k, l

    if (allocated(block%adv_batches)) then
      do l = 1, size(block%adv_batches)
        select case (block%adv_batches(l)%loc)
        case ('cell')
          call block%adv_batches(l)%accum_uv_cell( &
            block%state(itime)%u_lon             , &
            block%state(itime)%v_lat             )
          call block%adv_batches(l)%accum_mf_cell( &
            block%state(itime)%mf_lon_n          , &
            block%state(itime)%mf_lat_n          )
        case ('vtx')
          call block%adv_batches(l)%accum_uv_vtx ( &
            block%state(itime)%u_lat             , &
            block%state(itime)%v_lon             )
          call block%adv_batches(l)%accum_mf_vtx ( &
            block%state(itime)%mf_lat_t          , &
            block%state(itime)%mf_lon_t          )
        end select
      end do
    end if

  end subroutine adv_accum_wind

  subroutine adv_final()

    ntracer = 0
    if (allocated(batch_names      )) deallocate(batch_names      )
    if (allocated(tracer_names     )) deallocate(tracer_names     )
    if (allocated(tracer_long_names)) deallocate(tracer_long_names)
    if (allocated(tracer_units     )) deallocate(tracer_units     )
    if (allocated(tracer_dt        )) deallocate(tracer_dt        )

  end subroutine adv_final

end module adv_mod
