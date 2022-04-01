module adv_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use block_mod
  use process_mod
  use adv_batch_mod
  use ffsl_mod

  implicit none

  private

  public adv_init
  public adv_final
  public adv_add_tracer
  public adv_allocate_tracers
  public adv_accum_wind
  public adv_calc_mass_flux
  public adv_calc_tracer_flux

  interface
    subroutine calc_flux_interface(block, batch, m, mfx, mfy)
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
    end subroutine calc_flux_interface
  end interface

  procedure(calc_flux_interface), pointer :: adv_calc_mass_flux   => null()
  procedure(calc_flux_interface), pointer :: adv_calc_tracer_flux => null()

  integer ntracer
  character(30), allocatable :: batch_names(:)
  character(30), allocatable :: tracer_names(:)
  character(30), allocatable :: tracer_long_names(:)
  character(30), allocatable :: tracer_units(:)
  real(r8), allocatable :: tracer_dt(:)

contains

  subroutine adv_init()

    call adv_final()

    select case (adv_scheme)
    case ('ffsl')
      adv_calc_mass_flux   => ffsl_calc_mass_flux
      adv_calc_tracer_flux => ffsl_calc_tracer_flux
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

    integer nbatch, nbatch_tracer, i, j, k, l, is, ie, js, je, ks, ke
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
    call log_notice('There are ' // to_str(nbatch) // ' advection batches.')
    if (is_root_proc()) then
      do i = 1, nbatch
        write(*, *) '- ', trim(unique_batch_names(i)), real(unique_tracer_dt(i))
      end do
    end if

    ! Initialize advection batches.
    is = block%mesh%full_lon_lb; ie = block%mesh%full_lon_ub
    js = block%mesh%full_lat_lb; je = block%mesh%full_lat_ub
    ks = block%mesh%full_lev_lb; ke = block%mesh%full_lev_ub
    do i = 1, size(block%state)
      allocate(block%state(i)%adv_batches(nbatch))
      do j = 1, nbatch
        call block%state(i)%adv_batches(j)%init(block%mesh, unique_batch_names(j), unique_tracer_dt(j))
      end do
      allocate(block%state(i)%q      (is:ie,js:je,ks:ke,ntracer))
      allocate(block%state(i)%qmf_lon(is:ie,js:je,ks:ke,ntracer))
      allocate(block%state(i)%qmf_lat(is:ie,js:je,ks:ke,ntracer))
      allocate(block%state(i)%qmf_lev(is:ie,js:je,ks:ke,ntracer))
    end do

    do i = 1, size(block%state)
      do j = 1, nbatch
        nbatch_tracer = 0
        do k = 1, ntracer
          if (batch_names(k) == block%state(i)%adv_batches(j)%alert_key) then
            nbatch_tracer = nbatch_tracer + 1
          end if
        end do
        call block%state(i)%adv_batches(j)%allocate_tracers(nbatch_tracer)
        l = 0
        do k = 1, ntracer
          if (batch_names(k) == block%state(i)%adv_batches(j)%alert_key) then
            l = l + 1
            block%state(i)%adv_batches(j)%tracer_idx       (l) =                   k
            block%state(i)%adv_batches(j)%tracer_names     (l) = tracer_names     (k)
            block%state(i)%adv_batches(j)%tracer_long_names(l) = tracer_long_names(k)
            block%state(i)%adv_batches(j)%tracer_units     (l) = tracer_units     (k)
          end if
        end do
      end do
    end do

  end subroutine adv_allocate_tracers

  subroutine adv_accum_wind(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    integer i

    associate (state => block%state(itime))
    do i = 1, size(state%adv_batches)
      call state%adv_batches(i)%accum_wind(state%u, state%v)
    end do
    end associate

  end subroutine adv_accum_wind

  subroutine adv_final()

    if (allocated(batch_names      )) deallocate(batch_names      )
    if (allocated(tracer_names     )) deallocate(tracer_names     )
    if (allocated(tracer_long_names)) deallocate(tracer_long_names)
    if (allocated(tracer_units     )) deallocate(tracer_units     )
    if (allocated(tracer_dt        )) deallocate(tracer_dt        )

  end subroutine adv_final

end module adv_mod
