module adv_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use block_mod
  use process_mod
  use adv_batch_mod
  use tracer_mod
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
  public tracer_type

  interface
    subroutine calc_flux_interface(block, state, batch, tracer)
      import block_type, state_type, adv_batch_type, tracer_type
      type(block_type    ), intent(in   ) :: block
      type(state_type    ), intent(in   ) :: state
      type(adv_batch_type), intent(in   ) :: batch
      type(tracer_type   ), intent(inout) :: tracer
    end subroutine calc_flux_interface
  end interface

  procedure(calc_flux_interface), pointer :: adv_calc_mass_flux => null()
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

    integer nbatch, i, j, k
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

    do i = 1, size(block%state)
      allocate(block%state(i)%adv_batches(nbatch))
      do j = 1, nbatch
        call block%state(i)%adv_batches(j)%init(block%mesh, unique_batch_names(j), unique_tracer_dt(j))
      end do
    end do

    do i = 1, size(block%state)
      do j = 1, nbatch
        do k = 1, ntracer
          if (batch_names(k) == block%state(i)%adv_batches(j)%alert_key) then
            call block%state(i)%adv_batches(j)%add_tracer(tracer_names(k), tracer_long_names(k), tracer_units(k))
          end if
        end do
        call block%state(i)%adv_batches(j)%allocate_tracers()
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