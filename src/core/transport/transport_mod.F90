module transport_mod

  use flogger
  use string
  use const_mod
  use block_mod
  use process_mod

  implicit none

  integer ntracer
  character(30), allocatable :: batch_names(:)
  character(30), allocatable :: tracer_names(:)
  real(r8), allocatable :: dt_transport(:)

contains

  subroutine transport_init()

    call transport_final()

    ntracer = 0
    allocate(batch_names (1000))
    allocate(tracer_names(1000))
    allocate(dt_transport(1000))

  end subroutine transport_init

  subroutine transport_add_tracer(batch_name, tracer_name, dt)

    character(*), intent(in) :: batch_name
    character(*), intent(in) :: tracer_name
    real(r8), intent(in) :: dt

    ntracer = ntracer + 1
    if (ntracer > size(batch_names)) then
      call log_error('Insufficient character array!', __FILE__, __LINE__, pid=proc%id)
    end if
    batch_names (ntracer) = batch_name
    tracer_names(ntracer) = tracer_name
    dt_transport(ntracer) = dt

  end subroutine transport_add_tracer

  subroutine transport_allocate_tracers(block)

    type(block_type), intent(inout) :: block

    integer nbatch, i, j, k
    logical found
    character(30) unique_batch_names(100)
    real(r8) unique_dt_transport(100)

    nbatch = 0
    do i = 1, ntracer
      found = .false.
      do j = 1, i - 1
        if (i /= j .and. batch_names(i) == batch_names(j)) then
          found = .true.
          exit
        end if
      end do
      print *, i, j, found, nbatch
      if (.not. found) then
        nbatch = nbatch + 1
        unique_batch_names(nbatch ) = batch_names (i)
        unique_dt_transport(nbatch) = dt_transport(i)
      end if
    end do
    call log_notice('There are ' // to_str(nbatch) // ' transport batches.')
    if (is_root_proc()) then
      do i = 1, nbatch
        write(*, *) '- ', trim(unique_batch_names(i)), real(unique_dt_transport(i))
      end do
    end if

    do i = 1, size(block%state)
      allocate(block%state(i)%transport_batches(nbatch))
      do j = 1, nbatch
        call block%state(i)%transport_batches(j)%init(block%mesh, unique_batch_names(j), unique_dt_transport(j))
      end do
    end do

    do i = 1, size(block%state)
      do j = 1, nbatch
        do k = 1, ntracer
          if (batch_names(k) == block%state(i)%transport_batches(j)%alert_key) then
            call block%state(i)%transport_batches(j)%add_tracer(tracer_names(k))
          end if
        end do
      end do
    end do

  end subroutine transport_allocate_tracers

  subroutine transport_final()

    if (allocated(batch_names )) deallocate(batch_names )
    if (allocated(tracer_names)) deallocate(tracer_names)
    if (allocated(dt_transport)) deallocate(dt_transport)

  end subroutine transport_final

end module transport_mod