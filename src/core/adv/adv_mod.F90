module adv_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
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
      call ffsl_init()
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
    call log_notice('There are ' // to_str(nbatch) // ' advection batches.')
    if (is_root_proc()) then
      do i = 1, nbatch
        write(*, *) '- ', trim(unique_batch_names(i)), real(unique_tracer_dt(i))
      end do
    end if

    ! Initialize advection batches in block objects and allocate tracer arrays in state objects.
    is = block%mesh%full_lon_lb; ie = block%mesh%full_lon_ub
    js = block%mesh%full_lat_lb; je = block%mesh%full_lat_ub
    ks = block%mesh%full_lev_lb; ke = block%mesh%full_lev_ub
    allocate(block%adv_batches(nbatch))
    do i = 1, nbatch
      call block%adv_batches(i)%init(block%mesh, unique_batch_names(i), unique_tracer_dt(i))
    end do
    do i = 1, size(block%state)
      allocate(block%state(i)%q      (is:ie,js:je,ks:ke,ntracer))
      allocate(block%state(i)%qmf_lon(is:ie,js:je,ks:ke,ntracer))
      allocate(block%state(i)%qmf_lat(is:ie,js:je,ks:ke,ntracer))
      allocate(block%state(i)%qmf_lev(is:ie,js:je,ks:ke,ntracer))
    end do

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

    associate (mesh => block%mesh, u => block%state(itime)%u, v => block%state(itime)%v)
    do l = 1, size(block%adv_batches)
      associate (batch => block%adv_batches(l))
      if (batch%step == 0) then
        batch%u = u
        batch%v = v
        batch%step = batch%step + 1
      else if (batch%step == batch%nstep) then
        batch%u      = (batch%u + u) / (batch%nstep + 1)
        batch%v      = (batch%v + v) / (batch%nstep + 1)
        batch%step   = 0
        ! Calculate CFL numbers and divergence along each axis.
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              batch%cflx(i,j,k) = batch%dt * batch%u(i,j,k) / mesh%de_lon(j)
            end do
          end do
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              batch%cfly(i,j,k) = batch%dt * batch%v(i,j,k) / mesh%de_lat(j)
            end do
          end do
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              batch%divx(i,j,k) = (batch%u(i,j,k) - batch%u(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
              batch%divy(i,j,k) = (batch%v(i,j  ,k) * mesh%le_lat(j  ) - &
                                   batch%v(i,j-1,k) * mesh%le_lat(j-1)) / mesh%area_cell(j)
            end do
          end do
        end do
        if (mesh%has_south_pole()) then
          j = mesh%full_lat_ibeg
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              work(i,k) = batch%v(i,j,k)
            end do
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              batch%divy(i,j,k) = pole(k)
            end do
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%full_lat_iend
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              work(i,k) = -batch%v(i,j-1,k)
            end do
          end do
          call zonal_sum(proc%zonal_circle, work, pole)
          pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              batch%divy(i,j,k) = pole(k)
            end do
          end do
        end if
      else
        batch%u = batch%u + u
        batch%v = batch%v + v
        batch%step = batch%step + 1
      end if
      end associate
    end do
    end associate

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
