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
  use zonal_damp_mod

  implicit none

  private

  public adv_init
  public adv_run
  public adv_final
  public adv_add_tracer
  public adv_allocate_tracers
  public adv_accum_wind
  public adv_calc_mass_hflx_cell
  public adv_calc_mass_vflx_cell
  public adv_calc_mass_hflx_vtx
  public adv_calc_tracer_hflx_cell
  public adv_calc_tracer_vflx_cell
  public adv_calc_tracer_hflx_vtx
  public adv_calc_tracer_hval_vtx
  public adv_batch_type

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

  interface adv_allocate_tracers
    module procedure adv_allocate_tracers_1
    module procedure adv_allocate_tracers_2
  end interface adv_allocate_tracers

  procedure(calc_hflx_cell_interface), pointer :: adv_calc_mass_hflx_cell   => null()
  procedure(calc_vflx_cell_interface), pointer :: adv_calc_mass_vflx_cell   => null()
  procedure(calc_hflx_vtx_interface ), pointer :: adv_calc_mass_hflx_vtx    => null()
  procedure(calc_hflx_cell_interface), pointer :: adv_calc_tracer_hflx_cell => null()
  procedure(calc_vflx_cell_interface), pointer :: adv_calc_tracer_vflx_cell => null()
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
      adv_calc_mass_vflx_cell   => ffsl_calc_mass_vflx_cell
      adv_calc_mass_hflx_vtx    => ffsl_calc_mass_hflx_vtx
      adv_calc_tracer_hflx_cell => ffsl_calc_tracer_hflx_cell
      adv_calc_tracer_vflx_cell => ffsl_calc_tracer_vflx_cell
      adv_calc_tracer_hflx_vtx  => ffsl_calc_tracer_hflx_vtx
      adv_calc_tracer_hval_vtx  => ffsl_calc_tracer_hval_vtx
    case default
      call log_error('Invalid adv_scheme ' // trim(adv_scheme) // '!', pid=proc%id)
    end select

    ntracer = 0
    allocate(batch_names      (1000))
    allocate(tracer_names     (1000))
    allocate(tracer_long_names(1000))
    allocate(tracer_units     (1000))
    allocate(tracer_dt        (1000))

    call tvd_init()

  end subroutine adv_init

  subroutine adv_run(block, itime)

    type(block_type), intent(inout) :: block
    integer, intent(in) :: itime

    integer i, j, k, l, m
    real(r8) work(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend,block%mesh%num_full_lev)
    real(r8) pole(block%mesh%num_full_lev)

    call adv_accum_wind(block, itime)

    do m = 1, size(block%adv_batches)
      if (time_is_alerted(block%adv_batches(m)%alert_key) .and. time_step > 0) then
        do l = 1, size(block%adv_batches(m)%tracer_names)
          associate (mesh    => block%mesh                  , &
                     old     => block%adv_batches(m)%old    , &
                     new     => block%adv_batches(m)%new    , &
                     old_m   => block%adv_batches(m)%old_m  , &
                     q       => block%adv_batches(m)%q      , &
                     qmf_lon => block%adv_batches(m)%qmf_lon, &
                     qmf_lat => block%adv_batches(m)%qmf_lat, &
                     qmf_lev => block%adv_batches(m)%qmf_lev)
          ! Calculate tracer mass flux.
          call adv_calc_tracer_hflx_cell(block, block%adv_batches(m), q(:,:,:,l,old), qmf_lon, qmf_lat)
          call fill_halo(block, qmf_lon, full_lon=.false., full_lat=.true., full_lev=.true., &
                         south_halo=.false., north_halo=.false., east_halo=.false.)
          call fill_halo(block, qmf_lat, full_lon=.true., full_lat=.false., full_lev=.true., &
                         north_halo=.false.,  west_halo=.false., east_halo=.false.)
          ! Update tracer mixing ratio.
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                q(i,j,k,l,new) = old_m(i,j,k) * q(i,j,k,l,old) - ( &
                  (                                                &
                    qmf_lon(i  ,j,k) -                             &
                    qmf_lon(i-1,j,k)                               &
                  ) * mesh%le_lon(j) + (                           &
                    qmf_lat(i,j  ,k) * mesh%le_lat(j  ) -          &
                    qmf_lat(i,j-1,k) * mesh%le_lat(j-1)            &
                  )                                                &
                ) / mesh%area_cell(j) * dt_adv
              end do
            end do
          end do
          if (mesh%has_south_pole()) then
            j = mesh%full_lat_ibeg
            do k = mesh%full_lev_ibeg, mesh%full_lev_iend
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                work(i,k) = qmf_lat(i,j,k)
              end do
            end do
            call zonal_sum(proc%zonal_circle, work, pole)
            pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j) * dt_adv
            do k = mesh%full_lev_ibeg, mesh%full_lev_iend
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                q(i,j,k,l,new) = old_m(i,j,k) * q(i,j,k,l,old) - pole(k)
              end do
            end do
          end if
          if (mesh%has_north_pole()) then
            j = mesh%full_lat_iend
            do k = mesh%full_lev_ibeg, mesh%full_lev_iend
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                work(i,k) = qmf_lat(i,j-1,k)
              end do
            end do
            call zonal_sum(proc%zonal_circle, work, pole)
            pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j) * dt_adv
            do k = mesh%full_lev_ibeg, mesh%full_lev_iend
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                q(i,j,k,l,new) = old_m(i,j,k) * q(i,j,k,l,old) + pole(k)
              end do
            end do
          end if
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%full_lat_ibeg, mesh%full_lat_iend
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                q(i,j,k,l,new) = q(i,j,k,l,new) / block%state(itime)%m(i,j,k)
              end do
            end do
          end do
          ! Set upper and lower boundary conditions.
          do k = mesh%full_lev_lb, mesh%full_lev_ibeg - 1
            q(:,:,k,l,new) = q(:,:,mesh%full_lev_ibeg,l,new)
          end do
          do k = mesh%full_lev_iend + 1, mesh%full_lev_ub
            q(:,:,k,l,new) = q(:,:,mesh%full_lev_iend,l,new)
          end do
          call adv_calc_tracer_vflx_cell(block, block%adv_batches(m), q(:,:,:,l,new), qmf_lev)
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%full_lat_ibeg, mesh%full_lat_iend
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
                q(i,j,k,l,new) = q(i,j,k,l,new) - ((qmf_lev(i,j,k+1) - qmf_lev(i,j,k)) * dt_adv) / block%state(itime)%m(i,j,k)
              end do
            end do
          end do
          call fill_halo(block, q(:,:,:,l,new), full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
          do i = 1, 5
            call zonal_damp_on_cell(block, 2, dt_adv, q(:,:,:,l,new))
          end do
          call fill_halo(block, q(:,:,:,l,new), full_lon=.true., full_lat=.true., full_lev=.true.)
          end associate
        end do
        i = block%adv_batches(m)%old
        block%adv_batches(m)%old = block%adv_batches(m)%new
        block%adv_batches(m)%new = i
      end if
      call block%adv_batches(m)%copy_old_m(block%state(itime)%m)
    end do

  end subroutine adv_run

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

  subroutine adv_allocate_tracers_1(block)

    type(block_type), intent(inout) :: block

    integer nbatch, nbatch_tracer, i, j, k
    logical found
    character(30) unique_batch_names(100)
    real(r8) unique_tracer_dt(100)

    if (.not. advection) then
      call block%adv_batch_mass%init(block%mesh, 'cell', 'mass', dt_dyn)
      ! call block%adv_batch_pv  %init(block%mesh, 'vtx' , 'pv'  , dt_dyn)
    end if

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
    if (nbatch == 0) return
    if (allocated(block%adv_batches)) then
      call log_error('Advection batches have already been alllocated!', pid=proc%id)
    end if
    if (is_root_proc()) then
      call log_notice('There are ' // to_str(nbatch) // ' advection batches.')
      do i = 1, nbatch
        write(*, *) '- ', trim(unique_batch_names(i)), real(unique_tracer_dt(i))
      end do
    end if

    ! Initialize advection batches in block objects and allocate tracer arrays in state objects.
    allocate(block%adv_batches(nbatch))
    do i = 1, nbatch
      call block%adv_batches(i)%init(block%mesh, 'cell', unique_batch_names(i), unique_tracer_dt(i))
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
      associate (mesh => block%mesh)
      allocate(block%adv_batches(i)%q      (mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub,nbatch_tracer,2))
      allocate(block%adv_batches(i)%qmf_lon(mesh%half_lon_lb:mesh%half_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
      allocate(block%adv_batches(i)%qmf_lat(mesh%full_lon_lb:mesh%full_lon_ub,mesh%half_lat_lb:mesh%half_lat_ub,mesh%full_lev_lb:mesh%full_lev_ub))
      allocate(block%adv_batches(i)%qmf_lev(mesh%full_lon_lb:mesh%full_lon_ub,mesh%full_lat_lb:mesh%full_lat_ub,mesh%half_lev_lb:mesh%half_lev_ub))
      end associate
      k = 0
      do j = 1, ntracer
        if (batch_names(j) == block%adv_batches(i)%alert_key) then
          k = k + 1
          block%adv_batches(i)%tracer_names     (k) = tracer_names     (j)
          block%adv_batches(i)%tracer_long_names(k) = tracer_long_names(j)
          block%adv_batches(i)%tracer_units     (k) = tracer_units     (j)
        end if
      end do
    end do

  end subroutine adv_allocate_tracers_1

  subroutine adv_allocate_tracers_2(blocks)

    type(block_type), intent(inout) :: blocks(:)

    integer iblk

    do iblk = 1, size(blocks)
      call adv_allocate_tracers_1(blocks(iblk))
    end do

  end subroutine adv_allocate_tracers_2

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
          if (global_mesh%num_full_lev > 1) then
            call block%adv_batches(l)%accum_we_lev( &
              block%state(itime)%we_lev           , &
              block%state(itime)%m_lev            )
          end if
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
