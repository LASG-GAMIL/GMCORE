module adv_batch_mod

  use flogger
  use const_mod
  use namelist_mod
  use mesh_mod
  use allocator_mod
  use parallel_types_mod
  use parallel_zonal_mod

  implicit none

  private

  public adv_batch_type

  ! Different tracers can be combined into one batch, and adved in different frequencfly.
  type adv_batch_type
    type(mesh_type), pointer :: mesh => null()
    character(10) :: loc       = 'cell'
    character(30) :: alert_key = ''
    integer  :: nstep   = 0 ! Number of dynamic steps for one adv step
    integer  :: uv_step = 0 ! Step counter for u and v
    integer  :: mf_step = 0 ! Step counter for mass flux
    real(r8) :: dt          ! Advection time step size in seconds
    integer      , allocatable, dimension(:) :: tracer_idx
    character(10), allocatable, dimension(:) :: tracer_names
    character(30), allocatable, dimension(:) :: tracer_long_names
    character(10), allocatable, dimension(:) :: tracer_units
    real(r8), allocatable, dimension(:,:,:) :: mfx
    real(r8), allocatable, dimension(:,:,:) :: mfy
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:) :: cflx ! Fractional CFL number along x-axis
    real(r8), allocatable, dimension(:,:,:) :: cfly ! Fractional CFL number along y-axis
    real(r8), allocatable, dimension(:,:,:) :: divx ! Divergence along x-axis
    real(r8), allocatable, dimension(:,:,:) :: divy ! Divergence along y-axis
    real(r8), allocatable, dimension(:,:,:) :: qx   ! Tracer mixing ratio due to advective operator along x axis
    real(r8), allocatable, dimension(:,:,:) :: qy   ! Tracer mixing ratio due to advective operator along y axis
    real(r8), allocatable, dimension(:,:,:) :: qlx  ! Tracer mixing ratio at left cell edge along x axis
    real(r8), allocatable, dimension(:,:,:) :: qly  ! Tracer mixing ratio at left cell edge along y axis
    real(r8), allocatable, dimension(:,:,:) :: dqx  ! Tracer mixing ratio mismatch (or slope) at cell center along x axis
    real(r8), allocatable, dimension(:,:,:) :: dqy  ! Tracer mixing ratio mismatch (or slope) at cell center along y axis
    real(r8), allocatable, dimension(:,:,:) :: q6x  ! PPM mismatch at cell center along x axis
    real(r8), allocatable, dimension(:,:,:) :: q6y  ! PPM mismatch at cell center along y axis
  contains
    procedure :: init             => adv_batch_init
    procedure :: clear            => adv_batch_clear
    procedure :: allocate_tracers => adv_batch_allocate_tracers
    procedure :: accum_uv_cell    => adv_batch_accum_uv_cell
    procedure :: accum_mf_cell    => adv_batch_accum_mf_cell
    procedure :: accum_uv_vtx     => adv_batch_accum_uv_vtx
    procedure :: accum_mf_vtx     => adv_batch_accum_mf_vtx
    final :: adv_batch_final
  end type adv_batch_type

contains

  subroutine adv_batch_init(this, mesh, loc, alert_key, dt)

    class(adv_batch_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: loc
    character(*), intent(in) :: alert_key
    real(r8), intent(in) :: dt

    call this%clear()

    this%mesh      => mesh
    this%loc       = loc
    this%alert_key = alert_key
    this%dt        = dt
    this%nstep     = dt / dt_dyn
    this%uv_step   = 0
    this%mf_step   = 0

    select case (loc)
    case ('cell')
      call allocate_array(mesh, this%mfx , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%mfy , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%u   , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%v   , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%cflx, half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%cfly, full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%divx, full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%divy, full_lon=.true., full_lat=.true., full_lev=.true.)
      select case (adv_scheme)
      case ('ffsl')
        call allocate_array(mesh, this%qx, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%qy, full_lon=.true., full_lat=.true., full_lev=.true.)
        if (ffsl_flux_type == 'ppm') then
          call allocate_array(mesh, this%qlx, full_lon=.true., full_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%qly, full_lon=.true., full_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%dqx, full_lon=.true., full_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%dqy, full_lon=.true., full_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%q6x, full_lon=.true., full_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%q6y, full_lon=.true., full_lat=.true., full_lev=.true.)
        end if
      end select
    case ('vtx')
      call allocate_array(mesh, this%mfx , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%mfy , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%u   , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%v   , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%cflx, full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%cfly, half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%divx, half_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%divy, half_lon=.true., half_lat=.true., full_lev=.true.)
      select case (adv_scheme)
      case ('ffsl')
        call allocate_array(mesh, this%qx, half_lon=.true., half_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%qy, half_lon=.true., half_lat=.true., full_lev=.true.)
        if (ffsl_flux_type == 'ppm') then
          call allocate_array(mesh, this%qlx, half_lon=.true., half_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%qly, half_lon=.true., half_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%dqx, half_lon=.true., half_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%dqy, half_lon=.true., half_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%q6x, half_lon=.true., half_lat=.true., full_lev=.true.)
          call allocate_array(mesh, this%q6y, half_lon=.true., half_lat=.true., full_lev=.true.)
        end if
      end select
    case default
      call log_error('Invalid grid location ' // trim(loc) // '!', __FILE__, __LINE__)
    end select

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    this%mesh      => null()
    this%loc       = 'cell'
    this%alert_key = ''
    this%dt        = 0
    this%nstep     = 0
    this%uv_step   = 0
    this%mf_step   = 0

    if (allocated(this%tracer_idx       )) deallocate(this%tracer_idx       )
    if (allocated(this%tracer_names     )) deallocate(this%tracer_names     )
    if (allocated(this%tracer_long_names)) deallocate(this%tracer_long_names)
    if (allocated(this%tracer_units     )) deallocate(this%tracer_units     )

    if (allocated(this%mfx )) deallocate(this%mfx )
    if (allocated(this%mfy )) deallocate(this%mfy )
    if (allocated(this%u   )) deallocate(this%u   )
    if (allocated(this%v   )) deallocate(this%v   )
    if (allocated(this%cflx)) deallocate(this%cflx)
    if (allocated(this%cfly)) deallocate(this%cfly)
    if (allocated(this%divx)) deallocate(this%divx)
    if (allocated(this%divy)) deallocate(this%divy)
    if (allocated(this%qx  )) deallocate(this%qx  )
    if (allocated(this%qy  )) deallocate(this%qy  )
    if (allocated(this%qlx )) deallocate(this%qlx )
    if (allocated(this%qly )) deallocate(this%qly )
    if (allocated(this%dqx )) deallocate(this%dqx )
    if (allocated(this%dqy )) deallocate(this%dqy )
    if (allocated(this%q6x )) deallocate(this%q6x )
    if (allocated(this%q6y )) deallocate(this%q6y )

  end subroutine adv_batch_clear

  subroutine adv_batch_allocate_tracers(this, ntracer)

    class(adv_batch_type), intent(inout) :: this
    integer, intent(in) :: ntracer

    allocate(this%tracer_idx       (ntracer))
    allocate(this%tracer_names     (ntracer))
    allocate(this%tracer_long_names(ntracer))
    allocate(this%tracer_units     (ntracer))

  end subroutine adv_batch_allocate_tracers

  subroutine adv_batch_accum_uv_cell(this, u, v, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: u(this%mesh%half_lon_lb:this%mesh%half_lon_ub, &
                              this%mesh%full_lat_lb:this%mesh%full_lat_ub, &
                              this%mesh%full_lev_lb:this%mesh%full_lev_ub)
    real(r8), intent(in) :: v(this%mesh%full_lon_lb:this%mesh%full_lon_ub, &
                              this%mesh%half_lat_lb:this%mesh%half_lat_ub, &
                              this%mesh%full_lev_lb:this%mesh%full_lev_ub)
    real(r8), intent(in), optional :: dt

    real(r8) work(this%mesh%full_lon_ibeg:this%mesh%full_lon_iend,this%mesh%num_full_lev)
    real(r8) pole(this%mesh%num_full_lev)
    real(r8) dt_
    integer i, j, k

    dt_ = merge(dt, this%dt, present(dt))

    associate (mesh => this%mesh)
    if (this%uv_step == 0) then
      this%u = u
      this%v = v
      if (this%nstep > 1) this%uv_step = this%uv_step + 1
    else if (this%uv_step == this%nstep) then
      this%u = (this%u + u) / (this%nstep + 1)
      this%v = (this%v + v) / (this%nstep + 1)
      this%uv_step = 0
    else
      this%u = this%u + u
      this%v = this%v + v
      this%uv_step = this%uv_step + 1
    end if
    if (this%uv_step == 0) then
      ! Calculate CFL numbers and divergence along each axis.
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            this%cflx(i,j,k) = dt_ * this%u(i,j,k) / mesh%de_lon(j)
            if (abs(this%cflx(i,j,k)) > mesh%lon_halo_width) then
              call log_error('cflx exceeds mesh%lon_halo_width ' // to_str(mesh%lon_halo_width) // '!', __FILE__, __LINE__)
            end if
          end do
        end do
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            this%cfly(i,j,k) = dt_ * this%v(i,j,k) / mesh%de_lat(j)
            if (abs(this%cfly(i,j,k)) > 1) then
              call log_error('cfly exceeds 1!', __FILE__, __LINE__)
            end if
          end do
        end do
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            this%divx(i,j,k) = (this%u(i,j,k) - this%u(i-1,j,k)) * mesh%le_lon(j) / mesh%area_cell(j)
            this%divy(i,j,k) = (this%v(i,j  ,k) * mesh%le_lat(j  ) - &
                                this%v(i,j-1,k) * mesh%le_lat(j-1)) / mesh%area_cell(j)
          end do
        end do
      end do
      if (mesh%has_south_pole()) then
        j = mesh%full_lat_ibeg
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            work(i,k) = this%v(i,j,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            this%divy(i,j,k) = pole(k)
          end do
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%full_lat_iend
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            work(i,k) = -this%v(i,j-1,k)
          end do
        end do
        call zonal_sum(proc%zonal_circle, work, pole)
        pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            this%divy(i,j,k) = pole(k)
          end do
        end do
      end if
    end if
    end associate

  end subroutine adv_batch_accum_uv_cell

  subroutine adv_batch_accum_mf_cell(this, mfx, mfy)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: mfx(this%mesh%half_lon_lb:this%mesh%half_lon_ub, &
                                this%mesh%full_lat_lb:this%mesh%full_lat_ub, &
                                this%mesh%full_lev_lb:this%mesh%full_lev_ub)
    real(r8), intent(in) :: mfy(this%mesh%full_lon_lb:this%mesh%full_lon_ub, &
                                this%mesh%half_lat_lb:this%mesh%half_lat_ub, &
                                this%mesh%full_lev_lb:this%mesh%full_lev_ub)

    if (this%mf_step == 0) then
      this%mfx = mfx
      this%mfy = mfy
      if (this%nstep > 1) this%mf_step = this%mf_step + 1
    else if (this%mf_step == this%nstep) then
      this%mfx = this%mfx + mfx
      this%mfy = this%mfy + mfy
      this%mf_step = 0
    else
      this%mfx = this%mfx + mfx
      this%mfy = this%mfy + mfy
      this%mf_step = this%mf_step + 1
    end if

  end subroutine adv_batch_accum_mf_cell

  subroutine adv_batch_accum_uv_vtx(this, u, v, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: u(this%mesh%full_lon_lb:this%mesh%full_lon_ub, &
                              this%mesh%half_lat_lb:this%mesh%half_lat_ub, &
                              this%mesh%full_lev_lb:this%mesh%full_lev_ub)
    real(r8), intent(in) :: v(this%mesh%half_lon_lb:this%mesh%half_lon_ub, &
                              this%mesh%full_lat_lb:this%mesh%full_lat_ub, &
                              this%mesh%full_lev_lb:this%mesh%full_lev_ub)
    real(r8), intent(in), optional :: dt

    real(r8) dt_
    integer i, j, k

    dt_ = merge(dt, this%dt, present(dt))

    associate (mesh => this%mesh)
    if (this%uv_step == 0) then
      this%u = u
      this%v = v
      if (this%nstep > 1) this%uv_step = this%uv_step + 1
    else if (this%uv_step == this%nstep) then
      this%u = (this%u + u) / (this%nstep + 1)
      this%v = (this%v + v) / (this%nstep + 1)
      this%uv_step = 0
    else
      this%u = this%u + u
      this%v = this%v + v
      this%uv_step = this%uv_step + 1
    end if
    if (this%uv_step == 0) then
      ! Calculate CFL numbers and divergence along each axis.
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            this%cflx(i,j,k) = dt_ * this%u(i,j,k) / mesh%le_lat(j)
            if (abs(this%cflx(i,j,k)) > mesh%lon_halo_width) then
              call log_error('cflx exceeds mesh%lon_halo_width ' // to_str(mesh%lon_halo_width) // '!', __FILE__, __LINE__)
            end if
          end do
        end do
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            this%cfly(i,j,k) = dt_ * this%v(i,j,k) / mesh%le_lon(j)
            if (abs(this%cfly(i,j,k)) > 1) then
              call log_error('cfly exceeds 1!', __FILE__, __LINE__)
            end if
          end do
        end do
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            this%divx(i,j,k) = (this%u(i+1,j,k) - this%u(i,j,k)) * mesh%de_lat(j) / mesh%area_vtx(j)
            this%divy(i,j,k) = (this%v(i,j+1,k) * mesh%de_lon(j+1) - &
                                this%v(i,j  ,k) * mesh%de_lon(j  )) / mesh%area_vtx(j)
          end do
        end do
      end do
    end if
    end associate

  end subroutine adv_batch_accum_uv_vtx

  subroutine adv_batch_accum_mf_vtx(this, mfx, mfy)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: mfx(this%mesh%full_lon_lb:this%mesh%full_lon_ub, &
                                this%mesh%half_lat_lb:this%mesh%half_lat_ub, &
                                this%mesh%full_lev_lb:this%mesh%full_lev_ub)
    real(r8), intent(in) :: mfy(this%mesh%half_lon_lb:this%mesh%half_lon_ub, &
                                this%mesh%full_lat_lb:this%mesh%full_lat_ub, &
                                this%mesh%full_lev_lb:this%mesh%full_lev_ub)

    if (this%mf_step == 0) then
      this%mfx = mfx
      this%mfy = mfy
      if (this%nstep > 1) this%mf_step = this%mf_step + 1
    else if (this%mf_step == this%nstep) then
      this%mfx = this%mfx + mfx
      this%mfy = this%mfy + mfy
      this%mf_step = 0
    else
      this%mfx = this%mfx + mfx
      this%mfy = this%mfy + mfy
      this%mf_step = this%mf_step + 1
    end if

  end subroutine adv_batch_accum_mf_vtx

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

end module adv_batch_mod
