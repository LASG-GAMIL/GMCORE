module adv_batch_mod

  use flogger
  use const_mod
  use namelist_mod
  use mesh_mod
  use time_mod
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
    integer  :: we_step = 0 ! Step counter for we
    integer  :: mf_step = 0 ! Step counter for mass flux
    integer  :: old     = 1 ! Index for old time level
    integer  :: new     = 2 ! Index for new time level
    real(r8) :: dt          ! Advection time step size in seconds
    character(10), allocatable, dimension(:) :: tracer_names
    character(30), allocatable, dimension(:) :: tracer_long_names
    character(10), allocatable, dimension(:) :: tracer_units
    real(r8), allocatable, dimension(:,:,:) :: old_m ! Recorded old mass for converting mixing ratio
    real(r8), allocatable, dimension(:,:,:,:,:) :: q
    real(r8), allocatable, dimension(:,:,:) :: qmf_lon
    real(r8), allocatable, dimension(:,:,:) :: qmf_lat
    real(r8), allocatable, dimension(:,:,:) :: qmf_lev
    real(r8), allocatable, dimension(:,:,:) :: mfx
    real(r8), allocatable, dimension(:,:,:) :: mfy
    real(r8), allocatable, dimension(:,:,:) :: m_lev
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:) :: we
    real(r8), allocatable, dimension(:,:,:) :: cflx ! CFL number along x-axis
    real(r8), allocatable, dimension(:,:,:) :: cfly ! CFL number along y-axis
    real(r8), allocatable, dimension(:,:,:) :: cflz ! CFL number along z-axis
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
    procedure :: copy_old_m       => adv_batch_copy_old_m
    procedure :: accum_uv_cell    => adv_batch_accum_uv_cell
    procedure :: accum_mf_cell    => adv_batch_accum_mf_cell
    procedure :: accum_uv_vtx     => adv_batch_accum_uv_vtx
    procedure :: accum_mf_vtx     => adv_batch_accum_mf_vtx
    procedure :: accum_we_lev     => adv_batch_accum_we_lev
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
    this%we_step   = 0
    this%mf_step   = 0

    select case (loc)
    case ('cell')
      call allocate_array(mesh, this%old_m, full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%mfx , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%mfy , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%m_lev, full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%u   , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%v   , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%we  , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%cflx, half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%cfly, full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%cflz, full_lon=.true., full_lat=.true., half_lev=.true.)
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

    call time_add_alert(alert_key, seconds=dt)

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    this%mesh      => null()
    this%loc       = 'cell'
    this%alert_key = ''
    this%dt        = 0
    this%nstep     = 0
    this%uv_step   = 0
    this%we_step   = 0
    this%mf_step   = 0

    if (allocated(this%tracer_names     )) deallocate(this%tracer_names     )
    if (allocated(this%tracer_long_names)) deallocate(this%tracer_long_names)
    if (allocated(this%tracer_units     )) deallocate(this%tracer_units     )

    if (allocated(this%old_m  )) deallocate(this%old_m  )
    if (allocated(this%q      )) deallocate(this%q      )
    if (allocated(this%qmf_lon)) deallocate(this%qmf_lon)
    if (allocated(this%qmf_lat)) deallocate(this%qmf_lat)
    if (allocated(this%qmf_lev)) deallocate(this%qmf_lev)
    if (allocated(this%mfx )) deallocate(this%mfx )
    if (allocated(this%mfy )) deallocate(this%mfy )
    if (allocated(this%m_lev)) deallocate(this%m_lev)
    if (allocated(this%u   )) deallocate(this%u   )
    if (allocated(this%v   )) deallocate(this%v   )
    if (allocated(this%we  )) deallocate(this%we  )
    if (allocated(this%cflx)) deallocate(this%cflx)
    if (allocated(this%cfly)) deallocate(this%cfly)
    if (allocated(this%cflz)) deallocate(this%cflz)
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

    allocate(this%tracer_names     (ntracer))
    allocate(this%tracer_long_names(ntracer))
    allocate(this%tracer_units     (ntracer))

  end subroutine adv_batch_allocate_tracers

  subroutine adv_batch_copy_old_m(this, m)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: m(this%mesh%full_lon_lb:this%mesh%full_lon_ub, &
                              this%mesh%full_lat_lb:this%mesh%full_lat_ub, &
                              this%mesh%full_lev_lb:this%mesh%full_lev_ub)

    this%old_m = m

  end subroutine adv_batch_copy_old_m

  subroutine adv_batch_accum_uv_cell(this, u, v, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: u(this%mesh%half_lon_lb:this%mesh%half_lon_ub, &
                              this%mesh%full_lat_lb:this%mesh%full_lat_ub, &
                              this%mesh%full_lev_lb:this%mesh%full_lev_ub)
    real(r8), intent(in) :: v(this%mesh%full_lon_lb:this%mesh%full_lon_ub, &
                              this%mesh%half_lat_lb:this%mesh%half_lat_ub, &
                              this%mesh%full_lev_lb:this%mesh%full_lev_ub)
    real(8), intent(in), optional :: dt

    real(r8) work(this%mesh%full_lon_ibeg:this%mesh%full_lon_iend,this%mesh%num_full_lev)
    real(r8) pole(this%mesh%num_full_lev)
    real(8) dt_
    real(r8) dx, x0, x1, x2, x3, u1, u2, u3, u4
    real(r8) dy, y0, y1, y2, y3, v1, v2, v3, v4
    integer i, j, k, l

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
          dx = mesh%de_lon(j)
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            x0 = radius * mesh%half_lon(i) * mesh%full_cos_lat(j)
            ! Stage 1
            u1 = this%u(i,j,k)
            x1 = x0 + dt_ / 2 * u1
            ! Stage 2
            l = floor(i + (x1 - x0) / dx)
            u2 = (                                                                          &
              (radius * mesh%half_lon(l+1) * mesh%full_cos_lat(j) - x1) * this%u(l  ,j,k) + &
              (x1 - radius * mesh%half_lon(l  ) * mesh%full_cos_lat(j)) * this%u(l+1,j,k)   &
            ) / dx
            x2 = x0 + dt_ / 2 * u2
            ! Stage 3
            l = floor(i + (x2 - x0) / dx)
            u3 = (                                                                          &
              (radius * mesh%half_lon(l+1) * mesh%full_cos_lat(j) - x2) * this%u(l  ,j,k) + &
              (x2 - radius * mesh%half_lon(l  ) * mesh%full_cos_lat(j)) * this%u(l+1,j,k)   &
            ) / dx
            x3 = x0 + dt_ * u3
            ! Stage 4
            l = floor(i + (x3 - x0) / dx)
            u4 = (                                                                          &
              (radius * mesh%half_lon(l+1) * mesh%full_cos_lat(j) - x3) * this%u(l  ,j,k) + &
              (x3 - radius * mesh%half_lon(l  ) * mesh%full_cos_lat(j)) * this%u(l+1,j,k)   &
            ) / dx
            ! Final stage
            this%cflx(i,j,k) = (u1 + 2 * u2 + 2 * u3 + u4) / 6 * dt_ / dx
            if (abs(this%cflx(i,j,k)) > mesh%lon_halo_width) then
              call log_error('cflx exceeds mesh%lon_halo_width ' // &
                             to_str(mesh%lon_halo_width) // ' at j=' // to_str(j) // '!', __FILE__, __LINE__)
            end if
          end do
        end do
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          dy = mesh%de_lat(j)
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            y0 = radius * mesh%half_lat(j)
            ! Stage 1
            v1 = this%v(i,j,k)
            y1 = y0 + dt_ / 2 * v1
            ! Stage 2
            l = floor(j + (y1 - y0) / dy)
            if (l < 1 .or. l > global_mesh%num_half_lat - 1) then
              v2 = v1
            else
              v2 = (                                                   &
                (radius * mesh%half_lat(l+1) - y1) * this%v(i,l  ,k) + &
                (y1 - radius * mesh%half_lat(l  )) * this%v(i,l+1,k)   &
              ) / dy
            end if
            y2 = y0 + dt_ / 2 * v2
            ! Stage 3
            l = floor(j + (y2 - y0) / dy)
            if (l < 1 .or. l > global_mesh%num_half_lat - 1) then
              v3 = v1
            else
              v3 = (                                                   &
                (radius * mesh%half_lat(l+1) - y2) * this%v(i,l  ,k) + &
                (y2 - radius * mesh%half_lat(l  )) * this%v(i,l+1,k)   &
              ) / dy
            end if
            y3 = y0 + dt_ * v3
            ! Stage 4
            l = floor(j + (y3 - y0) / dy)
            if (l < 1 .or. l > global_mesh%num_half_lat - 1) then
              v4 = v1
            else
              v4 = (                                                   &
                (radius * mesh%half_lat(l+1) - y3) * this%v(i,l  ,k) + &
                (y3 - radius * mesh%half_lat(l  )) * this%v(i,l+1,k)   &
              ) / dy
            end if
            ! Final stage
            this%cfly(i,j,k) = (v1 + 2 * v2 + 2 * v3 + v4) / 6 * dt_ / dy
            if (abs(this%cfly(i,j,k)) > 1) then
              call log_error('cfly exceeds 1 at j=' // to_str(j) // '!', __FILE__, __LINE__)
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
    real(8), intent(in), optional :: dt

    real(8) dt_
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
              call log_error('cfly exceeds 1 at j=' // to_str(j) // '!', __FILE__, __LINE__)
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

  subroutine adv_batch_accum_we_lev(this, we, m_lev, dt)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: we   (this%mesh%full_lon_lb:this%mesh%full_lon_ub, &
                                  this%mesh%full_lat_lb:this%mesh%full_lat_ub, &
                                  this%mesh%half_lev_lb:this%mesh%half_lev_ub)
    real(r8), intent(in) :: m_lev(this%mesh%full_lon_lb:this%mesh%full_lon_ub, &
                                  this%mesh%full_lat_lb:this%mesh%full_lat_ub, &
                                  this%mesh%half_lev_lb:this%mesh%half_lev_ub)

    real(8), intent(in), optional :: dt

    real(8) dt_
    real(r8) z0, z1, z2, z3, w1, w2, w3, w4, deta
    integer i, j, k, l, ks, ke, s

    dt_ = merge(dt, this%dt, present(dt))

    associate (mesh => this%mesh)
    if (this%we_step == 0) then
      this%we    = we
      this%m_lev = m_lev
      if (this%nstep > 1) this%we_step = this%we_step + 1
    else if (this%we_step == this%nstep) then
      this%we    = (this%we    + we   ) / (this%nstep + 1)
      this%m_lev = (this%m_lev + m_lev) / (this%nstep + 1)
      this%we_step = 0
    else
      this%we    = this%we    + we
      this%m_lev = this%m_lev + m_lev
      this%we_step = this%we_step + 1
    end if
    if (this%we_step == 0) then
      do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            z0 = mesh%half_lev(k)
            ! Stage 1
            w1 = this%we(i,j,k) / this%m_lev(i,j,k) * mesh%half_dlev(k)
            z1 = z0 + dt_ / 2 * w1
            ! Stage 2
            if (w1 > 0) then
              do l = k, mesh%full_lev_iend
                if (mesh%half_lev(l) <= z1 .and. z1 <= mesh%half_lev(l+1)) exit
              end do
              if (l == mesh%full_lev_iend + 1) then
                call log_error('cflz exceeds range at i=' // to_str(i) // &
                                                   ', j=' // to_str(j) // &
                                                   ', k=' // to_str(k) // '!', __FILE__, __LINE__)
              end if
            else
              do l = k - 1, mesh%full_lev_ibeg, -1
                if (mesh%half_lev(l) <= z1 .and. z1 <= mesh%half_lev(l+1)) exit
              end do
              if (l == mesh%full_lev_ibeg - 1) then
                call log_error('cflz exceeds range at i=' // to_str(i) // &
                                                   ', j=' // to_str(j) // &
                                                   ', k=' // to_str(k) // '!', __FILE__, __LINE__)
              end if
            end if
            w2 = (                                                                                       &
              (mesh%half_lev(l+1) - z1) * this%we(i,j,l  ) / this%m_lev(i,j,l  ) * mesh%half_dlev(l  ) + &
              (z1 - mesh%half_lev(l  )) * this%we(i,j,l+1) / this%m_lev(i,j,l+1) * mesh%half_dlev(l+1)   &
            ) / mesh%full_dlev(l)
            z2 = z0 + dt_ / 2 * w2
            ! Stage 3
            if (w2 > 0) then
              do l = k, mesh%full_lev_iend
                if (mesh%half_lev(l) <= z2 .and. z2 <= mesh%half_lev(l+1)) exit
              end do
              if (l == mesh%full_lev_iend + 1) then
                call log_error('cflz exceeds range at i=' // to_str(i) // &
                                                   ', j=' // to_str(j) // &
                                                   ', k=' // to_str(k) // '!', __FILE__, __LINE__)
              end if
            else
              do l = k - 1, mesh%full_lev_ibeg, -1
                if (mesh%half_lev(l) <= z2 .and. z2 <= mesh%half_lev(l+1)) exit
              end do
              if (l == mesh%full_lev_ibeg - 1) then
                call log_error('cflz exceeds range at i=' // to_str(i) // &
                                                   ', j=' // to_str(j) // &
                                                   ', k=' // to_str(k) // '!', __FILE__, __LINE__)
              end if
            end if
            w3 = (                                                                                       &
              (mesh%half_lev(l+1) - z2) * this%we(i,j,l  ) / this%m_lev(i,j,l  ) * mesh%half_dlev(l  ) + &
              (z2 - mesh%half_lev(l  )) * this%we(i,j,l+1) / this%m_lev(i,j,l+1) * mesh%half_dlev(l+1)   &
            ) / mesh%full_dlev(l)
            z3 = z0 + dt_ * w3
            ! Stage 4
            if (w3 > 0) then
              do l = k, mesh%full_lev_iend
                if (mesh%half_lev(l) <= z3 .and. z3 <= mesh%half_lev(l+1)) exit
              end do
              if (l == mesh%full_lev_iend + 1) then
                call log_error('cflz exceeds range at i=' // to_str(i) // &
                                                   ', j=' // to_str(j) // &
                                                   ', k=' // to_str(k) // '!', __FILE__, __LINE__)
              end if
            else
              do l = k - 1, mesh%full_lev_ibeg, -1
                if (mesh%half_lev(l) <= z3 .and. z3 <= mesh%half_lev(l+1)) exit
              end do
              if (l == mesh%full_lev_ibeg - 1) then
                call log_error('cflz exceeds range at i=' // to_str(i) // &
                                                   ', j=' // to_str(j) // &
                                                   ', k=' // to_str(k) // '!', __FILE__, __LINE__)
              end if
            end if
            w4 = (                                                                                       &
              (mesh%half_lev(l+1) - z3) * this%we(i,j,l  ) / this%m_lev(i,j,l  ) * mesh%half_dlev(l  ) + &
              (z3 - mesh%half_lev(l  )) * this%we(i,j,l+1) / this%m_lev(i,j,l+1) * mesh%half_dlev(l+1)   &
            ) / mesh%full_dlev(l)
            ! Final stage
            deta = (w1 + 2 * w2 + 2 * w3 + w4) / 6 * dt_
            if (deta < 0) then
              ks = k - 1
              ke = mesh%full_lev_ibeg
              s = -1
            else
              ks = k
              ke = mesh%full_lev_iend
              s = 1
            end if
            deta = abs(deta)
            do l = ks, ke, s
              deta = deta - mesh%full_dlev(l)
              if (deta < 0) then
                this%cflz(i,j,k) = s * ((1 + deta / mesh%full_dlev(l)) + abs(l - ks))
                exit
              end if
            end do
          end do
        end do
      end do
    end if
    end associate

  end subroutine adv_batch_accum_we_lev

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

end module adv_batch_mod
