module adv_batch_mod

  use flogger
  use const_mod
  use namelist_mod
  use mesh_mod
  use allocator_mod

  private

  public adv_batch_type

  ! Different tracers can be combined into one batch, and adved in different frequencfly.
  type adv_batch_type
    type(mesh_type), pointer :: mesh => null()
    character(30) :: alert_key = ''
    integer  :: nstep  = 0 ! Number of dynamic steps for one adv step
    integer  :: step   = 0 ! Dynamic step counter
    real(r8) :: dt         ! Advection time step size in seconds
    integer      , allocatable, dimension(:) :: tracer_idx
    character(10), allocatable, dimension(:) :: tracer_names
    character(10), allocatable, dimension(:) :: tracer_long_names
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
    real(r8), allocatable, dimension(:,:,:) :: qxl  ! Tracer mixing ratio at left cell edge along x axis
    real(r8), allocatable, dimension(:,:,:) :: qyl  ! Tracer mixing ratio at left cell edge along y axis
    real(r8), allocatable, dimension(:,:,:) :: dqx  ! Tracer mixing ratio mismatch (or slope) at cell center along x axis
    real(r8), allocatable, dimension(:,:,:) :: dqy  ! Tracer mixing ratio mismatch (or slope) at cell center along y axis
    real(r8), allocatable, dimension(:,:,:) :: qx6  ! PPM mismatch at cell center along x axis
    real(r8), allocatable, dimension(:,:,:) :: qy6  ! PPM mismatch at cell center along y axis
  contains
    procedure :: init       => adv_batch_init
    procedure :: clear      => adv_batch_clear
    procedure :: accum_wind => adv_batch_accum_wind
    procedure :: allocate_tracers => adv_allocate_tracers
    final :: adv_batch_final
  end type adv_batch_type

contains

  subroutine adv_batch_init(this, mesh, alert_key, dt)

    class(adv_batch_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: alert_key
    real(r8), intent(in) :: dt

    call this%clear()

    this%mesh => mesh
    this%alert_key = alert_key
    this%dt        = dt
    this%nstep     = dt / dt_dyn
    this%step      = 0

    call allocate_array(mesh, this%mfx , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mfy , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u   , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v   , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%cflx, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%cfly, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%divx, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%divy, full_lon=.true., half_lat=.true., full_lev=.true.)
    select case (adv_scheme)
    case ('ffsl')
      call allocate_array(mesh, this%qx  , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%qy  , full_lon=.true., full_lat=.true., full_lev=.true.)
      if (ffsl_flux_type == 'ppm') then
        call allocate_array(mesh, this%qxl, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%qyl, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%dqx, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%dqy, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%qx6, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(mesh, this%qy6, full_lon=.true., full_lat=.true., full_lev=.true.)
      end if
    end select

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    this%alert_key = ''
    this%dt        = 0
    this%nstep     = 0
    this%step      = 0

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
    if (allocated(this%qxl )) deallocate(this%qxl )
    if (allocated(this%qyl )) deallocate(this%qyl )
    if (allocated(this%dqx )) deallocate(this%dqx )
    if (allocated(this%dqy )) deallocate(this%dqy )
    if (allocated(this%qx6 )) deallocate(this%qx6 )
    if (allocated(this%qy6 )) deallocate(this%qy6 )

  end subroutine adv_batch_clear

  subroutine adv_batch_accum_wind(this, u, v, mfx, mfy)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in)           :: u  (:,:,:)
    real(r8), intent(in)           :: v  (:,:,:)
    real(r8), intent(in), optional :: mfx(:,:,:)
    real(r8), intent(in), optional :: mfy(:,:,:)

    integer i, j, k

    if (this%step == 0) then
      this%u      = u
      this%v      = v
      if (present(mfx)) then
        this%mfx  = mfx
        this%mfy  = mfy
      end if
      this%step   = this%step + 1
    else if (this%step == this%nstep) then
      this%u      = this%u / this%nstep
      this%v      = this%v / this%nstep
      this%step   = 0
      ! Calculate CFL numbers and divergence along each axis.
      do k = this%mesh%full_lev_ibeg, this%mesh%full_lev_iend
        do j = this%mesh%full_lat_ibeg_no_pole, this%mesh%full_lat_iend_no_pole
          do i = this%mesh%half_lon_ibeg, this%mesh%half_lon_iend
            this%cflx(i,j,k) = this%dt * this%u(i,j,k) / this%mesh%de_lon(j)
            this%divx(i,j,k) = (this%u(i+1,j,k) - this%u(i,j,k)) / this%mesh%de_lon(j)
          end do
        end do
        do j = this%mesh%half_lat_ibeg, this%mesh%half_lat_iend
          do i = this%mesh%full_lon_ibeg, this%mesh%full_lon_iend
            this%cfly(i,j,k) = this%dt * this%v(i,j,k) / this%mesh%de_lat(j)
            this%divy(i,j,k) = (this%v(i,j+1,k) * this%mesh%half_cos_lat(j+1) - &
                                this%v(i,j  ,k) * this%mesh%half_cos_lat(j  )) / &
                               this%mesh%le_lon(j) / this%mesh%full_cos_lat(j)
          end do
        end do
      end do
    else
      this%u      = this%u   + u
      this%v      = this%v   + v
      if (present(mfx)) then
        this%mfx  = this%mfx + mfx
        this%mfy  = this%mfy + mfy
      end if
      this%step   = this%step + 1
    end if

  end subroutine adv_batch_accum_wind

  subroutine adv_allocate_tracers(this, ntracer)

    class(adv_batch_type), intent(inout) :: this
    integer, intent(in) :: ntracer

    allocate(this%tracer_idx       (ntracer))
    allocate(this%tracer_names     (ntracer))
    allocate(this%tracer_long_names(ntracer))
    allocate(this%tracer_units     (ntracer))

  end subroutine adv_allocate_tracers

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

end module adv_batch_mod