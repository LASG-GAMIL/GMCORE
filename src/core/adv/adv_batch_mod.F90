module adv_batch_mod

  use container
  use flogger
  use const_mod
  use namelist_mod
  use mesh_mod
  use allocator_mod
  use tracer_mod

  private

  public adv_batch_type
  public tracer_type

  ! Different tracers can be combined into one batch, and adved in different frequencfly.
  type adv_batch_type
    type(mesh_type), pointer :: mesh => null()
    character(30) :: alert_key = ''
    integer :: nstep  = 0 ! Number of dynamic steps for one adv step
    integer :: step   = 0 ! Dynamic step counter
    real(r8) :: dt        ! adv time step size in seconds
    real(r8), allocatable, dimension(:,:,:) :: mfx
    real(r8), allocatable, dimension(:,:,:) :: mfy
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:) :: cflx ! Fractional CFL number along x-axis
    real(r8), allocatable, dimension(:,:,:) :: cfly ! Fractional CFL number along y-axis
    real(r8), allocatable, dimension(:,:,:) :: divx ! Divergence along x-axis
    real(r8), allocatable, dimension(:,:,:) :: divy ! Divergence along y-axis
    type(hash_table_type) tracers
  contains
    procedure :: init       => adv_batch_init
    procedure :: clear      => adv_batch_clear
    procedure :: accum_wind => adv_batch_accum_wind
    procedure :: add_tracer => adv_add_tracer
    procedure :: get_tracer => adv_get_tracer
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

    call allocate_array(mesh, this%mfx, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mfy, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u  , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v  , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%cflx , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%cfly , full_lon=.true., half_lat=.true., full_lev=.true.)

    call this%tracers%init()

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    type(hash_table_iterator_type) it

    this%alert_key = ''
    this%dt        = 0
    this%nstep     = 0
    this%step      = 0
    if (allocated(this%mfx )) deallocate(this%mfx )
    if (allocated(this%mfy )) deallocate(this%mfy )
    if (allocated(this%u   )) deallocate(this%u   )
    if (allocated(this%v   )) deallocate(this%v   )
    if (allocated(this%cflx)) deallocate(this%cflx)
    if (allocated(this%cfly)) deallocate(this%cfly)
    if (allocated(this%divx)) deallocate(this%divx)
    if (allocated(this%divy)) deallocate(this%divy)

    it = hash_table_iterator(this%tracers)
    do while (.not. it%ended())
      select type (tracer => it%value)
      type is (tracer_type)
        call tracer%clear()
      end select
      call it%next()
    end do
    call this%tracers%clear()

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

  subroutine adv_add_tracer(this, name, long_name, units)

    class(adv_batch_type), intent(inout) :: this
    character(*), intent(in) :: name
    character(*), intent(in), optional :: long_name
    character(*), intent(in), optional :: units

    type(tracer_type) tracer

    call tracer%init(this%mesh, name, long_name, units)
    call this%tracers%insert(name, tracer)

  end subroutine adv_add_tracer

  function adv_get_tracer(this, name) result(res)

    class(adv_batch_type), intent(in) :: this
    character(*), intent(in) :: name

    type(tracer_type), pointer :: res

    select type (value => this%tracers%value(name))
    type is (tracer_type)
      res => value
    class default
      call log_error('Unknown tracer name ' // trim(name) // '!')
    end select

  end function adv_get_tracer

  subroutine adv_allocate_tracers(this)

    class(adv_batch_type), intent(inout) :: this

    type(hash_table_iterator_type) it

    it = hash_table_iterator(this%tracers)
    do while (.not. it%ended())
      select type (tracer => it%value)
      type is (tracer_type)
        call tracer%allocate_arrays()
      end select
      call it%next()
    end do

  end subroutine adv_allocate_tracers

  subroutine adv_batch_final(this)

    type(adv_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine adv_batch_final

end module adv_batch_mod