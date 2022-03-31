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

  ! Different tracers can be combined into one batch, and adved in different frequency.
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
    real(r8), allocatable, dimension(:,:,:) :: cx ! Fractional CFL number along x-axis
    real(r8), allocatable, dimension(:,:,:) :: cy ! Fractional CFL number along y-axis
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
    call allocate_array(mesh, this%cx , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%cy , full_lon=.true., half_lat=.true., full_lev=.true.)

    call this%tracers%init()

  end subroutine adv_batch_init

  subroutine adv_batch_clear(this)

    class(adv_batch_type), intent(inout) :: this

    type(hash_table_iterator_type) it

    this%alert_key = ''
    this%dt        = 0
    this%nstep     = 0
    this%step      = 0
    if (allocated(this%mfx)) deallocate(this%mfx)
    if (allocated(this%mfy)) deallocate(this%mfy)
    if (allocated(this%u  )) deallocate(this%u  )
    if (allocated(this%v  )) deallocate(this%v  )
    if (allocated(this%cx )) deallocate(this%cx )
    if (allocated(this%cy )) deallocate(this%cy )

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

  subroutine adv_batch_accum_wind(this, mfx, mfy, u, v)

    class(adv_batch_type), intent(inout) :: this
    real(r8), intent(in) :: mfx(:,:,:)
    real(r8), intent(in) :: mfy(:,:,:)
    real(r8), intent(in) :: u  (:,:,:)
    real(r8), intent(in) :: v  (:,:,:)

    if (this%step == 0) then
      this%mfx  = mfx
      this%mfy  = mfy
      this%u    = u
      this%v    = v
      this%step = this%step + 1
    else if (this%step == this%nstep) then
      this%u    = this%u / this%nstep
      this%v    = this%v / this%nstep
      this%step = 0
    else
      this%mfx  = this%mfx + mfx
      this%mfy  = this%mfy + mfy
      this%u    = this%u   + u
      this%v    = this%v   + v
      this%step = this%step + 1
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