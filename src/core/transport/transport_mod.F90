module transport_mod

  use mesh_type
  use allocator_mod

  implicit none

  type transport_batch_type
    integer :: nsteps = 0 ! Number of dynamic steps for one transport step
    integer :: step   = 0 ! Dynamic step counter
    real(r8) :: dt
    real(r8), allocatable, dimension(:,:,:) :: mfx
    real(r8), allocatable, dimension(:,:,:) :: mfy
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
  contains
    procedure :: init       => transport_batch_init
    procedure :: clear      => transport_batch_clear
    procedure :: accum_wind => transport_batch_accum_wind
    final :: transport_batch_final
  end type transport_batch_type

contains

  subroutine transport_batch_init(this, mesh, alert_key, dt)

    class(transport_batch_type), intent(inout) :: this
    type(mesh_type), intent(in) :: mesh
    character(*), intent(in) :: alert_key
    real(r8), intent(in) :: dt

    call this%clear()

    this%alert_key = alert_key
    this%dt = dt
    this%nsteps = dt / dt_in_seconds

    call allocate_array(mesh, this%mfx, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mfy, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u  , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v  , full_lon=.true., half_lat=.true., full_lev=.true.)

  end subroutine transport_batch_init

  subroutine transport_batch_clear(this)

    class(transport_batch_type), intent(inout) :: this

    if (allocated(this%mfx)) deallocate(this%mfx)
    if (allocated(this%mfy)) deallocate(this%mfy)
    if (allocated(this%u  )) deallocate(this%u  )
    if (allocated(this%v  )) deallocate(this%v  )

  end subroutine transport_batch_clear

  subroutine transport_batch_accum_wind(this, mfx, mfy, u, v)

    class(transport_batch_type), intent(inout) :: this
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
    else if (this%step == this%nsteps) then
      this%u    = this%u / this%nsteps
      this%v    = this%v / this%nsteps
      this%step = 0
    else
      this%mfx  = this%mfx + mfx
      this%mfy  = this%mfy + mfy
      this%u    = this%u   + u
      this%v    = this%v   + v
      this%step = this%step + 1
    end if

  end subroutine transport_batch_accum_wind

  subroutine transport_batch_final(this)

    type(transport_batch_type), intent(inout) :: this

    call this%clear()

  end subroutine transport_batch_final

end module transport_mod