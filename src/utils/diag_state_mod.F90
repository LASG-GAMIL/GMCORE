module diag_state_mod

  use datetime
  use const_mod
  use block_mod
  use namelist_mod
  use interp_mod
  use time_mod

  implicit none

  private

  public diag_state_init
  public diag_state_final
  public init_diag_state
  public diag_state

  integer, public, parameter :: height_levels   = 1
  integer, public, parameter :: pressure_levels = 2
  integer, public, parameter :: instance        = 1
  integer, public, parameter :: daily_mean      = 2
  integer, public, parameter :: monthly_mean    = 3
  integer, public, parameter :: yearly_mean     = 4

  type diag_state_type
    integer :: level_type = 0
    integer :: time_mean_type = 0
    integer :: time_step_counter = 0
    real(r8), allocatable, dimension(:    ) :: levels
    real(r8), allocatable, dimension(:,:,:) :: u
    real(r8), allocatable, dimension(:,:,:) :: v
    real(r8), allocatable, dimension(:,:,:) :: w
    real(r8), allocatable, dimension(:,:,:) :: z
    real(r8), allocatable, dimension(:,:,:) :: p
    real(r8), allocatable, dimension(:,:,:) :: pt
    real(r8), allocatable, dimension(:,:,:) :: t
    type(datetime_type) last_time
    real(r8), allocatable :: buf(:,:)
    real(r8), allocatable :: zm(:,:,:)
    real(r8), allocatable :: zm_lev(:,:,:)
  contains
    procedure :: is_init => diag_state_is_init
    procedure :: init_height_levels => diag_state_init_height_levels
    procedure :: init_pressure_levels => diag_state_init_pressure_levels
    procedure :: run => diag_state_run
    procedure :: clear => diag_state_clear
    procedure :: zeros => diag_state_zeros
    final :: finalize_diag_state
  end type diag_state_type

  interface
    subroutine init_diag_state_interface(blocks)
      import block_type
      type(block_type), intent(in) :: blocks(:)
    end subroutine init_diag_state_interface
  end interface

  type(diag_state_type), allocatable :: diag_state(:)
  procedure(init_diag_state_interface), pointer :: init_diag_state

contains

  subroutine diag_state_init(blocks)

    type(block_type), intent(in) :: blocks(:)

    allocate(diag_state(size(blocks)))

    if (associated(init_diag_state)) call init_diag_state(blocks)

  end subroutine diag_state_init

  subroutine diag_state_final()

    if (allocated(diag_state)) deallocate(diag_state)

  end subroutine diag_state_final

  logical function diag_state_is_init(this) result(res)

    class(diag_state_type), intent(in) :: this

    res = this%level_type /= 0

  end function diag_state_is_init

  subroutine diag_state_init_height_levels(this, block, z, time_mean_type)

    class(diag_state_type), intent(inout) :: this
    type(block_type), intent(in) :: block
    real(r8), intent(in) :: z(:) ! m
    integer, intent(in) :: time_mean_type

    integer is, ie, js, je

    is = block%mesh%full_lon_lb
    ie = block%mesh%full_lon_ub
    js = block%mesh%full_lat_lb
    je = block%mesh%full_lat_ub

    call this%clear()
    this%level_type = height_levels
    this%time_mean_type = time_mean_type
    allocate(this%levels(size(z)))
    allocate(this%u  (is:ie,js:je,size(z)))
    allocate(this%v  (is:ie,js:je,size(z)))
    allocate(this%pt (is:ie,js:je,size(z)))
    allocate(this%t  (is:ie,js:je,size(z)))
    allocate(this%p  (is:ie,js:je,size(z)))
    if (nonhydrostatic) then
      allocate(this%w(is:ie,js:je,size(z)))
    end if
    allocate(this%buf   (is:ie,js:je))
    allocate(this%zm    (is:ie,js:je,block%mesh%full_lev_lb:block%mesh%full_lev_ub))
    allocate(this%zm_lev(is:ie,js:je,block%mesh%half_lev_lb:block%mesh%half_lev_ub))

    this%levels = z
    this%last_time = curr_time

  end subroutine diag_state_init_height_levels

  subroutine diag_state_init_pressure_levels(this, block, p, time_mean_type)

    class(diag_state_type), intent(inout) :: this
    type(block_type), intent(in) :: block
    real(r8), intent(in) :: p(:) ! Pa
    integer, intent(in) :: time_mean_type

    integer is, ie, js, je

    is = block%mesh%full_lon_lb
    ie = block%mesh%full_lon_ub
    js = block%mesh%full_lat_lb
    je = block%mesh%full_lat_ub

    call this%clear()
    this%level_type = pressure_levels
    this%time_mean_type = time_mean_type
    allocate(this%levels(size(p)))
    allocate(this%u  (is:ie,js:je,size(p)))
    allocate(this%v  (is:ie,js:je,size(p)))
    allocate(this%pt (is:ie,js:je,size(p)))
    allocate(this%t  (is:ie,js:je,size(p)))
    allocate(this%z  (is:ie,js:je,size(p)))
    if (nonhydrostatic) then
      allocate(this%w(is:ie,js:je,size(p)))
    end if
    allocate(this%buf(is:ie,js:je))
    allocate(this%zm    (is:ie,js:je,block%mesh%full_lev_lb:block%mesh%full_lev_ub))
    allocate(this%zm_lev(is:ie,js:je,block%mesh%half_lev_lb:block%mesh%half_lev_ub))

    this%levels = p
    this%last_time = curr_time

  end subroutine diag_state_init_pressure_levels

  subroutine diag_state_clear(this)

    class(diag_state_type), intent(inout) :: this

    if (allocated(this%levels)) deallocate(this%levels)
    if (allocated(this%u     )) deallocate(this%u     )
    if (allocated(this%v     )) deallocate(this%v     )
    if (allocated(this%w     )) deallocate(this%w     )
    if (allocated(this%z     )) deallocate(this%z     )
    if (allocated(this%p     )) deallocate(this%p     )
    if (allocated(this%pt    )) deallocate(this%pt    )
    if (allocated(this%t     )) deallocate(this%t     )
    if (allocated(this%buf   )) deallocate(this%buf   )
    if (allocated(this%zm    )) deallocate(this%zm    )
    if (allocated(this%zm_lev)) deallocate(this%zm_lev)

  end subroutine diag_state_clear

  subroutine diag_state_zeros(this)

    class(diag_state_type), intent(inout) :: this

    if (allocated(this%u )) this%u  = 0
    if (allocated(this%v )) this%v  = 0
    if (allocated(this%w )) this%w  = 0
    if (allocated(this%pt)) this%pt = 0
    if (allocated(this%t )) this%t  = 0
    if (allocated(this%z )) this%z  = 0
    if (allocated(this%p )) this%p  = 0

  end subroutine diag_state_zeros

  subroutine finalize_diag_state(this)

    type(diag_state_type), intent(inout) :: this

    call this%clear()

  end subroutine finalize_diag_state

  subroutine diag_state_run(this, state)

    class(diag_state_type), intent(inout) :: this
    type(state_type), intent(in) :: state

    type(timedelta_type) dt
    integer, parameter :: calendar = datetime_noleap_calendar
    integer k

    select case (this%time_mean_type)
    case (instance)
      dt = create_timedelta(days=0)
    case (daily_mean)
      dt = create_timedelta(days=1)
    case (monthly_mean)
      dt = create_timedelta(days=days_of_month(this%last_time%year, this%last_time%month, calendar))
    case (yearly_mean)
      dt = create_timedelta(days=days_of_year(this%last_time%year, calendar))
    end select

    if (curr_time - this%last_time >= dt) then
      call this%zeros()
    end if

    select case (this%level_type)
    case (height_levels)
      this%zm = state%gz / g
      this%zm_lev = state%gz_lev / g
      do k = 1, size(this%levels)
        call interp_lon_edge_to_height_level(state%mesh, this%zm, state%u, this%levels(k), this%buf)
        this%u(:,:,k) = this%u(:,:,k) + this%buf
        call interp_lat_edge_to_height_level(state%mesh, this%zm, state%v, this%levels(k), this%buf)
        this%v(:,:,k) = this%v(:,:,k) + this%buf
        call interp_cell_to_height_level(state%mesh, this%zm, state%pt, this%levels(k), this%buf)
        this%pt(:,:,k) = this%pt(:,:,k) + this%buf
        call interp_cell_to_height_level(state%mesh, this%zm, state%t, this%levels(k), this%buf)
        this%t(:,:,k) = this%t(:,:,k) + this%buf
        call interp_lev_edge_to_height_level(state%mesh, this%zm_lev, state%p_lev, this%levels(k), this%buf)
        this%p(:,:,k) = this%p(:,:,k) + this%buf
        if (nonhydrostatic) then
          call interp_lev_edge_to_height_level(state%mesh, this%zm_lev, state%w, this%levels(k), this%buf)
          this%w(:,:,k) = this%w(:,:,k) + this%buf
        end if
      end do
    case (pressure_levels)
      this%zm_lev = state%gz_lev / g
      do k = 1, size(this%levels)
        call interp_lon_edge_to_pressure_level(state%mesh, state%p, state%u, this%levels(k), this%buf)
        this%u(:,:,k) = this%u(:,:,k) + this%buf
        call interp_lat_edge_to_pressure_level(state%mesh, state%p, state%v, this%levels(k), this%buf)
        this%v(:,:,k) = this%v(:,:,k) + this%buf
        call interp_cell_to_pressure_level(state%mesh, state%p, state%pt, this%levels(k), this%buf)
        this%pt(:,:,k) = this%pt(:,:,k) + this%buf
        call interp_cell_to_pressure_level(state%mesh, state%p, state%t, this%levels(k), this%buf)
        this%t(:,:,k) = this%t(:,:,k) + this%buf
        call interp_lev_edge_to_pressure_level(state%mesh, state%p, this%zm_lev, this%levels(k), this%buf)
        this%z(:,:,k) = this%z(:,:,k) + this%buf
        if (nonhydrostatic) then
          call interp_lev_edge_to_pressure_level(state%mesh, state%p_lev, state%w, this%levels(k), this%buf)
          this%w(:,:,k) = this%w(:,:,k) + this%buf
        end if
      end do
    end select

    this%time_step_counter = this%time_step_counter + 1
    if (curr_time - this%last_time >= dt) then
      this%u  = this%u  / this%time_step_counter
      this%v  = this%v  / this%time_step_counter
      this%pt = this%pt / this%time_step_counter
      this%t  = this%t  / this%time_step_counter
      select case (this%level_type)
      case (height_levels)
        this%p = this%p / this%time_step_counter
      case (pressure_levels)
        this%z = this%z / this%time_step_counter
      end select
      if (nonhydrostatic) then
        this%w = this%w / this%time_step_counter
      end if
      this%last_time = curr_time
      this%time_step_counter = 0
    end if

  end subroutine diag_state_run

end module diag_state_mod
