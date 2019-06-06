module state_mod

  use mesh_mod

  implicit none

  private

  public state_type
  public state
  public create_states

  type state_type
    real, allocatable :: u (:,:)
    real, allocatable :: v (:,:)
    real, allocatable :: gd(:,:)
    real, allocatable :: pv(:,:)
    real, allocatable :: dpvdlon(:,:)
    real, allocatable :: dpvdlat(:,:)
  contains
    procedure :: init => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

  type(state_type), allocatable :: state(:)

contains

  subroutine create_states()

    integer i

    if (.not. allocated(state)) then
      allocate(state(0:2))
      do i = lbound(state, 1), ubound(state, 1)
        call state(i)%init()
      end do
    end if

  end subroutine create_states

  subroutine state_init(this)

    class(state_type), intent(inout) :: this

    allocate(this%u      (1-mesh%halo_width:mesh%num_half_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_full_lat+mesh%halo_width))
    allocate(this%v      (1-mesh%halo_width:mesh%num_full_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_half_lat+mesh%halo_width))
    allocate(this%gd     (1-mesh%halo_width:mesh%num_full_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_full_lat+mesh%halo_width))
    allocate(this%pv     (1-mesh%halo_width:mesh%num_half_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_half_lat+mesh%halo_width))
    allocate(this%dpvdlon(1-mesh%halo_width:mesh%num_full_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_half_lat+mesh%halo_width))
    allocate(this%dpvdlat(1-mesh%halo_width:mesh%num_half_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_full_lat+mesh%halo_width))

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%u))       deallocate(this%u)
    if (allocated(this%v))       deallocate(this%v)
    if (allocated(this%gd))      deallocate(this%gd)
    if (allocated(this%pv))      deallocate(this%pv)
    if (allocated(this%dpvdlon)) deallocate(this%dpvdlon)
    if (allocated(this%dpvdlat)) deallocate(this%dpvdlat)

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
