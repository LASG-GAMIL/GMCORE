module block_mod

  use mpi
  use flogger
  use namelist_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use filter_types_mod
  use halo_mod
  use allocator_mod

  implicit none

  private

  public block_type
  public global_mesh
  public mesh_type
  public state_type
  public static_type
  public tend_type

  type block_type
    integer id
    type(mesh_type) mesh
    type(state_type), allocatable :: state(:)
    type(static_type) static
    type(tend_type), allocatable :: tend(:)
    type(filter_type) filter
    type(halo_type), allocatable :: halo(:)
  contains
    procedure :: init_stage_1 => block_init_stage_1
    procedure :: init_stage_2 => block_init_stage_2
    procedure :: clear => block_clear
    final :: block_final
  end type block_type

contains

  subroutine block_init_stage_1(this, id, lon_halo_width, lat_halo_width, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    class(block_type), intent(inout) :: this
    integer, intent(in) :: id
    integer, intent(in) :: lon_halo_width
    integer, intent(in) :: lat_halo_width
    integer, intent(in) :: lon_ibeg
    integer, intent(in) :: lon_iend
    integer, intent(in) :: lat_ibeg
    integer, intent(in) :: lat_iend

    this%id = id

    call this%mesh%init_from_parent(global_mesh, this%id, lon_halo_width, lat_halo_width, lon_ibeg, lon_iend, lat_ibeg, lat_iend)
    call this%filter%init(this%mesh)

  end subroutine block_init_stage_1

  subroutine block_init_stage_2(this, lon_halo_width)

    class(block_type), intent(inout) :: this
    integer, intent(in) :: lon_halo_width

    integer i

    call this%mesh%reinit(lon_halo_width)

    if (.not. allocated(this%state)) then
      select case (trim(time_scheme))
      case ('euler')
        allocate(this%state(2))
        allocate(this%tend (2))
      case ('pc2', 'wrfrk3')
        allocate(this%state(3))
        allocate(this%tend (3))
      case ('rk3', 'ssprk3')
        allocate(this%state(4))
        allocate(this%tend (4))
      case ('rk4')
        allocate(this%state(5))
        allocate(this%tend (5))
      case ('N/A')
        allocate(this%state(1))
      case default
        if (this%id == 0) call log_error('Unknown time scheme ' // trim(time_scheme))
      end select
      do i = 1, size(this%state)
        call this%state(i)%init(this%mesh)
      end do
      if (allocated(this%tend)) then
        do i = 1, size(this%tend)
          call this%tend(i)%init(this%mesh)
        end do
      end if
      call this%static%init(this%mesh)
    end if

  end subroutine block_init_stage_2

  subroutine block_clear(this)

    class(block_type), intent(inout) :: this

    integer i

    do i = 1, size(this%state)
      call this%state(i)%clear()
    end do
    do i = 1, size(this%tend)
      call this%tend(i)%clear()
    end do
    do i = 1, size(this%halo)
      call this%halo(i)%clear()
    end do

    if (allocated(this%halo )) deallocate(this%halo )
    if (allocated(this%state)) deallocate(this%state)
    if (allocated(this%tend )) deallocate(this%tend )
    if (allocated(this%halo )) deallocate(this%halo )

  end subroutine block_clear

  subroutine block_final(this)

    type(block_type), intent(inout) :: this

    call this%clear()

  end subroutine block_final

end module block_mod
