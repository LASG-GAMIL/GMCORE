module block_mod

  use mpi
  use flogger
  use namelist_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use reduced_types_mod
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

  !      |         |         |
  !  ____|_________|_________|_____
  !      |                   |
  !      |                   |
  !      |                   |
  !      |                   |
  !      |       BLOCK       |
  !      |                   |
  !      |                   |
  !      |                   |
  !  ____|___________________|_____
  !      |                   |
  !      |                   |

  type block_type
    integer id
    type(mesh_type) mesh
    type(state_type), allocatable :: state(:)
    type(static_type) static
    type(tend_type), allocatable :: tend(:)
    type(reduced_mesh_type), allocatable :: reduced_mesh(:)
    type(reduced_state_type), allocatable :: reduced_state(:)
    type(reduced_static_type), allocatable :: reduced_static(:)
    type(reduced_tend_type), allocatable :: reduced_tend(:)
    type(halo_type), allocatable :: halo(:)
  contains
    procedure :: init => block_init
    final :: block_final
  end type block_type

contains

  subroutine block_init(this, id, lon_halo_width, lat_halo_width, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    class(block_type), intent(inout) :: this
    integer, intent(in) :: id
    integer, intent(in) :: lon_halo_width
    integer, intent(in) :: lat_halo_width
    integer, intent(in) :: lon_ibeg
    integer, intent(in) :: lon_iend
    integer, intent(in) :: lat_ibeg
    integer, intent(in) :: lat_iend

    integer i

    this%id = id

    call this%mesh%init_from_parent(global_mesh, this%id, lon_halo_width, lat_halo_width, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    if (.not. allocated(this%state)) then
      select case (trim(time_scheme))
      case ('pc2', 'rk2')
        allocate(this%state(3))
        allocate(this%tend (3))
      case ('rk3', 'ssprk3')
        allocate(this%state(4))
        allocate(this%tend (4))
      case ('rk4')
        allocate(this%state(5))
        allocate(this%tend (5))
      case default
        if (id == 0) call log_error('Unknown time scheme ' // trim(time_scheme))
      end select
      do i = 1, size(this%state)
        call this%state(i)%init(this%mesh)
        call this%tend (i)%init(this%mesh)
      end do
      call this%static%init(this%mesh)
    end if

  end subroutine block_init

  subroutine block_final(this)

    type(block_type), intent(inout) :: this

    if (allocated(this%halo          )) deallocate(this%halo          )
    if (allocated(this%state         )) deallocate(this%state         )
    if (allocated(this%tend          )) deallocate(this%tend          )
    if (allocated(this%reduced_mesh  )) deallocate(this%reduced_mesh  )
    if (allocated(this%reduced_state )) deallocate(this%reduced_state )
    if (allocated(this%reduced_static)) deallocate(this%reduced_static)
    if (allocated(this%reduced_tend  )) deallocate(this%reduced_tend  )

  end subroutine block_final

end module block_mod
