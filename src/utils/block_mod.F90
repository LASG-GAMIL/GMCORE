module block_mod

  use mpi
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
    ! Work arrays
    real(r8), allocatable :: latlon_damp_cell_gx  (:,:,:)
    real(r8), allocatable :: latlon_damp_cell_gy  (:,:,:)
    real(r8), allocatable :: latlon_damp_cell_dfdx(:,:,:)
    real(r8), allocatable :: latlon_damp_cell_dfdy(:,:,:)
    real(r8), allocatable :: latlon_damp_lon_gx   (:,:,:)
    real(r8), allocatable :: latlon_damp_lon_gy   (:,:,:)
    real(r8), allocatable :: latlon_damp_lon_dfdx (:,:,:)
    real(r8), allocatable :: latlon_damp_lon_dfdy (:,:,:)
    real(r8), allocatable :: latlon_damp_lat_gx   (:,:,:)
    real(r8), allocatable :: latlon_damp_lat_gy   (:,:,:)
    real(r8), allocatable :: latlon_damp_lat_dfdx (:,:,:)
    real(r8), allocatable :: latlon_damp_lat_dfdy (:,:,:)
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
      case ('debug', 'pc2', 'pc2+fb', 'rk2')
        allocate(this%state(3))
        allocate(this%tend (3))
      case ('rk3')
        allocate(this%state(4))
        allocate(this%tend (4))
      case ('rk4')
        allocate(this%state(5))
        allocate(this%tend (5))
      end select
      do i = 1, size(this%state)
        call this%state(i)%init(this%mesh)
        call this%tend (i)%init(this%mesh)
      end do
      call this%static%init(this%mesh)
    end if

    call allocate_array(this%mesh, this%latlon_damp_cell_gx  , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_cell_gy  , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_cell_dfdy, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_cell_dfdx, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_lon_gx   , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_lon_gy   , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_lon_dfdy , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_lon_dfdx , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_lat_gx   , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_lat_gy   , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_lat_dfdy , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%latlon_damp_lat_dfdx , full_lon=.true., half_lat=.true., full_lev=.true.)

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

    if (allocated(this%latlon_damp_cell_gx  )) deallocate(this%latlon_damp_cell_gx  )
    if (allocated(this%latlon_damp_cell_gy  )) deallocate(this%latlon_damp_cell_gy  )
    if (allocated(this%latlon_damp_cell_dfdx)) deallocate(this%latlon_damp_cell_dfdx)
    if (allocated(this%latlon_damp_cell_dfdy)) deallocate(this%latlon_damp_cell_dfdy)
    if (allocated(this%latlon_damp_lon_gx   )) deallocate(this%latlon_damp_lon_gx   )
    if (allocated(this%latlon_damp_lon_gy   )) deallocate(this%latlon_damp_lon_gy   )
    if (allocated(this%latlon_damp_lon_dfdx )) deallocate(this%latlon_damp_lon_dfdx )
    if (allocated(this%latlon_damp_lon_dfdy )) deallocate(this%latlon_damp_lon_dfdy )
    if (allocated(this%latlon_damp_lat_gx   )) deallocate(this%latlon_damp_lat_gx   )
    if (allocated(this%latlon_damp_lat_gy   )) deallocate(this%latlon_damp_lat_gy   )
    if (allocated(this%latlon_damp_lat_dfdx )) deallocate(this%latlon_damp_lat_dfdx )
    if (allocated(this%latlon_damp_lat_dfdy )) deallocate(this%latlon_damp_lat_dfdy )

  end subroutine block_final

end module block_mod
