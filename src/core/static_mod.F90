module static_mod

  use const_mod
  use mesh_mod
  use allocator_mod
  use namelist_mod

  implicit none

  private

  public static_type

  type static_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:) :: landmask
    real(r8), allocatable, dimension(:,:) :: gzs
    real(r8), allocatable, dimension(:,:) :: zs_std
    real(r8), allocatable, dimension(:,:) :: dzsdlon
    real(r8), allocatable, dimension(:,:) :: dzsdlat
  contains
    procedure :: init => static_init
    procedure :: clear => static_clear
    final :: static_final
  end type static_type

contains

  subroutine static_init(this, mesh)

    class(static_type), intent(inout)         :: this
    type(mesh_type   ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%landmask, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%gzs     , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%zs_std  , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dzsdlon , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dzsdlat , full_lon=.true., half_lat=.true.)

  end subroutine static_init

  subroutine static_clear(this)

    class(static_type), intent(inout) :: this

    if (allocated(this%landmask)) deallocate(this%landmask)
    if (allocated(this%gzs     )) deallocate(this%gzs     )
    if (allocated(this%zs_std  )) deallocate(this%zs_std  )
    if (allocated(this%dzsdlon )) deallocate(this%dzsdlon )
    if (allocated(this%dzsdlat )) deallocate(this%dzsdlat )

  end subroutine static_clear

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

end module static_mod
