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
    real(r8), allocatable, dimension(:,:) :: gzs
    real(r8), allocatable, dimension(:,:) :: dzsdlon
    real(r8), allocatable, dimension(:,:) :: dzsdlat
  contains
    procedure :: init => static_init
    procedure :: clear => static_clear
    procedure :: prepare => static_prepare
    final :: static_final
  end type static_type

contains

  subroutine static_init(this, mesh)

    class(static_type), intent(inout)         :: this
    type(mesh_type   ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%gzs    , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dzsdlon, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dzsdlat, full_lon=.true., full_lat=.true.)

  end subroutine static_init

  subroutine static_clear(this)

    class(static_type), intent(inout) :: this

    if (allocated(this%gzs    )) deallocate(this%gzs    )
    if (allocated(this%dzsdlon)) deallocate(this%dzsdlon)
    if (allocated(this%dzsdlat)) deallocate(this%dzsdlat)

  end subroutine static_clear

  subroutine static_prepare(this)

    class(static_type), intent(inout) :: this

    integer i, j

    if (nonhydrostatic) then
      do j = this%mesh%full_lat_ibeg_no_pole, this%mesh%full_lat_iend_no_pole
        do i = this%mesh%full_lon_ibeg, this%mesh%full_lon_iend
          this%dzsdlon(i,j) = (this%gzs(i+1,j) - this%gzs(i-1,j)) / g / this%mesh%de_lon(j) * 0.5_r8
#ifdef V_POLE
          this%dzsdlat(i,j) = (this%gzs(i,j+1) - this%gzs(i,j-1)) / g / (this%mesh%de_lat(j+1) + this%mesh%de_lat(j))
#else
          this%dzsdlat(i,j) = (this%gzs(i,j+1) - this%gzs(i,j-1)) / g / (this%mesh%de_lat(j-1) + this%mesh%de_lat(j))
#endif
        end do
      end do
    end if

  end subroutine static_prepare

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    call this%clear()

  end subroutine static_final

end module static_mod
