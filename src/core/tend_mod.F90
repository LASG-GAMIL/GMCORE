module tend_mod

  use mesh_mod

  implicit none

  private

  public tend_type
  public tend
  public create_tends

  type tend_type
    real, allocatable :: du (:,:)
    real, allocatable :: dv (:,:)
    real, allocatable :: dgd(:,:)
  contains
    procedure :: init => tend_init
    procedure :: clear => tend_clear
    final :: tend_final
  end type tend_type

  type(tend_type), allocatable :: tend(:)

contains

  subroutine create_tends()

    integer i

    if (.not. allocated(tend)) then
      allocate(tend(0:2))
      do i = lbound(tend, 1), ubound(tend, 1)
        call tend(i)%init()
      end do
    end if

  end subroutine create_tends

  subroutine tend_init(this)

    class(tend_type), intent(inout) :: this

    allocate(this%du (1-mesh%halo_width:mesh%num_half_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_full_lat+mesh%halo_width))
    allocate(this%dv (1-mesh%halo_width:mesh%num_full_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_half_lat+mesh%halo_width))
    allocate(this%dgd(1-mesh%halo_width:mesh%num_full_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_full_lat+mesh%halo_width))

  end subroutine tend_init

  subroutine tend_clear(this)

    class(tend_type), intent(inout) :: this

    if (allocated(this%du))  deallocate(this%du)
    if (allocated(this%dv))  deallocate(this%dv)
    if (allocated(this%dgd)) deallocate(this%dgd)

  end subroutine tend_clear

  subroutine tend_final(this)

    type(tend_type), intent(inout) :: this

    call this%clear()

  end subroutine tend_final

end module tend_mod