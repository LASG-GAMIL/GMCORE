module static_mod

  use mesh_mod

  implicit none

  type static_type
    real, allocatable :: ghs(:,:)
  contains
    procedure :: init => static_init
    final :: static_final
  end type static_type

  type(static_type) static

contains

  subroutine create_static()

    call static%init()

  end subroutine create_static

  subroutine static_init(this)

    class(static_type), intent(inout) :: this

    if (allocated(this%ghs)) deallocate(this%ghs)

    allocate(this%ghs(1-mesh%halo_width:mesh%num_full_lon+mesh%halo_width,1-mesh%halo_width:mesh%num_full_lat+mesh%halo_width))

  end subroutine static_init

  subroutine static_final(this)

    type(static_type), intent(inout) :: this

    if (allocated(this%ghs)) deallocate(this%ghs)

  end subroutine static_final

end module static_mod