module refer_state_mod

  use flogger
  use const_mod
  use namelist_mod
  use allocator_mod
  use mesh_mod
  use static_mod
  use process_mod
  use parallel_mod
  use refer_state_wrf_mod
  use vert_coord_mod

  implicit none

  private

  type refer_profile_type
    integer num_lev
  end type refer_profile_type

  type refer_state_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:  ) :: phs
    real(r8), allocatable, dimension(:,:,:) :: ph
    real(r8), allocatable, dimension(:,:,:) :: gz
    real(r8), allocatable, dimension(:,:,:) :: t 
    real(r8), allocatable, dimension(:,:,:) :: pt
  contains
    procedure :: init  => refer_state_init
    procedure :: clear => refer_state_clear
    final :: refer_state_final
  end type refer_state_type

contains

  subroutine refer_state_init(this, mesh, static)

    class(refer_state_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh
    type(static_type), intent(in) :: static

    integer i, j, k

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%phs, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%ph , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%gz , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%t  , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pt , full_lon=.true., full_lat=.true., full_lev=.true.)

    select case (refer_state_scheme)
    case ('wrf')
      call refer_state_wrf_set_phs(mesh, static%gzs, this%phs)
    case default
      if (is_root_proc()) call log_error('Unknown refer_state_scheme ' // trim(refer_state_scheme) // '!')
    end select

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          this%ph(i,j,k) = vert_coord_calc_ph(k, this%phs(i,j))
        end do
      end do
    end do

  end subroutine refer_state_init

  subroutine refer_state_clear(this)

    class(refer_state_type), intent(inout) :: this

    if (allocated(this%phs)) deallocate(this%phs)
    if (allocated(this%ph )) deallocate(this%ph )
    if (allocated(this%gz )) deallocate(this%gz )
    if (allocated(this%t  )) deallocate(this%t  )
    if (allocated(this%pt )) deallocate(this%pt )

  end subroutine refer_state_clear

  subroutine refer_state_final(this)

    type(refer_state_type), intent(inout) :: this

    call this%clear()

  end subroutine refer_state_final

end module refer_state_mod
