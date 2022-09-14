module block_mod

  use mpi
  use flogger
  use namelist_mod
  use mesh_mod
  use state_mod
  use static_mod
  use tend_mod
  use adv_batch_mod
  use filter_types_mod
  use halo_mod
  use allocator_mod

  implicit none

  private

  public block_type
  public blocks
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
    type(adv_batch_type) adv_batch_mass ! For mass, potential temperature
    type(adv_batch_type) adv_batch_pv
    type(adv_batch_type), allocatable :: adv_batches(:)
    type(filter_type) big_filter
    type(filter_type) small_filter1
    type(filter_type) small_filter2
    type(halo_type), allocatable :: halo(:)
  contains
    procedure :: init_stage_1 => block_init_stage_1
    procedure :: init_stage_2 => block_init_stage_2
    procedure :: clear => block_clear
    final :: block_final
  end type block_type

  type(block_type), allocatable :: blocks(:)

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
    call this%big_filter%init(this%mesh, 'big_filter')
    call this%small_filter1%init(this%mesh, 'small_filter1')
    call this%small_filter2%init(this%mesh, 'small_filter2')

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
      case default
        if (this%id == 0) call log_error('Unknown time scheme ' // trim(time_scheme))
      end select
      do i = 1, size(this%state)
        call this%state(i)%init(this%mesh)
      end do
      do i = 1, size(this%tend)
        call this%tend(i)%init(this%mesh)
      end do
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
    call this%adv_batch_mass%clear()
    call this%adv_batch_pv  %clear()
    if (allocated(this%adv_batches)) then
      do i = 1, size(this%adv_batches)
        call this%adv_batches(i)%clear()
      end do
    end if
    do i = 1, size(this%halo)
      call this%halo(i)%clear()
    end do

    if (allocated(this%state)) deallocate(this%state)
    if (allocated(this%tend )) deallocate(this%tend )
    if (allocated(this%adv_batches)) deallocate(this%adv_batches)
    if (allocated(this%halo )) deallocate(this%halo )

  end subroutine block_clear

  subroutine block_final(this)

    type(block_type), intent(inout) :: this

    call this%clear()

  end subroutine block_final

end module block_mod
