module tracer_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use allocator_mod

  implicit none

  private

  public tracer_type

  type tracer_type
    type(mesh_type), pointer :: mesh => null()
    character(30) :: name = ''
    real(r8), allocatable, dimension(:,:,:) :: q   ! Tracer density or mixing ratio or any other quantity
    real(r8), allocatable, dimension(:,:,:) :: mfx ! Tracer mass flux along x axis
    real(r8), allocatable, dimension(:,:,:) :: mfy ! Tracer mass flux along y axis
    ! Work arrays
    real(r8), allocatable, dimension(:,:,:) :: qx    ! Tracer mixing ratio due to advective operator along x axis
    real(r8), allocatable, dimension(:,:,:) :: qy    ! Tracer mixing ratio due to advective operator along y axis
    real(r8), allocatable, dimension(:,:,:) :: qxl   ! Tracer mixing ratio at left cell edge along x axis
    real(r8), allocatable, dimension(:,:,:) :: qyl   ! Tracer mixing ratio at left cell edge along y axis
    real(r8), allocatable, dimension(:,:,:) :: dqx   ! Tracer mixing ratio mismatch (or slope) at cell center along x axis
    real(r8), allocatable, dimension(:,:,:) :: dqy   ! Tracer mixing ratio mismatch (or slope) at cell center along y axis
    real(r8), allocatable, dimension(:,:,:) :: qx6   ! PPM mismatch at cell center along x axis
    real(r8), allocatable, dimension(:,:,:) :: qy6   ! PPM mismatch at cell center along y axis
  contains
    procedure :: init            => tracer_init
    procedure :: allocate_arrays => tracer_allocate_arrays
    procedure :: clear           => tracer_clear
    final :: tracer_final
  end type tracer_type

contains

  subroutine tracer_init(this, mesh, name)

    class(tracer_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: name

    call this%clear()

    this%mesh => mesh
    this%name =  name

  end subroutine tracer_init

  subroutine tracer_allocate_arrays(this)

    class(tracer_type), intent(inout) :: this

    call allocate_array(this%mesh, this%q  , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%mfx, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%mfy, full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%qx , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(this%mesh, this%qy , full_lon=.true., full_lat=.true., full_lev=.true.)
    select case (adv_scheme)
    case ('ffsl')
      if (ffsl_flux_type == 'ppm') then
        call allocate_array(this%mesh, this%qxl, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%mesh, this%qyl, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%mesh, this%dqx, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%mesh, this%dqy, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%mesh, this%qx6, full_lon=.true., full_lat=.true., full_lev=.true.)
        call allocate_array(this%mesh, this%qy6, full_lon=.true., full_lat=.true., full_lev=.true.)
      end if
    end select

  end subroutine tracer_allocate_arrays

  subroutine tracer_clear(this)

    class(tracer_type), intent(inout) :: this

    nullify(this%mesh)
    this%name = ''
    if (allocated(this%q   )) deallocate(this%q   )
    if (allocated(this%mfx )) deallocate(this%mfx )
    if (allocated(this%mfy )) deallocate(this%mfy )
    if (allocated(this%qx  )) deallocate(this%qx  )
    if (allocated(this%qy  )) deallocate(this%qy  )
    if (allocated(this%qxl )) deallocate(this%qxl )
    if (allocated(this%qyl )) deallocate(this%qyl )
    if (allocated(this%dqx )) deallocate(this%dqx )
    if (allocated(this%dqy )) deallocate(this%dqy )
    if (allocated(this%qx6 )) deallocate(this%qx6 )
    if (allocated(this%qy6 )) deallocate(this%qy6 )

  end subroutine tracer_clear

  subroutine tracer_final(this)

    type(tracer_type), intent(inout) :: this

    call this%clear()

  end subroutine tracer_final

end module tracer_mod
