module state_mod

  use const_mod
  use mesh_mod
  use allocator_mod

  implicit none

  private

  public state_type
  public states
  public create_states

  type state_type
    type(mesh_type), pointer :: mesh => null()
    real(real_kind), allocatable, dimension(:,:) :: u
    real(real_kind), allocatable, dimension(:,:) :: v
    real(real_kind), allocatable, dimension(:,:) :: gd
    real(real_kind), allocatable, dimension(:,:) :: pv
    real(real_kind), allocatable, dimension(:,:) :: mass_vertex
    real(real_kind), allocatable, dimension(:,:) :: mass_lon
    real(real_kind), allocatable, dimension(:,:) :: mass_lat
    real(real_kind), allocatable, dimension(:,:) :: mass_flux_lon_n
    real(real_kind), allocatable, dimension(:,:) :: mass_flux_lat_n
    real(real_kind), allocatable, dimension(:,:) :: mass_flux_lon_t
    real(real_kind), allocatable, dimension(:,:) :: mass_flux_lat_t
    real(real_kind), allocatable, dimension(:,:) :: pv_lon
    real(real_kind), allocatable, dimension(:,:) :: pv_lat
    real(real_kind), allocatable, dimension(:,:) :: ke_cell
    real(real_kind) total_mass
    real(real_kind) total_energy
    real(real_kind) total_absolute_vorticity
    real(real_kind) total_potential_enstrophy
  contains
    procedure :: init => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

  type(state_type), allocatable, target :: states(:)

contains

  subroutine create_states()

    integer i

    if (.not. allocated(states)) then
      allocate(states(0:2))
      do i = lbound(states, 1), ubound(states, 1)
        call states(i)%init(mesh)
      end do
    end if

  end subroutine create_states

  subroutine state_init(this, mesh)

    class(state_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%u              , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%v              , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%gd             , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv             , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mass_vertex    , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mass_lon       , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%mass_lat       , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mass_flux_lon_n, half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%mass_flux_lat_n, full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mass_flux_lon_t, full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%mass_flux_lat_t, half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv_lon         , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv_lat         , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%ke_cell        , full_lon=.true., full_lat=.true.)

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%u              )) deallocate(this%u              )
    if (allocated(this%v              )) deallocate(this%v              )
    if (allocated(this%gd             )) deallocate(this%gd             )
    if (allocated(this%pv             )) deallocate(this%pv             )
    if (allocated(this%mass_vertex    )) deallocate(this%mass_vertex    )
    if (allocated(this%mass_lon       )) deallocate(this%mass_lon       )
    if (allocated(this%mass_lat       )) deallocate(this%mass_lat       )
    if (allocated(this%mass_flux_lon_n)) deallocate(this%mass_flux_lon_n)
    if (allocated(this%mass_flux_lat_n)) deallocate(this%mass_flux_lat_n)
    if (allocated(this%mass_flux_lon_t)) deallocate(this%mass_flux_lon_t)
    if (allocated(this%mass_flux_lat_t)) deallocate(this%mass_flux_lat_t)
    if (allocated(this%pv_lon         )) deallocate(this%pv_lon         )
    if (allocated(this%pv_lat         )) deallocate(this%pv_lat         )
    if (allocated(this%ke_cell        )) deallocate(this%ke_cell        )

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
