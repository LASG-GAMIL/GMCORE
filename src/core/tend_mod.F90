module tend_mod

  use const_mod
  use mesh_mod
  use allocator_mod

  implicit none

  private

  public tend_type
  public tend
  public create_tends

  type tend_type
    type(mesh_type), pointer :: mesh => null()
    real(real_kind), allocatable :: du (:,:)
    real(real_kind), allocatable :: dv (:,:)
    real(real_kind), allocatable :: dgd(:,:)
    ! Individual tendencies
    real(real_kind), allocatable :: coriolis_u(:,:)
    real(real_kind), allocatable :: coriolis_v(:,:)
    real(real_kind), allocatable :: grad_energy_u(:,:)
    real(real_kind), allocatable :: grad_energy_v(:,:)
    real(real_kind), allocatable :: div_mass_flux(:,:)
    ! Derived variables
    real(real_kind), allocatable :: pv_lon_edge(:,:)
    real(real_kind), allocatable :: pv_lat_edge(:,:)
    real(real_kind), allocatable :: dpvdlon(:,:)
    real(real_kind), allocatable :: dpvdlat(:,:)
    real(real_kind), allocatable :: ke_cell(:,:)
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
        call tend(i)%init(mesh)
      end do
    end if

  end subroutine create_tends

  subroutine tend_init(this, mesh)

    class(tend_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh

    call this%clear()

    this%mesh => mesh

#ifdef STAGGER_V_ON_POLE
    call allocate_array(mesh, this%du           , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dv           , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dgd          , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%coriolis_u   , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%coriolis_v   , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%grad_energy_u, half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%grad_energy_v, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%div_mass_flux, full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%pv_lon_edge  , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%pv_lat_edge  , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dpvdlon      , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dpvdlat      , half_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%ke_cell      , full_lon=.true., half_lat=.true.)
#else
    call allocate_array(mesh, this%du           , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%dv           , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dgd          , full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%coriolis_u   , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%coriolis_v   , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%grad_energy_u, half_lon=.true., full_lon=.true.)
    call allocate_array(mesh, this%grad_energy_v, full_lon=.true., half_lon=.true.)
    call allocate_array(mesh, this%div_mass_flux, full_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv_lon_edge  , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%pv_lat_edge  , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dpvdlon      , full_lon=.true., half_lat=.true.)
    call allocate_array(mesh, this%dpvdlat      , half_lon=.true., full_lat=.true.)
    call allocate_array(mesh, this%ke_cell      , full_lon=.true., full_lat=.true.)
#endif

  end subroutine tend_init

  subroutine tend_clear(this)

    class(tend_type), intent(inout) :: this

    if (allocated(this%du           )) deallocate(this%du           )
    if (allocated(this%dv           )) deallocate(this%dv           )
    if (allocated(this%dgd          )) deallocate(this%dgd          )
    if (allocated(this%coriolis_u   )) deallocate(this%coriolis_u   )
    if (allocated(this%coriolis_v   )) deallocate(this%coriolis_v   )
    if (allocated(this%grad_energy_u)) deallocate(this%grad_energy_u)
    if (allocated(this%grad_energy_v)) deallocate(this%grad_energy_v)
    if (allocated(this%div_mass_flux)) deallocate(this%div_mass_flux)
    if (allocated(this%pv_lon_edge  )) deallocate(this%pv_lon_edge  )
    if (allocated(this%pv_lat_edge  )) deallocate(this%pv_lat_edge  )
    if (allocated(this%dpvdlon      )) deallocate(this%dpvdlon      )
    if (allocated(this%dpvdlat      )) deallocate(this%dpvdlat      )
    if (allocated(this%ke_cell      )) deallocate(this%ke_cell      )

  end subroutine tend_clear

  subroutine tend_final(this)

    type(tend_type), intent(inout) :: this

    call this%clear()

  end subroutine tend_final

end module tend_mod