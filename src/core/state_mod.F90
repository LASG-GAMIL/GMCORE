module state_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use allocator_mod

  implicit none

  private

  public state_type

  ! NOTE:
  !   Variables with '_lon', '_lat' and '_lev' are on the half grids on the corresponding direction,
  !   and '_p' indicates that the variable is perturbed.
  type state_type
    type(mesh_type), pointer :: mesh => null()
    ! For nesting
    integer :: id = 0
    type(state_type), pointer :: parent => null()
    real(r8), allocatable, dimension(:,:,:) :: u_lon             ! Zonal wind speed (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: v_lon
    real(r8), allocatable, dimension(:,:,:) :: u_lat
    real(r8), allocatable, dimension(:,:,:) :: v_lat             ! Meridional wind speed (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: u_f               ! Zonal wind speed (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: v_f               ! Meridional wind speed (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: we_lev            ! Vertical coordinate speed multiplied by 𝛛π/𝛛η
    real(r8), allocatable, dimension(:,:,:) :: we_lev_lon        ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: we_lev_lat        ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: gz                ! Geopotential (m2 s-2)
    real(r8), allocatable, dimension(:,:,:) :: gz_f              ! Geopotential (m2 s-2)
    real(r8), allocatable, dimension(:,:,:) :: gz_lev            ! Geopotential height on half levels (m2 s-2)
    real(r8), allocatable, dimension(:,:,:) :: m                 ! Mass
    real(r8), allocatable, dimension(:,:,:) :: m_vtx             ! Mass on vertex
    real(r8), allocatable, dimension(:,:,:) :: m_lon             ! Mass on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: m_lat             ! Mass on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: m_lev             ! Mass on half levels
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_n          ! Normal mass flux on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_n          ! Normal mass flux on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_t          ! Tangient mass flux on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_t          ! Tangient mass flux on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: pv                ! Potential vorticity
    real(r8), allocatable, dimension(:,:,:) :: pv_lon            ! Potential vorticity on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: pv_lat            ! Potential vorticity on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: ke                ! Kinetic energy
    real(r8), allocatable, dimension(:,:,:) :: pt                ! Potential temperature
    real(r8), allocatable, dimension(:,:,:) :: pt_f              ! Potential temperature
    real(r8), allocatable, dimension(:,:,:) :: ptf_lon           ! Potential temperature on the zonal edge
    real(r8), allocatable, dimension(:,:,:) :: ptf_lat           ! Potential temperature on the merdional edge
    real(r8), allocatable, dimension(:,:,:) :: ptf_lev           ! Potential temperature on the vertical edge
    real(r8), allocatable, dimension(:,:,:) :: t                 ! Temperature
    real(r8), allocatable, dimension(:,:,:) :: ph                ! Hydrostatic pressure on full levels
    real(r8), allocatable, dimension(:,:,:) :: ph_lev            ! Hydrostatic pressure on half levels
    real(r8), allocatable, dimension(:,:  ) :: phs               ! Surface hydrostatic pressure
    real(r8), allocatable, dimension(:,:  ) :: phs_f             ! Surface hydrostatic pressure
    real(r8), allocatable, dimension(:,:,:) :: div               ! Divergence (s-1)
    real(r8), allocatable, dimension(:,:,:) :: div2              ! Laplacian of divergence (s-1)
    real(r8), allocatable, dimension(:,:,:) :: vor               ! Vorticity (s-1)
    ! Nonhydrostatic variables
    real(r8), allocatable, dimension(:,:,:) :: we
    real(r8), allocatable, dimension(:,:,:) :: w                 ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: w_lev             ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: w_lev_lon         ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: w_lev_lat         ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: gz_lev_lon        ! Geopotential
    real(r8), allocatable, dimension(:,:,:) :: gz_lev_lat        ! Geopotential
    real(r8), allocatable, dimension(:,:,:) :: rhod              ! Dry air density
    real(r8), allocatable, dimension(:,:,:) :: rhod_lon          ! Dry air density
    real(r8), allocatable, dimension(:,:,:) :: rhod_lat          ! Dry air density
    real(r8), pointer    , dimension(:,:,:) :: p                 ! Pressure on full levels
    real(r8), pointer    , dimension(:,:,:) :: p_lev             ! Pressure on half levels
    real(r8), allocatable, dimension(:,:,:) :: p_lev_lon         ! Pressure on half levels
    real(r8), allocatable, dimension(:,:,:) :: p_lev_lat         ! Pressure on half levels
    real(r8), allocatable, dimension(:,:,:) :: u_lev_lon
    real(r8), allocatable, dimension(:,:,:) :: v_lev_lat
    real(r8), allocatable, dimension(:,:,:) :: mf_lev_lon_n      ! Mass flux on zonal edge and half level
    real(r8), allocatable, dimension(:,:,:) :: mf_lev_lat_n      ! Mass flux on merdional edge and half level
    ! Smagorinsky damping variables
    real(r8), allocatable, dimension(:,:,:) :: tension_h         ! tension strain
    real(r8), allocatable, dimension(:,:,:) :: shear_h           ! shear strain on vertex
    real(r8), allocatable, dimension(:,:,:) :: kmh               ! nonlinear diffusion coef
    real(r8), allocatable, dimension(:,:,:) :: kmh_vtx           ! nonlinear diffusion coef on vertex
    real(r8), allocatable, dimension(:,:,:) :: kmh_lon           ! nonlinear diffusion coef on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: kmh_lat           ! nonlinear diffusion coef on meridional edge
    ! Tracers
    real(r8), allocatable, dimension(:,:,:,:) :: q               ! Tracer mass mixing ratio (kg kg-1)
    real(r8), allocatable, dimension(:,:,:,:) :: qmf_lon         ! Tracer mass flux along zonal direction
    real(r8), allocatable, dimension(:,:,:,:) :: qmf_lat         ! Tracer mass flux along meridional direction
    real(r8), allocatable, dimension(:,:,:,:) :: qmf_lev         ! Tracer mass flux along vertical direction
    ! Total diagnostics
    real(r8) tm
    real(r8) te
    real(r8) tpe
    real(r8) tav
  contains
    procedure :: init => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

contains

  subroutine state_init(this, mesh)

    class(state_type), intent(inout), target :: this
    type(mesh_type  ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%u_lon            , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v_lon            , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u_f              , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u_lat            , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v_lat            , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v_f              , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%we_lev           , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%we_lev_lon       , half_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%we_lev_lat       , full_lon=.true., half_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%gz               , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%gz_f             , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%gz_lev           , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%m                , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_vtx            , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_lon            , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_lat            , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_lev            , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%mf_lon_n         , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lon_t         , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lat_n         , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lat_t         , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv               , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv_lon           , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv_lat           , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ke               , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pt               , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pt_f             , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ptf_lon          , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ptf_lat          , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ptf_lev          , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%t                , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ph               , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ph_lev           , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%phs              , full_lon=.true., full_lat=.true.                 )
    call allocate_array(mesh, this%phs_f            , full_lon=.true., full_lat=.true.                 )
    call allocate_array(mesh, this%div              , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%vor              , half_lon=.true., half_lat=.true., full_lev=.true.)

    if (nonhydrostatic) then
      call allocate_array(mesh, this%we             , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%w              , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%w_lev          , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%w_lev_lon      , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%w_lev_lat      , full_lon=.true., half_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%gz_lev_lon     , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%gz_lev_lat     , full_lon=.true., half_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%rhod           , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%rhod_lon       , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%rhod_lat       , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%p              , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%p_lev          , full_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%p_lev_lon      , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%p_lev_lat      , full_lon=.true., half_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%u_lev_lon      , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%v_lev_lat      , full_lon=.true., half_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%mf_lev_lon_n   , half_lon=.true., full_lat=.true., half_lev=.true.)
      call allocate_array(mesh, this%mf_lev_lat_n   , full_lon=.true., half_lat=.true., half_lev=.true.)
    else
      this%p     => this%ph
      this%p_lev => this%ph_lev
    end if

    if (div_damp_order == 4) then
      call allocate_array(mesh, this%div2         , full_lon=.true., full_lat=.true., full_lev=.true.)
    end if

    if (use_smag_damp) then
      call allocate_array(mesh, this%tension_h    , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%shear_h      , half_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%kmh          , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%kmh_vtx      , half_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%kmh_lon      , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%kmh_lat      , full_lon=.true., half_lat=.true., full_lev=.true.)
    end if

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%u_lon            )) deallocate(this%u_lon            )
    if (allocated(this%v_lon            )) deallocate(this%v_lon            )
    if (allocated(this%u_lat            )) deallocate(this%u_lat            )
    if (allocated(this%v_lat            )) deallocate(this%v_lat            )
    if (allocated(this%u_f              )) deallocate(this%u_f              )
    if (allocated(this%v_f              )) deallocate(this%v_f              )
    if (allocated(this%we_lev           )) deallocate(this%we_lev           )
    if (allocated(this%we_lev_lon       )) deallocate(this%we_lev_lon       )
    if (allocated(this%we_lev_lat       )) deallocate(this%we_lev_lat       )
    if (allocated(this%gz               )) deallocate(this%gz               )
    if (allocated(this%gz_f             )) deallocate(this%gz_f             )
    if (allocated(this%gz_lev           )) deallocate(this%gz_lev           )
    if (allocated(this%m                )) deallocate(this%m                )
    if (allocated(this%m_vtx            )) deallocate(this%m_vtx            )
    if (allocated(this%m_lon            )) deallocate(this%m_lon            )
    if (allocated(this%m_lat            )) deallocate(this%m_lat            )
    if (allocated(this%mf_lon_n         )) deallocate(this%mf_lon_n         )
    if (allocated(this%mf_lat_n         )) deallocate(this%mf_lat_n         )
    if (allocated(this%mf_lat_t         )) deallocate(this%mf_lat_t         )
    if (allocated(this%mf_lon_t         )) deallocate(this%mf_lon_t         )
    if (allocated(this%pv               )) deallocate(this%pv               )
    if (allocated(this%pv_lon           )) deallocate(this%pv_lon           )
    if (allocated(this%pv_lat           )) deallocate(this%pv_lat           )
    if (allocated(this%ke               )) deallocate(this%ke               )
    if (allocated(this%pt               )) deallocate(this%pt               )
    if (allocated(this%pt_f             )) deallocate(this%pt_f             )
    if (allocated(this%ptf_lon          )) deallocate(this%ptf_lon          )
    if (allocated(this%ptf_lat          )) deallocate(this%ptf_lat          )
    if (allocated(this%ptf_lev          )) deallocate(this%ptf_lev          )
    if (allocated(this%t                )) deallocate(this%t                )
    if (allocated(this%ph               )) deallocate(this%ph               )
    if (allocated(this%ph_lev           )) deallocate(this%ph_lev           )
    if (allocated(this%phs              )) deallocate(this%phs              )
    if (allocated(this%phs_f            )) deallocate(this%phs_f            )
    if (allocated(this%div              )) deallocate(this%div              )
    if (allocated(this%div2             )) deallocate(this%div2             )
    if (allocated(this%vor              )) deallocate(this%vor              )

    if (allocated(this%m_lev            )) deallocate(this%m_lev            )
    if (allocated(this%we               )) deallocate(this%we               )
    if (allocated(this%w                )) deallocate(this%w                )
    if (allocated(this%w_lev            )) deallocate(this%w_lev            )
    if (allocated(this%w_lev_lon        )) deallocate(this%w_lev_lon        )
    if (allocated(this%w_lev_lat        )) deallocate(this%w_lev_lat        )
    if (allocated(this%gz_lev_lon       )) deallocate(this%gz_lev_lon       )
    if (allocated(this%gz_lev_lat       )) deallocate(this%gz_lev_lat       )
    if (allocated(this%rhod             )) deallocate(this%rhod             )
    if (allocated(this%rhod_lon         )) deallocate(this%rhod_lon         )
    if (allocated(this%rhod_lat         )) deallocate(this%rhod_lat         )
    if (allocated(this%p_lev_lon        )) deallocate(this%p_lev_lon        )
    if (allocated(this%p_lev_lat        )) deallocate(this%p_lev_lat        )
    if (allocated(this%u_lev_lon        )) deallocate(this%u_lev_lon        )
    if (allocated(this%v_lev_lat        )) deallocate(this%v_lev_lat        )
    if (allocated(this%mf_lev_lon_n     )) deallocate(this%mf_lev_lon_n     )
    if (allocated(this%mf_lev_lat_n     )) deallocate(this%mf_lev_lat_n     )

    if (nonhydrostatic) then
      if (associated(this%p             )) deallocate(this%p                )
      if (associated(this%p_lev         )) deallocate(this%p_lev            )
    end if

    if (allocated(this%tension_h        )) deallocate(this%tension_h        )
    if (allocated(this%shear_h          )) deallocate(this%shear_h          )
    if (allocated(this%kmh              )) deallocate(this%kmh              )
    if (allocated(this%kmh_vtx          )) deallocate(this%kmh_vtx          )
    if (allocated(this%kmh_lon          )) deallocate(this%kmh_lon          )
    if (allocated(this%kmh_lat          )) deallocate(this%kmh_lat          )

    if (allocated(this%q                )) deallocate(this%q                )
    if (allocated(this%qmf_lon          )) deallocate(this%qmf_lon          )
    if (allocated(this%qmf_lat          )) deallocate(this%qmf_lat          )
    if (allocated(this%qmf_lev          )) deallocate(this%qmf_lev          )

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
