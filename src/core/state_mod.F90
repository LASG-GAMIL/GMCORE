module state_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use allocator_mod
  use parallel_types_mod

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
    real(r8), allocatable, dimension(:,:,:) :: u             ! Zonal wind speed (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: u850          ! Zonal wind speed on 850hPa
    real(r8), allocatable, dimension(:,:,:) :: u700          ! Zonal wind speed on 700hPa
    real(r8), allocatable, dimension(:,:,:) :: v             ! Meridional wind speed (m s-1)
    real(r8), allocatable, dimension(:,:,:) :: v850          ! Meridional wind speed on 850hPa
    real(r8), allocatable, dimension(:,:,:) :: v700          ! Meridional wind speed on 700hPa
    real(r8), allocatable, dimension(:,:,:) :: w             ! Vertical wind speed
    real(r8), allocatable, dimension(:,:,:) :: wp            ! ω = dp / dt (Pa s-1)
    real(r8), allocatable, dimension(:,:,:) :: wedphdlev_lev     ! Vertical coordinate speed multiplied by 𝛛π/𝛛η
    real(r8), allocatable, dimension(:,:,:) :: wedphdlev_lev_lon ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: wedphdlev_lev_lat ! Vertical coordinate speed multiplied by 𝛛π/𝛛η on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: gz            ! Geopotential (m2 s-2)
    real(r8), allocatable, dimension(:,:,:) :: gz_p          ! Perturbed geopotential (m2 s-2)
    real(r8), allocatable, dimension(:,:,:) :: gz_lev        ! Geopotential height on half levels (m2 s-2)
    real(r8), allocatable, dimension(:,:,:) :: m             ! Mass
    real(r8), allocatable, dimension(:,:,:) :: m_vtx         ! Mass on vertex
    real(r8), allocatable, dimension(:,:,:) :: m_lon         ! Mass on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: m_lat         ! Mass on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_n      ! Normal mass flux on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_n      ! Normal mass flux on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: mf_lat_t      ! Tangient mass flux on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: mf_lon_t      ! Tangient mass flux on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: pv            ! Potential vorticity
    real(r8), allocatable, dimension(:,:,:) :: pv_lon        ! Potential vorticity on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: pv_lat        ! Potential vorticity on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_t     ! Meridional potential vorticity difference on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: dpv_lon_n     ! Zonal potential vorticity difference on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_t     ! Zonal potential vorticity difference on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: dpv_lat_n     ! Meridional potential vorticity difference on merdional edge
    real(r8), allocatable, dimension(:,:,:) :: ke            ! Kinetic energy
    real(r8), allocatable, dimension(:,:,:) :: pt            ! Potential temperature
    real(r8), allocatable, dimension(:,:,:) :: pt_lon        ! Potential temperature on the zonal edge
    real(r8), allocatable, dimension(:,:,:) :: pt_lat        ! Potential temperature on the merdional edge
    real(r8), allocatable, dimension(:,:,:) :: pt_lev        ! Potential temperature on the vertical edge
    real(r8), allocatable, dimension(:,:,:) :: pt_p          ! Perturbed potential temperature
    real(r8), allocatable, dimension(:,:,:) :: t             ! Temperature
    real(r8), allocatable, dimension(:,:  ) :: t850          ! Temperature on 850hPa
    real(r8), allocatable, dimension(:,:  ) :: t700          ! Temperature on 700hPa
    real(r8), allocatable, dimension(:,:,:) :: t_lnpop       ! T ln(p_{k+1/2} / p_{k-1/2})
    real(r8), allocatable, dimension(:,:,:) :: t_lnpop_lon   ! T ln(p_{k+1/2} / p_{k-1/2})
    real(r8), allocatable, dimension(:,:,:) :: t_lnpop_lat   ! T ln(p_{k+1/2} / p_{k-1/2})
    real(r8), allocatable, dimension(:,:,:) :: ak_t          ! ak T
    real(r8), allocatable, dimension(:,:,:) :: ak_t_lon      ! ak T
    real(r8), allocatable, dimension(:,:,:) :: ak_t_lat      ! ak T
    real(r8), allocatable, dimension(:,:,:) :: p             ! Pressure on full levels
    real(r8), allocatable, dimension(:,:,:) :: p_p           ! Perturbed pressure on full levels
    real(r8), allocatable, dimension(:,:,:) :: p_lev         ! Pressure on half levels
    real(r8), allocatable, dimension(:,:,:) :: ph            ! Hydrostatic pressure on full levels
    real(r8), allocatable, dimension(:,:,:) :: ph_lev        ! Hydrostatic pressure on half levels
    real(r8), allocatable, dimension(:,:  ) :: phs           ! Surface hydrostatic pressure
    real(r8), allocatable, dimension(:,:,:) :: ak            ! Coefficient for calculate gz on full levels
    real(r8), allocatable, dimension(:,:,:) :: div           ! Divergence (s-1)
    real(r8), allocatable, dimension(:,:,:) :: div2          ! Laplacian of divergence (s-1)
    real(r8), allocatable, dimension(:,:,:) :: vor           ! Vorticity (s-1)
    ! Smagorinsky damping variables
    real(r8), allocatable, dimension(:,:,:) :: tension_h     ! tension strain
    real(r8), allocatable, dimension(:,:,:) :: shear_h       ! shear strain on vertex
    real(r8), allocatable, dimension(:,:,:) :: kmh           ! nonlinear diffusion coef
    real(r8), allocatable, dimension(:,:,:) :: kmh_vtx       ! nonlinear diffusion coef on vertex
    real(r8), allocatable, dimension(:,:,:) :: kmh_lon       ! nonlinear diffusion coef on zonal edge
    real(r8), allocatable, dimension(:,:,:) :: kmh_lat       ! nonlinear diffusion coef on meridional edge
    real(r8) tm
    real(r8) te
    real(r8) tpe
    real(r8) tav
    type(async_type), allocatable :: async(:)
  contains
    procedure :: init => state_init
    procedure :: clear => state_clear
    final :: state_final
  end type state_type

contains

  subroutine state_init(this, mesh)

    class(state_type), intent(inout)         :: this
    type(mesh_type  ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%u            , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u850         , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%u700         , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v            , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v850         , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%v700         , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%wp           , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%wedphdlev_lev    , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%wedphdlev_lev_lon, half_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%wedphdlev_lev_lat, full_lon=.true., half_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%gz           , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%gz_lev       , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%m            , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_vtx        , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_lon        , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%m_lat        , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lon_n     , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lon_t     , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lat_n     , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%mf_lat_t     , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv           , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv_lon       , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pv_lat       , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpv_lon_t    , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpv_lon_n    , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpv_lat_t    , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpv_lat_n    , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ke           , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pt           , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pt_lon       , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pt_lat       , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pt_lev       , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%t            , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%t850         , full_lon=.true., full_lat=.true.                 )
    call allocate_array(mesh, this%t700         , full_lon=.true., full_lat=.true.                 )
    call allocate_array(mesh, this%p            , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%p_lev        , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%ph           , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ph_lev       , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%phs          , full_lon=.true., full_lat=.true.                 )
    call allocate_array(mesh, this%div          , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%vor          , half_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%ak           , full_lon=.true., full_lat=.true., full_lev=.true.)

    if (baroclinic .and. .not. hydrostatic) then
      call allocate_array(mesh, this%w            , full_lon=.true., full_lat=.true., half_lev=.true.)
    end if

    if (pgf_scheme == 'sb81') then
      call allocate_array(mesh, this%t_lnpop      , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%t_lnpop_lon  , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%t_lnpop_lat  , full_lon=.true., half_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%ak_t         , full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%ak_t_lon     , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%ak_t_lat     , full_lon=.true., half_lat=.true., full_lev=.true.)
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

    allocate(this%async(11))

  end subroutine state_init

  subroutine state_clear(this)

    class(state_type), intent(inout) :: this

    if (allocated(this%u                )) deallocate(this%u                )
    if (allocated(this%u850             )) deallocate(this%u850             )
    if (allocated(this%u700             )) deallocate(this%u700             )
    if (allocated(this%v                )) deallocate(this%v                )
    if (allocated(this%v850             )) deallocate(this%v850             )
    if (allocated(this%v700             )) deallocate(this%v700             )
    if (allocated(this%w                )) deallocate(this%w                )
    if (allocated(this%wp               )) deallocate(this%wp               )
    if (allocated(this%wedphdlev_lev    )) deallocate(this%wedphdlev_lev    )
    if (allocated(this%wedphdlev_lev_lon)) deallocate(this%wedphdlev_lev_lon)
    if (allocated(this%wedphdlev_lev_lat)) deallocate(this%wedphdlev_lev_lat)
    if (allocated(this%gz               )) deallocate(this%gz               )
    if (allocated(this%gz_p             )) deallocate(this%gz_p             )
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
    if (allocated(this%dpv_lon_t        )) deallocate(this%dpv_lon_t        )
    if (allocated(this%dpv_lon_n        )) deallocate(this%dpv_lon_n        )
    if (allocated(this%dpv_lat_t        )) deallocate(this%dpv_lat_t        )
    if (allocated(this%dpv_lat_n        )) deallocate(this%dpv_lat_n        )
    if (allocated(this%ke               )) deallocate(this%ke               )
    if (allocated(this%pt               )) deallocate(this%pt               )
    if (allocated(this%pt_lon           )) deallocate(this%pt_lon           )
    if (allocated(this%pt_lat           )) deallocate(this%pt_lat           )
    if (allocated(this%pt_lev           )) deallocate(this%pt_lev           )
    if (allocated(this%pt_p             )) deallocate(this%pt_p             )
    if (allocated(this%t                )) deallocate(this%t                )
    if (allocated(this%t850             )) deallocate(this%t850             )
    if (allocated(this%t700             )) deallocate(this%t700             )
    if (allocated(this%t_lnpop          )) deallocate(this%t_lnpop          )
    if (allocated(this%t_lnpop_lon      )) deallocate(this%t_lnpop_lon      )
    if (allocated(this%t_lnpop_lat      )) deallocate(this%t_lnpop_lat      )
    if (allocated(this%ak_t             )) deallocate(this%ak_t             )
    if (allocated(this%ak_t_lon         )) deallocate(this%ak_t_lon         )
    if (allocated(this%ak_t_lat         )) deallocate(this%ak_t_lat         )
    if (allocated(this%p                )) deallocate(this%p                )
    if (allocated(this%p_p              )) deallocate(this%p_p              )
    if (allocated(this%p_lev            )) deallocate(this%p_lev            )
    if (allocated(this%ph               )) deallocate(this%ph               )
    if (allocated(this%ph_lev           )) deallocate(this%ph_lev           )
    if (allocated(this%phs              )) deallocate(this%phs              )
    if (allocated(this%ak               )) deallocate(this%ak               )
    if (allocated(this%div              )) deallocate(this%div              )
    if (allocated(this%div2             )) deallocate(this%div2             )
    if (allocated(this%vor              )) deallocate(this%vor              )

    if (allocated(this%tension_h        )) deallocate(this%tension_h        )
    if (allocated(this%shear_h          )) deallocate(this%shear_h          )
    if (allocated(this%kmh              )) deallocate(this%kmh              )
    if (allocated(this%kmh_vtx          )) deallocate(this%kmh_vtx          )
    if (allocated(this%kmh_lon          )) deallocate(this%kmh_lon          )
    if (allocated(this%kmh_lat          )) deallocate(this%kmh_lat          )

    if (allocated(this%async            )) deallocate(this%async            )

  end subroutine state_clear

  subroutine state_final(this)

    type(state_type), intent(inout) :: this

    call this%clear()

  end subroutine state_final

end module state_mod
