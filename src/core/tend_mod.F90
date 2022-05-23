module tend_mod

  use const_mod
  use mesh_mod
  use namelist_mod
  use allocator_mod

  implicit none

  private

  public tend_type
  public operator(+)
  public operator(*)
  public operator(/)
  public assignment(=)

  type tend_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable, dimension(:,:,:) :: du
    real(r8), allocatable, dimension(:,:,:) :: dv
    real(r8), allocatable, dimension(:,:,:) :: dgz
    real(r8), allocatable, dimension(:,:,:) :: dpt
    real(r8), allocatable, dimension(:,:  ) :: dphs
    logical :: update_u   = .false.
    logical :: update_v   = .false.
    logical :: update_gz  = .false.
    logical :: update_pt  = .false.
    logical :: update_phs = .false.
    logical :: copy_gz    = .false.
    logical :: copy_pt    = .false.
    logical :: copy_phs   = .false.
    ! Individual tendencies
    real(r8), allocatable, dimension(:,:,:) :: qhv
    real(r8), allocatable, dimension(:,:,:) :: qhu
    real(r8), allocatable, dimension(:,:,:) :: dkedlon
    real(r8), allocatable, dimension(:,:,:) :: dkedlat
    real(r8), allocatable, dimension(:,:,:) :: dmfdlon
    real(r8), allocatable, dimension(:,:,:) :: dmfdlat
    real(r8), allocatable, dimension(:,:,:) :: dptfdlon ! Zonal potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: dptfdlat ! Meridional potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: dptfdlev ! Vertical potential temperature flux
    real(r8), allocatable, dimension(:,:,:) :: pgf_lon
    real(r8), allocatable, dimension(:,:,:) :: pgf_lat
    real(r8), allocatable, dimension(:,:,:) :: wedudlev
    real(r8), allocatable, dimension(:,:,:) :: wedvdlev
    real(r8), allocatable, dimension(:,:,:) :: smag_dptdt ! Smagorinsky damping potential temperature tendency
    real(r8), allocatable, dimension(:,:,:) :: smag_dudt
    real(r8), allocatable, dimension(:,:,:) :: smag_dvdt
    ! Nonhydrostatic tendencies
    real(r8), allocatable, dimension(:,:,:) :: adv_gz_lon ! Advection terms of geopotential
    real(r8), allocatable, dimension(:,:,:) :: adv_gz_lat ! Advection terms of geopotential
    real(r8), allocatable, dimension(:,:,:) :: adv_gz_lev ! Advection terms of geopotential
    real(r8), allocatable, dimension(:,:,:) :: adv_w_lon  ! Advection terms of vertical speed
    real(r8), allocatable, dimension(:,:,:) :: adv_w_lat  ! Advection terms of vertical speed
    real(r8), allocatable, dimension(:,:,:) :: adv_w_lev  ! Advection terms of vertical speed
  contains
    procedure :: init => tend_init
    procedure :: reset_flags => tend_reset_flags
    procedure :: clear => tend_clear
    final :: tend_final
  end type tend_type

  interface operator(+)
    module procedure add_tends
  end interface operator(+)

  interface operator(*)
    module procedure mult_scalar
  end interface operator(*)

  interface operator(/)
    module procedure div_scalar
  end interface operator(/)

  interface assignment(=)
    module procedure assign_tend
  end interface assignment(=)

contains

  subroutine tend_init(this, mesh)

    class(tend_type), intent(inout)         :: this
    type(mesh_type ), intent(in   ), target :: mesh

    call this%clear()

    this%mesh => mesh

    call allocate_array(mesh, this%du      , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dv      , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dgz     , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dpt     , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dphs    , full_lon=.true., full_lat=.true.                 )
    call allocate_array(mesh, this%qhv     , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%qhu     , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dkedlon , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dkedlat , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmfdlon , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dmfdlat , full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlon, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlat, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%dptfdlev, full_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pgf_lon , half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%pgf_lat , full_lon=.true., half_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%wedudlev, half_lon=.true., full_lat=.true., full_lev=.true.)
    call allocate_array(mesh, this%wedvdlev, full_lon=.true., half_lat=.true., full_lev=.true.)
    if (use_smag_damp) then
      call allocate_array(mesh, this%smag_dptdt, full_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%smag_dudt , half_lon=.true., full_lat=.true., full_lev=.true.)
      call allocate_array(mesh, this%smag_dvdt , full_lon=.true., half_lat=.true., full_lev=.true.)
    end if

    call allocate_array(mesh, this%adv_gz_lon, full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%adv_gz_lat, full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%adv_gz_lev, full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%adv_w_lon , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%adv_w_lat , full_lon=.true., full_lat=.true., half_lev=.true.)
    call allocate_array(mesh, this%adv_w_lev , full_lon=.true., full_lat=.true., half_lev=.true.)

  end subroutine tend_init

  subroutine tend_reset_flags(this)

    class(tend_type), intent(inout) :: this

    this%update_u   = .false.
    this%update_v   = .false.
    this%update_gz  = .false.
    this%update_pt  = .false.
    this%update_phs = .false.
    this%copy_gz    = .false.
    this%copy_pt    = .false.
    this%copy_phs   = .false.

  end subroutine tend_reset_flags

  subroutine tend_clear(this)

    class(tend_type), intent(inout) :: this

    if (allocated(this%du      )) deallocate(this%du      )
    if (allocated(this%dv      )) deallocate(this%dv      )
    if (allocated(this%dgz     )) deallocate(this%dgz     )
    if (allocated(this%dpt     )) deallocate(this%dpt     )
    if (allocated(this%dphs    )) deallocate(this%dphs    )
    if (allocated(this%qhv     )) deallocate(this%qhv     )
    if (allocated(this%qhu     )) deallocate(this%qhu     )
    if (allocated(this%dkedlon )) deallocate(this%dkedlon )
    if (allocated(this%dkedlat )) deallocate(this%dkedlat )
    if (allocated(this%dmfdlon )) deallocate(this%dmfdlon )
    if (allocated(this%dmfdlat )) deallocate(this%dmfdlat )
    if (allocated(this%dptfdlon)) deallocate(this%dptfdlon)
    if (allocated(this%dptfdlat)) deallocate(this%dptfdlat)
    if (allocated(this%dptfdlev)) deallocate(this%dptfdlev)
    if (allocated(this%pgf_lon )) deallocate(this%pgf_lon )
    if (allocated(this%pgf_lat )) deallocate(this%pgf_lat )
    if (allocated(this%wedudlev)) deallocate(this%wedudlev)
    if (allocated(this%wedvdlev)) deallocate(this%wedvdlev)

    if (allocated(this%smag_dptdt)) deallocate(this%smag_dptdt)
    if (allocated(this%smag_dudt )) deallocate(this%smag_dudt )
    if (allocated(this%smag_dvdt )) deallocate(this%smag_dvdt )

    if (allocated(this%adv_gz_lon)) deallocate(this%adv_gz_lon)
    if (allocated(this%adv_gz_lat)) deallocate(this%adv_gz_lat)
    if (allocated(this%adv_gz_lev)) deallocate(this%adv_gz_lev)
    if (allocated(this%adv_w_lon )) deallocate(this%adv_w_lon )
    if (allocated(this%adv_w_lat )) deallocate(this%adv_w_lat )
    if (allocated(this%adv_w_lev )) deallocate(this%adv_w_lev )

  end subroutine tend_clear

  subroutine tend_final(this)

    type(tend_type), intent(inout) :: this

    call this%clear()

  end subroutine tend_final

  function add_tends(x, y) result(res)

    type(tend_type), intent(in) :: x
    type(tend_type), intent(in) :: y

    type(tend_type) res

    if (x%update_u .and. y%update_u) then
      res%du = x%du + y%du
      res%update_u = .true.
    else
      res%update_u = .false.
    end if
    if (x%update_v .and. y%update_v) then
      res%dv = x%dv + y%dv
      res%update_v = .true.
    else
      res%update_v = .false.
    end if
    if (baroclinic) then
      if (x%update_phs .and. y%update_phs) then
        res%dphs = x%dphs + y%dphs
        res%update_phs = .true.
      else
        res%update_phs = .false.
      end if
      if (x%update_pt .and. y%update_pt) then
        res%dpt = x%dpt + y%dpt
        res%update_pt = .true.
      else
        res%update_pt = .false.
      end if
    else if (x%update_gz .and. y%update_gz) then
      res%dgz = x%dgz + y%dgz
      res%update_gz = .true.
    else
      res%update_gz = .false.
    end if

  end function add_tends

  function mult_scalar(s, x) result(res)

    real(r8), intent(in) :: s
    type(tend_type), intent(in) :: x

    type(tend_type) res

    if (x%update_u) then
      res%du = s * x%du
      res%update_u = .true.
    else
      res%update_u = .false.
    end if
    if (x%update_v) then
      res%dv = s * x%dv
      res%update_v = .true.
    else
      res%update_v = .false.
    end if
    if (baroclinic) then
      if (x%update_phs) then
        res%dphs = s * x%dphs
        res%update_phs = .true.
      else
        res%update_phs = .false.
      end if
      if (x%update_pt) then
        res%dpt = s * x%dpt
        res%update_pt = .true.
      else
        res%update_pt = .false.
      end if
    else if (x%update_gz) then
      res%dgz = s * x%dgz
      res%update_gz = .true.
    else
      res%update_gz = .false.
    end if

  end function mult_scalar

  function div_scalar(x, s) result(res)

    real(r8), intent(in) :: s
    type(tend_type), intent(in) :: x

    type(tend_type) res

    if (x%update_u) then
      res%du = x%du / s
      res%update_u = .true.
    else
      res%update_u = .false.
    end if
    if (x%update_v) then
      res%dv = x%dv / s
      res%update_v = .true.
    else
      res%update_v = .false.
    end if
    if (baroclinic) then
      if (x%update_phs) then
        res%dphs = x%dphs / s
        res%update_phs = .true.
      else
        res%update_phs = .false.
      end if
      if (x%update_pt) then
        res%dpt = x%dpt / s
        res%update_pt = .true.
      else
        res%update_pt = .false.
      end if
    else if (x%update_gz) then
      res%dgz = x%dgz / s
      res%update_gz = .true.
    else
      res%update_gz = .false.
    end if

  end function div_scalar

  subroutine assign_tend(x, y)

    type(tend_type), intent(inout) :: x
    type(tend_type), intent(in) :: y

    if (y%update_u) then
      x%du = y%du
      x%update_u = .true.
    else
      x%update_u = .false.
    end if
    if (y%update_v) then
      x%dv = y%dv
      x%update_v = .true.
    else
      x%update_v = .false.
    end if
    if (baroclinic) then
      if (y%update_phs) then
        x%dphs = y%dphs
        x%update_phs = .true.
      else
        x%update_phs = .false.
      end if
      if (y%update_pt) then
        x%dpt = y%dpt
        x%update_pt = .true.
      else
        x%update_pt = .false.
      end if
    else if (y%update_gz) then
      x%dgz = y%dgz
      x%update_gz = .true.
    else
      x%update_gz = .false.
    end if

  end subroutine assign_tend

end module tend_mod
