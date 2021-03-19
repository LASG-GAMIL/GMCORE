module mesh_mod

  use flogger
  use const_mod
  use namelist_mod
  use sphere_geometry_mod

  implicit none

  private

  public mesh_type
  public global_mesh

  type mesh_type
    ! For nesting
    integer :: id = 0
    type(mesh_type), pointer :: parent => null()
    integer lon_halo_width
    integer lat_halo_width
    integer num_full_lon
    integer num_half_lon
    integer num_full_lat
    integer num_half_lat
    integer num_full_lev
    integer num_half_lev
    integer full_lon_ibeg
    integer full_lon_iend
    integer full_lat_ibeg
    integer full_lat_iend
    integer full_lat_ibeg_no_pole
    integer full_lat_iend_no_pole
    integer full_lev_ibeg
    integer full_lev_iend
    integer half_lon_ibeg
    integer half_lon_iend
    integer half_lat_ibeg
    integer half_lat_iend
    integer half_lat_ibeg_no_pole
    integer half_lat_iend_no_pole
    integer half_lev_ibeg
    integer half_lev_iend
    integer full_lon_lb
    integer full_lon_ub
    integer half_lon_lb
    integer half_lon_ub
    integer full_lat_lb
    integer full_lat_ub
    integer half_lat_lb
    integer half_lat_ub
    integer full_lev_lb
    integer full_lev_ub
    integer half_lev_lb
    integer half_lev_ub
    real(r8) start_lon
    real(r8) end_lon
    real(r8) start_lat
    real(r8) end_lat
    real(r8) dlon
    real(r8), allocatable, dimension(:  ) :: dlat
    real(r8), allocatable, dimension(:  ) :: full_dlev
    real(r8), allocatable, dimension(:  ) :: half_dlev
    real(r8), allocatable, dimension(:  ) :: half_dlev_upper
    real(r8), allocatable, dimension(:  ) :: half_dlev_lower
    real(r8) total_area
    real(r8), allocatable, dimension(:  ) :: full_lon
    real(r8), allocatable, dimension(:  ) :: half_lon
    real(r8), allocatable, dimension(:  ) :: full_lat
    real(r8), allocatable, dimension(:  ) :: half_lat
    real(r8), allocatable, dimension(:  ) :: full_lev
    real(r8), allocatable, dimension(:  ) :: half_lev
    real(r8), allocatable, dimension(:  ) :: full_cos_lon
    real(r8), allocatable, dimension(:  ) :: half_cos_lon
    real(r8), allocatable, dimension(:  ) :: full_sin_lon
    real(r8), allocatable, dimension(:  ) :: half_sin_lon
    real(r8), allocatable, dimension(:  ) :: full_cos_lat
    real(r8), allocatable, dimension(:  ) :: half_cos_lat
    real(r8), allocatable, dimension(:  ) :: full_sin_lat
    real(r8), allocatable, dimension(:  ) :: half_sin_lat
    ! For output
    real(r8), allocatable, dimension(:  ) :: full_lon_deg
    real(r8), allocatable, dimension(:  ) :: half_lon_deg
    real(r8), allocatable, dimension(:  ) :: full_lat_deg
    real(r8), allocatable, dimension(:  ) :: half_lat_deg
    ! Area for weighting
    real(r8), allocatable, dimension(:  ) :: area_cell
    real(r8), allocatable, dimension(:  ) :: area_lon
    real(r8), allocatable, dimension(:  ) :: area_lon_west
    real(r8), allocatable, dimension(:  ) :: area_lon_east
    real(r8), allocatable, dimension(:  ) :: area_lon_north
    real(r8), allocatable, dimension(:  ) :: area_lon_south
    real(r8), allocatable, dimension(:  ) :: area_lat
    real(r8), allocatable, dimension(:  ) :: area_lat_west
    real(r8), allocatable, dimension(:  ) :: area_lat_east
    real(r8), allocatable, dimension(:  ) :: area_lat_north
    real(r8), allocatable, dimension(:  ) :: area_lat_south
    real(r8), allocatable, dimension(:  ) :: area_vtx
    real(r8), allocatable, dimension(:,:) :: area_subcell
    ! Edge length
    real(r8), allocatable, dimension(:  ) :: de_lon
    real(r8), allocatable, dimension(:  ) :: de_lat
    real(r8), allocatable, dimension(:  ) :: le_lat
    real(r8), allocatable, dimension(:  ) :: le_lon
    ! Coriolis parameters
    real(r8), allocatable, dimension(:  ) :: full_f
    real(r8), allocatable, dimension(:  ) :: half_f
    ! Weight for constructing tangential wind
    real(r8), allocatable, dimension(:,:) :: full_tangent_wgt
    real(r8), allocatable, dimension(:,:) :: half_tangent_wgt
  contains
    procedure :: init_global => mesh_init_global
    procedure :: init_from_parent => mesh_init_from_parent
    procedure :: common_init => mesh_common_init
    procedure :: has_south_pole => mesh_has_south_pole
    procedure :: has_north_pole => mesh_has_north_pole
    procedure :: is_south_pole => mesh_is_south_pole
    procedure :: is_north_pole => mesh_is_north_pole
    procedure :: is_pole => mesh_is_pole
    procedure :: is_full_lat_next_to_pole => mesh_is_full_lat_next_to_pole
    procedure :: is_half_lat_next_to_pole => mesh_is_half_lat_next_to_pole
    procedure :: is_inside_with_halo_full_lat => mesh_is_inside_with_halo_full_lat
    procedure :: is_inside_with_halo_half_lat => mesh_is_inside_with_halo_half_lat
    procedure :: is_outside_pole_full_lat => mesh_is_outside_pole_full_lat
    procedure :: is_outside_pole_half_lat => mesh_is_outside_pole_half_lat
    final :: mesh_final
  end type mesh_type

  type(mesh_type), target :: global_mesh

contains

  subroutine mesh_init_global(this, num_lon, num_lat, num_lev, id, lon_halo_width, lat_halo_width)

    class(mesh_type), intent(inout)           :: this
    integer         , intent(in   )           :: num_lon
    integer         , intent(in   )           :: num_lat
    integer         , intent(in   ), optional :: num_lev
    integer         , intent(in   ), optional :: id
    integer         , intent(in   ), optional :: lon_halo_width
    integer         , intent(in   ), optional :: lat_halo_width

    real(r8) dlat0
    real(16) x(3), y(3), z(3)
    integer i, j, j0, jj

    this%num_full_lon  = num_lon
    this%num_half_lon  = num_lon
    this%full_lon_ibeg = 1
    this%full_lon_iend = this%num_full_lon
    this%half_lon_ibeg = 1
    this%half_lon_iend = this%num_half_lon
    this%num_full_lat  = num_lat
    this%num_half_lat  = num_lat - 1
    this%full_lat_ibeg = 1
    this%full_lat_iend = this%num_full_lat
    this%half_lat_ibeg = 1
    this%half_lat_iend = this%num_half_lat
    this%num_full_lev  = merge(num_lev, 1, present(num_lev))
    this%num_half_lev  = this%num_full_lev + 1
    this%full_lev_ibeg = 1
    this%full_lev_iend = this%num_full_lev
    this%half_lev_ibeg = 1
    this%half_lev_iend = this%num_half_lev

    ! Here we set baroclinic according to levels.
    baroclinic = this%num_full_lev > 1
    if (.not. baroclinic) then
      hydrostatic = .false.
      nonhydrostatic = .false.
    end if

    this%id             = merge(id            , -1, present(id))
    this%lon_halo_width = merge(lon_halo_width,  1, present(lon_halo_width))
    this%lat_halo_width = merge(lat_halo_width,  1, present(lat_halo_width))
    this%start_lon      = 0.0_r8
    this%end_lon        =  pi2
    this%start_lat      = -pi05
    this%end_lat        =  pi05

    call this%common_init()

    this%dlon = (this%end_lon - this%start_lon) / this%num_full_lon
    do i = this%full_lon_lb, this%full_lon_ub
      this%full_lon(i) = this%start_lon + (i - 1) * this%dlon
      this%half_lon(i) = this%full_lon(i) + 0.5_r8 * this%dlon
      this%full_lon_deg(i) = this%full_lon(i) * deg
      this%half_lon_deg(i) = this%half_lon(i) * deg
    end do

    ! Set initial guess latitudes of full merdional grids.
    dlat0 = (this%end_lat - this%start_lat) / this%num_half_lat
    do j = 1, this%num_half_lat
      this%half_lat(j) = this%start_lat + (j - 0.5_r8) * dlat0
      if (abs(this%half_lat(j)) < 1.0e-12) this%half_lat(j) = 0.0_r8
    end do

    if (coarse_pole_mul /= 0) then
      ! Calculate real dlat which is large at polar region.
      dlat0 = this%dlon
      do j = 1, this%num_half_lat
        this%dlat(j) = dlat0 * (1 + (coarse_pole_mul - 1) * exp(-coarse_pole_decay * (abs(this%half_lat(j)) - pi05)**2))
      end do
      this%dlat(1:this%num_half_lat) = this%dlat(1:this%num_half_lat) / sum(this%dlat(1:this%num_half_lat)) * pi
    else
      this%dlat(1:this%num_half_lat) = dlat0
    end if

    ! Set latitudes of full merdional grids.
    this%full_lat(1) = this%start_lat
    this%full_lat_deg(1) = this%start_lat * deg
    do j = 2, this%num_full_lat - 1
      this%full_lat(j) = this%full_lat(j-1) + this%dlat(j-1)
      if (abs(this%full_lat(j)) < 1.0e-12) this%full_lat(j) = 0.0_r8
      this%full_lat_deg(j) = this%full_lat(j) * deg
    end do
    this%full_lat(this%num_full_lat) = this%end_lat
    this%full_lat_deg(this%num_full_lat) = this%end_lat * deg

    ! Set latitudes of half merdional grids.
    do j = 1, this%num_half_lat
      if (is_inf(this%full_lat(j)) .or. this%full_lat(j) == pi05) cycle
      this%half_lat(j) = this%full_lat(j) + 0.5_r8 * this%dlat(j)
      if (abs(this%half_lat(j)) < 1.0e-12) this%half_lat(j) = 0.0_r8
      this%half_lat_deg(j) = this%half_lat(j) * deg
    end do

    ! Ensure the grids are equatorial symmetry.
    do j = 1, this%num_full_lat
      if (this%full_lat(j) > 0) then
        this%full_lat(j) = -this%full_lat(this%num_full_lat-j+1)
        this%full_lat_deg(j) = -this%full_lat_deg(this%num_full_lat-j+1)
      end if
    end do
    do j = 1, this%num_half_lat
      if (this%half_lat(j) > 0) then
        this%half_lat(j) = -this%half_lat(this%num_half_lat-j+1)
        this%half_lat_deg(j) = -this%half_lat_deg(this%num_half_lat-j+1)
      end if
    end do

    do i = this%full_lon_lb, this%full_lon_ub
      this%full_cos_lon(i) = cos(this%full_lon(i))
      this%full_sin_lon(i) = sin(this%full_lon(i))
    end do

    do i = this%half_lon_lb, this%half_lon_ub
      this%half_cos_lon(i) = cos(this%half_lon(i))
      this%half_sin_lon(i) = sin(this%half_lon(i))
    end do

    do j = this%half_lat_lb, this%half_lat_ub
      if (this%half_lat(j) >= -pi05 .and. this%half_lat(j) <= pi05) then
        this%half_cos_lat(j) = cos(this%half_lat(j))
        this%half_sin_lat(j) = sin(this%half_lat(j))
      end if
    end do

    do j = this%full_lat_lb, this%full_lat_ub
      if (this%full_lat(j) >= -pi05 .and. this%full_lat(j) <= pi05) then
        this%full_cos_lat(j) = cos(this%full_lat(j))
        this%full_sin_lat(j) = sin(this%full_lat(j))
      end if
    end do

    ! Ensure the values of cos_lat and sin_lat are expected at the Poles.
    this%full_cos_lat(this%full_lat_ibeg) =  0.0_r8
    this%full_sin_lat(this%full_lat_ibeg) = -1.0_r8
    this%full_cos_lat(this%full_lat_iend) =  0.0_r8
    this%full_sin_lat(this%full_lat_iend) =  1.0_r8

    do j = this%full_lat_ibeg, this%full_lat_iend
      if (this%is_south_pole(j)) then
        this%area_cell(j) = radius**2 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
        this%area_subcell(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
      else if (this%is_north_pole(j)) then
        this%area_cell(j) = radius**2 * this%dlon * (1.0 - this%half_sin_lat(j-1))
        this%area_subcell(1,j) = radius**2 * 0.5d0 * this%dlon * (1.0d0 - this%half_sin_lat(j-1))
      else
        this%area_cell(j) = radius**2 * this%dlon * (this%half_sin_lat(j) - this%half_sin_lat(j-1))
        this%area_subcell(1,j) = radius**2 * 0.5d0 * this%dlon * (this%full_sin_lat(j) - this%half_sin_lat(j-1))
        this%area_subcell(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) - this%full_sin_lat(j))
        !
        !           1,j
        !           /|
        !          / |
        !         /  |
        !        /   |
        !    1,j \   |
        !         \  |
        !          \ |
        !           \|
        !          1,j-1
        !
        call cartesian_transform(this%full_lon(1), this%full_lat(j  ), x(1), y(1), z(1))
        call cartesian_transform(this%half_lon(1), this%half_lat(j-1), x(2), y(2), z(2))
        call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(3), y(3), z(3))
        this%area_lon_west(j) = calc_area(x, y, z)
        this%area_lon_east(j) = this%area_lon_west(j)
        this%area_lon(j) = this%area_lon_west(j) + this%area_lon_east(j)
        !
        !          1,j
        !           /\
        !          /  \
        !         /    \
        !        /______\
        !    1,j          2,j
        !
        call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(1), y(1), z(1))
        call cartesian_transform(this%full_lon(1), this%full_lat(j  ), x(2), y(2), z(2))
        call cartesian_transform(this%full_lon(2), this%full_lat(j  ), x(3), y(3), z(3))
        this%area_lon_north(j) = calc_area_with_last_small_arc(x, y, z)
        !
        !    1,j          2,j
        !        --------
        !        \      /
        !         \    /
        !          \  /
        !           \/
        !         1,j-1
        !
        call cartesian_transform(this%half_lon(1), this%half_lat(j-1), x(1), y(1), z(1))
        call cartesian_transform(this%full_lon(2), this%full_lat(j  ), x(2), y(2), z(2))
        call cartesian_transform(this%full_lon(1), this%full_lat(j  ), x(3), y(3), z(3))
        this%area_lon_south(j) = calc_area_with_last_small_arc(x, y, z)
      end if
    end do

    do j = this%half_lat_ibeg, this%half_lat_iend
      this%area_vtx(j) = radius**2 * this%dlon * (this%full_sin_lat(j+1) - this%full_sin_lat(j))
      !
      !          2,j+1
      !           /|
      !          / |
      !         /  |
      !        /   |
      !    1,j \   |
      !         \  |
      !          \ |
      !           \|
      !           2,j
      !
      call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(1), y(1), z(1))
      call cartesian_transform(this%full_lon(2), this%full_lat(j  ), x(2), y(2), z(2))
      call cartesian_transform(this%full_lon(2), this%full_lat(j+1), x(3), y(3), z(3))
      this%area_lat_west(j) = calc_area(x, y, z)
      this%area_lat_east(j) = this%area_lat_west(j)
      !
      !         2,j+1
      !           /\
      !          /  \
      !         /    \
      !        /______\
      !    1,j          2,j
      !
      call cartesian_transform(this%full_lon(2), this%full_lat(j+1), x(1), y(1), z(1))
      call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(2), y(2), z(2))
      call cartesian_transform(this%half_lon(2), this%half_lat(j  ), x(3), y(3), z(3))
      this%area_lat_north(j) = calc_area_with_last_small_arc(x, y, z)
      !
      !    1,j          2,j
      !        --------
      !        \      /
      !         \    /
      !          \  /
      !           \/
      !          2,j
      !
      call cartesian_transform(this%full_lon(2), this%full_lat(j), x(1), y(1), z(1))
      call cartesian_transform(this%half_lon(2), this%half_lat(j), x(2), y(2), z(2))
      call cartesian_transform(this%half_lon(1), this%half_lat(j), x(3), y(3), z(3))
      this%area_lat_south(j) = calc_area_with_last_small_arc(x, y, z)
      ! Reset up or down area to polar sector area.
      if (this%is_south_pole(j)) then
        this%area_lat_south(j) = this%area_cell(j)
      else if (this%is_north_pole(j+1)) then
        this%area_lat_north(j) = this%area_cell(j+1)
      end if
      this%area_lat(j) = this%area_lat_north(j) + this%area_lat_south(j)
    end do

    do j = this%full_lat_ibeg_no_pole, this%full_lat_iend_no_pole
      this%de_lon(j) = radius * this%full_cos_lat(j) * this%dlon
      this%le_lon(j) = 2.0d0 * this%area_lon(j) / this%de_lon(j)
    end do
    if (this%has_south_pole()) then
      this%le_lon(this%full_lat_ibeg) = 0.0_r8
      this%de_lon(this%full_lat_ibeg) = 0.0_r8
    end if
    if (this%has_north_pole()) then
      this%le_lon(this%full_lat_iend) = 0.0_r8
      this%de_lon(this%full_lat_iend) = 0.0_r8
    end if

    do j = this%half_lat_ibeg_no_pole, this%half_lat_iend_no_pole
      this%le_lat(j) = radius * this%half_cos_lat(j) * this%dlon
      this%de_lat(j) = 2.0d0 * this%area_lat(j) / this%le_lat(j)
    end do

    !  ____________________                 ____________________                  ____________________                  ____________________                 
    ! |          |         |               |          |         |                |          |         |                |          |         |                
    ! |          |         |               |          |         |                |          |         |                |          |         |                
    ! |          |         |               |          |         |                |          |         |                |          |         |                
    ! |          |         |               |          |         |                |          |         |                |          |         |                
    ! |_____o____|____o____|   j           |_____o____|____*____|   j            |_____*____|____o____|   j            |_____o____|____o____|   j             
    ! |          |////|////|               |          |////|    |                |/////|    |         |                |     |    |         |                
    ! |          |/3//|/2//|               |          |////|    |                |/////|    |         |                |     |    |         |                
    ! |          x---------|   j           |          x---------|   j            |-----|----x         |   j            |-----|----x         |   j            
    ! |          |    |/1//|               |          |    |    |                |/////|////|         |                |     |////|         |                
    ! |_____o____|____*____|   j - 1       |_____o____|____o____|   j - 1        |_____o____|____o____|   j - 1        |_____*____|____o____|   j - 1        
    !       i    i   i+1                         i    i   i+1                          i    i   i+1                          i
    !
    !
    !       [ 1    As_1 + As_2 + As_3]
    ! w = - [--- - ------------------]
    !       [ 2        A_{i+1,j}     ]
    !
    !

    select case (tangent_wgt_scheme)
    case ('classic')
      do j = this%full_lat_ibeg_no_pole, this%full_lat_iend_no_pole
        this%full_tangent_wgt(1,j) = this%le_lat(j-1) / this%de_lon(j) * 0.25d0
        this%full_tangent_wgt(2,j) = this%le_lat(j  ) / this%de_lon(j) * 0.25d0
      end do

      do j = this%half_lat_ibeg_no_pole, this%half_lat_iend_no_pole
        this%half_tangent_wgt(1,j) = this%le_lon(j  ) / this%de_lat(j) * 0.25d0
        this%half_tangent_wgt(2,j) = this%le_lon(j+1) / this%de_lat(j) * 0.25d0
      end do
    case ('thuburn09')
      do j = this%full_lat_ibeg_no_pole, this%full_lat_iend_no_pole
        this%full_tangent_wgt(1,j) = this%le_lat(j-1) / this%de_lon(j) * this%area_subcell(2,j  ) / this%area_cell(j  )
        this%full_tangent_wgt(2,j) = this%le_lat(j  ) / this%de_lon(j) * this%area_subcell(1,j  ) / this%area_cell(j  )
      end do

      do j = this%half_lat_ibeg_no_pole, this%half_lat_iend_no_pole
        this%half_tangent_wgt(1,j) = this%le_lon(j  ) / this%de_lat(j) * this%area_subcell(1,j  ) / this%area_cell(j  )
        this%half_tangent_wgt(2,j) = this%le_lon(j+1) / this%de_lat(j) * this%area_subcell(2,j+1) / this%area_cell(j+1)
      end do
    end select

    do j = this%full_lat_ibeg, this%full_lat_iend
      this%full_f(j) = 2.0_r8 * omega * this%full_sin_lat(j)
    end do
    do j = this%half_lat_ibeg, this%half_lat_iend
      this%half_f(j) = 2.0_r8 * omega * this%half_sin_lat(j)
    end do

  end subroutine mesh_init_global

  subroutine mesh_init_from_parent(this, parent, id,               &
                                   lon_halo_width, lat_halo_width, &
                                   lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    class(mesh_type), intent(inout) :: this
    class(mesh_type), intent(in) :: parent
    integer, intent(in) :: id
    integer, intent(in) :: lon_halo_width
    integer, intent(in) :: lat_halo_width
    integer, intent(in) :: lon_ibeg
    integer, intent(in) :: lon_iend
    integer, intent(in) :: lat_ibeg
    integer, intent(in) :: lat_iend

    integer i, j

    this%num_full_lon  = lon_iend - lon_ibeg + 1
    this%num_half_lon  = this%num_full_lon
    this%full_lon_ibeg = lon_ibeg
    this%full_lon_iend = lon_iend
    this%half_lon_ibeg = lon_ibeg
    this%half_lon_iend = lon_iend
    this%num_full_lat  = lat_iend - lat_ibeg + 1
    this%full_lat_ibeg = lat_ibeg
    this%full_lat_iend = lat_iend
    this%half_lat_ibeg = lat_ibeg
    this%half_lat_iend = merge(lat_iend - 1, lat_iend, this%has_north_pole())
    this%num_half_lat  = this%half_lat_iend - this%half_lat_ibeg + 1

    this%num_full_lev  = parent%num_full_lev
    this%num_half_lev  = parent%num_half_lev
    this%full_lev_lb   = parent%full_lev_lb
    this%full_lev_ub   = parent%full_lev_ub
    this%full_lev_ibeg = parent%full_lev_ibeg
    this%full_lev_iend = parent%full_lev_iend
    this%half_lev_lb   = parent%half_lev_lb
    this%half_lev_ub   = parent%half_lev_ub
    this%half_lev_ibeg = parent%half_lev_ibeg
    this%half_lev_iend = parent%half_lev_iend

    this%id             = id
    this%lon_halo_width = lon_halo_width
    this%lat_halo_width = lat_halo_width
    this%start_lon      = parent%full_lon(lon_ibeg)
    this%end_lon        = parent%full_lon(lon_iend) + parent%dlon
    this%start_lat      = merge(parent%half_lat(lat_ibeg-1), -pi05, .not. this%has_south_pole())
    this%end_lat        = merge(parent%half_lat(lat_iend  ),  pi05, .not. this%has_north_pole())

    call this%common_init()

    this%full_dlev = parent%full_dlev
    this%half_dlev = parent%half_dlev
    this%half_dlev_upper = parent%half_dlev_upper
    this%half_dlev_lower = parent%half_dlev_lower

    this%dlon = parent%dlon
    do i = this%full_lon_lb, this%full_lon_ub
      this%full_lon(i) = parent%full_lon(i)
      this%half_lon(i) = parent%half_lon(i)
      this%full_lon_deg(i) = parent%full_lon_deg(i)
      this%half_lon_deg(i) = parent%half_lon_deg(i)
      this%full_sin_lon(i) = parent%full_sin_lon(i)
      this%half_sin_lon(i) = parent%half_sin_lon(i)
      this%full_cos_lon(i) = parent%full_cos_lon(i)
      this%half_cos_lon(i) = parent%half_cos_lon(i)
    end do

    this%dlat = parent%dlat(lbound(this%dlat, 1):ubound(this%dlat, 1))
    do j = this%full_lat_lb, this%full_lat_ub
      this%full_lat(j) = parent%full_lat(j)
      this%full_lat_deg(j) = parent%full_lat_deg(j)
      this%full_sin_lat(j) = parent%full_sin_lat(j)
      this%full_cos_lat(j) = parent%full_cos_lat(j)
      this%area_cell(j) = parent%area_cell(j)
      this%area_subcell(:,j) = parent%area_subcell(:,j)
      this%area_lon_west(j) = parent%area_lon_west(j)
      this%area_lon_east(j) = parent%area_lon_east(j)
      this%area_lon_north(j) = parent%area_lon_north(j)
      this%area_lon_south(j) = parent%area_lon_south(j)
      this%area_lon(j) = parent%area_lon(j)
      this%le_lon(j) = parent%le_lon(j)
      this%de_lon(j) = parent%de_lon(j)
      this%full_tangent_wgt(:,j) = parent%full_tangent_wgt(:,j)
      this%full_f(j) = parent%full_f(j)
    end do
    do j = this%half_lat_lb, this%half_lat_ub
      this%half_lat(j) = parent%half_lat(j)
      this%half_lat_deg(j) = parent%half_lat_deg(j)
      this%half_sin_lat(j) = parent%half_sin_lat(j)
      this%half_cos_lat(j) = parent%half_cos_lat(j)
      this%area_vtx(j) = parent%area_vtx(j)
      this%area_lat_west(j) = parent%area_lat_west(j)
      this%area_lat_east(j) = parent%area_lat_east(j)
      this%area_lat_north(j) = parent%area_lat_north(j)
      this%area_lat_south(j) = parent%area_lat_south(j)
      this%area_lat(j) = parent%area_lat(j)
      this%le_lat(j) = parent%le_lat(j)
      this%de_lat(j) = parent%de_lat(j)
      this%half_tangent_wgt(:,j) = parent%half_tangent_wgt(:,j)
      this%half_f(j) = parent%half_f(j)
    end do

    this%full_lev = parent%full_lev
    this%half_lev = parent%half_lev

  end subroutine mesh_init_from_parent

  subroutine mesh_common_init(this)

    class(mesh_type), intent(inout) :: this

    this%total_area = radius**2 * (this%end_lon - this%start_lon) * (sin(this%end_lat) - sin(this%start_lat))

    this%full_lat_ibeg_no_pole = merge(this%full_lat_ibeg + 1, this%full_lat_ibeg, this%has_south_pole())
    this%full_lat_iend_no_pole = merge(this%full_lat_iend - 1, this%full_lat_iend, this%has_north_pole())
    this%half_lat_ibeg_no_pole = this%half_lat_ibeg
    this%half_lat_iend_no_pole = this%half_lat_iend

    this%full_lon_lb = this%full_lon_ibeg - this%lon_halo_width
    this%full_lon_ub = this%full_lon_iend + this%lon_halo_width
    this%full_lat_lb = this%full_lat_ibeg - this%lat_halo_width
    this%full_lat_ub = this%full_lat_iend + this%lat_halo_width
    this%half_lon_lb = this%half_lon_ibeg - this%lon_halo_width
    this%half_lon_ub = this%half_lon_iend + this%lon_halo_width
    this%half_lat_lb = this%half_lat_ibeg - this%lat_halo_width
    this%half_lat_ub = this%half_lat_iend + this%lat_halo_width
    this%full_lev_lb = this%full_lev_ibeg
    this%full_lev_ub = this%full_lev_iend
    this%half_lev_lb = this%half_lev_ibeg
    this%half_lev_ub = this%half_lev_iend

    allocate(this%dlat               (this%half_lat_lb:this%half_lat_ub)); this%dlat                = 0.0_r8
    allocate(this%full_dlev          (this%full_lev_lb:this%full_lev_ub)); this%full_dlev           = 0.0_r8
    allocate(this%half_dlev          (this%half_lev_lb:this%half_lev_ub)); this%half_dlev           = 0.0_r8
    allocate(this%half_dlev_upper    (this%half_lev_lb:this%half_lev_ub)); this%half_dlev_upper     = 0.0_r8
    allocate(this%half_dlev_lower    (this%half_lev_lb:this%half_lev_ub)); this%half_dlev_lower     = 0.0_r8
    allocate(this%full_lon           (this%full_lon_lb:this%full_lon_ub)); this%full_lon            = inf
    allocate(this%half_lon           (this%half_lon_lb:this%half_lon_ub)); this%half_lon            = inf
    allocate(this%full_lat           (this%full_lat_lb:this%full_lat_ub)); this%full_lat            = inf
    allocate(this%half_lat           (this%half_lat_lb:this%half_lat_ub)); this%half_lat            = inf
    allocate(this%full_lev           (this%full_lev_lb:this%full_lev_ub)); this%full_lev            = inf
    allocate(this%half_lev           (this%half_lev_lb:this%half_lev_ub)); this%half_lev            = inf
    allocate(this%full_cos_lon       (this%full_lon_lb:this%full_lon_ub)); this%full_cos_lon        = inf
    allocate(this%half_cos_lon       (this%half_lon_lb:this%half_lon_ub)); this%half_cos_lon        = inf
    allocate(this%full_sin_lon       (this%full_lon_lb:this%full_lon_ub)); this%full_sin_lon        = inf
    allocate(this%half_sin_lon       (this%half_lon_lb:this%half_lon_ub)); this%half_sin_lon        = inf
    allocate(this%full_cos_lat       (this%full_lat_lb:this%full_lat_ub)); this%full_cos_lat        = inf
    allocate(this%half_cos_lat       (this%half_lat_lb:this%half_lat_ub)); this%half_cos_lat        = inf
    allocate(this%full_sin_lat       (this%full_lat_lb:this%full_lat_ub)); this%full_sin_lat        = inf
    allocate(this%half_sin_lat       (this%half_lat_lb:this%half_lat_ub)); this%half_sin_lat        = inf
    allocate(this%full_lon_deg       (this%full_lon_lb:this%full_lon_ub)); this%full_lon_deg        = inf
    allocate(this%half_lon_deg       (this%half_lon_lb:this%half_lon_ub)); this%half_lon_deg        = inf
    allocate(this%full_lat_deg       (this%full_lat_lb:this%full_lat_ub)); this%full_lat_deg        = inf
    allocate(this%half_lat_deg       (this%half_lat_lb:this%half_lat_ub)); this%half_lat_deg        = inf
    allocate(this%area_cell          (this%full_lat_lb:this%full_lat_ub)); this%area_cell           = 0.0_r8
    allocate(this%area_lon           (this%full_lat_lb:this%full_lat_ub)); this%area_lon            = 0.0_r8
    allocate(this%area_lon_west      (this%full_lat_lb:this%full_lat_ub)); this%area_lon_west       = 0.0_r8
    allocate(this%area_lon_east      (this%full_lat_lb:this%full_lat_ub)); this%area_lon_east       = 0.0_r8
    allocate(this%area_lon_north     (this%full_lat_lb:this%full_lat_ub)); this%area_lon_north      = 0.0_r8
    allocate(this%area_lon_south     (this%full_lat_lb:this%full_lat_ub)); this%area_lon_south      = 0.0_r8
    allocate(this%area_lat           (this%half_lat_lb:this%half_lat_ub)); this%area_lat            = 0.0_r8
    allocate(this%area_lat_west      (this%half_lat_lb:this%half_lat_ub)); this%area_lat_west       = 0.0_r8
    allocate(this%area_lat_east      (this%half_lat_lb:this%half_lat_ub)); this%area_lat_east       = 0.0_r8
    allocate(this%area_lat_north     (this%half_lat_lb:this%half_lat_ub)); this%area_lat_north      = 0.0_r8
    allocate(this%area_lat_south     (this%half_lat_lb:this%half_lat_ub)); this%area_lat_south      = 0.0_r8
    allocate(this%area_vtx           (this%half_lat_lb:this%half_lat_ub)); this%area_vtx            = 0.0_r8
    allocate(this%area_subcell     (2,this%full_lat_lb:this%full_lat_ub)); this%area_subcell        = 0.0_r8
    allocate(this%de_lon             (this%full_lat_lb:this%full_lat_ub)); this%de_lon              = 0.0_r8
    allocate(this%de_lat             (this%half_lat_lb:this%half_lat_ub)); this%de_lat              = 0.0_r8
    allocate(this%le_lat             (this%half_lat_lb:this%half_lat_ub)); this%le_lat              = 0.0_r8
    allocate(this%le_lon             (this%full_lat_lb:this%full_lat_ub)); this%le_lon              = 0.0_r8
    allocate(this%full_f             (this%full_lat_lb:this%full_lat_ub)); this%full_f              = inf
    allocate(this%half_f             (this%half_lat_lb:this%half_lat_ub)); this%half_f              = inf
    allocate(this%full_tangent_wgt (2,this%full_lat_lb:this%full_lat_ub)); this%full_tangent_wgt    = inf
    allocate(this%half_tangent_wgt (2,this%half_lat_lb:this%half_lat_ub)); this%half_tangent_wgt    = inf

  end subroutine mesh_common_init

  logical function mesh_has_south_pole(this) result(res)

    class(mesh_type), intent(in) :: this

    res = this%full_lat_ibeg == 1

  end function mesh_has_south_pole

  logical function mesh_has_north_pole(this) result(res)

    class(mesh_type), intent(in) :: this

    res = this%full_lat_iend == global_mesh%num_full_lat

  end function mesh_has_north_pole

  logical function mesh_is_south_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    ! FIXME: has_south_pole should be removed.
    res = j == 1

  end function mesh_is_south_pole

  logical function mesh_is_north_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j == global_mesh%num_full_lat

  end function mesh_is_north_pole

  logical function mesh_is_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = this%is_south_pole(j) .or. this%is_north_pole(j)

  end function mesh_is_pole

  logical function mesh_is_full_lat_next_to_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j == 2 .or. j == global_mesh%num_full_lat - 1

  end function mesh_is_full_lat_next_to_pole

  logical function mesh_is_half_lat_next_to_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j == 1 .or. j == global_mesh%num_half_lat

  end function mesh_is_half_lat_next_to_pole

  logical function mesh_is_inside_with_halo_full_lat(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j >= this%full_lat_lb .and. j <= this%full_lat_ub

  end function mesh_is_inside_with_halo_full_lat

  logical function mesh_is_inside_with_halo_half_lat(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j >= this%half_lat_lb .and. j <= this%half_lat_ub

  end function mesh_is_inside_with_halo_half_lat

  logical function mesh_is_outside_pole_full_lat(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j < 1 .or. j > global_mesh%num_full_lat

  end function mesh_is_outside_pole_full_lat

  logical function mesh_is_outside_pole_half_lat(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = j < 1 .or. j > global_mesh%num_half_lat

  end function mesh_is_outside_pole_half_lat

  subroutine mesh_final(this)

    type(mesh_type), intent(inout) :: this

		if (allocated(this%dlat            )) deallocate(this%dlat            )
    if (allocated(this%full_dlev       )) deallocate(this%full_dlev       )
    if (allocated(this%half_dlev       )) deallocate(this%half_dlev       )
    if (allocated(this%half_dlev_upper )) deallocate(this%half_dlev_upper )
    if (allocated(this%half_dlev_lower )) deallocate(this%half_dlev_lower )
    if (allocated(this%full_lon        )) deallocate(this%full_lon        )
    if (allocated(this%full_lat        )) deallocate(this%full_lat        )
    if (allocated(this%full_lev        )) deallocate(this%full_lev        )
    if (allocated(this%half_lon        )) deallocate(this%half_lon        )
    if (allocated(this%half_lat        )) deallocate(this%half_lat        )
    if (allocated(this%half_lev        )) deallocate(this%half_lev        )
    if (allocated(this%full_cos_lon    )) deallocate(this%full_cos_lon    )
    if (allocated(this%half_cos_lon    )) deallocate(this%half_cos_lon    )
    if (allocated(this%full_sin_lon    )) deallocate(this%full_sin_lon    )
    if (allocated(this%half_sin_lon    )) deallocate(this%half_sin_lon    )
    if (allocated(this%full_cos_lat    )) deallocate(this%full_cos_lat    )
    if (allocated(this%half_cos_lat    )) deallocate(this%half_cos_lat    )
    if (allocated(this%full_sin_lat    )) deallocate(this%full_sin_lat    )
    if (allocated(this%half_sin_lat    )) deallocate(this%half_sin_lat    )
    if (allocated(this%full_lon_deg    )) deallocate(this%full_lon_deg    )
    if (allocated(this%half_lon_deg    )) deallocate(this%half_lon_deg    )
    if (allocated(this%full_lat_deg    )) deallocate(this%full_lat_deg    )
    if (allocated(this%half_lat_deg    )) deallocate(this%half_lat_deg    )
    if (allocated(this%area_cell       )) deallocate(this%area_cell       )
    if (allocated(this%area_lon        )) deallocate(this%area_lon        )
    if (allocated(this%area_lon_west   )) deallocate(this%area_lon_west   )
    if (allocated(this%area_lon_east   )) deallocate(this%area_lon_east   )
    if (allocated(this%area_lon_north  )) deallocate(this%area_lon_north  )
    if (allocated(this%area_lon_south  )) deallocate(this%area_lon_south  )
    if (allocated(this%area_lat        )) deallocate(this%area_lat        )
    if (allocated(this%area_lat_west   )) deallocate(this%area_lat_west   )
    if (allocated(this%area_lat_east   )) deallocate(this%area_lat_east   )
    if (allocated(this%area_lat_north  )) deallocate(this%area_lat_north  )
    if (allocated(this%area_lat_south  )) deallocate(this%area_lat_south  )
    if (allocated(this%area_vtx        )) deallocate(this%area_vtx        )
    if (allocated(this%area_subcell    )) deallocate(this%area_subcell    )
    if (allocated(this%de_lon          )) deallocate(this%de_lon          )
    if (allocated(this%de_lat          )) deallocate(this%de_lat          )
    if (allocated(this%le_lat          )) deallocate(this%le_lat          )
    if (allocated(this%le_lon          )) deallocate(this%le_lon          )
    if (allocated(this%full_f          )) deallocate(this%full_f          )
    if (allocated(this%half_f          )) deallocate(this%half_f          )
    if (allocated(this%full_tangent_wgt)) deallocate(this%full_tangent_wgt)
    if (allocated(this%half_tangent_wgt)) deallocate(this%half_tangent_wgt)

  end subroutine mesh_final

end module mesh_mod
