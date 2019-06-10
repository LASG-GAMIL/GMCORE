module mesh_mod

  use const_mod
  use sphere_geometry_mod
  use log_mod

  implicit none

  private

  public mesh_type
  public mesh
  public create_meshes
  public num_lon
  public num_lat

  type mesh_type
    integer :: id
    integer :: halo_width
    integer num_full_lon
    integer num_half_lon
    integer num_full_lat
    integer num_half_lat
    integer full_lon_start_idx
    integer full_lon_end_idx
    integer full_lat_start_idx
    integer full_lat_end_idx
    integer full_lat_start_idx_no_pole
    integer full_lat_end_idx_no_pole
    integer half_lon_start_idx
    integer half_lon_end_idx
    integer half_lat_start_idx
    integer half_lat_end_idx
    integer full_lon_lb
    integer full_lon_ub
    integer half_lon_lb
    integer half_lon_ub
    integer full_lat_lb
    integer full_lat_ub
    integer half_lat_lb
    integer half_lat_ub
    real(real_kind) start_lon
    real(real_kind) end_lon
    real(real_kind) start_lat
    real(real_kind) end_lat
    real(real_kind) dlon
    real(real_kind) dlat
    real(real_kind), allocatable :: full_lon(:)
    real(real_kind), allocatable :: half_lon(:)
    real(real_kind), allocatable :: full_lat(:)
    real(real_kind), allocatable :: half_lat(:)
    real(real_kind), allocatable :: full_cos_lon(:)
    real(real_kind), allocatable :: half_cos_lon(:)
    real(real_kind), allocatable :: full_sin_lon(:)
    real(real_kind), allocatable :: half_sin_lon(:)
    real(real_kind), allocatable :: full_cos_lat(:)
    real(real_kind), allocatable :: half_cos_lat(:)
    real(real_kind), allocatable :: full_sin_lat(:)
    real(real_kind), allocatable :: half_sin_lat(:)
    ! Area for weighting
    real(real_kind), allocatable :: cell_area(:)
    real(real_kind), allocatable :: lon_edge_area(:)
    real(real_kind), allocatable :: lon_edge_left_area(:)
    real(real_kind), allocatable :: lon_edge_right_area(:)
    real(real_kind), allocatable :: lat_edge_area(:)
    real(real_kind), allocatable :: lat_edge_up_area(:)
    real(real_kind), allocatable :: lat_edge_down_area(:)
    real(real_kind), allocatable :: vertex_area(:)
    real(real_kind), allocatable :: subcell_area(:,:)
    ! Edge length
    real(real_kind), allocatable :: cell_lon_distance(:)
    real(real_kind), allocatable :: cell_lat_distance(:)
    real(real_kind), allocatable :: vertex_lon_distance(:)
    real(real_kind), allocatable :: vertex_lat_distance(:)
    ! For output
    real(real_kind), allocatable :: full_lon_deg(:)
    real(real_kind), allocatable :: half_lon_deg(:)
    real(real_kind), allocatable :: full_lat_deg(:)
    real(real_kind), allocatable :: half_lat_deg(:)
  contains
    procedure :: init => mesh_init
    procedure :: has_south_pole => mesh_has_south_pole
    procedure :: has_north_pole => mesh_has_north_pole
    procedure :: is_south_pole => mesh_is_south_pole
    procedure :: is_north_pole => mesh_is_north_pole
    final :: mesh_final
  end type mesh_type

  integer num_lon
  integer num_lat

  type(mesh_type) mesh

contains

  subroutine mesh_init(this, num_lon, num_lat, id, halo_width, start_lon, end_lon, start_lat, end_lat)

    class(mesh_type), intent(inout) :: this
    integer, intent(in) :: num_lon
    integer, intent(in) :: num_lat
    integer, intent(in), optional :: id
    integer, intent(in), optional :: halo_width
    real(real_kind), intent(in), optional :: start_lon
    real(real_kind), intent(in), optional :: end_lon
    real(real_kind), intent(in), optional :: start_lat
    real(real_kind), intent(in), optional :: end_lat

    real(real_kind) x(3), y(3), z(3), total_area
    integer i, j

    this%num_full_lon = num_lon
    this%num_half_lon = num_lon
    this%num_full_lat = num_lat
    this%num_half_lat = num_lat - 1

    this%full_lon_start_idx = 1
    this%full_lon_end_idx = this%num_full_lon
    this%full_lat_start_idx = 1
    this%full_lat_end_idx = this%num_full_lat
    this%half_lon_start_idx = 1
    this%half_lon_end_idx = this%num_half_lon
    this%half_lat_start_idx = 1
    this%half_lat_end_idx = this%num_half_lat

    this%full_lat_start_idx_no_pole = merge(this%full_lat_start_idx, this%full_lat_start_idx + 1, this%has_south_pole())
    this%full_lat_end_idx_no_pole   = merge(this%full_lat_end_idx,   this%full_lat_end_idx   - 1, this%has_north_pole())

    this%full_lon_lb = this%full_lon_start_idx-this%halo_width
    this%full_lon_ub = this%full_lon_end_idx  +this%halo_width
    this%full_lat_lb = this%full_lat_start_idx-this%halo_width
    this%full_lat_ub = this%full_lat_end_idx  +this%halo_width
    this%half_lon_lb = this%half_lon_start_idx-this%halo_width
    this%half_lon_ub = this%half_lon_end_idx  +this%halo_width
    this%half_lat_lb = this%half_lat_start_idx-this%halo_width
    this%half_lat_ub = this%half_lat_end_idx  +this%halo_width

    this%id         = merge(id,         0,          present(id))
    this%halo_width = merge(halo_width, 1,          present(halo_width))
    this%start_lon  = merge(start_lon,  0.0d0,      present(start_lon))
    this%end_lon    = merge(end_lon,    2.0d0 * pi, present(end_lon))
    this%start_lat  = merge(start_lat, -0.5d0 * pi, present(start_lat))
    this%end_lat    = merge(end_lat,    0.5d0 * pi, present(end_lat))

    allocate(this%full_lon           (this%full_lon_lb:this%full_lon_ub))
    allocate(this%half_lon           (this%half_lon_lb:this%half_lon_ub))
    allocate(this%full_lat           (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_lat           (this%half_lat_lb:this%half_lat_ub))
    allocate(this%full_cos_lon       (this%full_lon_lb:this%full_lon_ub))
    allocate(this%half_cos_lon       (this%half_lon_lb:this%half_lon_ub))
    allocate(this%full_sin_lon       (this%full_lon_lb:this%full_lon_ub))
    allocate(this%half_sin_lon       (this%half_lon_lb:this%half_lon_ub))
    allocate(this%full_cos_lat       (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_cos_lat       (this%half_lat_lb:this%half_lat_ub))
    allocate(this%full_sin_lat       (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_sin_lat       (this%half_lat_lb:this%half_lat_ub))
    allocate(this%cell_area          (this%full_lat_lb:this%full_lat_ub))
    allocate(this%lon_edge_area      (this%full_lat_lb:this%full_lat_ub))
    allocate(this%lon_edge_left_area (this%full_lat_lb:this%full_lat_ub))
    allocate(this%lon_edge_right_area(this%full_lat_lb:this%full_lat_ub))
    allocate(this%lat_edge_area      (this%half_lat_lb:this%half_lat_ub))
    allocate(this%lat_edge_up_area   (this%half_lat_lb:this%half_lat_ub))
    allocate(this%lat_edge_down_area (this%half_lat_lb:this%half_lat_ub))
    allocate(this%vertex_area        (this%half_lat_lb:this%half_lat_ub))
    allocate(this%subcell_area     (2,this%full_lat_lb:this%full_lat_ub))
    allocate(this%cell_lon_distance  (this%full_lat_lb:this%full_lat_ub))
    allocate(this%cell_lat_distance  (this%half_lat_lb:this%half_lat_ub))
    allocate(this%vertex_lon_distance(this%half_lat_lb:this%half_lat_ub))
    allocate(this%vertex_lat_distance(this%full_lat_lb:this%full_lat_ub))
    allocate(this%full_lon_deg       (this%full_lon_lb:this%full_lon_ub))
    allocate(this%half_lon_deg       (this%half_lon_lb:this%half_lon_ub))
    allocate(this%full_lat_deg       (this%full_lat_lb:this%full_lat_ub))
    allocate(this%half_lat_deg       (this%half_lat_lb:this%half_lat_ub))

    this%dlon = (this%end_lon - this%start_lon) / this%num_full_lon
    do i = this%full_lon_lb, this%full_lon_ub
      this%full_lon(i) = this%start_lon + (i - 1) * this%dlon
      this%half_lon(i) = this%full_lon(i) + 0.5d0 * this%dlon
      this%full_lon_deg(i) = this%full_lon(i) * deg
      this%half_lon_deg(i) = this%half_lon(i) * deg
    end do

    this%dlat = (this%end_lat - this%start_lat) / this%num_half_lat
    do j = this%full_lat_lb, this%full_lat_ub
      this%full_lat(j) = this%start_lat + (j - 1) * this%dlat
      this%full_lat_deg(j) = this%full_lat(j) * deg
    end do
    this%full_lat(this%num_full_lat) = this%end_lat
    this%full_lat_deg(this%num_full_lat) = this%end_lat * deg

    do j = this%half_lat_lb, this%half_lat_ub
      this%half_lat(j) = this%full_lat(j) + 0.5d0 * this%dlat
      this%half_lat_deg(j) = this%half_lat(j) * deg
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
      this%half_cos_lat(j) = cos(this%half_lat(j))
      this%half_sin_lat(j) = sin(this%half_lat(j))
    end do

    do j = this%full_lat_lb, this%full_lat_ub
      this%full_cos_lat(j) = cos(this%full_lat(j))
      this%full_sin_lat(j) = sin(this%full_lat(j))
    end do

    do j = this%full_lat_start_idx, this%full_lat_end_idx
      if (this%is_south_pole(j)) then
        this%cell_area(j) = radius**2 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
        this%subcell_area(1,j) = 0.0d0
        this%subcell_area(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) + 1.0d0)
      else if (this%is_north_pole(j)) then
        this%cell_area(j) = radius**2 * this%dlon * (1.0 - this%half_sin_lat(j-1))
        this%subcell_area(1,j) = radius**2 * 0.5 * this%dlon * (1.0d0 - this%half_sin_lat(j-1))
        this%subcell_area(2,j) = 0.0d0
      else
        this%cell_area(j) = radius**2 * this%dlon * (this%half_sin_lat(j) - this%half_sin_lat(j-1))
        this%subcell_area(1,j) = radius**2 * 0.5d0 * this%dlon * (this%full_sin_lat(j) - this%half_sin_lat(j-1))
        this%subcell_area(2,j) = radius**2 * 0.5d0 * this%dlon * (this%half_sin_lat(j) - this%full_sin_lat(j))
        call cartesian_transform(this%full_lon(1), this%full_lat(j  ), x(1), y(1), z(1))
        call cartesian_transform(this%half_lon(1), this%half_lat(j-1), x(2), y(2), z(2))
        call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(3), y(3), z(3))
        this%lon_edge_left_area(j) = calc_area(x, y, z)
        this%lon_edge_right_area(j) = this%lon_edge_left_area(j)
        this%lon_edge_area(j) = this%lon_edge_left_area(j) + this%lon_edge_right_area(j)
      end if
      ! print *, j, this%lon_edge_left_area(j), this%lon_edge_right_area(j), this%lon_edge_area(j)
    end do

    do j = this%half_lat_start_idx, this%half_lat_end_idx
      this%vertex_area(j) = radius**2 * this%dlon * (this%full_sin_lat(j+1) - this%full_sin_lat(j))
      this%lat_edge_area(j) = this%vertex_area(j)
      call cartesian_transform(this%full_lon(2), this%full_lat(j+1), x(1), y(1), z(1))
      call cartesian_transform(this%half_lon(1), this%half_lat(j  ), x(2), y(2), z(2))
      call cartesian_transform(this%half_lon(2), this%half_lat(j  ), x(3), y(3), z(3))
      this%lat_edge_up_area(j) = calc_area_with_last_small_arc(x, y, z)
      call cartesian_transform(this%full_lon(2), this%full_lat(j), x(1), y(1), z(1))
      call cartesian_transform(this%half_lon(2), this%half_lat(j), x(2), y(2), z(2))
      call cartesian_transform(this%half_lon(1), this%half_lat(j), x(3), y(3), z(3))
      this%lat_edge_down_area(j) = calc_area_with_last_small_arc(x, y, z)
      this%lat_edge_area(j) = this%lat_edge_up_area(j) + this%lat_edge_down_area(j)
    end do

    total_area = 0.0d0
    do j = this%full_lat_start_idx, this%full_lat_end_idx
      total_area = total_area + this%cell_area(j) * this%num_full_lon
    end do
    if (abs((4 * pi * radius**2 - total_area) / (4 * pi * radius**2)) > 1.0d-12) then
      call log_error('Failed to calculate cell area!', __FILE__, __LINE__)
    end if

    total_area = 0.0d0
    do j = this%half_lat_start_idx, this%half_lat_end_idx
      total_area = total_area + this%vertex_area(j) * this%num_full_lon
    end do
    if (abs((4 * pi * radius**2 - total_area) / (4 * pi * radius**2)) > 1.0d-12) then
      call log_error('Failed to calculate vertex area!', __FILE__, __LINE__)
    end if

    total_area = 0.0d0
    do j = this%full_lat_start_idx, this%full_lat_end_idx
      total_area = total_area + this%lon_edge_area(j) * this%num_full_lon
    end do
    do j = this%half_lat_start_idx, this%half_lat_end_idx
      total_area = total_area + this%lat_edge_area(j) * this%num_full_lon
    end do
    if (abs((4 * pi * radius**2 - total_area) / (4 * pi * radius**2)) > 1.0d-12) then
      call log_error('Failed to calculate edge area!', __FILE__, __LINE__)
    end if

    do j = this%full_lat_start_idx_no_pole, this%full_lat_end_idx_no_pole
      this%vertex_lat_distance(j) = this%dlat * radius
      this%cell_lon_distance(j) = 2.0d0 * this%lon_edge_area(j) / this%vertex_lat_distance(j)
    end do

    do j = this%half_lat_start_idx, this%half_lat_end_idx
      this%vertex_lon_distance(j) = this%dlon * radius * this%half_cos_lat(j)
      this%cell_lat_distance(j) = 2.0d0 * this%lat_edge_area(j) / this%vertex_lon_distance(j)
    end do

  end subroutine mesh_init

  logical function mesh_has_south_pole(this) result(res)

    class(mesh_type), intent(in) :: this

    res = this%start_lat == -0.5 * pi

  end function mesh_has_south_pole

  logical function mesh_has_north_pole(this) result(res)

    class(mesh_type), intent(in) :: this

    res = this%end_lat == 0.5 * Pi

  end function mesh_has_north_pole

  logical function mesh_is_south_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = this%has_south_pole() .and. j == 1

  end function mesh_is_south_pole

  logical function mesh_is_north_pole(this, j) result(res)

    class(mesh_type), intent(in) :: this
    integer, intent(in) :: j

    res = this%has_north_pole() .and. j == this%num_full_lat

  end function mesh_is_north_pole

  subroutine mesh_final(this)

    type(mesh_type), intent(inout) :: this

    if (allocated(this%full_lon))            deallocate(this%full_lon)
    if (allocated(this%full_lat))            deallocate(this%full_lat)
    if (allocated(this%half_lon))            deallocate(this%half_lon)
    if (allocated(this%half_lat))            deallocate(this%half_lat)
    if (allocated(this%full_cos_lon))        deallocate(this%full_cos_lon)
    if (allocated(this%half_cos_lon))        deallocate(this%half_cos_lon)
    if (allocated(this%full_sin_lon))        deallocate(this%full_sin_lon)
    if (allocated(this%half_sin_lon))        deallocate(this%half_sin_lon)
    if (allocated(this%full_cos_lat))        deallocate(this%full_cos_lat)
    if (allocated(this%half_cos_lat))        deallocate(this%half_cos_lat)
    if (allocated(this%full_sin_lat))        deallocate(this%full_sin_lat)
    if (allocated(this%half_sin_lat))        deallocate(this%half_sin_lat)
    if (allocated(this%cell_area))           deallocate(this%cell_area)
    if (allocated(this%lon_edge_area))       deallocate(this%lon_edge_area)
    if (allocated(this%lon_edge_left_area))  deallocate(this%lon_edge_left_area)
    if (allocated(this%lon_edge_right_area)) deallocate(this%lon_edge_right_area)
    if (allocated(this%lat_edge_area))       deallocate(this%lat_edge_area)
    if (allocated(this%lat_edge_up_area))    deallocate(this%lat_edge_up_area)
    if (allocated(this%lat_edge_down_area))  deallocate(this%lat_edge_down_area)
    if (allocated(this%vertex_area))         deallocate(this%vertex_area)
    if (allocated(this%subcell_area))        deallocate(this%subcell_area)
    if (allocated(this%cell_lon_distance))   deallocate(this%cell_lon_distance)
    if (allocated(this%cell_lat_distance))   deallocate(this%cell_lat_distance)
    if (allocated(this%vertex_lon_distance)) deallocate(this%vertex_lon_distance)
    if (allocated(this%vertex_lat_distance)) deallocate(this%vertex_lat_distance)
    if (allocated(this%full_lon_deg))        deallocate(this%full_lon_deg)
    if (allocated(this%half_lon_deg))        deallocate(this%half_lon_deg)
    if (allocated(this%full_lat_deg))        deallocate(this%full_lat_deg)
    if (allocated(this%half_lat_deg))        deallocate(this%half_lat_deg)

  end subroutine mesh_final

  subroutine create_meshes()

    call mesh%init(num_lon, num_lat)

  end subroutine create_meshes

end module mesh_mod
