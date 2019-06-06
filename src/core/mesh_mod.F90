module mesh_mod

  use const_mod
  use sphere_geometry_mod
  use log_mod

  implicit none

  private

  public mesh_type
  public mesh
  public create_meshes
  public num_total_lon
  public num_total_lat

  type mesh_type
    integer :: id
    integer :: halo_width
    integer num_full_lon
    integer num_half_lon
    integer num_full_lat
    integer num_half_lat
    real start_lon
    real end_lon
    real start_lat
    real end_lat
    real dlon
    real dlat
    real, allocatable :: full_lon(:)
    real, allocatable :: half_lon(:)
    real, allocatable :: full_lat(:)
    real, allocatable :: half_lat(:)
    real, allocatable :: full_cos_lon(:)
    real, allocatable :: half_cos_lon(:)
    real, allocatable :: full_sin_lon(:)
    real, allocatable :: half_sin_lon(:)
    real, allocatable :: full_cos_lat(:)
    real, allocatable :: half_cos_lat(:)
    real, allocatable :: full_sin_lat(:)
    real, allocatable :: half_sin_lat(:)
    ! Area for weighting
    real, allocatable :: cell_area(:)
    real, allocatable :: lon_edge_area(:)
    real, allocatable :: lon_edge_left_area(:)
    real, allocatable :: lon_edge_right_area(:)
    real, allocatable :: lat_edge_area(:)
    real, allocatable :: lat_edge_up_area(:)
    real, allocatable :: lat_edge_down_area(:)
    real, allocatable :: vertex_area(:)
    real, allocatable :: subcell_area(:,:)
    ! For output
    real, allocatable :: full_lon_deg(:)
    real, allocatable :: half_lon_deg(:)
    real, allocatable :: full_lat_deg(:)
    real, allocatable :: half_lat_deg(:)
  contains
    procedure :: init => mesh_init
    procedure :: has_south_pole => mesh_has_south_pole
    procedure :: has_north_pole => mesh_has_north_pole
    final :: mesh_final
  end type mesh_type

  integer num_total_lon
  integer num_total_lat

  type(mesh_type) mesh

contains

  subroutine mesh_init(this, num_lon, num_lat, id, halo_width, start_lon, end_lon, start_lat, end_lat)

    class(mesh_type), intent(inout) :: this
    integer, intent(in) :: num_lon
    integer, intent(in) :: num_lat
    integer, intent(in), optional :: id
    integer, intent(in), optional :: halo_width
    real, intent(in), optional :: start_lon
    real, intent(in), optional :: end_lon
    real, intent(in), optional :: start_lat
    real, intent(in), optional :: end_lat

    real x(3), y(3), z(3), total_area
    integer i, j

    this%num_full_lon = num_lon
    this%num_half_lon = num_lon
    this%num_full_lat = num_lat
    this%num_half_lat = num_lat - 1


    this%id         = merge(id,         0,        present(id))
    this%halo_width = merge(halo_width, 1,        present(halo_width))
    this%start_lon  = merge(start_lon,  0.0,      present(start_lon))
    this%end_lon    = merge(end_lon,    2.0 * pi, present(end_lon))
    this%start_lat  = merge(start_lat, -0.5 * pi, present(start_lat))
    this%end_lat    = merge(end_lat,    0.5 * pi, present(end_lat))

    allocate(this%full_lon(this%num_full_lon))
    allocate(this%half_lon(this%num_half_lon))
    allocate(this%full_lat(1-this%halo_width:this%num_full_lat+this%halo_width))
    allocate(this%half_lat(1-this%halo_width:this%num_half_lat+this%halo_width))
    allocate(this%full_cos_lon(this%num_full_lon))
    allocate(this%half_cos_lon(this%num_half_lon))
    allocate(this%full_sin_lon(this%num_full_lon))
    allocate(this%half_sin_lon(this%num_half_lon))
    allocate(this%full_cos_lat(1-this%halo_width:this%num_full_lat+this%halo_width))
    allocate(this%half_cos_lat(1-this%halo_width:this%num_half_lat+this%halo_width))
    allocate(this%full_sin_lat(1-this%halo_width:this%num_full_lat+this%halo_width))
    allocate(this%half_sin_lat(1-this%halo_width:this%num_half_lat+this%halo_width))
    allocate(this%cell_area(this%num_full_lat))
    allocate(this%lon_edge_area(this%num_full_lat))
    allocate(this%lon_edge_left_area(this%num_full_lat))
    allocate(this%lon_edge_right_area(this%num_full_lat))
    allocate(this%lat_edge_area(this%num_half_lat))
    allocate(this%lat_edge_up_area(this%num_half_lat))
    allocate(this%lat_edge_down_area(this%num_half_lat))
    allocate(this%vertex_area(this%num_half_lat))
    allocate(this%subcell_area(2,this%num_full_lat))
    allocate(this%full_lon_deg(this%num_full_lon))
    allocate(this%half_lon_deg(this%num_half_lon))
    allocate(this%full_lat_deg(this%num_full_lat))
    allocate(this%half_lat_deg(this%num_half_lat))

    this%dlon = (this%end_lon - this%start_lon) / this%num_full_lon
    do i = 1, this%num_full_lon
      this%full_lon(i) = this%start_lon + (i - 1) * this%dlon
      this%half_lon(i) = this%full_lon(i) + 0.5 * this%dlon
      this%full_lon_deg(i) = this%full_lon(i) * deg
      this%half_lon_deg(i) = this%half_lon(i) * deg
    end do

    this%dlat = (this%end_lat - this%start_lat) / this%num_half_lat
    do j = 1 - this%halo_width, this%num_full_lat + this%halo_width
      this%full_lat(j) = this%start_lat + (j - 1) * this%dlat
      this%full_lat_deg(j) = this%full_lat(j) * deg
    end do
    this%full_lat(this%num_full_lat) = this%end_lat
    this%full_lat_deg(this%num_full_lat) = this%end_lat * deg

    do j = 1 - this%halo_width, this%num_half_lat + this%halo_width
      this%half_lat(j) = this%full_lat(j) + 0.5 * this%dlat
      this%half_lat_deg(j) = this%half_lat(j) * deg
    end do

    do i = 1, this%num_full_lon
      this%full_cos_lon(i) = cos(this%full_lon(i))
      this%full_sin_lon(i) = sin(this%full_lon(i))
    end do

    do i = 1, this%num_half_lon
      this%half_cos_lon(i) = cos(this%half_lon(i))
      this%half_sin_lon(i) = sin(this%half_lon(i))
    end do

    do j = 1 - this%halo_width, this%num_half_lat + this%halo_width
      this%half_cos_lat(j) = cos(this%half_lat(j))
      this%half_sin_lat(j) = sin(this%half_lat(j))
    end do

    do j = 1 - this%halo_width, this%num_full_lat + this%halo_width
      this%full_cos_lat(j) = cos(this%full_lat(j))
      this%full_sin_lat(j) = sin(this%full_lat(j))
    end do

    do j = 1, this%num_full_lat
      if (this%has_south_pole() .and. j == 1) then
        this%cell_area(j) = radius**2 * this%dlon * (this%half_sin_lat(j) + 1.0)
        this%subcell_area(1,j) = 0
        this%subcell_area(2,j) = radius**2 * 0.5 * this%dlon * (this%half_sin_lat(j) + 1.0)
      else if (this%has_north_pole() .and. j == this%num_full_lat) then
        this%cell_area(j) = radius**2 * this%dlon * (1.0 - this%half_sin_lat(j-1))
        this%subcell_area(1,j) = radius**2 * 0.5 * this%dlon * (1.0 - this%half_sin_lat(j-1))
        this%subcell_area(2,j) = 0
      else
        this%cell_area(j) = radius**2 * this%dlon * (this%half_sin_lat(j) - this%half_sin_lat(j-1))
        this%subcell_area(1,j) = radius**2 * 0.5 * this%dlon * (this%full_sin_lat(j) - this%half_sin_lat(j-1))
        this%subcell_area(2,j) = radius**2 * 0.5 * this%dlon * (this%half_sin_lat(j) - this%full_sin_lat(j))
        call cartesian_transform(mesh%full_lon(1), mesh%full_lat(j  ), x(1), y(1), z(1))
        call cartesian_transform(mesh%half_lon(1), mesh%half_lat(j-1), x(2), y(2), z(2))
        call cartesian_transform(mesh%half_lon(1), mesh%half_lat(j  ), x(3), y(3), z(3))
        this%lon_edge_left_area(j) = calc_area(x, y, z)
        this%lon_edge_right_area(j) = this%lon_edge_left_area(j)
        this%lon_edge_area(j) = this%lon_edge_left_area(j) + this%lon_edge_right_area(j)
      end if
    end do

    do j = 1, this%num_half_lat
      this%vertex_area(j) = radius**2 * this%dlon * (this%full_sin_lat(j+1) - this%full_sin_lat(j))
      this%lat_edge_area(j) = this%vertex_area(j)
      call cartesian_transform(mesh%full_lon(2), mesh%full_lat(j+1), x(1), y(1), z(1))
      call cartesian_transform(mesh%half_lon(1), mesh%half_lat(j  ), x(2), y(2), z(2))
      call cartesian_transform(mesh%half_lon(2), mesh%half_lat(j  ), x(3), y(3), z(3))
      this%lat_edge_up_area(j) = calc_area_with_last_small_arc(x, y, z)
      call cartesian_transform(mesh%full_lon(2), mesh%full_lat(j), x(1), y(1), z(1))
      call cartesian_transform(mesh%half_lon(2), mesh%half_lat(j), x(2), y(2), z(2))
      call cartesian_transform(mesh%half_lon(1), mesh%half_lat(j), x(3), y(3), z(3))
      this%lat_edge_down_area(j) = calc_area_with_last_small_arc(x, y, z)
      this%lat_edge_area(j) = this%lat_edge_up_area(j) + this%lat_edge_down_area(j)
    end do

    total_area = 0.0
    do j = 1, this%num_full_lat
      total_area = total_area + this%cell_area(j) * this%num_full_lon
    end do
    if (abs((4 * pi * radius**2 - total_area) / (4 * pi * radius**2)) > 1.0e-12) then
      call log_error('Failed to calculate cell area!', __FILE__, __LINE__)
    end if

    total_area = 0.0
    do j = 1, this%num_full_lat
      total_area = total_area + this%lon_edge_area(j) * this%num_full_lon
    end do
    do j = 1, this%num_half_lat
      total_area = total_area + this%lat_edge_area(j) * this%num_full_lon
    end do
    if (abs((4 * pi * radius**2 - total_area) / (4 * pi * radius**2)) > 1.0e-12) then
      call log_error('Failed to calculate edge area!', __FILE__, __LINE__)
    end if

  end subroutine mesh_init

  logical function mesh_has_south_pole(this) result(res)

    class(mesh_type), intent(in) :: this

    res = this%start_lat == -0.5 * pi

  end function mesh_has_south_pole

  logical function mesh_has_north_pole(this) result(res)

    class(mesh_type), intent(in) :: this

    res = this%end_lat == 0.5 * Pi

  end function mesh_has_north_pole

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
    if (allocated(this%full_lon_deg))        deallocate(this%full_lon_deg)
    if (allocated(this%half_lon_deg))        deallocate(this%half_lon_deg)
    if (allocated(this%full_lat_deg))        deallocate(this%full_lat_deg)
    if (allocated(this%half_lat_deg))        deallocate(this%half_lat_deg)

  end subroutine mesh_final

  subroutine create_meshes()

    call mesh%init(num_total_lon, num_total_lat)

  end subroutine create_meshes

end module mesh_mod
