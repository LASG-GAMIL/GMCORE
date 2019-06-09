module sphere_geometry_mod

  use const_mod
  use math_mod
  use log_mod

  implicit none

  private

  public euler_formula
  public cartesian_transform
  public inverse_cartesian_transform
  public rotation_transform
  public inverse_rotation_transform
  public calc_distance
  public calc_area
  public calc_area_with_last_small_arc
  public norm_vector
  public calc_sphere_angle
  public calc_arc_length
  public intersect
  public point_type

  integer, parameter :: ORIENT_LEFT  = 1
  integer, parameter :: ORIENT_RIGHT = 2
  integer, parameter :: ORIENT_ON    = 3

  type point_type
    real lon, lat
    real x, y, z
  contains
    procedure :: copy_coord => point_copy_coord
  end type point_type

  interface cartesian_transform
    module procedure cartesian_transform_1
    module procedure cartesian_transform_2
  end interface cartesian_transform

  interface inverse_cartesian_transform
    module procedure inverse_cartesian_transform_1
    module procedure inverse_cartesian_transform_2
  end interface inverse_cartesian_transform

  interface calc_sphere_angle
    module procedure calc_sphere_angle_1
    module procedure calc_sphere_angle_2
  end interface calc_sphere_angle

  interface calc_arc_length
    module procedure calc_arc_length_1
    module procedure calc_arc_length_2
    module procedure calc_arc_length_3
  end interface calc_arc_length

  interface orient
    module procedure orient1
    module procedure orient2
  end interface orient

contains

  integer function euler_formula(num_cell, num_vertex, num_edge) result(res)

    integer, intent(in), optional :: num_cell
    integer, intent(in), optional :: num_vertex
    integer, intent(in), optional :: num_edge

    if (present(num_cell) .and. present(num_vertex)) then
      res = num_cell + num_vertex - 2
    else if (present(num_cell) .and. present(num_edge)) then
      res = num_edge - num_cell + 2
    else if (present(num_vertex) .and. present(num_edge)) then
      res = num_edge - num_vertex + 2
    end if

  end function euler_formula

  subroutine cartesian_transform_1(lon, lat, x, y, z)

    real, intent(in)  :: lon, lat
    real, intent(out) :: x, y, z

    real cos_lat

    cos_lat = cos(lat)
    x = radius * cos_lat * cos(lon)
    y = radius * cos_lat * sin(lon)
    z = radius * sin(lat)

  end subroutine cartesian_transform_1

  subroutine cartesian_transform_2(point)

    class(point_type), intent(inout) :: point

    real cos_lat

    cos_lat = cos(point%lat)
    point%x = radius * cos_lat * cos(point%lon)
    point%y = radius * cos_lat * sin(point%lon)
    point%z = radius * sin(point%lat)

  end subroutine cartesian_transform_2

  subroutine inverse_cartesian_transform_1(lon, lat, x, y, z)

    real, intent(out) :: lon, lat
    real, intent(in)  :: x, y, z

    lon = atan2(y, x)
    lat = asin(z / radius)

    if (lon < 0.0d0) lon = lon + pi2

  end subroutine inverse_cartesian_transform_1

  subroutine inverse_cartesian_transform_2(point)

    class(point_type), intent(inout) :: point

    point%lon = atan2(point%y, point%x)
    point%lat = asin(point%z / radius)

    if (point%lon < 0.0d0) point%lon = point%lon + pi2

  end subroutine inverse_cartesian_transform_2

  ! ************************************************************************** !
  ! Rotation transform                                                         !
  ! Purpose:                                                                   !
  !   Calculate the rotating transformation and its inverse of the original    !
  !   coordinate system (lon_o,lat_o) to the rotated one (lon_r, lat_r) with   !
  !   the north pole (lon_p,lat_p) defined at the original coordinate system.  !
  ! ************************************************************************** !

  subroutine rotation_transform(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

    real, intent(in) :: lon_p, lat_p ! Rotated pole coordinate
    real, intent(in) :: lon_o, lat_o ! Original coordinate
    real, intent(out), optional :: lon_r, lat_r ! Rotated coordinate

    real tmp1, tmp2, tmp3, dlon

    dlon = lon_o - lon_p
    if (present(lon_r)) then
        tmp1 = cos(lat_o) * sin(dlon)
        tmp2 = cos(lat_o) * sin(lat_p) * cos(dlon) - cos(lat_p) * sin(lat_o)
        lon_r = atan2(tmp1, tmp2)
        if (lon_r < 0.0d0) lon_r = pi2 + lon_r
    end if
    if (present(lat_r)) then
        tmp1 = sin(lat_o) * sin(lat_p)
        tmp2 = cos(lat_o) * cos(lat_p) * cos(dlon)
        tmp3 = tmp1 + tmp2
        tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
        lat_r = asin(tmp3)
    end if

  end subroutine rotation_transform

  subroutine inverse_rotation_transform(lon_p, lat_p, lon_o, lat_o, lon_r, lat_r)

      real, intent(in)  :: lon_p, lat_p ! Rotated pole coordinate
      real, intent(out) :: lon_o, lat_o ! Original coordinate
      real, intent(in)  :: lon_r, lat_r ! Rotated coordinate

      real sin_lon_r, cos_lon_r, sin_lat_r, cos_lat_r, sin_lat_p, cos_lat_p
      real tmp1, tmp2, tmp3

      sin_lon_r = sin(lon_r)
      cos_lon_r = cos(lon_r)
      sin_lat_r = sin(lat_r)
      cos_lat_r = cos(lat_r)
      sin_lat_p = sin(lat_p)
      cos_lat_p = cos(lat_p)

      tmp1 = cos_lat_r * sin_lon_r
      tmp2 = sin_lat_r * cos_lat_p + cos_lat_r * cos_lon_r * sin_lat_p
      ! This trick is due to the inaccuracy of trigonometry calculation.
      if (abs(tmp2) < eps) tmp2 = 0.0d0
      lon_o = atan2(tmp1, tmp2)
      lon_o = lon_p + lon_o
      if (lon_o > pi2) lon_o = lon_o - pi2
      tmp1 = sin_lat_r * sin_lat_p
      tmp2 = cos_lat_r * cos_lat_p * cos_lon_r
      tmp3 = tmp1 - tmp2
      tmp3 = min(max(tmp3, -1.0d0), 1.0d0)
      lat_o = asin(tmp3)

  end subroutine inverse_rotation_transform

  real function calc_distance(lon1, lat1, lon2, lat2) result(res)

    real, intent(in) :: lon1
    real, intent(in) :: lat1
    real, intent(in) :: lon2
    real, intent(in) :: lat2

    res = radius * acos(min(1.0d0, max(-1.0d0, sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))))

  end function calc_distance

  real function calc_area(x, y, z) result(res)

    real, intent(in) :: x(:)
    real, intent(in) :: y(:)
    real, intent(in) :: z(:)

    integer n, im1, i, ip1
    real angle

    n = size(x)
#ifndef NDEBUG
    if (n < 3) then
      call log_error('Spherical polygon number is less than 3!', __FILE__, __LINE__)
    end if
#endif
    res = 0.0
    do i = 1, n
      im1 = merge(i - 1, n, i /= 1)
      ip1 = merge(i + 1, 1, i /= n)
      angle = calc_sphere_angle([x(im1),y(im1),z(im1)], [x(i),y(i),z(i)], [x(ip1),y(ip1),z(ip1)])
      res = res + angle
    end do
    res = radius**2 * (res - (n - 2) * pi)

  end function calc_area

  real function calc_area_with_last_small_arc(x, y, z) result(res)

    real, intent(in) :: x(:)
    real, intent(in) :: y(:)
    real, intent(in) :: z(:)

    integer n
    real xv(3), yv(3), zv(3)
    real lon0, lat0, lon1, lat1, lon2, lat2
    real dlon
    real area1, area2, area3

    if (size(x) /= 3) call log_error('Only support triangle with last edge as small arc!', __FILE__, __LINE__)

    res = calc_area(x, y, z)

    call inverse_cartesian_transform(lon0, lat0, x(1), y(1), z(1))
    call inverse_cartesian_transform(lon1, lat1, x(2), y(2), z(2))
    call inverse_cartesian_transform(lon2, lat2, x(3), y(3), z(3))
    if (lat1 /= lat2) call log_error('Small arc is not valid!', __FILE__, __LINE__)
    dlon = merge(lon2 - lon1, lon1 - lon2, lat0 > lat1)
    if (dlon < 0.0) dlon = dlon + pi2

    xv(1) = 0.0;  yv(1) = 0.0;
    if (lat0 * lat1 >= 0 .and. abs(lat0) > abs(lat1)) then
      ! Point 0 is at the side with the Pole.
      xv(2) = x(2); yv(2) = y(2); zv(2) = z(2)
      xv(3) = x(3); yv(3) = y(3); zv(3) = z(3)
    else
      xv(2) = x(3); yv(2) = y(3); zv(2) = z(3)
      xv(3) = x(2); yv(3) = y(2); zv(3) = z(2)
    end if
    if (lat1 >= 0.0) then
      ! Small arc is at the North Sphere.
      zv(1) = radius
      area1 = radius**2 * dlon * (1.0 - sin(lat1))
    else
      ! Small arc is at the South Sphere.
      zv(1) = -radius
      area1 = radius**2 * dlon * (sin(lat1) + 1.0)
    end if
    area2 = calc_area(xv, yv, zv)
    area3 = area1 - area2
    if (area3 < 0.0 .and. abs(area3) > 1.0e-10) call log_error('Lune area is negative!', __FILE__, __LINE__)

    if (lat0 * lat1 >= 0 .and. abs(lat0) > abs(lat1)) then
      res = res + area3
    else
      res = res - area3
    end if
    if (res < 0.0) call log_error('Failed to calculate area with small arc!', __FILE__, __LINE__)

  end function calc_area_with_last_small_arc

  function norm_vector(x) result(res)

    real, intent(in) :: x(:)
    real res(size(x))

    real n

    n = sqrt(sum(x * x))
    if (n /= 0) then
      res = x / n
    else
      res = x
    end if

  end function norm_vector

  ! Calculate the dihedra angle between plane AB and plane BC.

  real function calc_sphere_angle_1(a, b, c) result(res)

    real, intent(in) :: a(3)
    real, intent(in) :: b(3)
    real, intent(in) :: c(3)

    real nab(3) ! Normal vector of plane AB
    real nbc(3) ! Normal vector of plane BC

    nab = norm_vector(cross_product(a, b))
    nbc = norm_vector(cross_product(b, c))
    res = acos(- max(min(dot_product(nab, nbc), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(nab, nbc), a) < 0.0) res = pi2 - res

  end function calc_sphere_angle_1

  real function calc_sphere_angle_2(a, b, c) result(res)

    class(point_type), intent(in) :: a
    class(point_type), intent(in) :: b
    class(point_type), intent(in) :: c

    real xa(3), xb(3), xc(3)
    real nab(3) ! Normal vector of plane AB
    real nbc(3) ! Normal vector of plane BC

    xa = [a%x,a%y,a%z]
    xb = [b%x,b%y,b%z]
    xc = [c%x,c%y,c%z]
    nab = norm_vector(cross_product(xa, xb))
    nbc = norm_vector(cross_product(xb, xc))
    res = acos(-max(min(dot_product(nab, nbc), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(nab, nbc), xa) < 0.0) res = pi2 - res

  end function calc_sphere_angle_2

  ! Calculate the great circle arc length from A to B by assuming A and B are on the unit sphere surface.

  real function calc_arc_length_1(a, b) result(res)

    real, intent(in) :: a(3)
    real, intent(in) :: b(3)

    res = acos(max(min(dot_product(a, b), 1.0d0), -1.0d0))

  end function calc_arc_length_1

  real function calc_arc_length_2(a, b) result(res)

    class(point_type), intent(in) :: a
    class(point_type), intent(in) :: b

    res = acos(max(min(dot_product([a%x,a%y,a%z], [b%x,b%y,b%z]), 1.0d0), -1.0d0))

  end function calc_arc_length_2

  real function calc_arc_length_3(a, b) result(res)

    class(point_type), intent(in) :: a
    real, intent(in) :: b(3)

    res = acos(max(min(dot_product([a%x,a%y,a%z], b), 1.0d0), -1.0d0))

  end function calc_arc_length_3

  logical function intersect(a, b, c, d, e) result(res)

    class(point_type), intent(in) :: a
    class(point_type), intent(in) :: b
    class(point_type), intent(in) :: c
    class(point_type), intent(in) :: d
    class(point_type), intent(inout) :: e

    real n1(3), n2(3), v(3), r
    real lon1, lon2, lat1, lat2

    n1 = cross_product([a%x,a%y,a%z], [b%x,b%y,b%z])
    n2 = cross_product([c%x,c%y,c%z], [d%x,d%y,d%z])
    v  = cross_product(n1, n2)

    r = sqrt(sum(v * v))

    if (r > eps) then
      v = v / r
    else
      res = .false.
      return
    end if
    res = .true.

    lat1 = asin(v(3))
    lat2 = -lat1
    lon1 = atan2(v(2), v(1))
    lon2 = lon1 - pi

    if (lon1 < 0.0) lon1 = lon1 + pi2
    if (lon1 > pi2) lon1 = lon1 - pi2
    if (lon2 < 0.0) lon2 = lon2 + pi2
    if (lon2 > pi2) lon2 = lon2 - pi2

    e%lon = lon1
    e%lat = lat1
    call cartesian_transform(e)

    ! Check if e is the real intersection. If not use the other one.
    if (dot_product(cross_product([a%x,a%y,a%z], [e%x,e%y,e%z]), cross_product([b%x,b%y,b%z], [e%x,e%y,e%z])) >= 0 .or. &
        dot_product(cross_product([c%x,c%y,c%z], [e%x,e%y,e%z]), cross_product([d%x,d%y,d%z], [e%x,e%y,e%z])) >= 0) then
      e%lon = lon2
      e%lat = lat2
      call cartesian_transform(e)
    end if

  end function intersect

    integer function orient1(x1, y1, z1, x2, y2, z2, x0, y0, z0) result(res)

    real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x0, y0, z0

    real(8) det

    det = x0 * (y1 * z2 - y2 * z1) - y0 * (x1 * z2 - x2 * z1) + z0 * (x1 * y2 - x2 * y1)

    if (det > eps) then
      res = ORIENT_LEFT
    else if (-det > eps) then
      res = ORIENT_RIGHT
    else
      res = ORIENT_ON
    end if

  end function orient1

  integer function orient2(p1, p2, p0) result(res)

    class(point_type), intent(in) :: p1, p2, p0

    res = orient1(p1%x, p1%y, p1%z, p2%x, p2%y, p2%z, p0%x, p0%y, p0%z)

  end function orient2

  subroutine point_copy_coord(a, b)

    class(point_type), intent(inout) :: a
    class(point_type), intent(in)    :: b

    a%lon = b%lon
    a%lat = b%lat
    a%x   = b%x
    a%y   = b%y
    a%z   = b%z

  end subroutine point_copy_coord

end module sphere_geometry_mod
