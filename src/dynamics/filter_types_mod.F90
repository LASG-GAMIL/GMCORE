module filter_types_mod

  use const_mod
  use namelist_mod
  use math_mod
  use mesh_mod
  use state_mod

  implicit none

  private

  public filter_type

  type filter_type
    type(mesh_type), pointer :: mesh => null()
    real(r8), allocatable :: width_lon(:)
    integer , allocatable :: ngrid_lon(:)
    real(r8), allocatable :: wgt_lon(:,:)
    real(r8), allocatable :: width_lat(:)
    integer , allocatable :: ngrid_lat(:)
    real(r8), allocatable :: wgt_lat(:,:)
  contains
    procedure :: init  => filter_init
    procedure :: clear => filter_clear
    final :: filter_final
  end type filter_type

contains

  subroutine gaussian_weight(width, ngrid, w)

    real(r8), intent(in) :: width
    integer, intent(in) :: ngrid
    real(r8), intent(out) :: w(:)

    real(r8) s
    integer i, x

    s = filter_coef_b * width / 2.0_r8
    do i = 1, ngrid
      x = i - (ngrid + 1) / 2
      w(i) = exp(-x**2 / (2 * s**2)) / (s * sqrt(pi2))
    end do
    w = w / sum(w)

  end subroutine gaussian_weight

  subroutine filter_init(this, mesh, type)

    class(filter_type), intent(inout) :: this
    type(mesh_type), intent(in), target :: mesh
    character(*), intent(in) :: type

    real(r8) dx, dy, dt, cfl, w, lat0
    integer j, n

    call this%clear()

    this%mesh => mesh
    dt = dt_dyn
    allocate(this%width_lon(mesh%full_lat_lb:mesh%full_lat_ub)); this%width_lon = 0
    allocate(this%ngrid_lon(mesh%full_lat_lb:mesh%full_lat_ub)); this%ngrid_lon = 0
    allocate(this%width_lat(mesh%half_lat_lb:mesh%half_lat_ub)); this%width_lat = 0
    allocate(this%ngrid_lat(mesh%half_lat_lb:mesh%half_lat_ub)); this%ngrid_lat = 0

    if (max_wave_speed > 0) then
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        dx = mesh%de_lon(j)
        dy = mesh%le_lon(j)
        if (dx > 0) then
          cfl = max_wave_speed * dt / dx
          w = filter_coef_a * cfl / max_cfl * (filter_coef_c * (tanh(90 - abs(mesh%full_lat_deg(j))) - 1) + 1)
          n = ceiling(w) + 2; if (mod(n, 2) == 0) n = n + 1
          this%width_lon(j) = w
          this%ngrid_lon(j) = n
        end if
      end do
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        dx = mesh%le_lat(j)
        dy = mesh%de_lat(j)
        if (dx > 0) then
          cfl = max_wave_speed * dt / dx
          w = filter_coef_a * cfl / max_cfl * (filter_coef_c * (tanh(90 - abs(mesh%half_lat_deg(j))) - 1) + 1)
          n = ceiling(w) + 2; if (mod(n, 2) == 0) n = n + 1
          this%width_lat(j) = w
          this%ngrid_lat(j) = n
        end if
      end do
    end if

    allocate(this%wgt_lon(maxval(this%ngrid_lon),mesh%full_lat_lb:mesh%full_lat_ub)); this%wgt_lon = 0
    allocate(this%wgt_lat(maxval(this%ngrid_lat),mesh%half_lat_lb:mesh%half_lat_ub)); this%wgt_lat = 0
    select case (type)
    case ('big_filter')
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (this%ngrid_lon(j) > 1) then
          call gaussian_weight(this%width_lon(j), this%ngrid_lon(j), this%wgt_lon(:,j))
        end if
      end do
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        if (this%ngrid_lat(j) > 1) then
          call gaussian_weight(this%width_lat(j), this%ngrid_lat(j), this%wgt_lat(:,j))
        end if
      end do
    case ('small_filter1')
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (this%ngrid_lon(j) > 1) then
          w = filter_coef_d
          n = w * this%ngrid_lon(j); if (mod(n, 2) == 0) n = n + 1; this%ngrid_lon(j) = n
          this%width_lon(j) = w * this%width_lon(j)
          call gaussian_weight(this%width_lon(j), this%ngrid_lon(j), this%wgt_lon(:,j))
        end if
      end do
      lat0 = -global_mesh%half_lat_deg(1)
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        if (this%ngrid_lat(j) > 1) then
          w = filter_coef_d
          n = w * this%ngrid_lat(j); if (mod(n, 2) == 0) n = n + 1; this%ngrid_lat(j) = n
          this%width_lat(j) = w * this%width_lat(j)
          call gaussian_weight(this%width_lat(j), this%ngrid_lat(j), this%wgt_lat(:,j))
        end if
      end do
    case ('small_filter2')
      lat0 = -global_mesh%full_lat_deg(2)
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (this%ngrid_lon(j) > 1) then
          w = exp_two_values(filter_coef_e, filter_coef_d, lat0, 80.0_r8, abs(mesh%full_lat_deg(j)))
          n = w * this%ngrid_lon(j); if (mod(n, 2) == 0) n = n + 1; this%ngrid_lon(j) = n
          this%width_lon(j) = w * this%width_lon(j)
          call gaussian_weight(this%width_lon(j), this%ngrid_lon(j), this%wgt_lon(:,j))
        end if
      end do
      lat0 = -global_mesh%half_lat_deg(1)
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        if (this%ngrid_lat(j) > 1) then
          w = exp_two_values(filter_coef_e, filter_coef_d, lat0, 80.0_r8, abs(mesh%half_lat_deg(j)))
          n = w * this%ngrid_lat(j); if (mod(n, 2) == 0) n = n + 1; this%ngrid_lat(j) = n
          this%width_lat(j) = w * this%width_lat(j)
          call gaussian_weight(this%width_lat(j), this%ngrid_lat(j), this%wgt_lat(:,j))
        end if
      end do
    case default
      call log_error('Invalid filter type ' // trim(type) // '!')
    end select

  end subroutine filter_init

  subroutine filter_clear(this)

    class(filter_type), intent(inout) :: this

    if (allocated(this%width_lon)) deallocate(this%width_lon)
    if (allocated(this%ngrid_lon)) deallocate(this%ngrid_lon)
    if (allocated(this%wgt_lon  )) deallocate(this%wgt_lon  )
    if (allocated(this%width_lat)) deallocate(this%width_lat)
    if (allocated(this%ngrid_lat)) deallocate(this%ngrid_lat)
    if (allocated(this%wgt_lat  )) deallocate(this%wgt_lat  )

  end subroutine filter_clear

  subroutine filter_final(this)

    type(filter_type), intent(inout) :: this

    call this%clear()

  end subroutine filter_final

end module filter_types_mod
