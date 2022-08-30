module filter_types_mod

  use const_mod
  use namelist_mod
  use mesh_mod
  use state_mod

  implicit none

  private

  public filter_type

  type filter_type
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

  subroutine filter_init(this, mesh)

    class(filter_type), intent(inout) :: this
    type(mesh_type), intent(in) :: mesh

    real(r8) dx, dy, dt, cfl, w, s, x
    integer j, l, n

    call this%clear()

    dt = dt_dyn

    allocate(this%width_lon(mesh%full_lat_lb:mesh%full_lat_ub)); this%width_lon = 0
    allocate(this%ngrid_lon(mesh%full_lat_lb:mesh%full_lat_ub)); this%ngrid_lon = 0
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
    end if
    allocate(this%wgt_lon(maxval(this%ngrid_lon),mesh%full_lat_lb:mesh%full_lat_ub)); this%wgt_lon = 0
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (this%ngrid_lon(j) > 1) then
        s = filter_coef_b * this%width_lon(j) / 2.0_r8
        do l = 1, this%ngrid_lon(j)
          x = l - (this%ngrid_lon(j) + 1) / 2
          this%wgt_lon(l,j) = exp(-x**2 / (2 * s**2)) / (s * sqrt(pi2))
        end do
        this%wgt_lon(:,j) = this%wgt_lon(:,j) / sum(this%wgt_lon(:,j))
      end if
    end do

    allocate(this%width_lat(mesh%half_lat_lb:mesh%half_lat_ub)); this%width_lat = 0
    allocate(this%ngrid_lat(mesh%half_lat_lb:mesh%half_lat_ub)); this%ngrid_lat = 0
    if (max_wave_speed > 0) then
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
    allocate(this%wgt_lat(maxval(this%ngrid_lat),mesh%half_lat_lb:mesh%half_lat_ub)); this%wgt_lat = 0
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      if (this%ngrid_lat(j) > 1) then
        s = filter_coef_b * this%width_lat(j) / 2.0_r8
        do l = 1, this%ngrid_lat(j)
          x = l - (this%ngrid_lat(j) + 1) / 2
          this%wgt_lat(l,j) = exp(-x**2 / (2 * s**2)) / (s * sqrt(pi2))
        end do
        this%wgt_lat(:,j) = this%wgt_lat(:,j) / sum(this%wgt_lat(:,j))
      end if
    end do

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
