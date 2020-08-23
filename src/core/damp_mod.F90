module damp_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use operators_mod

  implicit none

  private

  public damp_init
  public damp_final
  public zonal_damp
  public latlon_damp_lon
  public latlon_damp_lat
  public latlon_damp_cell
  public div_damp

  integer, parameter :: diff_halo_width(2:8) = [1, 2, 2, 3, 3, 4, 4]

  real(r8) :: diff_weights(9,2:8) = reshape([  &
     1, -2,   1,   0,   0,   0,   0,   0,  0,  & ! 2
    -1,  3, - 3,   1,   0,   0,   0,   0,  0,  & ! 3
     1, -4,   6, - 4,   1,   0,   0,   0,  0,  & ! 4
    -1,  5, -10,  10, - 5,   1,   0,   0,  0,  & ! 5
     1, -6,  15, -20,  15, - 6,   1,   0,  0,  & ! 6
    -1,  7, -21,  35, -35,  21, - 7,   1,  0,  & ! 7
     1, -8,  28, -56,  70, -56,  28, - 8,  1   & ! 8
  ], [9, 7])

  real(r8), allocatable :: cx_full_lat(:,:,:)
  real(r8), allocatable :: cy_full_lat(:,:,:)
  real(r8), allocatable :: cx_half_lat(:,:,:)
  real(r8), allocatable :: cy_half_lat(:,:,:)
  real(r8), allocatable :: cdiv_full_lat(:,:)
  real(r8), allocatable :: cdiv_half_lat(:,:)

contains

  subroutine damp_init()

    integer j, k, o

    integer r

    call damp_final()

    allocate(cx_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cy_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cx_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cy_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))

    do o = 2, 4, 2
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cx_full_lat(j,k,o) = 0.5_r8**o * exp(- 100.0 * (pi05 - abs(global_mesh%full_lat(j)))**2) / dt_in_seconds
          cy_full_lat(j,k,o) = 0.5_r8**o * exp(- 100.0 * (pi05 - abs(global_mesh%full_lat(j)))**2) / dt_in_seconds
        end do
      end do
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cx_half_lat(j,k,o) = 0.5_r8**o * exp(- 100.0 * (pi05 - abs(global_mesh%half_lat(j)))**2) / dt_in_seconds
          cy_half_lat(j,k,o) = 0.5_r8**o * exp(- 100.0 * (pi05 - abs(global_mesh%half_lat(j)))**2) / dt_in_seconds
        end do
      end do
    end do

    allocate(cdiv_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend))
    allocate(cdiv_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend))

    select case (div_damp_order)
    case (2)
      r = 1

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cdiv_full_lat(j,k) = div_damp_coef2 * global_mesh%full_cos_lat(j)**r * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cdiv_half_lat(j,k) = div_damp_coef2 * global_mesh%half_cos_lat(j)**r * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do
    case (4)
      r = 2

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cdiv_full_lat(j,k) = div_damp_coef4 * global_mesh%full_cos_lat(j)**r * &
            radius**4 * global_mesh%dlat(j)**2 * global_mesh%dlon**2 / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cdiv_half_lat(j,k) = div_damp_coef4 * global_mesh%half_cos_lat(j)**r * &
            radius**4 * global_mesh%dlat(j)**2 * global_mesh%dlon**2 / dt_in_seconds
        end do
      end do
    end select

  end subroutine damp_init

  subroutine damp_final()

    if (allocated(cx_full_lat  )) deallocate(cx_full_lat  )
    if (allocated(cy_full_lat  )) deallocate(cy_full_lat  )
    if (allocated(cx_half_lat  )) deallocate(cx_half_lat  )
    if (allocated(cy_half_lat  )) deallocate(cy_half_lat  )
    if (allocated(cdiv_full_lat)) deallocate(cdiv_full_lat)
    if (allocated(cdiv_half_lat)) deallocate(cdiv_half_lat)

  end subroutine damp_final

  subroutine zonal_damp(block, order, dt, dx, lb, ub, n, f, async)

    type(block_type), intent(in) :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(in) :: dx
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    integer, intent(in) :: n
    real(r8), intent(inout) :: f(lb:ub)
    type(async_type), intent(inout), optional :: async

    integer i, ns
    real(r8) g(lb:ub)
    real(r8) df(lb:ub)
    real(r8) a, w(9)

    a  = (dx / 2.0_r8)**order / dt

    call fill_zonal_halo(block, 1 - lb, f)
    if (order == 2) then
      ns = diff_halo_width(order)
      w  = diff_weights(:,order)
      do i = 1, n
        g(i) = sum(f(i-ns:i+ns) * w(:2*ns+1))
      end do
      g = g * (-1)**(order / 2 + 1) * a / dx**order
      do i = 1, n
        f(i) = f(i) + dt * g(i)
      end do
    else
      ns = diff_halo_width(order - 1)
      w  = diff_weights(:,order - 1)
      do i = 1, n
        g (i) = sum(f(i-ns+1:i+ns) * w(:2*ns))
        df(i) = f(i+1) - f(i)
      end do
      g = g * (-1)**(order / 2) * a / dx**order
      do i = 1, n
        g(i) = g(i) * max(0.0_r8, sign(1.0_r8, -g(i) * df(i)))
      end do
      call fill_zonal_halo(block, 1 - lb, g, east_halo=.false.)
      do i = 1, n
        f(i) = f(i) - dt * (g(i) - g(i-1))
      end do
    end if
    if (present(async)) then
      call fill_zonal_halo(block, 1 - lb, f, async)
    else
      call fill_zonal_halo(block, 1 - lb, f)
    end if

  end subroutine zonal_damp

  subroutine latlon_damp_lon(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy, dfdx, dfdy
    integer i, j, k, ns, cyc
    real(r8) w(9)

    mesh => block%mesh

    if (order == 2) then
    else
      ns = diff_halo_width(order - 1)
      w  = diff_weights(:,order - 1)
      gx   => block%latlon_damp_lon_gx
      gy   => block%latlon_damp_lon_gy
      dfdx => block%latlon_damp_lon_dfdx
      dfdy => block%latlon_damp_lon_dfdy
      cycle_loop: do cyc = 1, damp_cycles
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          ! Calculate damping flux at interfaces.
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              !
              !                               i
              !   o      |      o      |      o      |      o      |      o
              !         i-2           i-1            i            i+1
              !                               ^
              !   o - cell   | - edge
              !
              gx  (i,j,k) = sum(f(i-ns:i+ns-1,j,k) * w(:2*ns))
              dfdx(i,j,k) = f(i,j,k) - f(i-1,j,k)
            end do
          end do
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
              !
              !                               j
              !   |      o      |      o      |      o      |      o      |
              !         j-2           j-1            j            j+1
              !                               ^
              !   o - cell   | - edge
              !
              gy  (i,j,k) = sum(f(i,j-ns:j+ns-1,k) * w(:2*ns))
              dfdy(i,j,k) = f(i,j,k) - f(i,j-1,k)
#else
              !
              !                               j
              !   |      o      |      o      |      o      |      o      |
              !         j-1            j            j+1           j+2
              !                               ^
              !   o - cell   | - edge
              !
              gy  (i,j,k) = sum(f(i,j-ns+1:j+ns,k) * w(:2*ns))
              dfdy(i,j,k) = f(i,j+1,k) - f(i,j,k)
#endif
            end do
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            gx(:,j,k) = gx(:,j,k) * (-1)**(order / 2)
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * dfdx(i,j,k)))
            end do
          end do
          call fill_halo(block, gx, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false.)
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            gy(:,j,k) = gy(:,j,k) * (-1)**(order / 2)
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * dfdy(i,j,k)))
            end do
          end do
#ifdef V_POLE
          call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false.)
#else
          call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true., north_halo=.false.)
#endif
          ! Damp physical variable at last.
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
              !
              !               - j+1
              !               |
              !               |
              !        o------|------o
              !        i      |     i+1
              !               |
              !               - j
              !
              f(i,j,k) = f(i,j,k) - dt / damp_cycles * (             &
                (gx(i+1,j,k) - gx(i,j,k)) * cx_full_lat(j,k,order) + &
                (gy(i,j+1,k) - gy(i,j,k)) * cy_full_lat(j,k,order)   &
              )
#else
              !
              !               - j
              !               |
              !               |
              !        o------|------o
              !        i      |     i+1
              !               |
              !               - j-1
              !
              f(i,j,k) = f(i,j,k) - dt / damp_cycles * (             &
                (gx(i+1,j,k) - gx(i,j,k)) * cx_full_lat(j,k,order) + &
                (gy(i,j,k) - gy(i,j-1,k)) * cy_full_lat(j,k,order)   &
              )
#endif
            end do
          end do
        end do
        call fill_halo(block, f, full_lon=.false., full_lat=.true., full_lev=.true.)
      end do cycle_loop
    end if

  end subroutine latlon_damp_lon

  subroutine latlon_damp_lat(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy, dfdx, dfdy
    integer i, j, k, ns, cyc
    real(r8) w(9)

    mesh => block%mesh

    if (order == 2) then
    else
      ns = diff_halo_width(order - 1)
      w  = diff_weights(:,order - 1)
      gx   => block%latlon_damp_lat_gx
      gy   => block%latlon_damp_lat_gy
      dfdx => block%latlon_damp_lat_dfdx
      dfdy => block%latlon_damp_lat_dfdy
      cycle_loop: do cyc = 1, damp_cycles
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          ! Calculate damping flux at interfaces.
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              !
              !                               i
              !   |      o      |      o      |      o      |      o      |
              !         i-1            i            i+1           i+2
              !                               ^
              !   o - cell   | - edge
              !
              gx  (i,j,k) = sum(f(i-ns+1:i+ns,j,k) * w(:2*ns))
              dfdx(i,j,k) = f(i+1,j,k) - f(i,j,k)
            end do
          end do
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              !
              !                               j
              !   o      |      o      |      o      |      o      |      o
              !         j-1            j            j+1           j+2
              !                               ^
              !   o - cell   | - edge
              !
              gy  (i,j,k) = sum(f(i,j-ns+1:j+ns,k) * w(:2*ns))
              dfdy(i,j,k) = f(i,j+1,k) - f(i,j,k)
#else
              !
              !                               j
              !   o      |      o      |      o      |      o      |      o
              !         j-2           j-1            j            j+1
              !                               ^
              !   o - cell   | - edge
              !
              gy  (i,j,k) = sum(f(i,j-ns:j+ns-1,k) * w(:2*ns))
              dfdy(i,j,k) = f(i,j,k) - f(i,j-1,k)
#endif
            end do
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            gx(:,j,k) = gx(:,j,k) * (-1)**(order / 2)
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * dfdx(i,j,k)))
            end do
          end do
          call fill_halo(block, gx, full_lon=.false., full_lat=.false., full_lev=.true., east_halo=.false.)
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            gy(:,j,k) = gy(:,j,k) * (-1)**(order / 2)
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * dfdy(i,j,k)))
            end do
          end do
#ifdef V_POLE
          call fill_halo(block, gy, full_lon=.true., full_lat=.true., full_lev=.true., north_halo=.false.)
#else
          call fill_halo(block, gy, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false.)
#endif
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              !
              !               o j
              !               |
              !               |
              !        |-------------|
              !       i-1     |      i
              !               |
              !               o j-1
              !
              f(i,j,k) = f(i,j,k) - dt / damp_cycles * (             &
                (gx(i,j,k) - gx(i-1,j,k)) * cx_half_lat(j,k,order) + &
                (gy(i,j,k) - gy(i,j-1,k)) * cy_half_lat(j,k,order)   &
              )
#else
              !
              !               o j+1
              !               |
              !               |
              !        |-------------|
              !       i-1     |      i
              !               |
              !               o j
              !
              f(i,j,k) = f(i,j,k) - dt / damp_cycles * (             &
                (gx(i,j,k) - gx(i-1,j,k)) * cx_half_lat(j,k,order) + &
                (gy(i,j+1,k) - gy(i,j,k)) * cy_half_lat(j,k,order)   &
              )
#endif
            end do
          end do
        end do
        call fill_halo(block, f, full_lon=.true., full_lat=.false., full_lev=.true.)
      end do cycle_loop
    end if

  end subroutine latlon_damp_lat

  subroutine latlon_damp_cell(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy, dfdx, dfdy
    integer i, j, k, ns, cyc
    real(r8) w(9)

    mesh => block%mesh

    if (order == 2) then
    else
      ns = diff_halo_width(order - 1)
      w  = diff_weights(:,order - 1)
      gx   => block%latlon_damp_cell_gx
      gy   => block%latlon_damp_cell_gy
      dfdx => block%latlon_damp_cell_dfdx
      dfdy => block%latlon_damp_cell_dfdy
      cycle_loop: do cyc = 1, damp_cycles
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          ! Calculate damping flux at interfaces.
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              !
              !                               i
              !   |      o      |      o      |      o      |      o      |
              !         i-1            i            i+1           i+2
              !                               ^
              !   o - cell   | - edge
              !
              gx  (i,j,k) = sum(f(i-ns+1:i+ns,j,k) * w(:2*ns))
              dfdx(i,j,k) = f(i+1,j,k) - f(i,j,k)
            end do
          end do
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              !
              !                               j
              !   |      o      |      o      |      o      |      o      |
              !         j-2           j-1            j            j+1
              !                               ^
              !   o - cell   | - edge
              !
              gy  (i,j,k) = sum(f(i,j-ns:j+ns-1,k) * w(:2*ns))
              dfdy(i,j,k) = f(i,j,k) - f(i,j-1,k)
#else
              !
              !                               j
              !   |      o      |      o      |      o      |      o      |
              !         j-1            j            j+1           j+2
              !                               ^
              !   o - cell   | - edge
              !
              gy  (i,j,k) = sum(f(i,j-ns+1:j+ns,k) * w(:2*ns))
              dfdy(i,j,k) = f(i,j+1,k) - f(i,j,k)
#endif
            end do
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            gx(:,j,k) = gx(:,j,k) * (-1)**(order / 2)
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * dfdx(i,j,k)))
            end do
          end do
          call fill_halo(block, gx, full_lon=.false., full_lat=.true., full_lev=.true., west_halo=.false.)
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            gy(:,j,k) = gy(:,j,k) * (-1)**(order / 2)
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * dfdy(i,j,k)))
            end do
          end do
#ifdef V_POLE
          call fill_halo(block, gy, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false.)
#else
          call fill_halo(block, gy, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)
#endif
          ! Damp physical variable at last.
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              !
              !               - j+1
              !               |
              !               |
              !        o------|------o
              !        i      |     i+1
              !               |
              !               - j
              !
              f(i,j,k) = f(i,j,k) - dt / damp_cycles * (             &
                (gx(i+1,j,k) - gx(i,j,k)) * cx_full_lat(j,k,order) + &
                (gy(i,j+1,k) - gy(i,j,k)) * cy_full_lat(j,k,order)   &
              )
#else
              !
              !               - j
              !               |
              !               |
              !        o------|------o
              !        i      |     i+1
              !               |
              !               - j-1
              !
              f(i,j,k) = f(i,j,k) - dt / damp_cycles * (             &
                (gx(i+1,j,k) - gx(i,j,k)) * cx_full_lat(j,k,order) + &
                (gy(i,j,k) - gy(i,j-1,k)) * cy_full_lat(j,k,order)   &
              )
#endif
            end do
          end do
        end do
        call fill_halo(block, f, full_lon=.true., full_lat=.true., full_lev=.true.)
      end do cycle_loop
    end if

  end subroutine latlon_damp_cell

  subroutine div_damp(block, old_state, new_state, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, cyc

    if (baroclinic) then
      mesh => old_state%mesh

      select case (div_damp_order)
      case (2)
        new_state%div = old_state%div
        cycle_loop: do cyc = 1, damp_cycles
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                new_state%u(i,j,k) = new_state%u(i,j,k) + dt / damp_cycles * cdiv_full_lat(j,k) * ( &
                  new_state%div(i+1,j,k) - new_state%div(i,j,k)) / mesh%de_lon(j)
              end do
            end do
          end do
          call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
                new_state%v(i,j,k) = new_state%v(i,j,k) + dt / damp_cycles * cdiv_half_lat(j,k) * ( &
                  new_state%div(i,j,k) - new_state%div(i,j-1,k)) / mesh%de_lat(j)
#else
                new_state%v(i,j,k) = new_state%v(i,j,k) + dt / damp_cycles * cdiv_half_lat(j,k) * ( &
                  new_state%div(i,j+1,k) - new_state%div(i,j,k)) / mesh%de_lat(j)
#endif
              end do
            end do
          end do
          call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)

          call calc_div(block, new_state)
        end do cycle_loop
      case (4)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              new_state%u(i,j,k) = new_state%u(i,j,k) - dt * cdiv_full_lat(j,k) * ( &
                old_state%div2(i+1,j,k) - old_state%div2(i,j,k)) / mesh%de_lon(j)
            end do
          end do
        end do
        call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              new_state%v(i,j,k) = new_state%v(i,j,k) - dt * cdiv_half_lat(j,k) * ( &
                old_state%div2(i,j,k) - old_state%div2(i,j-1,k)) / mesh%de_lat(j)
#else
              new_state%v(i,j,k) = new_state%v(i,j,k) - dt * cdiv_half_lat(j,k) * ( &
                old_state%div2(i,j+1,k) - old_state%div2(i,j,k)) / mesh%de_lat(j)
#endif
            end do
          end do
        end do
#ifndef V_POLE
        if (mesh%has_south_pole()) then
          j = mesh%half_lat_ibeg
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%v(i,j,k) = new_state%v(i,j,k) + dt * 3.5d4 * ( &
                old_state%div(i,j+1,k) - old_state%div(i,j,k)) / mesh%de_lat(j)
            end do
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%half_lat_iend
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%v(i,j,k) = new_state%v(i,j,k) + dt * 3.5d4 * ( &
                old_state%div(i,j+1,k) - old_state%div(i,j,k)) / mesh%de_lat(j)
            end do
          end do
        end if
#endif
        call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
      end select
    end if

  end subroutine div_damp

end module damp_mod
