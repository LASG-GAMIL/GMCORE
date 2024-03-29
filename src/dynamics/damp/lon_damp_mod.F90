module lon_damp_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public lon_damp_init
  public lon_damp_final
  public lon_damp_on_lon_edge
  public lon_damp_on_lat_edge
  public lon_damp_on_lev_edge
  public lon_damp_on_cell
  public lon_damp_on_vtx
  public lon_damp_1d

  real(r8), allocatable :: wgt_full_lat(:)
  real(r8), allocatable :: wgt_half_lat(:)

  integer, parameter :: diff_halo_width(2:8) = [1, 2, 2, 3, 3, 4, 4]

  real(r8), target :: diff_weights(9,2:8) = reshape([  &
     1, -2,   1,   0,   0,   0,   0,   0,  0,  & ! 2
    -1,  3, - 3,   1,   0,   0,   0,   0,  0,  & ! 3
     1, -4,   6, - 4,   1,   0,   0,   0,  0,  & ! 4
    -1,  5, -10,  10, - 5,   1,   0,   0,  0,  & ! 5
     1, -6,  15, -20,  15, - 6,   1,   0,  0,  & ! 6
    -1,  7, -21,  35, -35,  21, - 7,   1,  0,  & ! 7
     1, -8,  28, -56,  70, -56,  28, - 8,  1   & ! 8
  ], [9, 7])

  interface lon_damp_on_cell
    module procedure lon_damp_on_cell_2d
    module procedure lon_damp_on_cell_3d
  end interface lon_damp_on_cell

contains

  subroutine lon_damp_init()

    integer j
    real(r8) lat, lat0, lon_damp_lat0

    call lon_damp_final()

    allocate(wgt_full_lat(global_mesh%num_full_lat))
    allocate(wgt_half_lat(global_mesh%num_half_lat))

    lon_damp_lat0 = 80

    lat0 = -global_mesh%full_lat_deg(2)
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      lat = abs(global_mesh%full_lat_deg(j))
      wgt_full_lat(j) = exp((lat0 - lat)**2 * log(1.0e-3_r8) / (lat0 - lon_damp_lat0)**2)
      wgt_full_lat(j) = merge(0.0_r8, wgt_full_lat(j), wgt_full_lat(j) < 1.0e-6_r8)
    end do
    lat0 = -global_mesh%half_lat_deg(1)
    do j = global_mesh%half_lat_ibeg, global_mesh%half_lat_iend
      lat = abs(global_mesh%half_lat_deg(j))
      wgt_half_lat(j) = exp((lat0 - lat)**2 * log(1.0e-3_r8) / (lat0 - lon_damp_lat0)**2)
      wgt_half_lat(j) = merge(0.0_r8, wgt_half_lat(j), wgt_half_lat(j) < 1.0e-6_r8)
    end do

  end subroutine lon_damp_init

  subroutine lon_damp_final()

    if (allocated(wgt_full_lat)) deallocate(wgt_full_lat)
    if (allocated(wgt_half_lat)) deallocate(wgt_half_lat)

  end subroutine lon_damp_final

  subroutine lon_damp_on_lon_edge(block, order, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        call lon_damp_1d(block, order, mesh%half_lon_lb, mesh%half_lon_ub, mesh%lon_halo_width, f(:,j,k), wgt_full_lat(j))
      end do
    end do
    call fill_halo(block, f, full_lon=.false., full_lat=.true., full_lev=.true.)

  end subroutine lon_damp_on_lon_edge

  subroutine lon_damp_on_lat_edge(block, order, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        call lon_damp_1d(block, order, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j,k), wgt_half_lat(j))
      end do
    end do
    call fill_halo(block, f, full_lon=.true., full_lat=.false., full_lev=.true.)

  end subroutine lon_damp_on_lat_edge

  subroutine lon_damp_on_lev_edge(block, order, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%half_lev_lb:block%mesh%half_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        call lon_damp_1d(block, order, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j,k), wgt_full_lat(j))
      end do
    end do
    call fill_halo(block, f, full_lon=.true., full_lat=.true., full_lev=.false.)

  end subroutine lon_damp_on_lev_edge

  subroutine lon_damp_on_cell_2d(block, order, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      call lon_damp_1d(block, order, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j), wgt_full_lat(j))
    end do
    call fill_halo(block, f, full_lon=.true., full_lat=.true.)

  end subroutine lon_damp_on_cell_2d

  subroutine lon_damp_on_cell_3d(block, order, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        call lon_damp_1d(block, order, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j,k), wgt_full_lat(j))
      end do
    end do
    call fill_halo(block, f, full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine lon_damp_on_cell_3d

  subroutine lon_damp_on_vtx(block, order, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        call lon_damp_1d(block, order, mesh%half_lon_lb, mesh%half_lon_ub, mesh%lon_halo_width, f(:,j,k), wgt_half_lat(j))
      end do
    end do
    call fill_halo(block, f, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine lon_damp_on_vtx

  subroutine lon_damp_1d(block, order, lb, ub, hw, f, wgt)

    type(block_type), intent(in) :: block
    integer, intent(in) :: order
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    integer, intent(in) :: hw
    real(r8), intent(inout) :: f(lb:ub)
    real(r8), intent(in) :: wgt

    integer i, is, ie, ns
    real(r8) c, g(lb:ub)
    real(r8), pointer :: w(:)

    if (wgt == 0) return

    c = 0.5_r8**order
    is = lb + hw
    ie = ub - hw
    if (order > 2) then
      ns = diff_halo_width(order-1)
      w => diff_weights(:,order-1)
      ! Calculate damping flux at interfaces.
      do i = is, ie + 1
        g(i) = sum(f(i-ns:i+ns-1) * w(:2*ns))
      end do
      ! Limit damping flux to avoid upgradient (Xue 2000).
      g = g * (-1)**(order / 2)
      do i = is, ie + 1
        g(i) = g(i) * max(0.0_r8, sign(1.0_r8, -g(i) * (f(i) - f(i-1))))
      end do
      do i = is, ie
        f(i) = wgt * (f(i) - c * (g(i+1) - g(i))) + (1 - wgt) * f(i)
      end do
    else
      do i = is, ie
        g(i) = wgt * (f(i) + c * (f(i-1) - 2.0_r8 * f(i) + f(i+1))) + (1 - wgt) * f(i)
      end do
      f(is:ie) = g(is:ie)
    end if

  end subroutine lon_damp_1d

end module lon_damp_mod
