module zonal_damp_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public zonal_damp_init
  public zonal_damp_final
  public zonal_damp_on_lon_edge
  public zonal_damp_on_lat_edge
  public zonal_damp_on_lev_edge
  public zonal_damp_on_cell
  public zonal_damp_on_vtx
  public zonal_damp_1d

  logical , allocatable, target :: zonal_damp_on_full_lat(:)
  logical , allocatable, target :: zonal_damp_on_half_lat(:)

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

  interface zonal_damp_on_cell
    module procedure zonal_damp_on_cell_2d
    module procedure zonal_damp_on_cell_3d
  end interface zonal_damp_on_cell

contains

  subroutine zonal_damp_init()

    integer j

    call zonal_damp_final()

    allocate(zonal_damp_on_full_lat(global_mesh%num_full_lat)); zonal_damp_on_full_lat = .false.
    allocate(zonal_damp_on_half_lat(global_mesh%num_half_lat)); zonal_damp_on_half_lat = .false.

    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (abs(global_mesh%full_lat_deg(j)) >= polar_damp_lat0) then
        zonal_damp_on_full_lat(j) = .true.
      end if
    end do
    do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
#ifdef V_POLE
      zonal_damp_on_half_lat(j) = zonal_damp_on_full_lat(j) .or. zonal_damp_on_full_lat(j-1)
#else
      zonal_damp_on_half_lat(j) = zonal_damp_on_full_lat(j) .or. zonal_damp_on_full_lat(j+1)
#endif
    end do

  end subroutine zonal_damp_init

  subroutine zonal_damp_final()

    if (allocated(zonal_damp_on_full_lat)) deallocate(zonal_damp_on_full_lat)
    if (allocated(zonal_damp_on_half_lat)) deallocate(zonal_damp_on_half_lat)

  end subroutine zonal_damp_final

  subroutine zonal_damp_on_lon_edge(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (zonal_damp_on_full_lat(j)) then
          call zonal_damp_1d(block, order, dt, mesh%half_lon_lb, mesh%half_lon_ub, mesh%lon_halo_width, f(:,j,k))
        end if
      end do
    end do

  end subroutine zonal_damp_on_lon_edge

  subroutine zonal_damp_on_lat_edge(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        if (zonal_damp_on_half_lat(j)) then
          call zonal_damp_1d(block, order, dt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j,k))
        end if
      end do
    end do

  end subroutine zonal_damp_on_lat_edge

  subroutine zonal_damp_on_lev_edge(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%half_lev_lb:block%mesh%half_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (zonal_damp_on_full_lat(j)) then
          call zonal_damp_1d(block, order, dt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j,k))
        end if
      end do
    end do

  end subroutine zonal_damp_on_lev_edge

  subroutine zonal_damp_on_cell_2d(block, order, dt, f, lat0)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub)
    real(r8), intent(in), optional :: lat0

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    if (present(lat0)) then
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (abs(mesh%full_lat_deg(j)) >= lat0) then
          call zonal_damp_1d(block, order, dt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j))
        end if
      end do
    else
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (zonal_damp_on_full_lat(j)) then
          call zonal_damp_1d(block, order, dt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j))
        end if
      end do
    end if

  end subroutine zonal_damp_on_cell_2d

  subroutine zonal_damp_on_cell_3d(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        if (zonal_damp_on_full_lat(j)) then
          call zonal_damp_1d(block, order, dt, mesh%full_lon_lb, mesh%full_lon_ub, mesh%lon_halo_width, f(:,j,k))
        end if
      end do
    end do

  end subroutine zonal_damp_on_cell_3d

  subroutine zonal_damp_on_vtx(block, order, dt, f, lat0)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in), optional :: lat0

    type(mesh_type), pointer :: mesh
    integer j, k

    mesh => block%mesh
    if (present(lat0)) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          if (abs(mesh%half_lat_deg(j)) >= lat0) then
            call zonal_damp_1d(block, order, dt, mesh%half_lon_lb, mesh%half_lon_ub, mesh%lon_halo_width, f(:,j,k))
          end if
        end do
      end do
    else
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          if (zonal_damp_on_half_lat(j)) then
            call zonal_damp_1d(block, order, dt, mesh%half_lon_lb, mesh%half_lon_ub, mesh%lon_halo_width, f(:,j,k))
          end if
        end do
      end do
    end if

  end subroutine zonal_damp_on_vtx

  subroutine zonal_damp_1d(block, order, dt, lb, ub, hw, f)

    type(block_type), intent(in) :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    integer, intent(in) :: hw
    real(r8), intent(inout) :: f(lb:ub)

    integer i, is, ie, ns
    real(r8) c, g(lb:ub)
    real(r8), pointer :: w(:)

    c = 0.5_r8**order / dt
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
        f(i) = f(i) - c * dt * (g(i+1) - g(i))
      end do
    else
      do i = is, ie
        f(i) = f(i) + c * dt * (f(i-1) - 2.0_r8 * f(i) + f(i+1))
      end do
    end if
    call fill_zonal_halo(block, hw, f)

  end subroutine zonal_damp_1d

end module zonal_damp_mod
