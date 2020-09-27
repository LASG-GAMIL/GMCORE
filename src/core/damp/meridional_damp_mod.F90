module meridional_damp_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public meridional_damp_init
  public meridional_damp_final
  public meridional_damp_on_lon_edge
  public meridional_damp_on_lat_edge
  public meridional_damp_on_vtx

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

contains

  subroutine meridional_damp_init()

  end subroutine meridional_damp_init

  subroutine meridional_damp_final()

  end subroutine meridional_damp_final

  subroutine meridional_damp_on_lon_edge(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer i, io, j, k

    mesh => block%mesh

    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        io = i - global_mesh%num_half_lon / 2
        if (io > global_mesh%num_half_lon) io = io - global_mesh%num_half_lon
        f(i,j-1,:) = f(io,j+1,:)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        io = i - global_mesh%num_half_lon / 2
        if (io > global_mesh%num_half_lon) io = io - global_mesh%num_half_lon
        f(i,j+1,:) = f(io,j-1,:)
      end do
    end if

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        call meridional_damp_1d(block, order, dt, mesh%half_lat_lb, mesh%half_lat_ub, mesh%lat_halo_width, f(i,:,k))
      end do
    end do
    call fill_halo(block, f, full_lon=.true., full_lat=.false., full_lev=.true.)

  end subroutine meridional_damp_on_lon_edge

  subroutine meridional_damp_on_lat_edge(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer i, io, j, k

    mesh => block%mesh

    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        io = i - global_mesh%num_full_lon / 2
        if (io > global_mesh%num_full_lon) io = io - global_mesh%num_full_lon
        f(i,j-1,:) = f(io,j  ,:)
        f(i,j-2,:) = f(io,j+1,:)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        io = i - global_mesh%num_full_lon / 2
        if (io > global_mesh%num_full_lon) io = io - global_mesh%num_full_lon
        f(i,j+1,:) = f(io,j  ,:)
        f(i,j+2,:) = f(io,j-1,:)
      end do
    end if

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        call meridional_damp_1d(block, order, dt, mesh%half_lat_lb, mesh%half_lat_ub, mesh%lat_halo_width, f(i,:,k))
      end do
    end do
    call fill_halo(block, f, full_lon=.true., full_lat=.false., full_lev=.true.)

  end subroutine meridional_damp_on_lat_edge

  subroutine meridional_damp_on_vtx(block, order, dt, f)

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    integer i, io, j, k

    mesh => block%mesh

    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        io = i - global_mesh%num_half_lon / 2
        if (io > global_mesh%num_half_lon) io = io - global_mesh%num_half_lon
        f(i,j-1,:) = f(io,j  ,:)
        f(i,j-2,:) = f(io,j+1,:)
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        io = i - global_mesh%num_half_lon / 2
        if (io > global_mesh%num_half_lon) io = io - global_mesh%num_half_lon
        f(i,j+1,:) = f(io,j  ,:)
        f(i,j+2,:) = f(io,j-1,:)
      end do
    end if

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        call meridional_damp_1d(block, order, dt, mesh%half_lat_lb, mesh%half_lat_ub, mesh%lat_halo_width, f(i,:,k))
      end do
    end do
    call fill_halo(block, f, full_lon=.false., full_lat=.false., full_lev=.true.)

  end subroutine meridional_damp_on_vtx

  subroutine meridional_damp_1d(block, order, dt, lb, ub, hw, f)

    type(block_type), intent(in) :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    integer, intent(in) :: hw
    real(r8), intent(inout) :: f(lb:ub)

    integer j, js, je, ns
    real(r8) c, g(lb:ub)
    real(r8), pointer :: w(:)

    c = 0.5_r8**order / dt
    js = lb + hw
    je = ub - hw
    if (order > 2) then
      ns = diff_halo_width(order-1)
      w => diff_weights(:,order-1)
      ! Calculate damping flux at interfaces.
      do j = js, je + 1
        g(j) = sum(f(j-ns:j+ns-1) * w(:2*ns))
      end do
      ! Limit damping flux to avoid upgradient (Xue 2000).
      g = g * (-1)**(order / 2)
      do j = js, je + 1
        g(j) = g(j) * max(0.0_r8, sign(1.0_r8, -g(j) * (f(j) - f(j-1))))
      end do
      do j = js, je
        f(j) = f(j) - c * dt * (g(j+1) - g(j))
      end do
    else
      do j = js, je
        f(j) = f(j) + c * dt * (f(j-1) - 2.0_r8 * f(j) + f(j+1))
      end do
    end if

  end subroutine meridional_damp_1d

end module meridional_damp_mod
