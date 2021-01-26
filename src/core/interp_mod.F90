module interp_mod

  use const_mod
  use namelist_mod
  use block_mod
  use process_mod
  use parallel_mod

  implicit none

  private

  !                              / lev_lat_edge
  !               o-------------o------------o lev_vtx
  !              /|            /            /|
  !             / |                        / |
  !            /  |        |              /  |
  !           o   |        o lev_edge   -o- lev_lon_edge
  !          /    |        |            /    |
  !         /     o vtx                /     o vtx
  !        /      |                   /      |
  !       o-------+-----o------------o       |
  !       |       |                  |       |
  ! lon_edge -o-  |        o cell    |  -o- lon_edge
  !       |       |                  |       |
  !       |       o------------------+-------o
  !       |      /       /           |      /
  !       o vtx /       o lat_edge   o vtx /
  !       |    /       /             |    /
  !       |   o                      |   o
  !       |  /                       |  /
  !       | /                        | /
  !       |/                         |/
  !       o-------------o------------o
  !

  public interp_init
  public interp_final
  public interp_cell_to_lon_edge
  public interp_cell_to_lat_edge
  public interp_cell_to_lev_edge
  public interp_cell_to_vtx
  public interp_lon_edge_to_cell
  public interp_lat_edge_to_cell
  public interp_lon_edge_to_lev_lon_edge
  public interp_lat_edge_to_lev_lat_edge
  public interp_lev_edge_to_cell
  public interp_lev_edge_to_lev_lon_edge
  public interp_lev_edge_to_lev_lat_edge
  public interp_cell_to_isobaric_level
  public interp_lon_edge_to_isobaric_level
  public interp_lat_edge_to_isobaric_level

  real(r8), allocatable :: pole_wgt(:)

contains

  subroutine interp_init()

    integer j

    call interp_final()

    allocate(pole_wgt(global_mesh%num_full_lat))

    do j = 1, global_mesh%num_full_lat
      if (abs(global_mesh%full_lat_deg(j)) > 85) then
        pole_wgt(j) = 1 / upwind_wgt_pt
      else
        pole_wgt(j) = 1 + (1 / upwind_wgt_pt - 1) * &
          exp(-0.1 * (abs(global_mesh%full_lat_deg(j)) - 85)**2)
      end if
    end do

  end subroutine interp_init

  subroutine interp_final()

    if (allocated(pole_wgt)) deallocate(pole_wgt)

  end subroutine interp_final

  subroutine interp_cell_to_lon_edge(mesh, x, x_lon, reversed_area, u, upwind_wgt_, enhance_pole)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_lon(mesh%half_lon_lb:mesh%half_lon_ub, &
                                     mesh%full_lat_lb:mesh%full_lat_ub, &
                                     mesh%full_lev_lb:mesh%full_lev_ub)
    logical, intent(in), optional :: reversed_area
    real(r8), intent(in), optional :: u(mesh%half_lon_lb:mesh%half_lon_ub, &
                                        mesh%full_lat_lb:mesh%full_lat_ub, &
                                        mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in), optional :: upwind_wgt_
    logical, intent(in), optional :: enhance_pole

    real(r8), parameter :: c11 =  0.5_r8
    real(r8), parameter :: c12 = -0.5_r8
    real(r8), parameter :: c31 =  7.0_r8 / 12.0_r8
    real(r8), parameter :: c32 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c33 =  1.0_r8 / 12.0_r8
    real(r8) beta(mesh%full_lat_ibeg:mesh%full_lat_iend)
    integer i, j, k

    if (present(u)) then
      beta = merge(upwind_wgt_, upwind_wgt, present(upwind_wgt_))
      if (present(enhance_pole)) then
        if (enhance_pole) beta = beta * pole_wgt(mesh%full_lat_ibeg:mesh%full_lat_iend)
      end if
      ! Upwind-biased interpolation
      select case (upwind_order)
      case (1)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lon(i,j,k) = c11 * (x(i+1,j,k) + x(i,j,k)) + &
                             c12 * (x(i+1,j,k) - x(i,j,k)) * &
                         beta(j) * sign(1.0_r8, u(i,j,k))
            end do
          end do
        end do
        return
      case (3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lon(i,j,k) = c31 * (x(i+1,j,k) + x(i  ,j,k))  + &
                             c32 * (x(i+2,j,k) + x(i-1,j,k))  + &
                             c33 * (x(i+2,j,k) - x(i-1,j,k)   - &
                          3.0_r8 * (x(i+1,j,k) - x(i  ,j,k))) * &
                         beta(j) * sign(1.0_r8, u(i,j,k))
            end do
          end do
        end do
        return
      end select
    end if

    if (merge(reversed_area, .false., present(reversed_area))) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            x_lon(i,j,k) = (mesh%area_lon_east(j) * x(i  ,j,k) + &
                            mesh%area_lon_west(j) * x(i+1,j,k)   &
                           ) / mesh%area_lon(j)
          end do
        end do
      end do
    else ! reversed_area == .false.
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            x_lon(i,j,k) = (mesh%area_lon_west(j) * x(i  ,j,k) + &
                            mesh%area_lon_east(j) * x(i+1,j,k)   &
                           ) / mesh%area_lon(j)
          end do
        end do
      end do
    end if

  end subroutine interp_cell_to_lon_edge

  subroutine interp_cell_to_lat_edge(mesh, x, x_lat, reversed_area, handle_pole, v, upwind_wgt_, enhance_pole)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                     mesh%half_lat_lb:mesh%half_lat_ub, &
                                     mesh%full_lev_lb:mesh%full_lev_ub)
    logical, intent(in), optional :: reversed_area
    logical, intent(in), optional :: handle_pole
    real(r8), intent(in), optional :: v(mesh%full_lon_lb:mesh%full_lon_ub, &
                                        mesh%half_lat_lb:mesh%half_lat_ub, &
                                        mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in), optional :: upwind_wgt_
    logical, intent(in), optional :: enhance_pole

    real(r8), parameter :: c11 =  0.5_r8
    real(r8), parameter :: c12 = -0.5_r8
    real(r8), parameter :: c31 =  7.0_r8 / 12.0_r8
    real(r8), parameter :: c32 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c33 =  1.0_r8 / 12.0_r8
    real(r8) beta(mesh%full_lat_ibeg:mesh%full_lat_iend)
    integer i, j, k, jm1, jp2, io

    if (present(v)) then
      beta = merge(upwind_wgt_, upwind_wgt, present(upwind_wgt_))
      if (present(enhance_pole)) then
        if (enhance_pole) beta = beta * pole_wgt(mesh%full_lat_ibeg:mesh%full_lat_iend)
      end if
      ! Upwind-biased interpolation
      select case (upwind_order)
      case (1)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
#ifdef V_POLE
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              x_lat(i,j,k) = c11 * (x(i,j,k) + x(i,j-1,k)) + &
                             c12 * (x(i,j,k) - x(i,j-1,k)) * &
                         beta(j) * sign(1.0_r8, v(i,j,k))
            end do
#else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              x_lat(i,j,k) = c11 * (x(i,j+1,k) + x(i,j,k)) + &
                             c12 * (x(i,j+1,k) - x(i,j,k)) * &
                         beta(j) * sign(1.0_r8, v(i,j,k))
            end do
#endif
          end do
        end do
        return
      case(3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            if (mesh%is_half_lat_next_to_pole(j)) then
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
                x_lat(i,j,k) = c11 * (x(i,j,k) + x(i,j-1,k)) + &
                               c12 * (x(i,j,k) - x(i,j-1,k)) * &
                           beta(j) * sign(1.0_r8, v(i,j,k))
#else
                x_lat(i,j,k) = c11 * (x(i,j+1,k) + x(i,j,k)) + &
                               c12 * (x(i,j+1,k) - x(i,j,k)) * &
                           beta(j) * sign(1.0_r8, v(i,j,k))
#endif
              end do
            else
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
                x_lat(i,j,k) = c31 * (x(i,j  ,k) + x(i,j-1,k))  + &
                               c32 * (x(i,j+1,k) + x(i,j-2,k))  + &
                               c33 * (x(i,j+1,k) - x(i,j-2,k)   - &
                            3.0_r8 * (x(i,j  ,k) - x(i,j-1,k))) * &
                           beta(j) * sign(1.0_r8, v(i,j,k))
#else
                x_lat(i,j,k) = c31 * (x(i,j+1,k) + x(i,j  ,k))  + &
                               c32 * (x(i,j+2,k) + x(i,j-1,k))  + &
                               c33 * (x(i,j+2,k) - x(i,j-1,k)   - &
                            3.0_r8 * (x(i,j+1,k) - x(i,j  ,k))) * &
                           beta(j) * sign(1.0_r8, v(i,j,k))
              end do
#endif
            end if
          end do
        end do
        return
      end select
    end if

    if (merge(reversed_area, .false., present(reversed_area))) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            x_lat(i,j,k) = (mesh%area_lat_south(j) * x(i,j  ,k) + &
                            mesh%area_lat_north(j) * x(i,j-1,k)   &
                           ) / mesh%area_lat(j)
#else
            x_lat(i,j,k) = (mesh%area_lat_south(j) * x(i,j+1,k) + &
                            mesh%area_lat_north(j) * x(i,j  ,k)   &
                           ) / mesh%area_lat(j)
#endif
          end do
        end do
      end do
    else ! reversed_area == .false.
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            x_lat(i,j,k) = (mesh%area_lat_north(j) * x(i,j  ,k) + &
                            mesh%area_lat_south(j) * x(i,j-1,k)   &
                           ) / mesh%area_lat(j)
#else
            x_lat(i,j,k) = (mesh%area_lat_north(j) * x(i,j+1,k) + &
                            mesh%area_lat_south(j) * x(i,j  ,k)   &
                           ) / mesh%area_lat(j)
#endif
          end do
        end do
      end do
    end if
#ifndef V_POLE
    if (merge(handle_pole, .false., present(handle_pole))) then
      if (mesh%has_south_pole()) then
        j = mesh%half_lat_ibeg
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lat(i,j,:) = 2 * x(i,j+1,:) - x(i,j+2,:)
        end do
      end if
      if (mesh%has_north_pole()) then
        j = mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lat(i,j,:) = 2 * x(i,j-1,:) - x(i,j-2,:)
        end do
      end if
    end if
#endif

  end subroutine interp_cell_to_lat_edge

  subroutine interp_lev_edge_to_cell(mesh, x_lev, x)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lev(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(inout) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                                 mesh%full_lat_lb:mesh%full_lat_ub, &
                                 mesh%full_lev_lb:mesh%full_lev_ub)

    integer i, j, k
    real(r8) a, b

    ! =======
    !
    ! ---o--- k
    !
    ! ===?=== k
    !
    ! ---o--- k+1
    !
    ! =======
    ! NOTE: a and b should be 1/2.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      a = mesh%half_dlev_upper(k+1) / mesh%full_dlev(k)
      b = mesh%half_dlev_lower(k  ) / mesh%full_dlev(k)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x(i,j,k) = a * x_lev(i,j,k) + b * x_lev(i,j,k+1)
        end do
      end do
    end do

  end subroutine interp_lev_edge_to_cell

  subroutine interp_lev_edge_to_lev_lon_edge(mesh, x_lev, x_lev_lon, u, upwind_wgt_)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lev(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(inout) :: x_lev_lon(mesh%half_lon_lb:mesh%half_lon_ub, &
                                         mesh%full_lat_lb:mesh%full_lat_ub, &
                                         mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(in), optional :: u(mesh%half_lon_lb:mesh%half_lon_ub, &
                                        mesh%full_lat_lb:mesh%full_lat_ub, &
                                        mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(in), optional :: upwind_wgt_

    real(r8), parameter :: c11 =  0.5_r8
    real(r8), parameter :: c12 = -0.5_r8
    real(r8), parameter :: c31 =  7.0_r8 / 12.0_r8
    real(r8), parameter :: c32 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c33 =  1.0_r8 / 12.0_r8
    real(r8) beta
    integer i, j, k

    if (present(u)) then
      beta = merge(upwind_wgt_, upwind_wgt, present(upwind_wgt_))
      ! Upwind-biased interpolation
      select case (upwind_order)
      case (1)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lev_lon(i,j,k) = c11 * (x_lev(i+1,j,k) + x_lev(i,j,k)) + &
                                 c12 * (x_lev(i+1,j,k) - x_lev(i,j,k)) * &
                                beta * sign(1.0_r8, u(i,j,k))
            end do
          end do
        end do
        return
      case (3)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lev_lon(i,j,k) = c31 * (x_lev(i+1,j,k) + x_lev(i  ,j,k))  + &
                                 c32 * (x_lev(i+2,j,k) + x_lev(i-1,j,k))  + &
                                 c33 * (x_lev(i+2,j,k) - x_lev(i-1,j,k)   - &
                              3.0_r8 * (x_lev(i+1,j,k) - x_lev(i  ,j,k))) * &
                                beta * sign(1.0_r8, u(i,j,k))
            end do
          end do
        end do
        return
      end select
    end if
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_lev_lon(i,j,k) = (mesh%area_lon_west(j) * x_lev(i  ,j,k) + &
                              mesh%area_lon_east(j) * x_lev(i+1,j,k)   &
                             ) / mesh%area_lon(j)
        end do
      end do
    end do

  end subroutine interp_lev_edge_to_lev_lon_edge

  subroutine interp_lev_edge_to_lev_lat_edge(mesh, x_lev, x_lev_lat, v, upwind_wgt_)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lev(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(inout) :: x_lev_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                         mesh%half_lat_lb:mesh%half_lat_ub, &
                                         mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(in), optional :: v(mesh%full_lon_lb:mesh%full_lon_ub, &
                                        mesh%half_lat_lb:mesh%half_lat_ub, &
                                        mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(in), optional :: upwind_wgt_

    real(r8), parameter :: c11 =  0.5_r8
    real(r8), parameter :: c12 = -0.5_r8
    real(r8), parameter :: c31 =  7.0_r8 / 12.0_r8
    real(r8), parameter :: c32 = -1.0_r8 / 12.0_r8
    real(r8), parameter :: c33 =  1.0_r8 / 12.0_r8
    real(r8) beta
    integer i, j, k

    if (present(v)) then
      beta = merge(upwind_wgt_, upwind_wgt, present(upwind_wgt_))
      ! Upwind-biased interpolation
      select case (upwind_order)
      case (1)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              x_lev_lat(i,j,k) = c11 * (x_lev(i,j,k) + x_lev(i,j-1,k)) + &
                                 c12 * (x_lev(i,j,k) - x_lev(i,j-1,k)) * &
                                beta * sign(1.0_r8, v(i,j,k))
#else
              x_lev_lat(i,j,k) = c11 * (x_lev(i,j+1,k) + x_lev(i,j,k)) + &
                                 c12 * (x_lev(i,j+1,k) - x_lev(i,j,k)) * &
                                beta * sign(1.0_r8, v(i,j,k))
#endif
            end do
          end do
        end do
        return
      case (3)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            if (mesh%is_half_lat_next_to_pole(j)) then
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
                x_lev_lat(i,j,k) = c11 * (x_lev(i,j,k) + x_lev(i,j-1,k)) + &
                                   c12 * (x_lev(i,j,k) - x_lev(i,j-1,k)) * &
                                  beta * sign(1.0_r8, v(i,j,k))
#else
                x_lev_lat(i,j,k) = c11 * (x_lev(i,j+1,k) + x_lev(i,j,k)) + &
                                   c12 * (x_lev(i,j+1,k) - x_lev(i,j,k)) * &
                                  beta * sign(1.0_r8, v(i,j,k))
#endif
              end do
            else
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
                x_lev_lat(i,j,k) = c31 * (x_lev(i,j  ,k) + x_lev(i,j-1,k))  + &
                                   c32 * (x_lev(i,j+1,k) + x_lev(i,j-2,k))  + &
                                   c33 * (x_lev(i,j+1,k) - x_lev(i,j-2,k)   - &
                                3.0_r8 * (x_lev(i,j  ,k) - x_lev(i,j-1,k))) * &
                                  beta * sign(1.0_r8, v(i,j,k))
#else
                x_lev_lat(i,j,k) = c31 * (x_lev(i,j+1,k) + x_lev(i,j  ,k))  + &
                                   c32 * (x_lev(i,j+2,k) + x_lev(i,j-1,k))  + &
                                   c33 * (x_lev(i,j+2,k) - x_lev(i,j-1,k)   - &
                                3.0_r8 * (x_lev(i,j+1,k) - x_lev(i,j  ,k))) * &
                                  beta * sign(1.0_r8, v(i,j,k))
#endif
              end do
            end if
          end do
        end do
        return
      end select
    end if
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          x_lev_lat(i,j,k) = (mesh%area_lat_north(j) * x_lev(i,j  ,k) + &
                              mesh%area_lat_south(j) * x_lev(i,j-1,k)   &
                             ) / mesh%area_lat(j)
#else
          x_lev_lat(i,j,k) = (mesh%area_lat_north(j) * x_lev(i,j+1,k) + &
                              mesh%area_lat_south(j) * x_lev(i,j  ,k)   &
                             ) / mesh%area_lat(j)
#endif
        end do
      end do
    end do

  end subroutine interp_lev_edge_to_lev_lat_edge

  subroutine interp_cell_to_vtx(mesh, x, x_vtx)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_vtx(mesh%half_lon_lb:mesh%half_lon_ub, &
                                     mesh%half_lat_lb:mesh%half_lat_ub, &
                                     mesh%full_lev_lb:mesh%full_lev_ub)

    integer i, j, k
    real(r8) work(mesh%full_lon_ibeg:mesh%full_lon_iend,mesh%num_full_lev)
    real(r8) pole(mesh%num_full_lev)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          x_vtx(i,j,k) = (                                           &
            (x(i,j-1,k) + x(i+1,j-1,k)) * mesh%area_subcell(2,j-1) + &
            (x(i,j  ,k) + x(i+1,j  ,k)) * mesh%area_subcell(1,j  )   &
          ) / mesh%area_vtx(j)
#else
          x_vtx(i,j,k) = (                                           &
            (x(i,j  ,k) + x(i+1,j  ,k)) * mesh%area_subcell(2,j  ) + &
            (x(i,j+1,k) + x(i+1,j+1,k)) * mesh%area_subcell(1,j+1)   &
          ) / mesh%area_vtx(j)
#endif
        end do
      end do
    end do
#ifdef V_POLE
    if (mesh%has_south_pole()) then
      j = mesh%half_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = x(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_comm, work, pole)
      pole = pole / global_mesh%num_half_lon
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_vtx(i,j,k) = pole(k)
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%half_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = x(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_comm, work, pole)
      pole = pole / global_mesh%num_half_lon
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_vtx(i,j,k) = pole(k)
        end do
      end do
    end if
#endif

  end subroutine interp_cell_to_vtx

  subroutine interp_cell_to_lev_edge(mesh, x, x_lev, handle_top_bottom)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_lev(mesh%full_lon_lb:mesh%full_lon_ub, &
                                     mesh%full_lat_lb:mesh%full_lat_ub, &
                                     mesh%half_lev_lb:mesh%half_lev_ub)
    logical, intent(in), optional :: handle_top_bottom

    integer i, j, k
    real(r8) x1, x2, a, b

    ! -------
    !
    ! ===o=== k-1
    !
    ! ---?--- k
    !
    ! ===o=== k
    !
    ! -------
    do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
      a = mesh%full_dlev(k  ) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      b = mesh%full_dlev(k-1) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev(i,j,k) = a * x(i,j,k-1) + b * x(i,j,k)
        end do
      end do
    end do

    if (merge(handle_top_bottom, .false., present(handle_top_bottom))) then
      k = mesh%half_lev_ibeg
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      x1 = mesh%full_lev(k  ) - mesh%half_lev(k)
      x2 = mesh%full_lev(k+1) - mesh%half_lev(k)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev(i,j,k) = a * x(i,j,k) + b * x(i,j,k+1)
        end do
      end do
      k = mesh%half_lev_iend
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = mesh%half_lev(k) - mesh%full_lev(k-1)
      x2 = mesh%half_lev(k) - mesh%full_lev(k-2)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev(i,j,k) = a * x(i,j,k-1) + b * x(i,j,k-2)
        end do
      end do
    end if

  end subroutine interp_cell_to_lev_edge

  subroutine interp_cell_to_isobaric_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: po
    real(r8), intent(inout) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                                 mesh%full_lat_lb:mesh%full_lat_ub)
    logical, intent(in), optional :: logp

    logical logp_
    real(r8) p0, dp1, dp2
    integer i, j, k

    logp_ = merge(logp, .false., present(logp))

    p0 = merge(log(po), po, logp_)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%full_lev_iend, mesh%full_lev_ibeg + 1, -1
          if (p(i,j,k) >= po .and. p(i,j,k-1) <= po) then
            if (logp_) then
              dp1 = p0 - log(p(i,j,k-1))
              dp2 = log(p(i,j,k)) - p0
            else
              dp1 = p0 - p(i,j,k-1)
              dp2 = p(i,j,k) - p0
            end if
            y(i,j) = (dp2 * x(i,j,k-1) + dp1 * x(i,j,k)) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_cell_to_isobaric_level

  subroutine interp_lon_edge_to_isobaric_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: x(mesh%half_lon_lb:mesh%half_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: po
    real(r8), intent(inout) :: y(mesh%half_lon_lb:mesh%half_lon_ub, &
                                 mesh%full_lat_lb:mesh%full_lat_ub)
    logical, intent(in), optional :: logp

    logical logp_
    real(r8) p0, p1_lon, p2_lon, dp1, dp2
    integer i, j, k

    logp_ = merge(logp, .false., present(logp))

    p0 = merge(log(po), po, logp_)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        do k = mesh%full_lev_iend, mesh%full_lev_ibeg + 1, -1
          p1_lon = 0.5_r8 * (p(i,j,k-1) + p(i+1,j,k-1))
          p2_lon = 0.5_r8 * (p(i,j,k  ) + p(i+1,j,k  ))
          if (p2_lon >= po .and. p1_lon <= po) then
            if (logp_) then
              dp1 = p0 - log(p1_lon)
              dp2 = log(p2_lon) - p0
            else
              dp1 = p0 - p1_lon
              dp2 = p2_lon - p0
            end if
            y(i,j) = (dp2 * x(i,j,k-1) + dp1 * x(i,j,k)) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lon_edge_to_isobaric_level

  subroutine interp_lat_edge_to_isobaric_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%half_lat_lb:mesh%half_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: po
    real(r8), intent(inout) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                                 mesh%half_lat_lb:mesh%half_lat_ub)
    logical, intent(in), optional :: logp

    logical logp_
    real(r8) p0, p1_lat, p2_lat, dp1, dp2
    integer i, j, k

    logp_ = merge(logp, .false., present(logp))

    p0 = merge(log(po), po, logp_)

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%full_lev_iend, mesh%full_lev_ibeg + 1, -1
#ifdef V_POLE
          p1_lat = 0.5_r8 * (p(i,j-1,k-1) + p(i,j,k-1))
          p2_lat = 0.5_r8 * (p(i,j-1,k  ) + p(i,j,k  ))
#else
          p1_lat = 0.5_r8 * (p(i,j,k-1) + p(i,j+1,k-1))
          p2_lat = 0.5_r8 * (p(i,j,k  ) + p(i,j+1,k  ))
#endif
          if (p2_lat >= po .and. p1_lat <= po) then
            if (logp_) then
              dp1 = p0 - log(p1_lat)
              dp2 = log(p2_lat) - p0
            else
              dp1 = p0 - p1_lat
              dp2 = p2_lat - p0
            end if
            y(i,j) = (dp2 * x(i,j,k-1) + dp1 * x(i,j,k)) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lat_edge_to_isobaric_level

  subroutine interp_lon_edge_to_cell(mesh, x_lon, x)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lon(mesh%half_lon_lb:mesh%half_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(out) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                               mesh%full_lat_lb:mesh%full_lat_ub, &
                               mesh%full_lev_lb:mesh%full_lev_ub)

    integer i, j, k

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x(i,j,k) = (mesh%area_lon_east(j) * x_lon(i-1,j,k) + &
                      mesh%area_lon_west(j) * x_lon(i  ,j,k)   &
                     ) / mesh%area_lon(j)
        end do
      end do
    end do

  end subroutine interp_lon_edge_to_cell

  subroutine interp_lat_edge_to_cell(mesh, x_lat, x)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%half_lat_lb:mesh%half_lat_ub, &
                                  mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(out) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                               mesh%full_lat_lb:mesh%full_lat_ub, &
                               mesh%full_lev_lb:mesh%full_lev_ub)

    integer i, j, k

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          x(i,j,k) = (mesh%area_lat_north(j  ) * x_lat(i,j  ,k) + &
                      mesh%area_lat_south(j+1) * x_lat(i,j+1,k)   &
                     ) / (mesh%area_lat_north(j) + mesh%area_lat_south(j+1))
#else
          x(i,j,k) = (mesh%area_lat_south(j  ) * x_lat(i,j  ,k) + &
                      mesh%area_lat_north(j-1) * x_lat(i,j-1,k)   &
                     ) / (mesh%area_lat_south(j) + mesh%area_lat_north(j-1))
#endif
        end do
      end do
    end do

  end subroutine interp_lat_edge_to_cell

  subroutine interp_lon_edge_to_lev_lon_edge(mesh, x_lon, x_lev_lon, handle_top_bottom)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lon(mesh%half_lon_lb:mesh%half_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(out) :: x_lev_lon(mesh%half_lon_lb:mesh%half_lon_ub, &
                                       mesh%full_lat_lb:mesh%full_lat_ub, &
                                       mesh%half_lev_lb:mesh%half_lev_ub)
    logical, intent(in), optional :: handle_top_bottom

    integer i, j, k
    real(r8) x1, x2, a, b

    ! -------
    !
    ! ===o=== k-1
    !
    ! ---?--- k
    !
    ! ===o=== k
    !
    ! ----o--
    do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
      a = mesh%full_dlev(k-1) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      b = mesh%full_dlev(k  ) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_lev_lon(i,j,k) = a * x_lon(i,j,k) + b * x_lon(i,j,k-1)
        end do
      end do
    end do

    if (merge(handle_top_bottom, .false., present(handle_top_bottom))) then
      k = mesh%half_lev_ibeg
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      x1 = mesh%full_lev(k  ) - mesh%half_lev(k)
      x2 = mesh%full_lev(k+1) - mesh%half_lev(k)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev_lon(i,j,k) = a * x_lon(i,j,k) + b * x_lon(i,j,k+1)
        end do
      end do
      k = mesh%half_lev_iend
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = mesh%half_lev(k) - mesh%full_lev(k-1)
      x2 = mesh%half_lev(k) - mesh%full_lev(k-2)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev_lon(i,j,k) = a * x_lon(i,j,k-1) + b * x_lon(i,j,k-2)
        end do
      end do
    end if

  end subroutine interp_lon_edge_to_lev_lon_edge

  subroutine interp_lat_edge_to_lev_lat_edge(mesh, x_lat, x_lev_lat, handle_top_bottom)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%half_lat_lb:mesh%half_lat_ub, &
                                  mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(out) :: x_lev_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                       mesh%half_lat_lb:mesh%half_lat_ub, &
                                       mesh%half_lev_lb:mesh%half_lev_ub)
    logical, intent(in), optional :: handle_top_bottom

    integer i, j, k
    real(r8) x1, x2, a, b

    ! -------
    !
    ! ===o=== k-1
    !
    ! ---?--- k
    !
    ! ===o=== k
    !
    ! -------
    do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
      a = mesh%full_dlev(k-1) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      b = mesh%full_dlev(k  ) / (mesh%full_dlev(k-1) + mesh%full_dlev(k))
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev_lat(i,j,k) = a * x_lat(i,j,k) + b * x_lat(i,j,k-1)
        end do
      end do
    end do

    if (merge(handle_top_bottom, .false., present(handle_top_bottom))) then
      k = mesh%half_lev_ibeg
      ! ---?--- 1
      !
      ! ===o=== 1   x1
      !
      ! -------
      !
      ! ===o=== 2   x2
      x1 = mesh%full_lev(k  ) - mesh%half_lev(k)
      x2 = mesh%full_lev(k+1) - mesh%half_lev(k)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev_lat(i,j,k) = a * x_lat(i,j,k) + b * x_lat(i,j,k+1)
        end do
      end do
      k = mesh%half_lev_iend
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = mesh%half_lev(k) - mesh%full_lev(k-1)
      x2 = mesh%half_lev(k) - mesh%full_lev(k-2)
      a =  x2 / (x2 - x1)
      b = -x1 / (x2 - x1)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev_lat(i,j,k) = a * x_lat(i,j,k-1) + b * x_lat(i,j,k-2)
        end do
      end do
    end if

  end subroutine interp_lat_edge_to_lev_lat_edge

end module interp_mod
