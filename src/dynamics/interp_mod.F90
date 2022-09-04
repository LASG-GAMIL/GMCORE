module interp_mod

  use const_mod
  use namelist_mod
  use block_mod
  use process_mod
  use parallel_mod
  use upwind_mod
  use weno_mod

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
  public average_cell_to_lon_edge
  public average_cell_to_lat_edge
  public interp_cell_to_vtx
  public interp_lon_edge_to_cell
  public interp_lat_edge_to_cell
  public interp_lon_edge_to_lev_lon_edge
  public interp_lat_edge_to_lev_lat_edge
  public interp_lev_edge_to_cell
  public interp_lev_edge_to_lev_lon_edge
  public interp_lev_edge_to_lev_lat_edge
  public interp_cell_to_height_level
  public interp_lon_edge_to_height_level
  public interp_lat_edge_to_height_level
  public interp_lev_edge_to_height_level
  public interp_cell_to_pressure_level
  public interp_lon_edge_to_pressure_level
  public interp_lat_edge_to_pressure_level
  public interp_lev_edge_to_pressure_level

contains

  subroutine interp_init()
  
    call interp_final()

  end subroutine interp_init

  subroutine interp_final()

  end subroutine interp_final

  subroutine interp_cell_to_lon_edge(mesh, x, x_lon, reversed_area, u, upwind_wgt_)

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

    real(r8) beta
    integer i, j, k

    if (present(u)) then
      ! WENO interpolation
      select case (weno_order)
      case (3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lon(i,j,k) = weno3(sign(1.0_r8, u(i,j,k)), x(i-1:i+2,j,k))
            end do
          end do
        end do
        return
      end select
      ! Upwind-biased interpolation
      beta = merge(upwind_wgt_, upwind_wgt, present(upwind_wgt_))
      select case (upwind_order)
      case (1)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lon(i,j,k) = upwind1(sign(1.0_r8, u(i,j,k)), beta, x(i:i+1,j,k))
            end do
          end do
        end do
        return
      case (3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lon(i,j,k) = upwind3(sign(1.0_r8, u(i,j,k)), beta, x(i-1:i+2,j,k))
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

  subroutine interp_cell_to_lat_edge(mesh, x, x_lat, reversed_area, v, upwind_wgt_)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                     mesh%half_lat_lb:mesh%half_lat_ub, &
                                     mesh%full_lev_lb:mesh%full_lev_ub)
    logical, intent(in), optional :: reversed_area
    real(r8), intent(in), optional :: v(mesh%full_lon_lb:mesh%full_lon_ub, &
                                        mesh%half_lat_lb:mesh%half_lat_ub, &
                                        mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in), optional :: upwind_wgt_

    real(r8) beta
    integer i, j, k

    if (present(v)) then
      beta = merge(upwind_wgt_, upwind_wgt, present(upwind_wgt_))
      ! WENO interpolation
      select case (weno_order)
      case (3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              x_lat(i,j,k) = weno3(sign(1.0_r8, v(i,j,k)), x(i,j-1:j+2,k))
            end do
          end do
        end do
        return
      end select
      ! Upwind-biased interpolation
      select case (upwind_order)
      case (1)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              x_lat(i,j,k) = upwind1(sign(1.0_r8, v(i,j,k)), beta, x(i,j:j+1,k))
            end do
          end do
        end do
        return
      case(3)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              x_lat(i,j,k) = upwind3(sign(1.0_r8, v(i,j,k)), beta, x(i,j-1:j+2,k))
            end do
          end do
        end do
        return
      end select
    end if

    if (merge(reversed_area, .false., present(reversed_area))) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            x_lat(i,j,k) = (mesh%area_lat_south(j) * x(i,j+1,k) + &
                            mesh%area_lat_north(j) * x(i,j  ,k)   &
                           ) / mesh%area_lat(j)
          end do
        end do
      end do
    else ! reversed_area == .false.
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            x_lat(i,j,k) = (mesh%area_lat_north(j) * x(i,j+1,k) + &
                            mesh%area_lat_south(j) * x(i,j  ,k)   &
                           ) / mesh%area_lat(j)
          end do
        end do
      end do
    end if

  end subroutine interp_cell_to_lat_edge

  subroutine average_cell_to_lon_edge(mesh, x, x_lon)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in ) :: x    (mesh%full_lon_lb:mesh%full_lon_ub, &
                                   mesh%full_lat_lb:mesh%full_lat_ub, &
                                   mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(out) :: x_lon(mesh%half_lon_lb:mesh%half_lon_ub, &
                                   mesh%full_lat_lb:mesh%full_lat_ub, &
                                   mesh%full_lev_lb:mesh%full_lev_ub)

    integer i, j, k

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_lon(i,j,k) = (x(i,j,k) + x(i+1,j,k)) * 0.5_r8
        end do
      end do
    end do

  end subroutine average_cell_to_lon_edge

  subroutine average_cell_to_lat_edge(mesh, x, x_lat)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in ) :: x    (mesh%full_lon_lb:mesh%full_lon_ub, &
                                   mesh%full_lat_lb:mesh%full_lat_ub, &
                                   mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(out) :: x_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                   mesh%half_lat_lb:mesh%half_lat_ub, &
                                   mesh%full_lev_lb:mesh%full_lev_ub)

    integer i, j, k

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lat(i,j,k) = (x(i,j,k) + x(i,j+1,k)) * 0.5_r8
        end do
      end do
    end do

  end subroutine average_cell_to_lat_edge

  subroutine interp_lev_edge_to_cell(mesh, x_lev, x)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x_lev(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(inout) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                                 mesh%full_lat_lb:mesh%full_lat_ub, &
                                 mesh%full_lev_lb:mesh%full_lev_ub)

    integer i, j, k

    ! =======
    !
    ! ---o--- k
    !
    ! ===?=== k
    !
    ! ---o--- k+1
    !
    ! =======
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x(i,j,k) = 0.5_r8 * (x_lev(i,j,k) + x_lev(i,j,k+1))
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

    real(r8) beta
    integer i, j, k

    if (present(u)) then
      ! WENO interpolation
      select case (weno_order)
      case (3)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lev_lon(i,j,k) = weno3(sign(1.0_r8, u(i,j,k)), x_lev(i-1:i+2,j,k))
            end do
          end do
        end do
        return
      end select
      ! Upwind-biased interpolation
      beta = merge(upwind_wgt_, upwind_wgt, present(upwind_wgt_))
      select case (upwind_order)
      case (1)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lev_lon(i,j,k) = upwind1(sign(1.0_r8, u(i,j,k)), beta, x_lev(i:i+1,j,k))
            end do
          end do
        end do
        return
      case (3)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              x_lev_lon(i,j,k) = upwind3(sign(1.0_r8, u(i,j,k)), beta, x_lev(i-1:i+2,j,k))
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

    real(r8) beta
    integer i, j, k

    if (present(v)) then
      beta = merge(upwind_wgt_, upwind_wgt, present(upwind_wgt_))
      ! WENO interpolation
      select case (weno_order)
      case (3)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              x_lev_lat(i,j,k) = weno3(sign(1.0_r8, v(i,j,k)), x_lev(i,j-1:j+2,k))
            end do
          end do
        end do
        return
      end select
      ! Upwind-biased interpolation
      select case (upwind_order)
      case (1)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              x_lev_lat(i,j,k) = upwind1(sign(1.0_r8, v(i,j,k)), beta, x_lev(i,j:j+1,k))
            end do
          end do
        end do
        return
      case (3)
        do k = mesh%half_lev_ibeg, mesh%half_lev_iend
          do j = mesh%half_lat_ibeg, mesh%half_lat_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              x_lev_lat(i,j,k) = upwind3(sign(1.0_r8, v(i,j,k)), beta, x_lev(i,j-1:j+2,k))
            end do
          end do
        end do
        return
      end select
    end if
    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev_lat(i,j,k) = (mesh%area_lat_north(j) * x_lev(i,j+1,k) + &
                              mesh%area_lat_south(j) * x_lev(i,j  ,k)   &
                             ) / mesh%area_lat(j)
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
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_vtx(i,j,k) = (                                           &
            (x(i,j  ,k) + x(i+1,j  ,k)) * mesh%area_subcell(2,j  ) + &
            (x(i,j+1,k) + x(i+1,j+1,k)) * mesh%area_subcell(1,j+1)   &
          ) / mesh%area_vtx(j)
        end do
      end do
    end do

  end subroutine interp_cell_to_vtx

  subroutine interp_cell_to_lev_edge(mesh, x, x_lev)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in   ) :: x    (mesh%full_lon_lb:mesh%full_lon_ub, &
                                     mesh%full_lat_lb:mesh%full_lat_ub, &
                                     mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_lev(mesh%full_lon_lb:mesh%full_lon_ub, &
                                     mesh%full_lat_lb:mesh%full_lat_ub, &
                                     mesh%half_lev_lb:mesh%half_lev_ub)

    integer i, j, k
    real(r8) x1, x2, a, b

    if (mesh%num_full_lev == 1) return

    ! --------------------------------------------------------------------------
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
      a = mesh%full_dlev(k-1) / (2 * mesh%half_dlev(k))
      b = mesh%full_dlev(k  ) / (2 * mesh%half_dlev(k))
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev(i,j,k) = a * x(i,j,k-1) + b * x(i,j,k)
        end do
      end do
    end do

    k = mesh%half_lev_ibeg
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        x_lev(i,j,k) = x(i,j,k)
      end do
    end do
    k = mesh%half_lev_iend
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        x_lev(i,j,k) = x(i,j,k-1)
      end do
    end do

  end subroutine interp_cell_to_lev_edge

  subroutine interp_cell_to_height_level(mesh, z, x, zo, y)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: z(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: zo
    real(r8), intent(inout) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                                 mesh%full_lat_lb:mesh%full_lat_ub)

    real(r8) dz1, dz2, z1, z2, a, b
    integer i, j, k

    ! --o-- z(k-1)
    !
    ! --?-- zo
    !
    ! --o-- z(k)
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%full_lev_iend, mesh%full_lev_ibeg + 1, -1
          if (z(i,j,k) <= zo .and. zo <= z(i,j,k-1)) then
            dz1 = z(i,j,k-1) - zo
            dz2 = zo - z(i,j,k)
            y(i,j) = (dz2 * x(i,j,k-1) + dz1 * x(i,j,k)) / (dz1 + dz2)
            exit
          else if (zo < z(i,j,k) .and. k == mesh%full_lev_iend) then
            z1 = z(i,j,k  ) - zo
            z2 = z(i,j,k-1) - zo
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            y(i,j) = a * x(i,j,k) + b * x(i,j,k-1)
            exit
          else if (zo > z(i,j,k-1) .and. k == mesh%full_lev_ibeg + 1) then
            z1 = zo - z(i,j,k-1)
            z2 = zo - z(i,j,k  )
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            y(i,j) = a * x(i,j,k-1) + b * x(i,j,k)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_cell_to_height_level

  subroutine interp_lon_edge_to_height_level(mesh, z, x, zo, y)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: z(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: x(mesh%half_lon_lb:mesh%half_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: zo
    real(r8), intent(out) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                               mesh%full_lat_lb:mesh%full_lat_ub)

    real(r8) dz1, dz2, x1, x2, z1, z2, a, b
    integer i, j, k

    !               x1
    ! x(i-1,k-1) o--x--o x(i,k-1)
    !
    !               ? zo
    !
    ! x(i-1,k  ) o--x--o x(i,k  )
    !               x2
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%full_lev_iend, mesh%full_lev_ibeg + 1, -1
          if (z(i,j,k) <= zo .and. zo <= z(i,j,k-1)) then
            dz1 = z(i,j,k-1) - zo
            dz2 = zo - z(i,j,k)
            x1 = 0.5_r8 * (x(i-1,j,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i-1,j,k  ) + x(i,j,k  ))
            y(i,j) = (dz2 * x1 + dz1 * x2) / (dz1 + dz2)
            exit
          else if (zo < z(i,j,k) .and. k == mesh%full_lev_iend) then
            z1 = z(i,j,k  ) - zo
            z2 = z(i,j,k-1) - zo
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            x1 = 0.5_r8 * (x(i-1,j,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i-1,j,k  ) + x(i,j,k  ))
            y(i,j) = a * x1 + b * x2
            exit
          else if (zo > z(i,j,k-1) .and. k == mesh%full_lev_ibeg + 1) then
            z1 = zo - z(i,j,k-1)
            z2 = zo - z(i,j,k  )
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            x1 = 0.5_r8 * (x(i-1,j,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i-1,j,k  ) + x(i,j,k  ))
            y(i,j) = a * x1 + b * x2
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lon_edge_to_height_level

  subroutine interp_lat_edge_to_height_level(mesh, z, x, zo, y)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: z(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%half_lat_lb:mesh%half_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: zo
    real(r8), intent(out) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                               mesh%full_lat_lb:mesh%full_lat_ub)

    real(r8) dz1, dz2, x1, x2, z1, z2, a, b
    integer i, j, k

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%full_lev_iend, mesh%full_lev_ibeg + 1, -1
          if (z(i,j,k) <= zo .and. zo <= z(i,j,k-1)) then
            dz1 = z(i,j,k-1) - zo
            dz2 = zo - z(i,j,k)
            x1 = 0.5_r8 * (x(i,j-1,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i,j-1,k  ) + x(i,j,k  ))
            y(i,j) = (dz2 * x1 + dz1 * x2) / (dz1 + dz2)
            exit
          else if (zo < z(i,j,k) .and. k == mesh%full_lev_iend) then
            z1 = z(i,j,k  ) - zo
            z2 = z(i,j,k-1) - zo
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            x1 = 0.5_r8 * (x(i,j-1,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i,j-1,k  ) + x(i,j,k  ))
            y(i,j) = a * x1 + b * x2
            exit
          else if (zo > z(i,j,k-1) .and. k == mesh%full_lev_ibeg + 1) then
            z1 = zo - z(i,j,k-1)
            z2 = zo - z(i,j,k  )
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            x1 = 0.5_r8 * (x(i,j-1,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i,j-1,k  ) + x(i,j,k  ))
            y(i,j) = a * x1 + b * x2
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lat_edge_to_height_level

  subroutine interp_lev_edge_to_height_level(mesh, z, x, zo, y)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: z(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(in) :: zo
    real(r8), intent(out) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                               mesh%full_lat_lb:mesh%full_lat_ub)

    real(r8) dz1, dz2, x1, x2, z1, z2, a, b
    integer i, j, k

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%half_lev_iend, mesh%half_lev_ibeg + 1, -1
          if (z(i,j,k) <= zo .and. zo <= z(i,j,k-1)) then
            dz1 = z(i,j,k-1) - zo
            dz2 = zo - z(i,j,k)
            y(i,j) = (dz2 * x(i,j,k-1) + dz1 * x(i,j,k)) / (dz1 + dz2)
            exit
          else if (zo < z(i,j,k) .and. k == mesh%full_lev_iend) then
            z1 = z(i,j,k  ) - zo
            z2 = z(i,j,k-1) - zo
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            y(i,j) = a * x(i,j,k) + b * x(i,j,k-1)
            exit
          else if (zo > z(i,j,k-1) .and. k == mesh%full_lev_ibeg + 1) then
            z1 = zo - z(i,j,k-1)
            z2 = zo - z(i,j,k  )
            a  =  z2 / (z2 - z1)
            b  = -z1 / (z2 - z1)
            y(i,j) = a * x(i,j,k-1) + b * x(i,j,k)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lev_edge_to_height_level

  subroutine interp_cell_to_pressure_level(mesh, p, x, po, y, logp)

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
          if (p(i,j,k-1) <= po .and. po <= p(i,j,k)) then
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

  end subroutine interp_cell_to_pressure_level

  subroutine interp_lon_edge_to_pressure_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: x(mesh%half_lon_lb:mesh%half_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: po
    real(r8), intent(out) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                               mesh%full_lat_lb:mesh%full_lat_ub)
    logical, intent(in), optional :: logp

    logical logp_
    real(r8) p0, dp1, dp2, x1, x2
    integer i, j, k

    logp_ = merge(logp, .false., present(logp))

    p0 = merge(log(po), po, logp_)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%full_lev_iend, mesh%full_lev_ibeg + 1, -1
          if (p(i,j,k-1) <= po .and. po <= p(i,j,k)) then
            if (logp_) then
              dp1 = p0 - log(p(i,j,k-1))
              dp2 = log(p(i,j,k)) - p0
            else
              dp1 = p0 - p(i,j,k-1)
              dp2 = p(i,j,k) - p0
            end if
            x1 = 0.5_r8 * (x(i-1,j,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i-1,j,k  ) + x(i,j,k  ))
            y(i,j) = (dp2 * x1 + dp1 * x2) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lon_edge_to_pressure_level

  subroutine interp_lat_edge_to_pressure_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%half_lat_lb:mesh%half_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(in) :: po
    real(r8), intent(out) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                               mesh%full_lat_lb:mesh%full_lat_ub)
    logical, intent(in), optional :: logp

    logical logp_
    real(r8) p0, dp1, dp2, x1, x2
    integer i, j, k

    logp_ = merge(logp, .false., present(logp))

    p0 = merge(log(po), po, logp_)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%full_lev_iend, mesh%full_lev_ibeg + 1, -1
          if (p(i,j,k-1) <= po .and. po <= p(i,j,k)) then
            if (logp_) then
              dp1 = p0 - log(p(i,j,k-1))
              dp2 = log(p(i,j,k)) - p0
            else
              dp1 = p0 - p(i,j,k-1)
              dp2 = p(i,j,k) - p0
            end if
            x1 = 0.5_r8 * (x(i,j-1,k-1) + x(i,j,k-1))
            x2 = 0.5_r8 * (x(i,j-1,k  ) + x(i,j,k  ))
            y(i,j) = (dp2 * x1 + dp1 * x2) / (dp1 + dp2)
            exit
          end if
        end do
      end do
    end do

  end subroutine interp_lat_edge_to_pressure_level

  subroutine interp_lev_edge_to_pressure_level(mesh, p, x, po, y, logp)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: p(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(in) :: po
    real(r8), intent(out) :: y(mesh%full_lon_lb:mesh%full_lon_ub, &
                               mesh%full_lat_lb:mesh%full_lat_ub)
    logical, intent(in), optional :: logp

    logical logp_
    real(r8) p0, dp1, dp2, x1, x2
    integer i, j, k

    logp_ = merge(logp, .false., present(logp))

    p0 = merge(log(po), po, logp_)

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        do k = mesh%half_lev_iend, mesh%half_lev_ibeg + 1, -1
          if (p(i,j,k-1) <= po .and. po <= p(i,j,k)) then
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

  end subroutine interp_lev_edge_to_pressure_level

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
          x(i,j,k) = (mesh%area_lat_south(j  ) * x_lat(i,j  ,k) + &
                      mesh%area_lat_north(j-1) * x_lat(i,j-1,k)   &
                     ) / (mesh%area_lat_south(j) + mesh%area_lat_north(j-1))
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
    real(r8) x1, x2, x3, a, b, c

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
      !
      ! -------
      !
      ! ===o=== 3   x3
      x1 = mesh%full_lev(k  ) - mesh%half_lev(k)
      x2 = mesh%full_lev(k+1) - mesh%half_lev(k)
      x3 = mesh%full_lev(k+2) - mesh%half_lev(k)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_lev_lon(i,j,k) = a * x_lon(i,j,k) + b * x_lon(i,j,k+1) + c * x_lon(i,j,k+2)
        end do
      end do
      k = mesh%half_lev_iend
      ! ===o=== NLEV - 2  x3
      !
      ! -------
      !
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = mesh%half_lev(k) - mesh%full_lev(k-1)
      x2 = mesh%half_lev(k) - mesh%full_lev(k-2)
      x3 = mesh%half_lev(k) - mesh%full_lev(k-3)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_lev_lon(i,j,k) = a * x_lon(i,j,k-1) + b * x_lon(i,j,k-2) + c * x_lon(i,j,k-3)
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
    real(r8) x1, x2, x3, a, b, c

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
      !
      ! -------
      !
      ! ===o=== 3   x3
      x1 = mesh%full_lev(k  ) - mesh%half_lev(k)
      x2 = mesh%full_lev(k+1) - mesh%half_lev(k)
      x3 = mesh%full_lev(k+2) - mesh%half_lev(k)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev_lat(i,j,k) = a * x_lat(i,j,k) + b * x_lat(i,j,k+1) + c * x_lat(i,j,k+2)
        end do
      end do
      k = mesh%half_lev_iend
      ! ===o=== NLEV - 2  x3
      !
      ! -------
      !
      ! ===o=== NLEV - 1  x2
      !
      ! -------
      !
      ! ===o=== NLEV      x1
      !
      ! ---?--- NLEV + 1
      x1 = mesh%half_lev(k) - mesh%full_lev(k-1)
      x2 = mesh%half_lev(k) - mesh%full_lev(k-2)
      x3 = mesh%half_lev(k) - mesh%full_lev(k-3)
      a =  x2 * x3 / (x1**2 - x1 * x2 - x1 * x3 + x2 * x3)
      b =  x1 * x3 / (x2**2 - x2 * x1 - x2 * x3 + x1 * x3)
      c =  x1 * x2 / (x3**2 - x3 * x1 - x3 * x2 + x1 * x2)
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev_lat(i,j,k) = a * x_lat(i,j,k-1) + b * x_lat(i,j,k-2) + c * x_lat(i,j,k-3)
        end do
      end do
    end if

  end subroutine interp_lat_edge_to_lev_lat_edge

end module interp_mod
