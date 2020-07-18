module interp_mod

  use const_mod
  use block_mod

  implicit none

  private

  public interp_cell_to_edge_on_full_level
  public interp_cell_to_edge_on_half_level
  public interp_full_level_to_half_level_on_cell
  public interp_cell_to_isobaric_level

contains

  subroutine interp_cell_to_edge_on_full_level(mesh, x, x_lon, x_lat)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_lon(mesh%half_lon_lb:mesh%half_lon_ub, &
                                     mesh%full_lat_lb:mesh%full_lat_ub, &
                                     mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                     mesh%half_lat_lb:mesh%half_lat_ub, &
                                     mesh%full_lev_lb:mesh%full_lev_ub)

    integer i, j, k

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_lon(i,j,k) = (mesh%area_lon_west(j) * x(i,j,k) + mesh%area_lon_east(j) * x(i+1,j,k)) / mesh%area_lon(j)
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          x_lat(i,j,k) = (mesh%area_lat_north(j) * x(i,j,k) + mesh%area_lat_south(j) * x(i,j-1,k)) / mesh%area_lat(j)
#else
          x_lat(i,j,k) = (mesh%area_lat_north(j) * x(i,j+1,k) + mesh%area_lat_south(j) * x(i,j,k)) / mesh%area_lat(j)
#endif
        end do
      end do
    end do

  end subroutine interp_cell_to_edge_on_full_level

  subroutine interp_cell_to_edge_on_half_level(mesh, x, x_lon, x_lat)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(inout) :: x_lon(mesh%half_lon_lb:mesh%half_lon_ub, &
                                     mesh%full_lat_lb:mesh%full_lat_ub, &
                                     mesh%half_lev_lb:mesh%half_lev_ub)
    real(r8), intent(inout) :: x_lat(mesh%full_lon_lb:mesh%full_lon_ub, &
                                     mesh%half_lat_lb:mesh%half_lat_ub, &
                                     mesh%half_lev_lb:mesh%half_lev_ub)

    integer i, j, k

    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          x_lon(i,j,k) = (mesh%area_lon_west(j) * x(i,j,k) + mesh%area_lon_east(j) * x(i+1,j,k)) / mesh%area_lon(j)
        end do
      end do
    end do

    do k = mesh%half_lev_ibeg, mesh%half_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          x_lat(i,j,k) = (mesh%area_lat_north(j) * x(i,j,k) + mesh%area_lat_south(j) * x(i,j-1,k)) / mesh%area_lat(j)
#else
          x_lat(i,j,k) = (mesh%area_lat_north(j) * x(i,j+1,k) + mesh%area_lat_south(j) * x(i,j,k)) / mesh%area_lat(j)
#endif
        end do
      end do
    end do

  end subroutine interp_cell_to_edge_on_half_level

  subroutine interp_full_level_to_half_level_on_cell(mesh, x, x_lev)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: x(mesh%full_lon_lb:mesh%full_lon_ub, &
                              mesh%full_lat_lb:mesh%full_lat_ub, &
                              mesh%full_lev_lb:mesh%full_lev_ub)
    real(r8), intent(inout) :: x_lev(mesh%full_lon_lb:mesh%full_lon_ub, &
                                     mesh%full_lat_lb:mesh%full_lat_ub, &
                                     mesh%half_lev_lb:mesh%half_lev_ub)

    integer i, j, k
    real(r8) deta1, deta2, deta3

    do k = mesh%half_lev_ibeg + 1, mesh%half_lev_iend - 1
      deta1 = mesh%half_lev(k) - mesh%half_lev(k-1)
      deta2 = mesh%half_lev(k+1) - mesh%half_lev(k)
      deta3 = deta1 + deta2
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          x_lev(i,j,k) = (x(i,j,k) * deta1 + x(i,j,k-1) * deta2) / deta3
        end do
      end do
    end do

  end subroutine interp_full_level_to_half_level_on_cell

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

end module interp_mod
