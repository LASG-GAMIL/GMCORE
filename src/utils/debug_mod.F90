module debug_mod

  use, intrinsic :: ieee_arithmetic
  use flogger
  use const_mod
  use namelist_mod
  use mesh_mod
  use static_mod
  use state_mod
  use tend_mod
  use block_mod
  use parallel_mod

  private

  public debug_check_areas
  public debug_check_space_operators
  public debug_print_min_max
  public debug_print_min_max_cell
  public debug_print_min_max_lev_edge
  public debug_print_min_max_lon_edge
  public debug_print_min_max_lat_edge
  public debug_is_inf

  interface debug_print_min_max
    module procedure debug_print_min_max_1d
    module procedure debug_print_min_max_2d
    module procedure debug_print_min_max_3d
  end interface debug_print_min_max

contains

  subroutine debug_check_areas()

    type(mesh_type), pointer :: mesh
    real(r8) total_area
    integer j

    mesh => global_mesh

    total_area = 0.0_r8
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      total_area = total_area + mesh%area_cell(j) * mesh%num_full_lon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-12) then
      call log_error('Failed to calculate cell area!', __FILE__, __LINE__)
    end if

    total_area = 0.0_r8
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      total_area = total_area + mesh%area_vtx(j) * mesh%num_half_lon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-12) then
      call log_error('Failed to calculate vertex area!', __FILE__, __LINE__)
    end if

    total_area = 0.0_r8
    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      total_area = total_area + sum(mesh%area_subcell(:,j)) * mesh%num_full_lon * 2
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-12) then
      call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
    end if

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      if (abs(mesh%area_cell(j) - 2.0d0 * sum(mesh%area_subcell(:,j))) / mesh%area_cell(j) > 1.0d-12) then
        call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
      end if
    end do

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      if (abs(mesh%area_vtx(j) - 2.0_r8 * (mesh%area_subcell(2,j) + mesh%area_subcell(1,j+1))) / mesh%area_vtx(j) > 1.0d-12) then
        call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
      end if
    end do

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      if (abs(mesh%area_lon_north(j) + mesh%area_lon_south(j) - mesh%area_lon(j)) / mesh%area_lon(j) > 1.0d-12) then
        call log_error('Failed to calculate north and south subcell on lon grids!', __FILE__, __LINE__)
      end if
    end do

    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      if (abs(mesh%area_lat_west(j) + mesh%area_lat_east(j) - mesh%area_lat(j)) / mesh%area_lat(j) > 1.0d-11) then
        call log_error('Failed to calculate west and east subcell on lat grids!', __FILE__, __LINE__)
      end if
    end do

    total_area = 0.0_r8
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      total_area = total_area + mesh%area_lon(j) * mesh%num_full_lon
    end do
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      total_area = total_area + mesh%area_lat(j) * mesh%num_full_lon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-9) then
      call log_error('Failed to calculate edge area!', __FILE__, __LINE__)
    end if

    total_area = 0.0_r8
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      total_area = total_area + (mesh%area_lon_north(j) + mesh%area_lon_south(j)) * mesh%num_full_lon
    end do
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      total_area = total_area + (mesh%area_lat_west(j) + mesh%area_lat_east(j)) * mesh%num_full_lon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-9) then
      call log_error('Failed to calculate edge area!', __FILE__, __LINE__)
    end if

  end subroutine debug_check_areas

  subroutine debug_check_space_operators(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type ), intent(in) :: state
    type(tend_type  ), intent(in) :: tend

    integer i, j, k
    type(mesh_type), pointer :: mesh
    real(r8) ip_cf
    real(r8) ip_ke, ip_ke_h, ip_ke_v
    real(r8) ip_pe
    real(r8) ip_mf
    real(r8) ip_ptf

    mesh => state%mesh
    ip_cf  = 0.0_r8
    ip_ke  = 0.0_r8; ip_ke_h = 0.0_r8; ip_ke_v = 0.0_r8
    ip_pe  = 0.0_r8
    ip_mf  = 0.0_r8
    ip_ptf = 0.0_r8

    if (baroclinic .and. hydrostatic) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            ip_cf   = ip_cf   + tend%qhv     (i,j,k) * state%mfx_lon(i,j,k) * mesh%area_lon(j)
            ip_ke_h = ip_ke_h + tend%dkedlon (i,j,k) * state%mfx_lon(i,j,k) * mesh%area_lon(j) * 2
            ip_ke_v = ip_ke_v + ( &
              tend%wedudlev(i,j,k) * state%mfx_lon(i,j,k) + &
              0.5_r8 * state%u_lon(i,j,k)**2 * (state%we_lev_lon(i,j,k+1) - state%we_lev_lon(i,j,k)) &
            ) * mesh%area_lon(j)
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ip_cf   = ip_cf   - tend%qhu     (i,j,k) * state%mfy_lat(i,j,k) * mesh%area_lat(j)
            ip_ke_h = ip_ke_h + tend%dkedlat (i,j,k) * state%mfy_lat(i,j,k) * mesh%area_lat(j) * 2
            ip_ke_v = ip_ke_v + ( &
              tend%wedvdlev(i,j,k) * state%mfy_lat(i,j,k) + &
              0.5_r8 * state%v_lat(i,j,k)**2 * (state%we_lev_lat(i,j,k+1) - state%we_lev_lat(i,j,k)) &
            ) * mesh%area_lat(j)
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ip_ke_h = ip_ke_h + (tend%dmfdlon(i,j,k) + tend%dmfdlat(i,j,k)) * state%ke(i,j,k) * mesh%area_cell(j)
            ip_ptf  = ip_ptf + (tend%dptfdlon(i,j,k) + tend%dptfdlat(i,j,k) + tend%dptfdlev(i,j,k)) * mesh%area_cell(j)
          end do
        end do
      end do

      call global_sum(proc%comm, ip_cf)
      call global_sum(proc%comm, ip_ke_h)
      call global_sum(proc%comm, ip_ke_v)
      call global_sum(proc%comm, ip_ptf)

      if (proc%id == 0) then
        write(6, *) &
          ip_cf   / (4 * pi * radius**2), &
          ip_ke_h / (4 * pi * radius**2), &
          ip_ke_v / (4 * pi * radius**2), &
          ip_ptf  / (4 * pi * radius**2)
      end if
    else
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            ip_cf = ip_cf  + tend%qhv    (i,j,k) * state%mfx_lon(i,j,k) * mesh%area_lon(j)
            ip_ke = ip_ke + tend%dkedlon(i,j,k) * state%mfx_lon(i,j,k) * mesh%area_lon(j) * 2
            ip_pe = ip_pe + tend%pgf_lon(i,j,k) * state%mfx_lon(i,j,k) * mesh%area_lon(j) * 2
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ip_cf = ip_cf - tend%qhu    (i,j,k) * state%mfy_lat(i,j,k) * mesh%area_lat(j)
            ip_ke = ip_ke + tend%dkedlat(i,j,k) * state%mfy_lat(i,j,k) * mesh%area_lat(j) * 2
            ip_pe = ip_pe + tend%pgf_lat(i,j,k) * state%mfy_lat(i,j,k) * mesh%area_lat(j) * 2
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            ip_ke = ip_ke + (tend%dmfdlon(i,j,k) + tend%dmfdlat(i,j,k)) * state%ke(i,j,k) * mesh%area_cell(j)
            ip_pe = ip_pe + (tend%dmfdlon(i,j,k) + tend%dmfdlat(i,j,k)) * state%gz(i,j,k) * mesh%area_cell(j)
            ip_mf = ip_mf + (tend%dmfdlon(i,j,k) + tend%dmfdlat(i,j,k)) * mesh%area_cell(j)
          end do
        end do
      end do

      call global_sum(proc%comm, ip_cf)
      call global_sum(proc%comm, ip_ke)
      call global_sum(proc%comm, ip_pe)

      if (proc%id == 0) then
        write(6, *) &
          ip_cf / (4 * pi * radius**2), &
          ip_ke / (4 * pi * radius**2), &
          ip_pe / (4 * pi * radius**2)
      end if
    end if

  end subroutine debug_check_space_operators

  subroutine debug_print_min_max_1d(array, label)

    real(r8), intent(in) :: array(:)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array), maxval(array)

  end subroutine debug_print_min_max_1d

  subroutine debug_print_min_max_2d(array, label)

    real(r8), intent(in) :: array(:,:)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array), maxval(array)

  end subroutine debug_print_min_max_2d

  subroutine debug_print_min_max_3d(array, label)

    real(r8), intent(in) :: array(:,:,:)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array), maxval(array)

  end subroutine debug_print_min_max_3d

  subroutine debug_print_min_max_cell(mesh, array, label)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: array(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%full_lev_lb:mesh%full_lev_ub)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array(mesh%full_lon_ibeg:mesh%full_lon_iend,   &
                                          mesh%full_lat_ibeg:mesh%full_lat_iend,   &
                                          mesh%full_lev_ibeg:mesh%full_lev_iend)), &
                             maxval(array(mesh%full_lon_ibeg:mesh%full_lon_iend,   &
                                          mesh%full_lat_ibeg:mesh%full_lat_iend,   &
                                          mesh%full_lev_ibeg:mesh%full_lev_iend))

  end subroutine debug_print_min_max_cell

  subroutine debug_print_min_max_lev_edge(mesh, array, label)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: array(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%half_lev_lb:mesh%half_lev_ub)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array(mesh%full_lon_ibeg:mesh%full_lon_iend,   &
                                          mesh%full_lat_ibeg:mesh%full_lat_iend,   &
                                          mesh%half_lev_ibeg:mesh%half_lev_iend)), &
                             maxval(array(mesh%full_lon_ibeg:mesh%full_lon_iend,   &
                                          mesh%full_lat_ibeg:mesh%full_lat_iend,   &
                                          mesh%half_lev_ibeg:mesh%half_lev_iend))

  end subroutine debug_print_min_max_lev_edge

  subroutine debug_print_min_max_lon_edge(mesh, array, label)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: array(mesh%half_lon_lb:mesh%half_lon_ub, &
                                  mesh%full_lat_lb:mesh%full_lat_ub, &
                                  mesh%full_lev_lb:mesh%full_lev_ub)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array(mesh%half_lon_ibeg:mesh%half_lon_iend,   &
                                          mesh%full_lat_ibeg:mesh%full_lat_iend,   &
                                          mesh%full_lev_ibeg:mesh%full_lev_iend)), &
                             maxval(array(mesh%half_lon_ibeg:mesh%half_lon_iend,   &
                                          mesh%full_lat_ibeg:mesh%full_lat_iend,   &
                                          mesh%full_lev_ibeg:mesh%full_lev_iend))

  end subroutine debug_print_min_max_lon_edge

  subroutine debug_print_min_max_lat_edge(mesh, array, label)

    type(mesh_type), intent(in) :: mesh
    real(r8), intent(in) :: array(mesh%full_lon_lb:mesh%full_lon_ub, &
                                  mesh%half_lat_lb:mesh%half_lat_ub, &
                                  mesh%full_lev_lb:mesh%full_lev_ub)
    character(*), intent(in) :: label

    write(6, *) trim(label), minval(array(mesh%full_lon_ibeg:mesh%full_lon_iend,   &
                                          mesh%half_lat_ibeg:mesh%half_lat_iend,   &
                                          mesh%full_lev_ibeg:mesh%full_lev_iend)), &
                             maxval(array(mesh%full_lon_ibeg:mesh%full_lon_iend,   &
                                          mesh%half_lat_ibeg:mesh%half_lat_iend,   &
                                          mesh%full_lev_ibeg:mesh%full_lev_iend))

  end subroutine debug_print_min_max_lat_edge

  logical function debug_is_inf(x) result(res)

    real(r8), intent(in) :: x

    res = .not. ieee_is_finite(x)

  end function debug_is_inf

end module debug_mod
