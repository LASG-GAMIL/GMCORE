module debug_mod

  use flogger
  use const_mod
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

    total_area = 0.0d0
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

#ifdef V_POLE
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      if (mesh%is_south_pole(j)) then
        if (abs(mesh%area_vtx(j) - 2.0d0 * mesh%area_subcell(1,j)) / mesh%area_vtx(j) > 1.0d-12) then
          call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
        end if
      else if (mesh%is_north_pole(j)) then
        if (abs(mesh%area_vtx(j) - 2.0d0 * mesh%area_subcell(2,j-1)) / mesh%area_vtx(j) > 1.0d-12) then
          call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
        end if
      else
        if (abs(mesh%area_vtx(j) - 2.0d0 * (mesh%area_subcell(2,j-1) + mesh%area_subcell(1,j))) / mesh%area_vtx(j) > 1.0d-12) then
          call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
        end if
      end if
    end do
#else
    do j = mesh%half_lat_ibeg, mesh%half_lat_iend
      if (abs(mesh%area_vtx(j) - 2.0_r8 * (mesh%area_subcell(2,j) + mesh%area_subcell(1,j+1))) / mesh%area_vtx(j) > 1.0d-12) then
        call log_error('Failed to calculate subcell area!', __FILE__, __LINE__)
      end if
    end do
#endif

    total_area = 0.0d0
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      total_area = total_area + mesh%area_lon(j) * mesh%num_full_lon
    end do
    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      total_area = total_area + mesh%area_lat(j) * mesh%num_full_lon
    end do
    if (abs(global_mesh%total_area - total_area) / global_mesh%total_area > 1.0d-9) then
      call log_error('Failed to calculate edge area!', __FILE__, __LINE__)
    end if

  end subroutine debug_check_areas

  subroutine debug_check_space_operators(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type ), intent(in) :: state
    type(tend_type  ), intent(in) :: tend

    integer i, j
    type(mesh_type), pointer :: mesh
    real(r8) ip_cf
    real(r8) ip_dke
    real(r8) ip_dpe
    real(r8) ip_dmf

    mesh => state%mesh
    ip_cf  = 0.0_r8
    ip_dke = 0.0_r8
    ip_dpe = 0.0_r8
    ip_dmf = 0.0_r8

    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        ip_cf = ip_cf + tend%qhv(i,j) * state%mf_lon_n(i,j) * mesh%area_lon(j)
        ip_dke = ip_dke + tend%dkedlon(i,j) * state%mf_lon_n(i,j) * mesh%area_lon(j) * 2
        ip_dpe = ip_dpe + tend%dpedlon(i,j) * state%mf_lon_n(i,j) * mesh%area_lon(j) * 2
      end do
    end do

    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        ip_cf = ip_cf - tend%qhu(i,j) * state%mf_lat_n(i,j) * mesh%area_lat(j)
        ip_dke = ip_dke + tend%dkedlat(i,j) * state%mf_lat_n(i,j) * mesh%area_lat(j) * 2
        ip_dpe = ip_dpe + tend%dpedlat(i,j) * state%mf_lat_n(i,j) * mesh%area_lat(j) * 2
      end do
    end do

    do j = mesh%full_lat_ibeg, mesh%full_lat_iend
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        ip_dke = ip_dke + (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * state%ke(i,j) * mesh%area_cell(j)
        ip_dpe = ip_dpe + (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * state%gz(i,j,1) * mesh%area_cell(j)
        ip_dmf = ip_dmf + (tend%dmfdlon(i,j) + tend%dmfdlat(i,j)) * mesh%area_cell(j)
      end do
    end do

    call log_add_diag('cf' , ip_cf  / radius**2)
    call log_add_diag('dke', ip_dke / radius**2)
    call log_add_diag('dpe', ip_dpe / radius**2)
    ! call log_add_diag('dmf', ip_dmf / radius**2)

  end subroutine debug_check_space_operators

  subroutine debug_print_min_max(array)

    real(r8), intent(in) :: array(:,:)

    print *, minval(array), maxval(array)

  end subroutine debug_print_min_max

end module debug_mod
