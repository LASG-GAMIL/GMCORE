module parallel_mod

  use mesh_mod

  implicit none

  interface fill_halo
    module procedure fill_halo_1d_r8_1
    module procedure fill_halo_1d_r8_2
    module procedure fill_halo_2d_r8
  end interface fill_halo

  interface zero_halo
    module procedure zero_halo_1d_r8
  end interface zero_halo

  interface zonal_sum
    module procedure zonal_sum_0d_r8
  end interface zonal_sum

contains

  subroutine fill_halo_1d_r8_1(halo_width, array, left_halo, right_halo)

    integer, intent(in   )           :: halo_width
    real(8), intent(inout)           :: array(:)
    logical, intent(in   ), optional :: left_halo
    logical, intent(in   ), optional :: right_halo

    integer i, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * halo_width
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - halo_width
      n = lbound(array, 1) + halo_width - 1
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine fill_halo_1d_r8_1

  subroutine fill_halo_1d_r8_2(mesh, array, left_halo, right_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(:)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo

    integer i, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * mesh%lon_halo_width
      do i = 1, mesh%lon_halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - mesh%lon_halo_width
      n = lbound(array, 1) + mesh%lon_halo_width - 1
      do i = 1, mesh%lon_halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine fill_halo_1d_r8_2

  subroutine fill_halo_2d_r8(mesh, array, left_halo, right_halo, top_halo, bottom_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(:,:)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo
    logical        , intent(in   ), optional :: top_halo
    logical        , intent(in   ), optional :: bottom_halo

    integer i, j, m, n

    if (merge(left_halo, .true., present(left_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * mesh%lon_halo_width
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, mesh%lon_halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

    if (merge(right_halo, .true., present(right_halo))) then
      m = ubound(array, 1) - mesh%lon_halo_width
      n = lbound(array, 1) + mesh%lon_halo_width - 1
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, mesh%lon_halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

  end subroutine fill_halo_2d_r8

  subroutine zero_halo_1d_r8(mesh, array, left_halo, right_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(mesh%full_lon_lb:mesh%full_lon_ub)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo

    if (merge(left_halo, .false., present(left_halo))) then
      array(mesh%full_lon_lb:mesh%full_lon_start_idx-1) = 0.0d0
    end if

    if (merge(right_halo, .false., present(right_halo))) then
      array(mesh%full_lon_end_idx+1:mesh%full_lon_ub) = 0.0d0
    end if

  end subroutine zero_halo_1d_r8

  subroutine overlay_inner_halo(mesh, array, left_halo, right_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(mesh%full_lon_lb:mesh%full_lon_ub)
    logical        , intent(in   ), optional :: left_halo
    logical        , intent(in   ), optional :: right_halo

    integer i

    if (merge(left_halo, .false., present(left_halo))) then
      do i = mesh%full_lon_start_idx, mesh%full_lon_start_idx + mesh%lon_halo_width - 2
        array(i) = array(i) + array(mesh%full_lon_end_idx+i-mesh%full_lon_start_idx+1)
      end do
    end if

  end subroutine overlay_inner_halo

  subroutine zonal_sum_0d_r8(value)

    real(8), intent(inout) :: value

  end subroutine zonal_sum_0d_r8

end module parallel_mod