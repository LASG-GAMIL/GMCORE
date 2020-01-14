module parallel_mod

  use mpi
  use mesh_mod
  use block_mod

  implicit none

  private

  public fill_halo
  public zero_halo
  public zonal_sum
  public global_sum
  public overlay_inner_halo

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

  interface global_sum
    module procedure global_sum_0d_r8
  end interface global_sum

contains

  subroutine fill_halo_1d_r8_1(halo_width, array, west_halo, east_halo)

    integer, intent(in   )           :: halo_width
    real(8), intent(inout)           :: array(:)
    logical, intent(in   ), optional :: west_halo
    logical, intent(in   ), optional :: east_halo

    integer i, m, n

    if (merge(west_halo, .true., present(west_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * halo_width
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      m = ubound(array, 1) - halo_width
      n = lbound(array, 1) + halo_width - 1
      do i = 1, halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine fill_halo_1d_r8_1

  subroutine fill_halo_1d_r8_2(block, array, west_halo, east_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo

    integer i, m, n

    if (merge(west_halo, .true., present(west_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * block%mesh%lon_halo_width
      do i = 1, block%mesh%lon_halo_width
        array(m+i) = array(n+i)
      end do
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      m = ubound(array, 1) - block%mesh%lon_halo_width
      n = lbound(array, 1) + block%mesh%lon_halo_width - 1
      do i = 1, block%mesh%lon_halo_width
        array(m+i) = array(n+i)
      end do
    end if

  end subroutine fill_halo_1d_r8_2

  subroutine fill_halo_2d_r8(block, array, west_halo, east_halo, top_halo, bottom_halo)

    type(block_type), intent(in) :: block
    real(8), intent(inout) :: array(:,:)
    logical, intent(in), optional :: west_halo
    logical, intent(in), optional :: east_halo
    logical, intent(in), optional :: top_halo
    logical, intent(in), optional :: bottom_halo

    integer i, j, m, n

    if (merge(west_halo, .true., present(west_halo))) then
      m = lbound(array, 1) - 1
      n = ubound(array, 1) - 2 * block%mesh%lon_halo_width
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, block%mesh%lon_halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

    if (merge(east_halo, .true., present(east_halo))) then
      m = ubound(array, 1) - block%mesh%lon_halo_width
      n = lbound(array, 1) + block%mesh%lon_halo_width - 1
      do j = lbound(array, 2), ubound(array, 2)
        do i = 1, block%mesh%lon_halo_width
          array(m+i,j) = array(n+i,j)
        end do
      end do
    end if

  end subroutine fill_halo_2d_r8

  subroutine zero_halo_1d_r8(mesh, array, west_halo, east_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(mesh%full_lon_lb:mesh%full_lon_ub)
    logical        , intent(in   ), optional :: west_halo
    logical        , intent(in   ), optional :: east_halo

    if (merge(west_halo, .false., present(west_halo))) then
      array(mesh%full_lon_lb:mesh%full_lon_ibeg-1) = 0.0d0
    end if

    if (merge(east_halo, .false., present(east_halo))) then
      array(mesh%full_lon_iend+1:mesh%full_lon_ub) = 0.0d0
    end if

  end subroutine zero_halo_1d_r8

  subroutine overlay_inner_halo(mesh, array, west_halo, east_halo)

    type(mesh_type), intent(in   )           :: mesh
    real(8)        , intent(inout)           :: array(mesh%full_lon_lb:mesh%full_lon_ub)
    logical        , intent(in   ), optional :: west_halo
    logical        , intent(in   ), optional :: east_halo

    integer i

    if (merge(west_halo, .false., present(west_halo))) then
      do i = mesh%full_lon_ibeg, mesh%full_lon_ibeg + mesh%lon_halo_width - 2
        array(i) = array(i) + array(mesh%full_lon_iend+i-mesh%full_lon_ibeg+1)
      end do
    end if

  end subroutine overlay_inner_halo

  subroutine zonal_sum_0d_r8(value)

    real(8), intent(inout) :: value

  end subroutine zonal_sum_0d_r8

  subroutine global_sum_0d_r8(comm, value)

    integer, intent(in) :: comm
    real(8), intent(inout) :: value

    integer ierr
    real(8) res

    call MPI_ALLREDUCE(value, res, 1, MPI_DOUBLE, MPI_SUM, comm, ierr)
    value = res

  end subroutine global_sum_0d_r8

end module parallel_mod