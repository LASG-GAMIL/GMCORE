module lat_damp_mod

  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod

  implicit none

  private

  public lat_damp_init
  public lat_damp_final
  public lat_damp_on_cell

  interface lat_damp_on_cell
    module procedure lat_damp_on_cell_2d
  end interface lat_damp_on_cell

  real(r8), allocatable :: wgt_full_lat(:)
  real(r8), allocatable :: wgt_half_lat(:)

  real(r8), parameter :: coef = 0.1_r8

contains

  subroutine lat_damp_init()

    integer j
    real(r8) lat0

    allocate(wgt_full_lat(global_mesh%num_full_lat)); wgt_full_lat = 0
    lat0 = abs(global_mesh%full_lat_deg(2))
    do j = 2, global_mesh%num_full_lat - 1
      wgt_full_lat(j) = exp((abs(global_mesh%full_lat_deg(j)) - lat0)**2 / (60 - lat0)**2 * log(1.0e-3_r8))
    end do
    allocate(wgt_half_lat(global_mesh%num_half_lat)); wgt_half_lat = 0
    lat0 = abs(global_mesh%half_lat_deg(2))
    do j = 1, global_mesh%num_half_lat
      wgt_half_lat(j) = exp((abs(global_mesh%half_lat_deg(j)) - lat0)**2 / (60 - lat0)**2 * log(1.0e-3_r8))
    end do

  end subroutine lat_damp_init

  subroutine lat_damp_final()

    if (allocated(wgt_full_lat)) deallocate(wgt_full_lat)
    if (allocated(wgt_half_lat)) deallocate(wgt_half_lat)

  end subroutine lat_damp_final

  subroutine lat_damp_on_cell_2d(block, f)

    type(block_type), intent(in), target :: block
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub)

    type(mesh_type), pointer :: mesh
    real(r8) tmp(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                 block%mesh%full_lat_lb:block%mesh%full_lat_ub)
    integer i, j

    mesh => block%mesh
    tmp = f
    do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
      do i = mesh%full_lon_ibeg, mesh%full_lon_iend
        f(i,j) = (1 - wgt_full_lat(j)) * tmp(i,j) + wgt_full_lat(j) * ( &
          (1 - 2 * coef) * tmp(i,j) + coef * (tmp(i,j-1) + tmp(i,j+1))  &
        )
      end do
    end do
    call fill_halo(block, f, full_lon=.true., full_lat=.true.)

  end subroutine lat_damp_on_cell_2d

end module lat_damp_mod
