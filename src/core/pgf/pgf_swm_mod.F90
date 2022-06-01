module pgf_swm_mod

  use const_mod
  use parallel_mod
  use block_mod

  implicit none

contains

  subroutine pgf_swm_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

  end subroutine pgf_swm_prepare

  subroutine pgf_swm_run(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j, k

    associate (mesh    => block%mesh  , &
               gz      => state%gz_f  , & ! in
               pgf_lon => tend%pgf_lon, & ! out
               pgf_lat => tend%pgf_lat)   ! out
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          pgf_lon(i,j,k) = (gz(i+1,j,k) - gz(i,j,k)) / mesh%de_lon(j)
        end do
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          pgf_lat(i,j,k) = (gz(i,j+1,k) - gz(i,j,k)) / mesh%de_lat(j)
        end do
      end do
    end do
    end associate

  end subroutine pgf_swm_run

end module pgf_swm_mod
