module pgf_lin97_mod

  use const_mod
  use namelist_mod
  use block_mod

  implicit none

contains

  subroutine pgf_lin97_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

  end subroutine pgf_lin97_prepare

  subroutine pgf_lin97_run(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh

    if (baroclinic .and. hydrostatic) then
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            
          end do
        end do
      end do
    end if

  end subroutine pgf_lin97_run

end module pgf_lin97_mod
