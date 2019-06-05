module parallel_mod

  use mpi
  use array_mod

  implicit none

  private

  type halo_grid_type
    integer, allocatable :: 
  end type halo_grid_type

  type neighbor_domain_type
    type(halo_map_type) halo_map
    type(distributed_domain_type), pointer :: domain
  end type neighbor_domain_type

  type distributed_domain_type
    integer full_i_start
    integer full_i_end
    integer full_j_start
    integer full_j_end
    type(array_mod) neighbor_domains
  end type distributed_domain_type

  type(distributed_domain_type) domain

contains

  subroutine parallel_init()


  end subroutine parallel_init

  subroutine parallel_fill_halo_r8(array)

    real(8), intent(inout) :: array(

  end subroutine parallel_fill_halo_r8

end module parallel_mod
