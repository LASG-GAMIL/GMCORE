module parallel_mod

  use mpi
  use mesh_mod

  implicit none

  private

  public parallel_init
  public parallel_final
  public parallel_fill_halo

  type parallel_info_type
    integer comm
    integer num_proc
    integer rank
    character(30) proc_name
  end type parallel_info_type

  type(parallel_info_type) parallel_info

  interface parallel_fill_halo
    module procedure parallel_fill_halo_2d_r8
  end interface parallel_fill_halo

contains

  subroutine parallel_init()

    integer ierr
    integer n

    parallel_info%comm = MPI_COMM_WORLD

    call MPI_INIT(ierr)
    call MPI_CART_CREATE(MPI_COMM_WORLD, 2, [2, 2], [.true., .false.], .true., parallel_info%comm, ierr)
    call MPI_COMM_SIZE(parallel_info%comm, parallel_info%num_proc, ierr)
    call MPI_COMM_RANK(parallel_info%comm, parallel_info%rank, ierr)
    call MPI_GET_PROCESSOR_NAME(parallel_info%proc_name, n, ierr)

  end subroutine parallel_init

  subroutine parallel_final()

    integer ierr

    call MPI_FINALIZE(ierr)

  end subroutine parallel_final

  subroutine parallel_fill_halo_2d_r8(array, all_halo)

    real(8), intent(inout) :: array(:,:)
    logical, intent(in), optional :: all_halo

  end subroutine parallel_fill_halo_2d_r8

end module parallel_mod
