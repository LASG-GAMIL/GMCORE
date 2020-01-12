module process_mod

  use mpi
  use namelist_mod
  use mesh_mod
  use block_mod

  implicit none

  private

  public process_init
  public process_final
  public proc

  type process_type
    integer comm
    integer id
    type(block_type), allocatable :: blocks(:)
  end type process_type

  type(process_type) proc

contains

  subroutine process_init()

    integer ierr
    integer nproc

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

    if (.not. allocated(proc%blocks)) allocate(proc%blocks(1))

    call proc%blocks(1)%init(       &
      max(global_mesh%lon_halo_width, maxval(reduce_factors)), &
      global_mesh%lat_halo_width,      &
      global_mesh%full_lon_ibeg,  &
      global_mesh%full_lon_iend,    &
      global_mesh%full_lat_ibeg,  &
      global_mesh%full_lat_iend)

  end subroutine process_init

  subroutine process_final()

    integer ierr

    if (allocated(proc%blocks)) deallocate(proc%blocks)

    call MPI_FINALIZE(ierr)

  end subroutine process_final

end module process_mod
