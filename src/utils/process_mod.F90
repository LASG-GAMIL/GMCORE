module process_mod

  use mpi
  use namelist_mod
  use mesh_mod
  use block_mod

  implicit none

  private

  public process_init
  public process_final
  public process

  type process_type
    integer comm
    integer id
    type(block_type), allocatable :: blocks(:)
  end type process_type

  type(process_type) process

contains

  subroutine process_init()

    integer ierr
    integer nproc

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

    if (.not. allocated(process%blocks)) allocate(process%blocks(1))

    call process%blocks(1)%init(       &
      max(global_mesh%lon_halo_width, maxval(reduce_factors)), &
      global_mesh%lat_halo_width,      &
      global_mesh%full_lon_start_idx,  &
      global_mesh%full_lon_end_idx,    &
      global_mesh%full_lat_start_idx,  &
      global_mesh%full_lat_end_idx)

  end subroutine process_init

  subroutine process_final()

    integer ierr

    if (allocated(process%blocks)) deallocate(process%blocks)

    call MPI_FINALIZE(ierr)

  end subroutine process_final

end module process_mod
