module parallel_types_mod

  use mpi

  implicit none

  private

  public async_type

  integer, public, parameter :: async_v         = 1
  integer, public, parameter :: async_u         = 2
  integer, public, parameter :: async_gz        = 3
  integer, public, parameter :: async_pv        = 4
  integer, public, parameter :: async_pv_lon    = 5
  integer, public, parameter :: async_pv_lat    = 6
  integer, public, parameter :: async_ke        = 7
  integer, public, parameter :: async_mf_lon_n  = 8
  integer, public, parameter :: async_mf_lat_n  = 9
  integer, public, parameter :: async_dpv_lon_t = 10
  integer, public, parameter :: async_dpv_lat_t = 11

  type async_type
    integer, allocatable :: send_req(:)
    integer, allocatable :: recv_req(:)
  contains
    procedure :: init => async_init
    procedure :: wait => async_wait
    final :: async_final
  end type async_type

contains

  subroutine async_init(this, num_ngb)

    class(async_type), intent(inout) :: this
    integer, intent(in) :: num_ngb

    if (allocated(this%send_req)) deallocate(this%send_req)
    if (allocated(this%recv_req)) deallocate(this%recv_req)

    allocate(this%send_req(num_ngb)); this%send_req = MPI_REQUEST_NULL
    allocate(this%recv_req(num_ngb)); this%recv_req = MPI_REQUEST_NULL

  end subroutine async_init

  subroutine async_wait(this)

    class(async_type), intent(inout) :: this

    integer i, ierr

    do i = 1, size(this%send_req)
      if (this%send_req(i) /= MPI_REQUEST_NULL) then
        call MPI_WAIT(this%send_req(i), MPI_STATUS_IGNORE, ierr)
      end if
      if (this%recv_req(i) /= MPI_REQUEST_NULL) then
        call MPI_WAIT(this%recv_req(i), MPI_STATUS_IGNORE, ierr)
      end if
    end do

  end subroutine async_wait

  subroutine async_final(this)

    type(async_type), intent(inout) :: this

    if (allocated(this%send_req)) deallocate(this%send_req)
    if (allocated(this%recv_req)) deallocate(this%recv_req)

  end subroutine async_final

end module parallel_types_mod
