module parallel_types_mod

  use mpi

  implicit none

  integer, parameter :: async_v         = 1
  integer, parameter :: async_u         = 2
  integer, parameter :: async_gd        = 3
  integer, parameter :: async_pv        = 4
  integer, parameter :: async_pv_lon    = 5
  integer, parameter :: async_pv_lat    = 6
  integer, parameter :: async_ke        = 7
  integer, parameter :: async_mf_lon_n  = 8
  integer, parameter :: async_mf_lat_n  = 9
  integer, parameter :: async_dpv_lon_t = 10
  integer, parameter :: async_dpv_lat_t = 11

  type async_type
    integer :: west_send_req  = MPI_REQUEST_NULL
    integer :: west_recv_req  = MPI_REQUEST_NULL
    integer :: east_send_req  = MPI_REQUEST_NULL
    integer :: east_recv_req  = MPI_REQUEST_NULL
    integer :: south_send_req = MPI_REQUEST_NULL
    integer :: south_recv_req = MPI_REQUEST_NULL
    integer :: north_send_req = MPI_REQUEST_NULL
    integer :: north_recv_req = MPI_REQUEST_NULL
    integer :: j = 0
    integer, allocatable :: zonal_west_send_req(:)
    integer, allocatable :: zonal_west_recv_req(:)
    integer, allocatable :: zonal_east_send_req(:)
    integer, allocatable :: zonal_east_recv_req(:)
  contains
    procedure :: init => async_init
    ! FIXME: gfortran <= 9 has sigfault bug for this finalizer when async_type array is used in other types.
    ! final :: async_final
  end type async_type

contains

  subroutine async_init(this, num_async_lat)

    class(async_type), intent(inout) :: this
    integer, intent(in) :: num_async_lat

    if (allocated(this%zonal_west_send_req)) deallocate(this%zonal_west_send_req)
    if (allocated(this%zonal_west_recv_req)) deallocate(this%zonal_west_recv_req)
    if (allocated(this%zonal_east_send_req)) deallocate(this%zonal_east_send_req)
    if (allocated(this%zonal_east_recv_req)) deallocate(this%zonal_east_recv_req)

    allocate(this%zonal_west_send_req(num_async_lat))
    allocate(this%zonal_west_recv_req(num_async_lat))
    allocate(this%zonal_east_send_req(num_async_lat))
    allocate(this%zonal_east_recv_req(num_async_lat))

  end subroutine async_init

  subroutine async_final(this)

    type(async_type), intent(inout) :: this

    if (allocated(this%zonal_west_send_req)) deallocate(this%zonal_west_send_req)
    if (allocated(this%zonal_west_recv_req)) deallocate(this%zonal_west_recv_req)
    if (allocated(this%zonal_east_send_req)) deallocate(this%zonal_east_send_req)
    if (allocated(this%zonal_east_recv_req)) deallocate(this%zonal_east_recv_req)

  end subroutine async_final

end module parallel_types_mod