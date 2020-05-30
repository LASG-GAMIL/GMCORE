module parallel_types_mod

  use mpi

  implicit none

  integer, parameter :: async_v         = 1
  integer, parameter :: async_u         = 2
  integer, parameter :: async_gz        = 3
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
  end type async_type

end module parallel_types_mod
