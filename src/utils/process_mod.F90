module process_mod

  use mpi
  use flogger
  use string
  use namelist_mod
  use mesh_mod
  use block_mod

  implicit none

  private

  public process_init
  public process_stop
  public process_final
  public proc
  public is_root_proc

  type process_type
    integer comm
    integer group
    integer :: zonal_comm = MPI_COMM_NULL
    integer :: zonal_group = MPI_GROUP_NULL
    integer cart_dims(2)
    integer cart_coords(2)
    integer id                              ! MPI process ID
    integer idom                            ! Nest domain index (root domain is 1)
    integer lon_ibeg                        ! Index of beginning longitude grid
    integer lon_iend                        ! Index of ending longitude grid
    integer lat_ibeg                        ! Index of beginning latitude grid
    integer lat_iend                        ! Index of ending latitude grid
    integer, allocatable :: parent_ngb(:)   ! Neighbor  parent domain indices
    integer, allocatable ::  child_ngb(:)   ! Neighbor   child domain indices
    integer :: ngb(4) = MPI_PROC_NULL       ! Neighbor sibling domain indices
    type(block_type), allocatable :: blocks(:)
  end type process_type

  type(process_type) proc

contains

  subroutine process_init()

    integer ierr
    integer num_total_lon, num_total_lat

    call setup_mpi()
    call decompose_domains(proc%idom, proc%cart_dims, proc%cart_coords, num_total_lon, num_total_lat, &
                           proc%lon_ibeg, proc%lon_iend, proc%lat_ibeg, proc%lat_iend)
    call setup_zonal_comm_for_reduce(num_total_lat, proc%lat_ibeg, proc%lat_iend)
    call connect_parent() ! <-- FIXME: Needs implementation.

    if (.not. allocated(proc%blocks)) allocate(proc%blocks(1))

    call proc%blocks(1)%init(proc%id, global_mesh%lon_halo_width, global_mesh%lat_halo_width, &
                             proc%lon_ibeg, proc%lon_iend, proc%lat_ibeg, proc%lat_iend)

    ! Setup halos (only normal halos for the time being).
    allocate(proc%blocks(1)%halo(4))
    call proc%blocks(1)%halo(1)%init_normal(proc%blocks(1)%mesh, ngb_proc_id=proc%ngb(1), &
                                            west_lat_ibeg=proc%lat_ibeg, west_lat_iend=proc%lat_iend)
    call proc%blocks(1)%halo(2)%init_normal(proc%blocks(1)%mesh, ngb_proc_id=proc%ngb(2), &
                                            east_lat_ibeg=proc%lat_ibeg, east_lat_iend=proc%lat_iend)
    call proc%blocks(1)%halo(3)%init_normal(proc%blocks(1)%mesh, ngb_proc_id=proc%ngb(3), &
                                            south_lon_ibeg=proc%lon_ibeg, south_lon_iend=proc%lon_iend)
    call proc%blocks(1)%halo(4)%init_normal(proc%blocks(1)%mesh, ngb_proc_id=proc%ngb(4), &
                                            north_lon_ibeg=proc%lon_ibeg, north_lon_iend=proc%lon_iend)

  end subroutine process_init

  subroutine process_stop(code)

    integer, intent(in) :: code

    integer ierr

    call MPI_BARRIER(proc%comm, ierr)
    call MPI_ABORT(proc%comm, code, ierr)

  end subroutine process_stop

  subroutine process_final()

    integer ierr

    if (allocated(proc%blocks)) deallocate(proc%blocks)

    call MPI_FINALIZE(ierr)

  end subroutine process_final

  pure logical function is_root_proc()

    is_root_proc = proc%id == 0

  end function is_root_proc

  subroutine setup_mpi()

    integer ierr, np, np2, pid, tmp_comm, i
    logical periods(2)

    call MPI_INIT(ierr)
    ! Get basic MPI information.
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, pid, ierr)

    if (num_proc_lon(1) /= 0 .and. num_proc_lat(1) /= 0) then
      ! Check if process topology in namelist is compatible with MPI runtime.
      np2 = 0
      do i = 1, nest_max_dom
        np2 = np2 + num_proc_lon(i) * num_proc_lat(i)
      end do
      if (np /= np2 .and. is_root_proc()) then
        call log_error('Namelist num_proc_lon and num_proc_lat are not compatible with MPI runtime!')
      end if
      ! Set the process topology into proc object.
      np2 = 0
      do i = 1, nest_max_dom
        np2 = np2 + num_proc_lon(i) * num_proc_lat(i)
        if (pid + 1 <= np2) then
          proc%cart_dims(1) = num_proc_lon(i)
          proc%cart_dims(2) = num_proc_lat(i)
          proc%idom = i
          exit
        end if
      end do
    else
      proc%cart_dims = 0
      call MPI_DIMS_CREATE(np, 2, proc%cart_dims, ierr)
    end if
    periods = [.true.,.false.]
    if (proc%idom > 1 .and. (nest_lon_beg(proc%idom) /= 0 .or. nest_lon_end(proc%idom) /= 360)) then
      periods(1) = .false.
    end if
    ! Set MPI process topology.
    ! - Create new communicator for each domain.
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, proc%idom, pid, tmp_comm, ierr)
    call MPI_CART_CREATE(tmp_comm, 2, proc%cart_dims, periods, .true., proc%comm, ierr)
    call MPI_COMM_GROUP(proc%comm, proc%group, ierr)
    ! - Get the basic information and neighborhood.
    call MPI_COMM_RANK(proc%comm, proc%id, ierr)
    call MPI_CART_COORDS(proc%comm, proc%id, 2, proc%cart_coords, ierr)
    call MPI_CART_SHIFT(proc%comm, 0, 1, proc%ngb(1), proc%ngb(2), ierr)
    call MPI_CART_SHIFT(proc%comm, 1, 1, proc%ngb(3), proc%ngb(4), ierr)
    call MPI_COMM_FREE(tmp_comm, ierr)

  end subroutine setup_mpi

  subroutine decompose_domains(idom, cart_dims, cart_coords, num_total_lon, num_total_lat, lon_ibeg, lon_iend, lat_ibeg, lat_iend)

    integer, intent(in ) :: idom
    integer, intent(in ) :: cart_dims(2)
    integer, intent(in ) :: cart_coords(2)
    integer, intent(out) :: num_total_lon
    integer, intent(out) :: num_total_lat
    integer, intent(out) :: lon_ibeg
    integer, intent(out) :: lon_iend
    integer, intent(out) :: lat_ibeg
    integer, intent(out) :: lat_iend

    integer ierr, i, j
    integer num_lon, num_lat, res_num, half_num

    if (idom > 1) then
      ! Get the start and end indices according to nest domain range.
      ! Zonal direction
      lon_ibeg = 0; lon_iend = 0
      do i = 1, global_mesh%num_full_lon
        if (global_mesh%full_lon_deg(i) >= nest_lon_beg(idom-1)) then
          lon_ibeg = i
          exit
        end if
      end do
      do i = 1, global_mesh%num_full_lon
        if (global_mesh%full_lon_deg(i) >= nest_lon_end(idom-1)) then
          lon_iend = i
          exit
        end if
      end do
      if (lon_iend == 0) lon_iend = global_mesh%num_full_lon
      num_total_lon = lon_iend - lon_ibeg + 1
      ! Meridional direction
      lat_ibeg = 0; lat_iend = 0
#ifdef V_POLE
      do j = 1, global_mesh%num_half_lat
        if (global_mesh%half_lat_deg(j) >= nest_lat_beg(idom-1)) then
          lat_ibeg = j
          exit
        end if
      end do
      do j = 1, global_mesh%num_half_lat
        if (global_mesh%half_lat_deg(j) >= nest_lat_end(idom-1)) then
          lat_iend = j
          exit
        end if
      end do
      if (lat_iend == 0) lat_iend = global_mesh%num_full_lat
#else
      do j = 1, global_mesh%num_full_lat
        if (global_mesh%full_lat_deg(j) >= nest_lat_beg(idom-1)) then
          lat_ibeg = j
          exit
        end if
      end do
      do j = 1, global_mesh%num_full_lat
        if (global_mesh%full_lat_deg(j) >= nest_lat_end(idom-1)) then
          lat_iend = j
          exit
        end if
      end do
      if (lat_iend == 0) lat_iend = global_mesh%num_full_lat
#endif
      num_total_lat = lat_iend - lat_ibeg + 1
    else
      lon_ibeg = 1
      lat_ibeg = 1
      num_total_lon = global_mesh%num_full_lon
#ifdef V_POLE
      num_total_lat = global_mesh%num_half_lat
#else
      num_total_lat = global_mesh%num_full_lat
#endif
    end if

    res_num = mod(num_total_lon, cart_dims(1))
    do i = 0, cart_coords(1) - 1
      if (res_num /= 0 .and. i < res_num - 1) then
        num_lon = num_total_lon / cart_dims(1) + 1
      else
        num_lon = num_total_lon / cart_dims(1)
      end if
      lon_ibeg = lon_ibeg + num_lon
    end do
    if (res_num /= 0 .and. cart_coords(1) < res_num) then
      num_lon = num_total_lon / cart_dims(1) + 1
    else
      num_lon = num_total_lon / cart_dims(1)
    end if
    lon_iend = lon_ibeg + num_lon - 1

    res_num = mod(num_total_lat, cart_dims(2))
    half_num = cart_dims(2) - cart_dims(2) / 2
    do j = 0, cart_coords(2) - 1
      if (res_num /= 0 .and. j >= half_num .and. j < half_num + res_num) then
        num_lat = num_total_lat / cart_dims(2) + 1
      else
        num_lat = num_total_lat / cart_dims(2)
      end if
      lat_ibeg = lat_ibeg + num_lat
    end do
    if (res_num /= 0 .and. cart_coords(2) >= half_num .and. cart_coords(2) < half_num + res_num) then
      num_lat = num_total_lat / cart_dims(2) + 1
    else
      num_lat = num_total_lat / cart_dims(2)
    end if
    lat_iend = lat_ibeg + num_lat - 1

  end subroutine decompose_domains

  subroutine setup_zonal_comm_for_reduce(num_total_lat, lat_ibeg, lat_iend)

    integer, intent(in) :: num_total_lat
    integer, intent(in) :: lat_ibeg
    integer, intent(in) :: lat_iend

    integer ierr, i, j, jr
    integer, allocatable :: zonal_proc_id(:)

    ! Create zonal communicator for reduce algorithm.
    if (proc%idom == 1) then ! Only root domain has reduce region.
      jr = 0
      do j = 1, size(reduce_factors)
        if (reduce_factors(j) > 0) then
          jr = j
        else if (jr /= 0) then
          exit
        end if
      end do
      allocate(zonal_proc_id(proc%cart_dims(1)))
      if (global_mesh%is_south_pole(lat_ibeg) .or. global_mesh%is_north_pole(lat_iend) .or. &
          lat_ibeg <= jr .or. lat_iend > num_total_lat - jr) then
        call log_notice('Create zonal communicator on process ' // to_string(proc%id) // '.')
        do i = 1, proc%cart_dims(1)
          call MPI_CART_RANK(proc%comm, [i-1,proc%cart_coords(2)], zonal_proc_id(i), ierr)
        end do
        call MPI_GROUP_INCL(proc%group, size(zonal_proc_id), zonal_proc_id, proc%zonal_group, ierr)
        call MPI_COMM_CREATE_GROUP(proc%comm, proc%zonal_group, sum(zonal_proc_id), proc%zonal_comm, ierr)
      end if
      deallocate(zonal_proc_id)
    end if

  end subroutine setup_zonal_comm_for_reduce

  subroutine connect_parent()

    integer parent_idom, parent_cart_dims(2)
    integer lon_halo_width, lat_halo_width
    integer num_total_lon, num_total_lat
    integer parent_lon_ibeg, parent_lon_iend
    integer parent_lat_ibeg, parent_lat_iend
    integer i, j

    !  ____________________________________________________________________________
    ! |          |          |          |          |          |          |          |
    ! |          |          |          |          |          |          |          |
    ! |          |          |          |          |          |          |          |
    ! |          |          |          |          |          |          |          |
    ! |__________|__________|__________|__________|__________|__________|__________|
    ! |          |      ____|\\\\\\\\\\|____      |          |          |          |
    ! |          |     |    |\|\\\\\\|\|    |     |          |          |          |
    ! |          |     |    |\|\\\\\\|\|    |     |          |          |          |
    ! |          |     |____|\|\\\\\\|\|____|     |          |          |          |
    ! |__________|_____|____|\|//////|\|____|_____|__________|__________|__________|
    ! |          |     |    |\|//////|\|    |     |          |          |          |
    ! |          |     |____|\|//////|\|____|     |          |          |          |
    ! |          |          |\\\\\\\\\\|          |          |          |          |
    ! |          |          |\\\\\\\\\\|          |          |          |          |
    ! |__________|__________|\\\\\\\\\\|__________|__________|__________|__________|
    !
    ! //////                \\\\\\
    ! ////// - child domain \\\\\\ - parent domain
    ! //////                \\\\\\

    if (proc%idom > 1) then ! Only child processes needs to care about.
      parent_idom = proc%idom - 1
      parent_cart_dims = [num_proc_lon(parent_idom),num_proc_lat(parent_idom)]
      lon_halo_width = 1
      lat_halo_width = 1
      do j = 0, num_proc_lat(parent_idom) - 1
        do i = 0, num_proc_lon(parent_idom) - 1
          call decompose_domains(parent_idom, parent_cart_dims, [i,j], num_total_lon, num_total_lat, &
                                 parent_lon_ibeg, parent_lon_iend, parent_lat_ibeg, parent_lat_iend)
          if (((proc%lon_ibeg - lon_halo_width <= parent_lon_ibeg .and. parent_lon_ibeg <= proc%lon_iend + lon_halo_width)   .or. &
               (proc%lon_ibeg - lon_halo_width <= parent_lon_iend .and. parent_lon_iend <= proc%lon_iend + lon_halo_width)) .and. &
              ((proc%lat_ibeg - lat_halo_width <= parent_lat_ibeg .and. parent_lat_ibeg <= proc%lat_iend + lat_halo_width)   .or. &
               (proc%lat_ibeg - lat_halo_width <= parent_lat_iend .and. parent_lat_iend <= proc%lat_iend + lat_halo_width))) then
            print *, proc%id, proc%cart_coords, i, j, parent_lon_ibeg, parent_lon_iend, proc%lon_ibeg - lon_halo_width, proc%lon_iend + lon_halo_width, &
                                                      parent_lat_ibeg, parent_lat_iend, proc%lat_ibeg - lat_halo_width, proc%lat_iend + lat_halo_width
          end if
        end do
      end do
    end if

  end subroutine connect_parent

end module process_mod
