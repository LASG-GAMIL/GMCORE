module vor_damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use process_mod
  use block_mod
#ifdef USE_HZD
  use tridiag_hzd_mod
#else
  use tridiag_mkl_mod
#endif

  implicit none

  private

  public vor_damp_init
  public vor_damp_final
  public vor_damp_run

  real(r8), allocatable :: cv_full_lat(:,:)
  real(r8), allocatable :: cv_half_lat(:,:)

  logical, allocatable :: use_implicit_solver(:)
  real(r8), parameter :: beta = 0.5_r8
#ifdef USE_HZD
  type(hzd_tridiag_solver_type), allocatable :: zonal_solver(:,:)
#else
  type(mkl_tridiag_solver_type), allocatable :: zonal_solver(:,:)
#endif

contains

  subroutine vor_damp_init()

    integer j, j0, jr, k
    integer iblk, js, je
    real(r8) a, b

    call vor_damp_final()

    ! Only do vorticity damping in reduced regions.
    ! First, find the interface when reduce starts.
    j0 = 0
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat_deg(j) >= -vor_damp_lat0) then
        j0 = j
        exit
      end if
    end do
    if (is_root_proc()) then
      call log_notice('Vorticity damping control latitude index is ' // to_str(j0) // '.')
    end if

    allocate(cv_full_lat(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(cv_half_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    select case (vor_damp_order)
    case (2)
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          if (global_mesh%full_lat(j) <= 0) then
            jr = j - global_mesh%full_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%full_lat_iend_no_pole - j + 1
          end if
          cv_full_lat(j,k) = (vor_damp_coef2 * exp(jr**2 * log(vor_damp_decay) / j0**2) + vor_damp_bkg_coef) * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          if (global_mesh%half_lat(j) <= 0) then
            jr = j - global_mesh%half_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%half_lat_iend_no_pole - j + 1
          end if
          cv_half_lat(j,k) = (vor_damp_coef2 * exp(jr**2 * log(vor_damp_decay) / j0**2) + vor_damp_bkg_coef) * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do
    case default
      call log_error('Unsupported vor_damp_order ' // trim(to_str(vor_damp_order)) // '!')
    end select

    ! Initialize cyclic tridiagonal solvers on each zonal circles if need implicit integration.
    js =  1e8
    je = -1e8
    do iblk = 1, size(proc%blocks)
      js = min(proc%blocks(iblk)%mesh%half_lat_ibeg_no_pole, js)
      je = max(proc%blocks(iblk)%mesh%half_lat_iend_no_pole, je)
    end do
    allocate(use_implicit_solver(js:je))
    use_implicit_solver = .false.
    allocate(zonal_solver(js:je,global_mesh%num_full_lev))
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = js, je
        if (abs(global_mesh%half_lat_deg(j)) > vor_damp_imp_lat0) then
          use_implicit_solver(j) = .true.
          if (k > 1) then
            if (cv_half_lat(j,k) == cv_half_lat(j,k-1)) then
              call zonal_solver(j,k)%clone(zonal_solver(j,k-1))
              cycle
            end if
          end if
          b = -cv_half_lat(j,k) * (1 - beta) * dt_in_seconds / global_mesh%le_lat(j)**2
          a = 2 * (-b) + 1
          call zonal_solver(j,k)%init_sym_const(global_mesh%num_full_lon / num_proc_lon(1), a, b)
        end if
      end do
    end do

  end subroutine vor_damp_init

  subroutine vor_damp_final()

    if (allocated(cv_full_lat)) deallocate(cv_full_lat)
    if (allocated(cv_half_lat)) deallocate(cv_half_lat)

    if (allocated(use_implicit_solver)) deallocate(use_implicit_solver)
    if (allocated(zonal_solver)) deallocate(zonal_solver)

  end subroutine vor_damp_final

  subroutine vor_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) rhs(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)
    integer i, j, k

    mesh => state%mesh

    select case (vor_damp_order)
    case (2)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            state%u(i,j,k) = state%u(i,j,k) - dt * cv_full_lat(j,k) * ( &
              state%vor(i,j,k) - state%vor(i,j-1,k)) / mesh%le_lon(j)
          end do
        end do
      end do
      call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          if (use_implicit_solver(j)) then
            ! Set right hand side.
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              rhs(i) = state%v(i,j,k) + &
                       cv_half_lat(j,k) * beta * dt / mesh%le_lat(j) * (       &
                         state%vor(i,j,k) - state%vor(i-1,j,k)                 &
                       ) +                                                     &
                       cv_half_lat(j,k) * (1 - beta) * dt / mesh%le_lat(j) * ( &
                         state%u(i-1,j+1,k) - state%u(i-1,j,k) -               &
                         state%u(i  ,j+1,k) + state%u(i  ,j,k)                 &
                       ) / mesh%de_lat(j)
            end do
            call zonal_solver(j,k)%solve(rhs, state%v(mesh%full_lon_ibeg:mesh%full_lon_iend,j,k))
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              state%v(i,j,k) = state%v(i,j,k) + dt * cv_half_lat(j,k) * ( &
                state%vor(i,j,k) - state%vor(i-1,j,k)) / mesh%le_lat(j)
            end do
          end if
        end do
      end do
      call fill_halo(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
    case (4)
    end select

  end subroutine vor_damp_run

end module vor_damp_mod
