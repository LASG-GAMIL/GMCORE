module vor_damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use process_mod
  use block_mod
  use tridiag_mod
  use operators_mod

  implicit none

  private

  public vor_damp_init
  public vor_damp_final
  public vor_damp_run

  real(r8), allocatable :: c_lon(:,:)
  real(r8), allocatable :: c_lat(:,:)

  logical, allocatable :: use_implicit_solver(:)
  real(r8), parameter :: beta = 0.5_r8
  type(tridiag_solver_type), allocatable :: zonal_solver(:,:)

contains

  subroutine vor_damp_init(blocks)

    type(block_type), intent(in) :: blocks(:)

    integer j, j0, jr, k
    integer iblk, js, je
    real(r8) a, b

    if (.not. use_vor_damp) return

    call vor_damp_final()

    j0 = 0
    do j = global_mesh%full_lat_ibeg, global_mesh%full_lat_iend
      if (global_mesh%full_lat(j) < 0 .and. -global_mesh%full_lat_deg(j) >= vor_damp_lat0) j0 = j
    end do

    allocate(c_lon(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(c_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    select case (vor_damp_order)
    case (2)
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          jr = merge(j - global_mesh%full_lat_ibeg_no_pole + 1, global_mesh%full_lat_iend_no_pole - j + 1, global_mesh%full_lat(j) < 0)
          if (j0 == 0) then
            c_lon(j,k) = vor_damp_coef2 * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
          else
            c_lon(j,k) = vor_damp_coef2 * &
              exp(jr**2 * log(1.0e-10) / j0**2) * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
          end if
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg, global_mesh%half_lat_iend
          jr = merge(j - global_mesh%half_lat_ibeg + 1, global_mesh%half_lat_iend - j + 1, global_mesh%half_lat(j) < 0)
          if (j0 == 0) then
            c_lat(j,k) = vor_damp_coef2 * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
          else
            c_lat(j,k) = vor_damp_coef2 * &
              exp(jr**2 * log(1.0e-10) / j0**2) * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
          end if
        end do
      end do
    case default
      call log_error('Unsupported vor_damp_order ' // trim(to_str(vor_damp_order)) // '!')
    end select

    ! Initialize cyclic tridiagonal solvers on each zonal circles if need implicit integration.
    js =  1e8
    je = -1e8
    do iblk = 1, size(blocks)
      js = min(blocks(iblk)%mesh%half_lat_ibeg, js)
      je = max(blocks(iblk)%mesh%half_lat_iend, je)
    end do

    allocate(use_implicit_solver(js:je))
    use_implicit_solver = .false.
    allocate(zonal_solver(js:je,global_mesh%num_full_lev))
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = js, je
        if (abs(global_mesh%half_lat_deg(j)) > vor_damp_imp_lat0) then
          use_implicit_solver(j) = .true.
          if (k > 1) then
            if (c_lat(j,k) == c_lat(j,k-1)) then
              call zonal_solver(j,k)%clone(zonal_solver(j,k-1))
              cycle
            end if
          end if
          b = -c_lat(j,k) * (1 - beta) * dt_dyn / global_mesh%le_lat(j)**2
          a = 2 * (-b) + 1
          call zonal_solver(j,k)%init_sym_const(blocks(1)%mesh%num_full_lon, a, b, zonal_tridiag_solver)
        end if
      end do
    end do

  end subroutine vor_damp_init

  subroutine vor_damp_final()

    if (allocated(c_lon)) deallocate(c_lon)
    if (allocated(c_lat)) deallocate(c_lat)

    if (allocated(use_implicit_solver)) deallocate(use_implicit_solver)
    if (allocated(zonal_solver)) deallocate(zonal_solver)

  end subroutine vor_damp_final

  subroutine vor_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    real(r8) rhs(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend)
    integer i, j, k

    associate (mesh => block%mesh , &
               vor  => state%vor  , &
               u    => state%u_lon, &
               v    => state%v_lat)
    call calc_vor(block, state, u, v)
    select case (vor_damp_order)
    case (2)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            u(i,j,k) = u(i,j,k) - dt * c_lon(j,k) * (vor(i,j,k) - vor(i,j-1,k)) / mesh%le_lon(j)
          end do
        end do
      end do
      call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          if (use_implicit_solver(j)) then
            ! Set right hand side.
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              rhs(i) = v(i,j,k) +                                 &
                c_lat(j,k) * beta * dt / mesh%le_lat(j) * (       &
                  vor(i,j,k) - vor(i-1,j,k)                       &
                ) +                                               &
                c_lat(j,k) * (1 - beta) * dt / mesh%le_lat(j) * ( &
                  u(i-1,j+1,k) - u(i-1,j,k) -                     &
                  u(i  ,j+1,k) + u(i  ,j,k)                       &
                ) / mesh%de_lat(j)
            end do
            call zonal_solver(j,k)%solve(rhs, v(mesh%full_lon_ibeg:mesh%full_lon_iend,j,k))
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              v(i,j,k) = v(i,j,k) + dt * c_lat(j,k) * ( &
                vor(i,j,k) - vor(i-1,j,k)) / mesh%le_lat(j)
            end do
          end if
        end do
      end do
      call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    case (4)
    end select
    end associate

  end subroutine vor_damp_run

end module vor_damp_mod
