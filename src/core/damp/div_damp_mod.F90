module div_damp_mod

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

  public div_damp_init
  public div_damp_final
  public div_damp_run

  real(r8), allocatable :: c_lon(:,:)
  real(r8), allocatable :: c_lat(:,:)

  logical, allocatable :: use_implicit_solver(:)
  real(r8), parameter :: beta = 0.5_r8
  type(tridiag_solver_type), allocatable :: zonal_solver(:,:)

contains

  subroutine div_damp_init(blocks)

    type(block_type), intent(in) :: blocks(:)

    integer j, k, r, j0, jr
    integer iblk, js, je
    real(r8) a, b

    if (.not. use_div_damp) return

    call div_damp_final()

    j0 = 0
    do j = global_mesh%full_lat_ibeg, global_mesh%full_lat_iend
      if (global_mesh%full_lat(j) < 0 .and. -global_mesh%full_lat_deg(j) >= div_damp_lat0) j0 = j
    end do

    allocate(c_lon(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(c_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    select case (div_damp_order)
    case (2)
      r = 1

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          jr = merge(j - global_mesh%full_lat_ibeg_no_pole + 1, global_mesh%full_lat_iend_no_pole - j + 1, global_mesh%full_lat(j) < 0)
          if (baroclinic) then
            if (j0 == 0) then
              c_lon(j,k) = div_damp_coef2 * &
                (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
                global_mesh%full_cos_lat(j)**r * &
                radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
            else
              c_lon(j,k) = div_damp_coef2 * &
                (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
                (global_mesh%full_cos_lat(j)**r + div_damp_pole * exp(jr**2 * log(0.01_r8) / j0**2)) * &
                radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
            end if
          else
            c_lon(j,k) = div_damp_coef2 * &
              global_mesh%full_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
          end if
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg, global_mesh%half_lat_iend
          jr = merge(j - global_mesh%half_lat_ibeg + 1, global_mesh%half_lat_iend - j + 1, global_mesh%half_lat(j) < 0)
          if (baroclinic) then
            if (j0 == 0) then
              c_lat(j,k) = div_damp_coef2 * &
                (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
                global_mesh%half_cos_lat(j)**r * &
                radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
            else
              c_lat(j,k) = div_damp_coef2 * &
                (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
                (global_mesh%half_cos_lat(j)**r + div_damp_pole * exp(jr**2 * log(0.01_r8) / j0**2)) * &
                radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
            end if
          else
            c_lat(j,k) = div_damp_coef2 * &
              global_mesh%half_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_dyn
          end if
        end do
      end do
    case default
      call log_error('Unsupported div_damp_order ' // trim(to_str(div_damp_order)) // '!')
    end select

    ! Initialize cyclic tridiagonal solvers on each zonal circles if need implicit integration.
    js =  1e8
    je = -1e8
    do iblk = 1, size(proc%blocks)
      js = min(proc%blocks(iblk)%mesh%full_lat_ibeg_no_pole, js)
      je = max(proc%blocks(iblk)%mesh%full_lat_iend_no_pole, je)
    end do

    allocate(use_implicit_solver(js:je))
    use_implicit_solver = .false.
    allocate(zonal_solver(js:je,global_mesh%num_full_lev))
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = js, je
        if (abs(global_mesh%full_lat_deg(j)) > div_damp_imp_lat0) then
          use_implicit_solver(j) = .true.
          if (k > 1) then
            if (c_lon(j,k) == c_lon(j,k-1)) then
              call zonal_solver(j,k)%clone(zonal_solver(j,k-1))
              cycle
            end if
          end if
          b = -c_lon(j,k) * (1 - beta) * dt_dyn / global_mesh%de_lon(j)**2
          a = 2 * (-b) + 1
          call zonal_solver(j,k)%init_sym_const(blocks(1)%mesh%num_half_lon, a, b, zonal_tridiag_solver)
        end if
      end do
    end do

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(c_lon)) deallocate(c_lon)
    if (allocated(c_lat)) deallocate(c_lat)

    if (allocated(use_implicit_solver)) deallocate(use_implicit_solver)
    if (allocated(zonal_solver)) deallocate(zonal_solver)

  end subroutine div_damp_final

  subroutine div_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    real(r8) rhs(block%mesh%half_lon_ibeg:block%mesh%half_lon_iend)
    integer i, j, k

    call calc_div(block, state)

    associate (mesh => block%mesh , &
               div  => state%div  , &
               u    => state%u_lon, &
               v    => state%v_lat)
    select case (div_damp_order)
    case (2)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg, mesh%half_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            v(i,j,k) = v(i,j,k) + dt * c_lat(j,k) * ( &
              div(i,j+1,k) - div(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (use_implicit_solver(j)) then
            ! Set right hand side.
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              rhs(i) = u(i,j,k) +                                        &
                       c_lon(j,k) * beta * dt / mesh%de_lon(j) * (       &
                         div(i+1,j,k) - div(i,j,k)                       &
                       ) +                                               &
                       c_lon(j,k) * (1 - beta) * dt / mesh%de_lon(j) * ( &
                         v(i+1,j,k) - v(i+1,j-1,k) -                     &
                         v(i  ,j,k) + v(i  ,j-1,k)                       &
                       ) / mesh%le_lon(j)
            end do
            call zonal_solver(j,k)%solve(rhs, u(mesh%half_lon_ibeg:mesh%half_lon_iend,j,k))
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              u(i,j,k) = u(i,j,k) + dt * c_lon(j,k) * ( &
                div(i+1,j,k) - div(i,j,k)) / mesh%de_lon(j)
            end do
          end if
        end do
      end do
      call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    case (4)
    end select
    end associate

  end subroutine div_damp_run

end module div_damp_mod
