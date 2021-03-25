module div_damp_mod

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

  public div_damp_init
  public div_damp_final
  public div_damp_run

  real(r8), allocatable :: cd_full_lat(:,:)
  real(r8), allocatable :: cd_half_lat(:,:)

  logical, allocatable :: use_implicit_solver(:)
  real(r8), parameter :: beta = 0.5_r8
#ifdef USE_HZD
  type(hzd_tridiag_solver_type), allocatable :: zonal_solver(:,:)
#else
  type(mkl_tridiag_solver_type), allocatable :: zonal_solver(:,:)
#endif

contains

  subroutine div_damp_init()

    integer j, k, r, jr, j0
    integer iblk, js, je
    real(r8) a, b

    call div_damp_final()

    j0 = 0
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
        if (jr > size(reduce_factors)) exit
        if (reduce_factors(jr) > 1) j0 = jr
      end if
    end do
    j0 = max(div_damp_j0, j0)

    allocate(cd_full_lat(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(cd_half_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    select case (div_damp_order)
    case (2)
      r = 1

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          if (global_mesh%full_lat(j) <= 0) then
            jr = j - global_mesh%full_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%full_lat_iend_no_pole - j + 1
          end if
          if (baroclinic) then
            cd_full_lat(j,k) = div_damp_coef2 * &
              (1.0_r8 + div_damp_upper * exp(k**2 * log(0.2_r8) / div_damp_k0**2)) * &
              (global_mesh%full_cos_lat(j)**r + div_damp_polar * exp(jr**2 * log(div_damp_exp) / j0**2)) * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          else
            cd_full_lat(j,k) = div_damp_coef2 * &
              global_mesh%full_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          end if
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          if (global_mesh%half_lat(j) <= 0) then
            jr = j - global_mesh%half_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%half_lat_iend_no_pole - j + 1
          end if
          if (baroclinic) then
            cd_half_lat(j,k) = div_damp_coef2 * &
              (1.0_r8 + div_damp_upper * exp(k**2 * log(0.2_r8) / div_damp_k0**2)) * &
              (global_mesh%half_cos_lat(j)**r + div_damp_polar * exp(jr**2 * log(div_damp_exp) / j0**2)) * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
          else
            cd_half_lat(j,k) = div_damp_coef2 * &
              global_mesh%half_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
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
            if (cd_full_lat(j,k) == cd_full_lat(j,k-1)) then
              call zonal_solver(j,k)%clone(zonal_solver(j,k-1))
              cycle
            end if
          end if
          b = -cd_full_lat(j,k) * (1 - beta) * dt_in_seconds / global_mesh%de_lon(j)**2
          a = 2 * (-b) + 1
          call zonal_solver(j,k)%init_sym_const(global_mesh%num_half_lon / num_proc_lon(1), a, b)
        end if
      end do
    end do

    ! j = 2
    ! k = 8
    ! if (proc%id == 0) then
    !   print *, "proc = 0"
    !   print *, "matrix"
    !   print *, zonal_solver(j,k)%A(1)%matrix
    !   print *, "w"
    !   print *, zonal_solver(j,k)%A(1)%w
    !   print *, "w_next"
    !   print *, zonal_solver(j,k)%A(1)%w_next
    !   print *, "v"
    !   print *, zonal_solver(j,k)%A(1)%v
    !   print *, "v_prev"
    !   print *, zonal_solver(j,k)%A(1)%v_prev
    ! end if
    ! call barrier()
    ! if (proc%id == 0 + num_proc_lat(1)) then 
    !   print *, "proc = ", 0 + num_proc_lat(1)
    !   print *, "matrix"
    !   print *, zonal_solver(j,k)%A(1)%matrix
    !   print *, "w"
    !   print *, zonal_solver(j,k)%A(1)%w
    !   print *, "w_next"
    !   print *, zonal_solver(j,k)%A(1)%w_next
    !   print *, "v"
    !   print *, zonal_solver(j,k)%A(1)%v
    !   print *, "v_prev"
    !   print *, zonal_solver(j,k)%A(1)%v_prev
    ! end if

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(cd_full_lat)) deallocate(cd_full_lat)
    if (allocated(cd_half_lat)) deallocate(cd_half_lat)

    if (allocated(use_implicit_solver)) deallocate(use_implicit_solver)
    if (allocated(zonal_solver)) deallocate(zonal_solver)

  end subroutine div_damp_final

  subroutine div_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    real(r8) rhs(block%mesh%half_lon_ibeg:block%mesh%half_lon_iend)
    real(r8) rhs_copy(block%mesh%half_lon_ibeg:block%mesh%half_lon_iend)
    integer i, j, k

    mesh => state%mesh

    select case (div_damp_order)
    case (2)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            state%v(i,j,k) = state%v(i,j,k) + dt * cd_half_lat(j,k) * ( &
              state%div(i,j+1,k) - state%div(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      call fill_halo(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (use_implicit_solver(j)) then
            ! Set right hand side.
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              rhs(i) = state%u(i,j,k) + &
                       cd_full_lat(j,k) * beta * dt / mesh%de_lon(j) * (       &
                         state%div(i+1,j,k) - state%div(i,j,k)                 &
                       ) +                                                     &
                       cd_full_lat(j,k) * (1 - beta) * dt / mesh%de_lon(j) * ( &
                         state%v(i+1,j,k) - state%v(i+1,j-1,k) -               &
                         state%v(i  ,j,k) + state%v(i  ,j-1,k)                 &
                       ) / mesh%le_lon(j)
            end do
            if (j == 2 .and. k == 8) rhs_copy = rhs
            call zonal_solver(j,k)%solve(rhs, state%u(mesh%half_lon_ibeg:mesh%half_lon_iend,j,k))
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              state%u(i,j,k) = state%u(i,j,k) + dt * cd_full_lat(j,k) * ( &
                state%div(i+1,j,k) - state%div(i,j,k)) / mesh%de_lon(j)
            end do
          end if
        end do
      end do
      call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)
    case (4)
    end select

    ! if (proc%id == 0) then
    !   print *, "proc = 0, j = 2, k = 8"
    !   print *, "rhs"
    !   print *, rhs_copy(mesh%half_lon_ibeg:mesh%half_lon_iend)
    !   print *, "x"
    !   print *, state%u(mesh%half_lon_ibeg:mesh%half_lon_iend,2,8)
    ! end if
    ! call barrier()
    ! if (proc%id == 0 + num_proc_lat(1)) then
    !   print *, "proc = ", 0+num_proc_lat(1), ", j = 2, k = 8"
    !   print *, "rhs"
    !   print *, rhs_copy(mesh%half_lon_ibeg:mesh%half_lon_iend)
    !   print *, "x"
    !   print *, state%u(mesh%half_lon_ibeg:mesh%half_lon_iend,2,8)
    ! end if

  end subroutine div_damp_run

end module div_damp_mod
