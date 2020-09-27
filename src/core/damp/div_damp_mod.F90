module div_damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use operators_mod

  implicit none

  private

  public div_damp_init
  public div_damp_final
  public div_damp

  real(r8), allocatable :: cd_full_lat(:,:)
  real(r8), allocatable :: cd_half_lat(:,:)

contains

  subroutine div_damp_init()

    integer j, k, r

    call div_damp_final()

    allocate(cd_full_lat(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(cd_half_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    select case (div_damp_order)
    case (2)
      r = 1

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cd_full_lat(j,k) = div_damp_coef2 * global_mesh%full_cos_lat(j)**r * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cd_half_lat(j,k) = div_damp_coef2 * global_mesh%half_cos_lat(j)**r * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do
    case (4)
      r = 2

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cd_full_lat(j,k) = div_damp_coef4 * global_mesh%full_cos_lat(j)**r * &
            radius**4 * global_mesh%dlat(j)**2 * global_mesh%dlon**2 / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cd_half_lat(j,k) = div_damp_coef4 * global_mesh%half_cos_lat(j)**r * &
            radius**4 * global_mesh%dlat(j)**2 * global_mesh%dlon**2 / dt_in_seconds
        end do
      end do
    case default
      call log_error('Unsupported div_damp_order ' // trim(to_string(div_damp_order)) // '!')
    end select

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(cd_full_lat)) deallocate(cd_full_lat)
    if (allocated(cd_half_lat)) deallocate(cd_half_lat)

  end subroutine div_damp_final

  subroutine div_damp(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k, cyc

    mesh => state%mesh

    select case (div_damp_order)
    case (2)
      cycle_loop: do cyc = 1, div_damp_cycles
        call calc_div(block, state)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              state%u(i,j,k) = state%u(i,j,k) + dt / div_damp_cycles * cd_full_lat(j,k) * ( &
                state%div(i+1,j,k) - state%div(i,j,k)) / mesh%de_lon(j)
            end do
          end do
        end do
        call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              state%v(i,j,k) = state%v(i,j,k) + dt / div_damp_cycles * cd_half_lat(j,k) * ( &
                state%div(i,j,k) - state%div(i,j-1,k)) / mesh%de_lat(j)
#else
              state%v(i,j,k) = state%v(i,j,k) + dt / div_damp_cycles * cd_half_lat(j,k) * ( &
                state%div(i,j+1,k) - state%div(i,j,k)) / mesh%de_lat(j)
#endif
            end do
          end do
        end do
        call fill_halo(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
      end do cycle_loop
    case (4)
    end select

  end subroutine div_damp

end module div_damp_mod
