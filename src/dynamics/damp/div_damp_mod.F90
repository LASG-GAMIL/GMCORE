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

contains

  subroutine div_damp_init(blocks)

    type(block_type), intent(in) :: blocks(:)

    integer j, k, r
    integer iblk

    if (.not. use_div_damp) return

    call div_damp_final()

    allocate(c_lon(global_mesh%num_full_lat,global_mesh%num_full_lev))
    allocate(c_lat(global_mesh%num_half_lat,global_mesh%num_full_lev))

    select case (div_damp_order)
    case (2)
      r = 1
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          if (baroclinic) then
            c_lon(j,k) = div_damp_coef2 * &
              (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
              global_mesh%full_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon
          else
            c_lon(j,k) = div_damp_coef2 * &
              global_mesh%full_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon
          end if
        end do
      end do
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg, global_mesh%half_lat_iend
          if (baroclinic) then
            c_lat(j,k) = div_damp_coef2 * &
              (1.0_r8 + div_damp_top * exp((k-1)**2 * log(0.2_r8) / (div_damp_k0-1)**2)) * &
              global_mesh%half_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon
          else
            c_lat(j,k) = div_damp_coef2 * &
              global_mesh%half_cos_lat(j)**r * &
              radius**2 * global_mesh%dlat(j) * global_mesh%dlon
          end if
        end do
      end do
    case default
      call log_error('Unsupported div_damp_order ' // trim(to_str(div_damp_order)) // '!')
    end select

  end subroutine div_damp_init

  subroutine div_damp_final()

    if (allocated(c_lon)) deallocate(c_lon)
    if (allocated(c_lat)) deallocate(c_lat)

  end subroutine div_damp_final

  subroutine div_damp_run(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

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
            v(i,j,k) = v(i,j,k) + c_lat(j,k) * (div(i,j+1,k) - div(i,j,k)) / mesh%de_lat(j)
          end do
        end do
      end do
      call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            u(i,j,k) = u(i,j,k) + c_lon(j,k) * (div(i+1,j,k) - div(i,j,k)) / mesh%de_lon(j)
          end do
        end do
      end do
      call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)
    case (4)
    end select
    end associate

  end subroutine div_damp_run

end module div_damp_mod
