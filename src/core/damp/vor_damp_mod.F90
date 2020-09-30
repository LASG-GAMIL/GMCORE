module vor_damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use pv_mod
  use reduce_mod

  implicit none

  private

  public vor_damp_init
  public vor_damp_final
  public vor_damp

  real(r8), allocatable :: cv_full_lat(:,:)
  real(r8), allocatable :: cv_half_lat(:,:)

contains

  subroutine vor_damp_init()

    integer j, jr0, jr, k

    call vor_damp_final()

    ! Only do vorticity damping in reduced regions.
    ! First, find the interface when reduce starts.
    jr0 = 0
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
        if (reduce_factors(jr) > 1) jr0 = jr
      end if
    end do

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
          cv_full_lat(j,k) = vor_damp_coef2 * exp(jr**2 * log(0.5_r8) / jr0**2) * &
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
          cv_half_lat(j,k) = vor_damp_coef2 * exp(jr**2 * log(0.5_r8) / jr0**2) * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do
    case default
      call log_error('Unsupported vor_damp_order ' // trim(to_string(vor_damp_order)) // '!')
    end select

  end subroutine vor_damp_init

  subroutine vor_damp_final()

    if (allocated(cv_full_lat)) deallocate(cv_full_lat)
    if (allocated(cv_half_lat)) deallocate(cv_half_lat)

  end subroutine vor_damp_final

  subroutine vor_damp(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    type(mesh_type), pointer :: mesh
    integer i, j, k, jr, jb, move

    mesh => state%mesh

    select case (vor_damp_order)
    case (2)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            tend%dvordlat(i,j,k) = cv_full_lat(j,k) * (state%vor(i,j+1,k) - state%vor(i,j,k)) / mesh%le_lon(j)
#else
            tend%dvordlat(i,j,k) = cv_full_lat(j,k) * (state%vor(i,j,k) - state%vor(i,j-1,k)) / mesh%le_lon(j)
#endif
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
#ifdef V_POLE
          if (block%reduced_mesh(j)%reduce_factor > block%reduced_mesh(j+1)%reduce_factor) then
            jr = j - 1; jb =  1
          else
            jr = j    ; jb =  0
          end if
#else
          if (block%reduced_mesh(j)%reduce_factor > block%reduced_mesh(j+1)%reduce_factor) then
            jr = j    ; jb =  0
          else
            jr = j + 1; jb = -1
          end if
#endif
          if (block%reduced_mesh(jr)%reduce_factor > 1) then
            tend%dvordlon(:,j,k) = 0.0_r8
            do move = 1, block%reduced_mesh(jr)%reduce_factor
              do i = block%reduced_mesh(jr)%full_lon_ibeg, block%reduced_mesh(jr)%full_lon_iend
                block%reduced_tend(jr)%dvordlon(i,k) = cv_half_lat(j,k) * ( &
                  block%reduced_state(jr)%vor(k,i  ,jb,move) -              &
                  block%reduced_state(jr)%vor(k,i-1,jb,move)                &
                ) / block%reduced_mesh(jr)%le_lat(jb)
              end do
              call reduce_append_array(move, block%reduced_mesh(jr), block%reduced_tend(jr)%dvordlon(:,k), mesh, tend%dvordlon(:,j,k))
            end do
            call overlay_inner_halo(block, tend%dvordlon(:,j,k), west_halo=.true.)
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              tend%dvordlon(i,j,k) = cv_half_lat(j,k) * (state%vor(i,j,k) - state%vor(i-1,j,k)) / mesh%le_lat(j)
            end do
          end if
        end do
      end do
    case (4)
    end select

  end subroutine vor_damp

end module vor_damp_mod
