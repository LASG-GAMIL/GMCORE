module smag_damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public smag_damp_run

contains

  subroutine smag_damp_run(block, dt, tend, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(tend_type), intent(inout) :: tend
    type(state_type), intent(inout) :: state

    integer i, j, k
    real(r8) ls2

    associate (mesh      => block%mesh    , &
               smag_t    => state%smag_t  , & ! working array
               smag_s    => state%smag_s  , & ! working array
               kmh_lon   => state%kmh_lon , & ! working array
               kmh_lat   => state%kmh_lat , & ! working array
               dudt      => tend%smag_dudt, & ! working array
               dvdt      => tend%smag_dvdt, & ! working array
               u         => state%u_lon   , & ! inout
               v         => state%v_lat   , & ! inout
               pt        => state%pt      )   ! not used yet
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          smag_t(i,j,k) = (                       &
            u(i,j,k) - u(i-1,j,k)                 &
          ) / mesh%de_lon(j) - (                  &
            v(i,j  ,k) * mesh%half_cos_lat(j  ) - &
            v(i,j-1,k) * mesh%half_cos_lat(j-1)   &
          ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
        end do
      end do
    end do
    call fill_halo(block, smag_t, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false., south_halo=.false.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          smag_s(i,j,k) = (                       &
            v(i+1,j,k) - v(i,j,k)                 &
          ) / mesh%le_lat(j) + (                  &
            u(i,j+1,k) * mesh%full_cos_lat(j+1) - &
            u(i,j  ,k) * mesh%full_cos_lat(j  )   &
          ) / mesh%de_lat(j) / mesh%half_cos_lat(j)
        end do
      end do
    end do
    call fill_halo(block, smag_s, full_lon=.false., full_lat=.false., full_lev=.true., east_halo=.false., north_halo=.false.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        ls2 = smag_damp_coef / (1 / mesh%de_lon(j)**2 + 1 / mesh%le_lon(j)**2)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          kmh_lon(i,j,k) = ls2 * sqrt(                           &
            0.5_r8 * (smag_t(i,j,k)**2 + smag_t(i+1,j  ,k)**2) + &
            0.5_r8 * (smag_s(i,j,k)**2 + smag_s(i  ,j-1,k)**2)   &
          )
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        ls2 = smag_damp_coef / (1 / mesh%le_lat(j)**2 + 1 / mesh%de_lat(j)**2)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          kmh_lat(i,j,k) = ls2 * sqrt(                           &
            0.5_r8 * (smag_t(i,j,k)**2 + smag_t(i  ,j+1,k)**2) + &
            0.5_r8 * (smag_s(i,j,k)**2 + smag_s(i-1,j  ,k)**2)   &
          )
        end do
      end do
    end do

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          dudt(i,j,k) = kmh_lon(i,j,k) * (                                           &
            (u(i-1,j,k) - 2 * u(i,j,k) + u(i+1,j,k)) / mesh%de_lon(j)**2 +           &
            ((u(i,j+1,k) - u(i,j  ,k)) / mesh%de_lat(j  ) * mesh%half_cos_lat(j  ) - &
             (u(i,j  ,k) - u(i,j-1,k)) / mesh%de_lat(j-1) * mesh%half_cos_lat(j-1)   &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)                                &
          )
        end do
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          u(i,j,k) = u(i,j,k) + dt * dudt(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        if (j == global_mesh%half_lat_ibeg .or. j == global_mesh%half_lat_iend) then
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dvdt(i,j,k) = kmh_lat(i,j,k) * (                                           &
              (v(i-1,j,k) - 2 * v(i,j,k) + v(i+1,j,k)) / mesh%le_lat(j)**2             &
            )
          end do
        else
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dvdt(i,j,k) = kmh_lat(i,j,k) * (                                           &
              (v(i-1,j,k) - 2 * v(i,j,k) + v(i+1,j,k)) / mesh%le_lat(j)**2 +           &
              ((v(i,j+1,k) - v(i,j  ,k)) / mesh%le_lon(j+1) * mesh%full_cos_lat(j+1) - &
               (v(i,j  ,k) - v(i,j-1,k)) / mesh%le_lon(j  ) * mesh%full_cos_lat(j  )   &
              ) / mesh%de_lat(j) / mesh%half_cos_lat(j)                                &
            )
          end do
        end if
      end do
    end do
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          v(i,j,k) = v(i,j,k) + dt * dvdt(i,j,k)
        end do
      end do
    end do
    call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)
    end associate

  end subroutine smag_damp_run

end module smag_damp_mod
