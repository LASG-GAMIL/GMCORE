module smag_damp_mod

  use flogger
  use string
  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public smag_damp_init
  public smag_damp_run
  public smag_damp_final

  real(r8), allocatable, dimension(:), target :: decay_from_top

contains

  subroutine smag_damp_init()

    integer k, k0

    call smag_damp_final()

    allocate(decay_from_top(global_mesh%num_full_lev))

    k0 = 8
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      decay_from_top(k) = exp((k - 1)**2 * log(0.01d0) / k0**2) + 1
    end do
    ! FIXME: Disable the decay for the time being.
    decay_from_top = 1

  end subroutine smag_damp_init

  subroutine smag_damp_final()

    if (allocated(decay_from_top)) deallocate(decay_from_top)

  end subroutine smag_damp_final

  subroutine smag_damp_run(block, dt, tend, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(tend_type), intent(inout) :: tend
    type(state_type), intent(inout) :: state

    integer i, j, k
    real(r8) ls2

    associate (mesh      => block%mesh     , &
               smag_t    => state%smag_t   , & ! working array
               smag_s    => state%smag_s   , & ! working array
               kmh_lon   => state%kmh_lon  , & ! working array
               kmh_lat   => state%kmh_lat  , & ! working array
               kmh       => state%kmh      , & ! working array
               dudt      => tend%smag_dudt , & ! working array
               dvdt      => tend%smag_dvdt , & ! working array
               dptdt     => tend%smag_dptdt, & ! working array
               u         => state%u_lon    , & ! inout
               v         => state%v_lat    , & ! inout
               pt        => state%pt       )   ! inout
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
        ls2 = smag_damp_coef / (1 / mesh%de_lon(j)**2 + 1 / mesh%le_lon(j)**2) * decay_from_top(k)
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
        ls2 = smag_damp_coef / (1 / mesh%le_lat(j)**2 + 1 / mesh%de_lat(j)**2) * decay_from_top(k)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          kmh_lat(i,j,k) = ls2 * sqrt(                           &
            0.5_r8 * (smag_t(i,j,k)**2 + smag_t(i  ,j+1,k)**2) + &
            0.5_r8 * (smag_s(i,j,k)**2 + smag_s(i-1,j  ,k)**2)   &
          )
        end do
      end do
    end do

    ! do k = mesh%full_lev_ibeg, mesh%full_lev_iend
    !   do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
    !     ls2 = smag_damp_coef / (1 / mesh%de_lon(j)**2 + 1 / mesh%le_lon(j)**2) * decay_from_top(k)
    !     do i = mesh%full_lon_ibeg, mesh%full_lon_iend
    !       kmh(i,j,k) = ls2 * sqrt(                        &
    !         smag_t(i,j,k)**2 + 0.25_r8 * (                &
    !           smag_s(i-1,j-1,k)**2 + smag_s(i-1,j,k)**2 + &
    !           smag_s(i  ,j-1,k)**2 + smag_s(i  ,j,k)**2   &
    !         )                                             &
    !       )
    !     end do
    !   end do
    ! end do

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

    ! do k = mesh%full_lev_ibeg, mesh%full_lev_iend
    !   do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
    !     do i = mesh%full_lon_ibeg, mesh%full_lon_iend
    !       dptdt(i,j,k) = kmh(i,j,k) * (                                                &
    !         (pt(i-1,j,k) - 2 * pt(i,j,k) + pt(i+1,j,k)) / mesh%de_lon(j)**2 +          &
    !         ((pt(i,j+1,k) - pt(i,j  ,k)) / mesh%de_lat(j  ) * mesh%half_cos_lat(j  ) - &
    !          (pt(i,j  ,k) - pt(i,j-1,k)) / mesh%de_lat(j-1) * mesh%half_cos_lat(j-1)   &
    !         ) / mesh%le_lon(j) / mesh%full_cos_lat(j)                                  &
    !       )
    !     end do
    !   end do
    ! end do
    ! do k = mesh%full_lev_ibeg, mesh%full_lev_iend
    !   do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
    !     do i = mesh%full_lon_ibeg, mesh%full_lon_iend
    !       pt(i,j,k) = pt(i,j,k) + dt * dptdt(i,j,k)
    !     end do
    !   end do
    ! end do
    ! call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine smag_damp_run

end module smag_damp_mod
