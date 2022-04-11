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
    real(r8) ls

    associate (mesh      => block%mesh     , &
               tension_h => state%tension_h, &
               shear_h   => state%shear_h  , &
               kmh       => state%kmh      , &
               kmh_lon   => state%kmh_lon  , &
               kmh_lat   => state%kmh_lat  , &
               kmh_vtx   => state%kmh_vtx  , &
               dpt       => tend%smag_dpt  , &
               u         => state%u_lon    , &
               v         => state%v_lat    , &
               pt        => state%pt)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            tension_h(i,j,k) = (                    &
              u(i,j,k) - u(i-1,j,k)                 &
            ) / mesh%de_lon(j) - (                  &
              v(i,j  ,k) * mesh%half_cos_lat(j  ) - &
              v(i,j-1,k) * mesh%half_cos_lat(j-1)   &
            ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
          end do
        end do
      end do
      call fill_halo(block, tension_h, full_lon=.true., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            shear_h(i,j,k) = (                      &
              v(i+1,j,k) - v(i,j,k)                 &
            ) / mesh%le_lat(j) + (                  &
              u(i,j+1,k) * mesh%full_cos_lat(j+1) - &
              u(i,j  ,k) * mesh%full_cos_lat(j  )   &
            ) / mesh%de_lat(j) / mesh%half_cos_lat(j)
          end do
        end do
      end do
      call fill_halo(block, shear_h, full_lon=.false., full_lat=.false., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          ls = radius * mesh%full_cos_lat(j) * mesh%dlon
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            kmh(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt( &
              tension_h(i,j,k)**2 +                             &
              0.25_r8 * (                                       &
                shear_h(i-1,j  ,k)**2 + shear_h(i,j  ,k)**2 +   &
                shear_h(i-1,j-1,k)**2 + shear_h(i,j-1,k)**2     &
              )                                                 &
            )
          end do
        end do
      end do
      call fill_halo(block, kmh, full_lon=.true., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          ls = radius * mesh%half_cos_lat(j) * mesh%dlon
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            kmh_vtx(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt( &
              shear_h(i,j,k)**2 +                                   &
              0.25_r8 * (                                           &
                tension_h(i,j  ,k)**2 + tension_h(i+1,j  ,k)**2 +   &
                tension_h(i,j+1,k)**2 + tension_h(i+1,j+1,k)**2     &
              )                                                     &
            )
          end do
        end do
      end do
      call fill_halo(block, kmh_vtx, full_lon=.false., full_lat=.false., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          ls = radius * mesh%full_cos_lat(j) * mesh%dlon
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            kmh_lon(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(      &
              0.5_r8 * (tension_h(i,j,k)**2 + tension_h(i+1,j  ,k)**2) + &
              0.5_r8 * (shear_h  (i,j,k)**2 + shear_h  (i  ,j-1,k)**2)   &
            )
          end do
        end do
      end do
      call fill_halo(block, kmh_lon, full_lon=.false., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          ls = radius * mesh%half_cos_lat(j) * mesh%dlon
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            kmh_lat(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(      &
              0.5_r8 * (tension_h(i,j,k)**2 + tension_h(i  ,j+1,k)**2) + &
              0.5_r8 * (  shear_h(i,j,k)**2 +   shear_h(i-1,j  ,k)**2)   &
            )
          end do
        end do
      end do
      call fill_halo(block, kmh_lat, full_lon=.true., full_lat=.false., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            u(i,j,k) = u(i,j,k) + dt * ((           &
              kmh(i+1,j,k) * tension_h(i+1,j,k) -   &
              kmh(i  ,j,k) * tension_h(i  ,j,k)     &
            ) / mesh%de_lon(j) + (                  &
              kmh_vtx(i,j  ,k) * shear_h(i,j  ,k) - &
              kmh_vtx(i,j-1,k) * shear_h(i,j-1,k)   &
            ) / mesh%le_lon(j))
          end do
        end do
      end do
      call fill_halo(block, u, full_lon=.false., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            v(i,j,k) = v(i,j,k) + dt * ((           &
              kmh_vtx(i  ,j,k) * shear_h(i  ,j,k) - &
              kmh_vtx(i-1,j,k) * shear_h(i-1,j,k)   &
            ) / mesh%le_lat(j) - (                  &
              kmh(i,j+1,k) * tension_h(i,j+1,k) -   &
              kmh(i,j  ,k) * tension_h(i,j  ,k)     &
            ) / mesh%de_lat(j))
          end do
        end do
      end do
      call fill_halo(block, v, full_lon=.true., full_lat=.false., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            dpt(i,j,k) = (                                     &
              kmh_lon(i  ,j,k) * (pt(i+1,j,k) - pt(i  ,j,k)) - &
              kmh_lon(i-1,j,k) * (pt(i  ,j,k) - pt(i-1,j,k))   &
            ) / mesh%de_lon(j)**2 + (                          &
              kmh_lat(i,j  ,k) * mesh%half_cos_lat(j  ) * (    &
                pt(i,j+1,k) - pt(i,j  ,k)                      &
              ) / mesh%de_lat(j  ) -                           &
              kmh_lat(i,j-1,k) * mesh%half_cos_lat(j-1) * (    &
                pt(i,j  ,k) - pt(i,j-1,k)                      &
              ) / mesh%de_lat(j-1)                             &
            ) / mesh%full_cos_lat(j) / mesh%le_lon(j)
          end do
        end do
      end do
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            pt(i,j,k) = pt(i,j,k) + dt * dpt(i,j,k)
          end do
        end do
      end do
      call fill_halo(block, pt, full_lon=.true., full_lat=.true., full_lev=.true.)
    end associate

  end subroutine smag_damp_run

end module smag_damp_mod
