module pgf_dflx_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use reduce_mod

  implicit none

contains

  subroutine pgf_dflx_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

  end subroutine pgf_dflx_prepare

  subroutine pgf_dflx_run(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j, k, move

    associate (mesh       => block%mesh      , &
               gz         => state%gz        , &
               gz_lev_lon => state%gz_lev_lon, &
               gz_lev_lat => state%gz_lev_lat, &
               p_lev      => state%p_lev     , &
               m          => state%m         , &
               m_lon      => state%m_lon     , &
               m_lat      => state%m_lat     , &
               pgf_lon    => tend%pgf_lon    , &
               pgf_lat    => tend%pgf_lat)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            pgf_lon(i,j,k) = -(                                                &
              gz(i  ,j,k) * (p_lev(i  ,j,k+1) - p_lev(i  ,j,k)) / m(i  ,j,k) - &
              gz(i-1,j,k) * (p_lev(i-1,j,k+1) - p_lev(i-1,j,k)) / m(i-1,j,k)   &
            ) / mesh%de_lon(j) + (                                             &
              gz_lev_lon(i,j,k+1) * (p_lev(i,j,k+1) - p_lev(i-1,j,k+1)) -      &
              gz_lev_lon(i,j,k  ) * (p_lev(i,j,k  ) - p_lev(i-1,j,k  ))        &
            ) / mesh%de_lon(j) / m_lon(i,j,k)
          end do
        end do

        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            pgf_lat(i,j,k) = -(                                                &
              gz(i,j  ,k) * (p_lev(i,j  ,k+1) - p_lev(i,j  ,k)) / m(i,j  ,k) - &
              gz(i,j-1,k) * (p_lev(i,j-1,k+1) - p_lev(i,j-1,k)) / m(i,j-1,k)   &
            ) / mesh%de_lat(j) + (                                             &
              gz_lev_lat(i,j,k+1) * (p_lev(i,j,k+1) - p_lev(i,j-1,k+1)) -      &
              gz_lev_lat(i,j,k  ) * (p_lev(i,j,k  ) - p_lev(i,j-1,k  ))        &
            ) / mesh%de_lat(j) / m_lat(i,j,k)
#else
            pgf_lat(i,j,k) = -(                                                &
              gz(i,j+1,k) * (p_lev(i,j+1,k+1) - p_lev(i,j+1,k)) / m(i,j+1,k) - &
              gz(i,j  ,k) * (p_lev(i,j  ,k+1) - p_lev(i,j  ,k)) / m(i,j  ,k)   &
            ) / mesh%de_lat(j) + (                                             &
              gz_lev_lat(i,j,k+1) * (p_lev(i,j+1,k+1) - p_lev(i,j,k+1)) -      &
              gz_lev_lat(i,j,k  ) * (p_lev(i,j+1,k  ) - p_lev(i,j,k  ))        &
            ) / mesh%de_lat(j) / m_lat(i,j,k)
#endif
          end do
        end do
      end do
    end associate

  end subroutine pgf_dflx_run

end module pgf_dflx_mod
