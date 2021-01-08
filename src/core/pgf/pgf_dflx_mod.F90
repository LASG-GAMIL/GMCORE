module pgf_dflx_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use interp_mod
  use reduce_mod
  use debug_mod

  implicit none

contains

  subroutine pgf_dflx_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    call interp_lev_edge_to_cell(block%mesh, state%gz_lev, state%gz)
    call interp_lev_edge_to_lev_lon_edge(block%mesh, state%gz_lev, state%gz_lev_lon)
    call interp_lev_edge_to_lev_lat_edge(block%mesh, state%gz_lev, state%gz_lev_lat)
    call fill_halo(block, state%gz, full_lon=.true., full_lat=.true., full_lev=.true.)
    call fill_halo(block, state%gz_lev_lon, full_lon=.false., full_lat=.true. , full_lev=.false., south_halo=.false., north_halo=.false.)
    call fill_halo(block, state%gz_lev_lat, full_lon=.true. , full_lat=.false., full_lev=.false., west_halo=.false., east_halo=.false.)

  end subroutine pgf_dflx_prepare

  subroutine pgf_dflx_run(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j, k, move

    !
    ! pgf_lon:
    !          _                    _
    !      1  |  ∂(ϕ𝞭p)         ∂p   |
    !     --- | -------- - 𝞭(ϕ ----) |
    !     𝞭π  |_   ∂x           ∂x  _|
    !
    ! Same for pgf_lat.
    !
    associate (mesh          => block%mesh         , &
               reduced_mesh  => block%reduced_mesh , &
               reduced_state => block%reduced_state, &
               reduced_tend  => block%reduced_tend , &
               gz            => state%gz           , &
               gz_lev_lon    => state%gz_lev_lon   , &
               gz_lev_lat    => state%gz_lev_lat   , &
               p_lev         => state%p_lev        , &
               m_lon         => state%m_lon        , &
               m_lat         => state%m_lat        , &
               pgf_lon       => tend%pgf_lon       , &
               pgf_lat       => tend%pgf_lat)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        !
        !      4        3,4        3
        !    i,j,k             i+1,j,k
        !      o---------x---------o
        !      |                   |
        !      |                   |
        !      |                   |
        !      |                   |
        ! 1,4  x       i,j,k       x  2,3
        !      |                   |
        !      |                   |
        !      |                   |
        !      |                   |
        !      o---------x---------o
        !    i,j,k+1           i+1,j,k+1  --> east
        !      1        1,2       2
        !
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (reduced_mesh(j)%reduce_factor > 1) then
            pgf_lon(:,j,k) = 0.0_r8
            do move = 1, reduced_mesh(j)%reduce_factor
              do i = reduced_mesh(j)%half_lon_ibeg, reduced_mesh(j)%half_lon_iend
                reduced_tend(j)%pgf_lon(i,k) = ((               &
                  reduced_state(j)%gz(k,i+1,0,move) * (         &
                    reduced_state(j)%p_lev(k+1,i+1,0,move) -    &
                    reduced_state(j)%p_lev(k  ,i+1,0,move)      &
                  ) -                                           &
                  reduced_state(j)%gz(k,i  ,0,move) * (         &
                    reduced_state(j)%p_lev(k+1,i  ,0,move) -    &
                    reduced_state(j)%p_lev(k  ,i  ,0,move)      &
                  )                                             &
                ) - (                                           &
                  reduced_state(j)%gz_lev_lon(k+1,i,0,move) * ( &
                    reduced_state(j)%p_lev(k+1,i+1,0,move) -    &
                    reduced_state(j)%p_lev(k+1,i  ,0,move)      &
                  ) -                                           &
                  reduced_state(j)%gz_lev_lon(k  ,i,0,move) * ( &
                    reduced_state(j)%p_lev(k  ,i+1,0,move) -    &
                    reduced_state(j)%p_lev(k  ,i  ,0,move)      &
                  )                                             &
                )) / reduced_mesh(j)%de_lon(0)                  &
                   / reduced_state(j)%m_lon(k,i,0,move)
              end do
              call reduce_append_array(move, reduced_mesh(j), reduced_tend(j)%pgf_lon(:,k), mesh, pgf_lon(:,j,k))
            end do
            call overlay_inner_halo(block, pgf_lon(:,j,k), west_halo=.true.)
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              pgf_lon(i,j,k) = ((                                           &
                gz(i+1,j,k) * (p_lev(i+1,j,k+1) - p_lev(i+1,j,k)) -         &
                gz(i  ,j,k) * (p_lev(i  ,j,k+1) - p_lev(i  ,j,k))           &
              ) - (                                                         &
                gz_lev_lon(i,j,k+1) * (p_lev(i+1,j,k+1) - p_lev(i,j,k+1)) - &
                gz_lev_lon(i,j,k  ) * (p_lev(i+1,j,k  ) - p_lev(i,j,k  ))   &
              )) / mesh%de_lon(j) / m_lon(i,j,k)
            end do
          end if
        end do

        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            pgf_lat(i,j,k) = ((                                           &
              gz(i,j  ,k) * (p_lev(i,j  ,k+1) - p_lev(i,j  ,k)) -         &
              gz(i,j-1,k) * (p_lev(i,j-1,k+1) - p_lev(i,j-1,k))           &
            ) - (                                                         &
              gz_lev_lat(i,j,k+1) * (p_lev(i,j,k+1) - p_lev(i,j-1,k+1)) - &
              gz_lev_lat(i,j,k  ) * (p_lev(i,j,k  ) - p_lev(i,j-1,k  ))   &
            )) / mesh%de_lat(j) / m_lat(i,j,k)
#else
            pgf_lat(i,j,k) = ((                                           &
              gz(i,j+1,k) * (p_lev(i,j+1,k+1) - p_lev(i,j+1,k)) -         &
              gz(i,j  ,k) * (p_lev(i,j  ,k+1) - p_lev(i,j  ,k))           &
            ) - (                                                         &
              gz_lev_lat(i,j,k+1) * (p_lev(i,j+1,k+1) - p_lev(i,j,k+1)) - &
              gz_lev_lat(i,j,k  ) * (p_lev(i,j+1,k  ) - p_lev(i,j,k  ))   &
            )) / mesh%de_lat(j) / m_lat(i,j,k)
#endif
          end do
        end do
      end do
    end associate

  end subroutine pgf_dflx_run

end module pgf_dflx_mod
