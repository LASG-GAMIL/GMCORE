module pgf_sb81_mod

  use namelist_mod
  use parallel_mod
  use block_mod
  use interp_mod
  use reduce_mod

  implicit none

contains

  subroutine pgf_sb81_prepare(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k
    real(r8), parameter :: ln2 = log(2.0_r8)

    associate (mesh        => block%mesh       , &
               ph_lev      => state%ph_lev     , & ! in
               m           => state%m          , & ! in
               t           => state%t          , & ! in
               gz_lev      => state%gz_lev     , & ! in
               ak          => state%ak         , & ! out
               ak_t        => state%ak_t       , & ! out
               ak_t_lon    => state%ak_t_lon   , & ! out
               ak_t_lat    => state%ak_t_lat   , & ! out
               gz          => state%gz         , & ! out
               t_lnpop     => state%t_lnpop    , & ! out
               t_lnpop_lon => state%t_lnpop_lon, & ! out
               t_lnpop_lat => state%t_lnpop_lat)   ! out
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            if (k == 1) then
              ak(i,j,k) = ln2
            else
              ak(i,j,k) = 1.0_r8 - ph_lev(i,j,k) / m(i,j,k) * log(ph_lev(i,j,k+1) / ph_lev(i,j,k))
            end if
          end do
        end do
      end do
  
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            gz(i,j,k) = gz_lev(i,j,k+1) + ak(i,j,k) * Rd * t(i,j,k) ! Simmons and Burridge (1981)
          end do
        end do
      end do
      call fill_halo(block, gz, full_lon=.true., full_lat=.true., full_lev=.true.)

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg, mesh%full_lat_iend
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            t_lnpop(i,j,k) = t(i,j,k) * log(ph_lev(i,j,k+1) / ph_lev(i,j,k))
            ak_t(i,j,k) = ak(i,j,k) * t(i,j,k)
          end do
        end do
      end do
      call fill_halo(block, t_lnpop, full_lon=.true., full_lat=.true., full_lev=.true.)
      call fill_halo(block, ak_t   , full_lon=.true., full_lat=.true., full_lev=.true.)

      call interp_cell_to_lon_edge(mesh, t_lnpop, t_lnpop_lon)
      call interp_cell_to_lat_edge(mesh, t_lnpop, t_lnpop_lat)
      call interp_cell_to_lon_edge(mesh, ak_t, ak_t_lon)
      call interp_cell_to_lat_edge(mesh, ak_t, ak_t_lat)
    end associate

  end subroutine pgf_sb81_prepare

  subroutine pgf_sb81_run(block, state, tend)

    type(block_type), intent(inout) :: block
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    type(mesh_type), pointer :: mesh
    integer i, j, k, move

    if (baroclinic) then
      mesh => state%mesh

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (block%reduced_mesh(j)%reduce_factor > 1) then
            tend%pgf_lon(:,j,k) = 0.0_r8
            do move = 1, block%reduced_mesh(j)%reduce_factor
              do i = block%reduced_mesh(j)%half_lon_ibeg, block%reduced_mesh(j)%half_lon_iend
                block%reduced_tend(j)%pgf_lon(i,k) = Rd /               &
                  block%reduced_mesh(j)%de_lon(0) /                     &
                  block%reduced_state(j)%m_lon(k,i,0,move) * (          &
                    block%reduced_state(j)%t_lnpop_lon(k,i,0,move) * (  &
                      block%reduced_state(j)%ph_lev(k,i+1,0,move) -     &
                      block%reduced_state(j)%ph_lev(k,i  ,0,move)       &
                    ) + block%reduced_state(j)%ak_t_lon(k,i,0,move) * ( &
                      block%reduced_state(j)%m(k,i+1,0,move) -          &
                      block%reduced_state(j)%m(k,i  ,0,move)            &
                    )                                                   &
                )
              end do
              call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%pgf_lon(:,k), mesh, tend%pgf_lon(:,j,k))
            end do
            call overlay_inner_halo(block, tend%pgf_lon(:,j,k), west_halo=.true.)
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%pgf_lon(i,j,k) = Rd / mesh%de_lon(j) / state%m_lon(i,j,k) * (            &
                state%t_lnpop_lon(i,j,k) * (state%ph_lev(i+1,j,k) - state%ph_lev(i,j,k)) + &
                state%ak_t_lon(i,j,k) * (state%m(i+1,j,k) - state%m(i,j,k))                &
              )
            end do
          end if
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            tend%pgf_lat(i,j,k) = Rd / mesh%de_lat(j) / state%m_lat(i,j,k) * (           &
              state%t_lnpop_lat(i,j,k) * (state%ph_lev(i,j,k) - state%ph_lev(i,j-1,k)) + &
              state%ak_t_lat(i,j,k) * (state%m(i,j,k) - state%m(i,j-1,k))                &
            )
#else
            tend%pgf_lat(i,j,k) = Rd / mesh%de_lat(j) / state%m_lat(i,j,k) * (           &
              state%t_lnpop_lat(i,j,k) * (state%ph_lev(i,j+1,k) - state%ph_lev(i,j,k)) + &
              state%ak_t_lat(i,j,k) * (state%m(i,j+1,k) - state%m(i,j,k))                &
            )
#endif
          end do
        end do
      end do

      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          if (block%reduced_mesh(j)%reduce_factor > 1) then
            call zero_halo(block, tend%pgf_lon(:,j,k), east_halo=.true.)
            do move = 1, block%reduced_mesh(j)%reduce_factor
              do i = block%reduced_mesh(j)%full_lon_ibeg, block%reduced_mesh(j)%full_lon_iend
                block%reduced_tend(j)%pgf_lon(i,k) = (                                            &
                  block%reduced_state(j)%gz(k,i+1,0,move) - block%reduced_state(j)%gz(k,i,0,move) &
                ) / block%reduced_mesh(j)%de_lon(0)
              end do
              call reduce_append_array(move, block%reduced_mesh(j), block%reduced_tend(j)%pgf_lon(:,k), mesh, tend%pgf_lon(:,j,k))
            end do
            call overlay_inner_halo(block, tend%pgf_lon(:,j,k), west_halo=.true.)
          else
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              tend%pgf_lon(i,j,k) = tend%pgf_lon(i,j,k) + (state%gz(i+1,j,k) - state%gz(i,j,k)) / mesh%de_lon(j)
            end do
          end if
        end do
      end do
  
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            tend%pgf_lat(i,j,k) = tend%pgf_lat(i,j,k) + (state%gz(i,j  ,k) - state%gz(i,j-1,k)) / mesh%de_lat(j)
#else
            tend%pgf_lat(i,j,k) = tend%pgf_lat(i,j,k) + (state%gz(i,j+1,k) - state%gz(i,j  ,k)) / mesh%de_lat(j)
#endif
          end do
        end do
      end do
    end if

  end subroutine pgf_sb81_run

end module pgf_sb81_mod
