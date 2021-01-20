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

  subroutine smag_damp_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k
    real(r8) ls

    mesh => state%mesh
    
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          state%tension_h(i,j,k) = (                    &
            state%u(i,j,k) - state%u(i-1,j,k)           &
          ) / mesh%de_lon(j) - (                        &
            state%v(i,j+1,k) * mesh%half_cos_lat(j+1) - &
            state%v(i,j  ,k) * mesh%half_cos_lat(j  )   &
          ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
#else
          state%tension_h(i,j,k) = (                    &
            state%u(i,j,k) - state%u(i-1,j,k)           &
          ) / mesh%de_lon(j) - (                        &
            state%v(i,j  ,k) * mesh%half_cos_lat(j  ) - &
            state%v(i,j-1,k) * mesh%half_cos_lat(j-1)   &
          ) / mesh%le_lon(j) / mesh%full_cos_lat(j)
#endif
        end do
      end do
    end do
    call fill_halo(block, state%tension_h, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%shear_h(i,j,k) = (                      &
            state%v(i+1,j,k) - state%v(i,j,k)           &
          ) / mesh%le_lat(j) + (                        &
            state%u(i,j  ,k) * mesh%full_cos_lat(j  ) - &
            state%u(i,j-1,k) * mesh%full_cos_lat(j-1)   &
          ) / mesh%de_lat(j) / mesh%half_cos_lat(j)
#else
          state%shear_h(i,j,k) = (                      &
            state%v(i+1,j,k) - state%v(i,j,k)           &
          ) / mesh%le_lat(j) + (                        &
            state%u(i,j+1,k) * mesh%full_cos_lat(j+1) - &
            state%u(i,j  ,k) * mesh%full_cos_lat(j  )   &
          ) / mesh%de_lat(j) / mesh%half_cos_lat(j)
#endif
        end do
      end do
    end do
    call fill_halo(block, state%shear_h, full_lon=.false., full_lat=.false., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        ls = radius * mesh%full_cos_lat(j) * mesh%dlon
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          state%kmh(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(     &
            state%tension_h(i,j,k)**2 +                                 &
            0.25_r8 * (                                                 &
              state%shear_h(i-1,j  ,k)**2 + state%shear_h(i,j  ,k)**2 + &
              state%shear_h(i-1,j+1,k)**2 + state%shear_h(i,j+1,k)**2   &
            )                                                           &
          )
#else
          state%kmh(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(     &
            state%tension_h(i,j,k)**2 +                                 &
            0.25_r8 * (                                                 &
              state%shear_h(i-1,j  ,k)**2 + state%shear_h(i,j  ,k)**2 + &
              state%shear_h(i-1,j-1,k)**2 + state%shear_h(i,j-1,k)**2   &
            )                                                           &
          )
#endif
        end do
      end do
    end do
    call fill_halo(block, state%kmh, full_lon=.true., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        ls = radius * mesh%half_cos_lat(j) * mesh%dlon
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%kmh_vtx(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(     &
            state%shear_h(i,j,k)**2 +                                       &
            0.25_r8 * (                                                     &
              state%tension_h(i,j  ,k)**2 + state%tension_h(i+1,j  ,k)**2 + &
              state%tension_h(i,j-1,k)**2 + state%tension_h(i+1,j-1,k)**2   &
            )
          )
#else
          state%kmh_vtx(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(     &
            state%shear_h(i,j,k)**2 +                                       &
            0.25_r8 * (                                                     &
              state%tension_h(i,j  ,k)**2 + state%tension_h(i+1,j  ,k)**2 + &
              state%tension_h(i,j+1,k)**2 + state%tension_h(i+1,j+1,k)**2   &
            )                                                               &
          )
#endif
        end do
      end do
    end do
    call fill_halo(block, state%kmh_vtx, full_lon=.false., full_lat=.false., full_lev=.true.)
    
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        ls = radius * mesh%full_cos_lat(j) * mesh%dlon
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%kmh_lon(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(            &
            0.5_r8 * (state%tension_h(i,j,k)**2 + state%tension_h(i+1,j  ,k)**2) + &
            0.5_r8 * (state%shear_h  (i,j,k)**2 + state%shear_h  (i  ,j+1,k)**2)   &
          )
#else
          state%kmh_lon(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(            &
            0.5_r8 * (state%tension_h(i,j,k)**2 + state%tension_h(i+1,j  ,k)**2) + &
            0.5_r8 * (state%shear_h  (i,j,k)**2 + state%shear_h  (i  ,j-1,k)**2)   &
          )
#endif
        end do
      end do
    end do
    call fill_halo(block, state%kmh_lon, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        ls = radius * mesh%half_cos_lat(j) * mesh%dlon
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          state%kmh_lat(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(            &
            0.5_r8 * (state%tension_h(i,j,k)**2 + state%tension_h(i  ,j-1,k)**2) + &
            0.5_r8 * (  state%shear_h(i,j,k)**2 +   state%shear_h(i-1,j  ,k)**2)   &
          )
#else 
          state%kmh_lat(i,j,k) = 2.0 * (smag_damp_coef * ls)**2 * sqrt(            &
            0.5_r8 * (state%tension_h(i,j,k)**2 + state%tension_h(i  ,j+1,k)**2) + &
            0.5_r8 * (  state%shear_h(i,j,k)**2 +   state%shear_h(i-1,j  ,k)**2)   &
          )
#endif
        end do
      end do
    end do
    call fill_halo(block, state%kmh_lat, full_lon=.true., full_lat=.false., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          state%u(i,j,k) = state%u(i,j,k) + dt * ((           &
            state%kmh(i+1,j,k) * state%tension_h(i+1,j,k) -   &
            state%kmh(i  ,j,k) * state%tension_h(i  ,j,k)     &
          ) / mesh%de_lon(j) + (                              &
            state%kmh_vtx(i,j+1,k) * state%shear_h(i,j+1,k) - &
            state%kmh_vtx(i,j  ,k) * state%shear_h(i,j  ,k)   &
          ) / mesh%le_lon(j))
#else
          state%u(i,j,k) = state%u(i,j,k) + dt * ((           &
            state%kmh(i+1,j,k) * state%tension_h(i+1,j,k) -   &
            state%kmh(i  ,j,k) * state%tension_h(i  ,j,k)     &
          ) / mesh%de_lon(j) + (                              &
            state%kmh_vtx(i,j  ,k) * state%shear_h(i,j  ,k) - &
            state%kmh_vtx(i,j-1,k) * state%shear_h(i,j-1,k)   &
          ) / mesh%le_lon(j))
#endif
        end do
      end do
    end do
    call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)    

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          state%v(i,j,k) = state%v(i,j,k) + dt * ((           &
            state%kmh_vtx(i  ,j,k) * state%shear_h(i  ,j,k) - &
            state%kmh_vtx(i-1,j,k) * state%shear_h(i-1,j,k)   &
          ) / mesh%le_lat(j) - (                              &
            state%kmh(i,j  ,k) * state%tension_h(i,j  ,k) -   &
            state%kmh(i,j-1,k) * state%tension_h(i,j-1,k)     &
          ) / mesh%de_lat(j))
#else
          state%v(i,j,k) = state%v(i,j,k) + dt * ((           &
            state%kmh_vtx(i  ,j,k) * state%shear_h(i  ,j,k) - &
            state%kmh_vtx(i-1,j,k) * state%shear_h(i-1,j,k)   &
          ) / mesh%le_lat(j) - (                              &
            state%kmh(i,j+1,k) * state%tension_h(i,j+1,k) -   &
            state%kmh(i,j  ,k) * state%tension_h(i,j  ,k)     &
          ) / mesh%de_lat(j))
#endif
        end do
      end do
    end do 
    call fill_halo(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
    
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          if (j == mesh%full_lat_ibeg_no_pole) then
            state%pt(i,j,k) = state%pt(i,j,k) + dt * ((                          &
              state%kmh_lon(i  ,j,k) * (state%pt(i+1,j,k) - state%pt(i  ,j,k)) - &
              state%kmh_lon(i-1,j,k) * (state%pt(i  ,j,k) - state%pt(i-1,j,k))   &
            ) / mesh%de_lon(j)**2 + (                                            &
              state%kmh_lat(i,j+1,k) * mesh%half_cos_lat(j+1) * (                &
                state%pt(i,j+1,k) - state%pt(i,j  ,k)                            &
              ) / mesh%de_lat(j+1)                                               &
            ) / mesh%full_cos_lat(j) / mesh%le_lon(j))
          else if (j == mesh%full_lat_iend_no_pole) then
            state%pt(i,j,k) = state%pt(i,j,k) + dt * ((                          &
              state%kmh_lon(i  ,j,k) * (state%pt(i+1,j,k) - state%pt(i  ,j,k)) - &
              state%kmh_lon(i-1,j,k) * (state%pt(i  ,j,k) - state%pt(i-1,j,k))   &
            ) / mesh%de_lon(j)**2 - (                                            &
              state%kmh_lat(i,j  ,k) * mesh%half_cos_lat(j  ) * (                &
                state%pt(i,j  ,k) - state%pt(i,j-1,k)                            &
              ) / mesh%de_lat(j  )                                               &
            ) / mesh%full_cos_lat(j) / mesh%le_lon(j))
          else
            state%pt(i,j,k) = state%pt(i,j,k) + dt * ((                          &
              state%kmh_lon(i  ,j,k) * (state%pt(i+1,j,k) - state%pt(i  ,j,k)) - &
              state%kmh_lon(i-1,j,k) * (state%pt(i  ,j,k) - state%pt(i-1,j,k))   &
            ) / mesh%de_lon(j)**2 + (                                            &
              state%kmh_lat(i,j+1,k) * mesh%half_cos_lat(j+1) * (                &
                state%pt(i,j+1,k) - state%pt(i,j  ,k)                            &
              ) / mesh%de_lat(j+1) -                                             &
              state%kmh_lat(i,j  ,k) * mesh%half_cos_lat(j  ) * (                &
                state%pt(i,j  ,k) - state%pt(i,j-1,k)                            &
              ) / mesh%de_lat(j  )                                               &
            ) / mesh%full_cos_lat(j) / mesh%le_lon(j))
          end if
#else
          state%pt(i,j,k) = state%pt(i,j,k) + dt * ((                          &
            state%kmh_lon(i  ,j,k) * (state%pt(i+1,j,k) - state%pt(i  ,j,k)) - &
            state%kmh_lon(i-1,j,k) * (state%pt(i  ,j,k) - state%pt(i-1,j,k))   &
          ) / mesh%de_lon(j)**2 + (                                            &
            state%kmh_lat(i,j  ,k) * mesh%half_cos_lat(j  ) * (                &
              state%pt(i,j+1,k) - state%pt(i,j  ,k)                            &
            ) / mesh%de_lat(j  ) -                                             &
            state%kmh_lat(i,j-1,k) * mesh%half_cos_lat(j-1) * (                &
              state%pt(i,j  ,k) - state%pt(i,j-1,k)                            &
            ) / mesh%de_lat(j-1)                                               &
          ) / mesh%full_cos_lat(j) / mesh%le_lon(j))
#endif
        end do
      end do
    end do
    call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine smag_damp_run

end module smag_damp_mod
