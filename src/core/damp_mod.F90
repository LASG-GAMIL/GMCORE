module damp_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use operators_mod

  implicit none

  private

  public damp_init
  public damp_final
  public zonal_damp
  public latlon_damp_vtx
  public div_damp

  real(r8), allocatable :: cdiv_full_lat(:,:)
  real(r8), allocatable :: cdiv_half_lat(:,:)

contains

  subroutine damp_init()

    integer j, k

    integer r

    allocate(cdiv_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend))
    allocate(cdiv_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend))

    select case (div_damp_order)
    case (2)
      r = 1

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cdiv_full_lat(j,k) = div_damp_coef2 * global_mesh%full_cos_lat(j)**r * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cdiv_half_lat(j,k) = div_damp_coef2 * global_mesh%half_cos_lat(j)**r * &
            radius**2 * global_mesh%dlat(j) * global_mesh%dlon / dt_in_seconds
        end do
      end do
    case (4)
      r = 2

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          cdiv_full_lat(j,k) = div_damp_coef4 * global_mesh%full_cos_lat(j)**r * &
            radius**4 * global_mesh%dlat(j)**2 * global_mesh%dlon**2 / dt_in_seconds
        end do
      end do

      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          cdiv_half_lat(j,k) = div_damp_coef4 * global_mesh%half_cos_lat(j)**r * &
            radius**4 * global_mesh%dlat(j)**2 * global_mesh%dlon**2 / dt_in_seconds
        end do
      end do
    end select

  end subroutine damp_init

  subroutine damp_final()

    if (allocated(cdiv_full_lat)) deallocate(cdiv_full_lat)
    if (allocated(cdiv_half_lat)) deallocate(cdiv_half_lat)

  end subroutine damp_final

  subroutine zonal_damp(block, order, dt, dx, lb, ub, n, f, async)

    type(block_type), intent(in) :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(in) :: dx
    integer, intent(in) :: lb
    integer, intent(in) :: ub
    integer, intent(in) :: n
    real(r8), intent(inout) :: f(lb:ub)
    type(async_type), intent(inout), optional :: async

    integer, parameter :: diff_halo_width(2:8) = [1, 2, 2, 3, 3, 4, 4]

    real(r8) :: diff_weights(9,2:8) = reshape([  &
       1, -2,   1,   0,   0,   0,   0,   0,  0,  & ! 2
      -1,  3, - 3,   1,   0,   0,   0,   0,  0,  & ! 3
       1, -4,   6, - 4,   1,   0,   0,   0,  0,  & ! 4
      -1,  5, -10,  10, - 5,   1,   0,   0,  0,  & ! 5
       1, -6,  15, -20,  15, - 6,   1,   0,  0,  & ! 6
      -1,  7, -21,  35, -35,  21, - 7,   1,  0,  & ! 7
       1, -8,  28, -56,  70, -56,  28, - 8,  1   & ! 8
    ], [9, 7])

    integer i, ns
    real(r8) g(lb:ub)
    real(r8) df(lb:ub)
    real(r8) a, w(9)

    a  = (dx / 2.0_r8)**order / dt

    call fill_zonal_halo(block, 1 - lb, f)
    if (order == 2) then
      ns = diff_halo_width(order)
      w  = diff_weights(:,order)
      do i = 1, n
        g(i) = sum(f(i-ns:i+ns) * w(:2*ns+1))
      end do
      g = g * (-1)**(order / 2 + 1) * a / dx**order
      do i = 1, n
        f(i) = f(i) + dt * g(i)
      end do
    else
      ns = diff_halo_width(order - 1)
      w  = diff_weights(:,order - 1)
      do i = 1, n
        g (i) = sum(f(i-ns+1:i+ns) * w(:2*ns))
        df(i) = f(i+1) - f(i)
      end do
      g = g * (-1)**(order / 2) * a / dx**order
      do i = 1, n
        g(i) = g(i) * max(0.0_r8, sign(1.0_r8, -g(i) * df(i)))
      end do
      call fill_zonal_halo(block, 1 - lb, g, east_halo=.false.)
      do i = 1, n
        f(i) = f(i) - dt * (g(i) - g(i-1))
      end do
    end if
    if (present(async)) then
      call fill_zonal_halo(block, 1 - lb, f, async)
    else
      call fill_zonal_halo(block, 1 - lb, f)
    end if

  end subroutine zonal_damp

  subroutine latlon_damp_vtx(block, order, dt, f)

    ! Scalar diffusion:
    !
    ! 2nd order:
    !
    !   ∂ F         1      ∂² F        1      ∂          ∂ F
    !   --- = 𝞶 ---------- ---- + 𝞶 --------- --- cos(φ) ---
    !   ∂ t     a² cos²(φ) ∂ λ²     a² cos(φ) ∂ φ        ∂ φ

    type(block_type), intent(in), target :: block
    integer, intent(in) :: order
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub)

    type(mesh_type), pointer :: mesh
    real(r8), dimension(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                        block%mesh%half_lat_lb:block%mesh%half_lat_ub) :: f0, df
    real(r8), parameter :: latlon_damp_coef = 1.0d-12
    integer i, j, o
    integer sign

    mesh => block%mesh

    df = 0
    f0 = f
    do o = 1, order / 2
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          df(i,j) = (f0(i+1,j) - 2 * f0(i,j) + f0(i-1,j)) / mesh%le_lat(j)**2 + ( &
#ifdef V_POLE
                    (f0(i,j+1) - f0(i,j  )) * mesh%full_cos_lat(j  ) - &
                    (f0(i,j  ) - f0(i,j-1)) * mesh%full_cos_lat(j-1)   &
#else
                    (f0(i,j+1) - f0(i,j  )) * mesh%full_cos_lat(j+1) - &
                    (f0(i,j  ) - f0(i,j-1)) * mesh%full_cos_lat(j  )   &
#endif
                    ) / mesh%de_lat(j)**2 / mesh%half_cos_lat(j)
        end do
      end do
      if (o /= order / 2) then
        call fill_halo(block, df, full_lon=.false., full_lat=.false.)
        f0 = df
      end if
    end do

    sign = (-1)**(order / 2 + 1)

    do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
      do i = mesh%half_lon_ibeg, mesh%half_lon_iend
        f(i,j) = f(i,j) + sign * dt * latlon_damp_coef * df(i,j)
      end do
    end do

    call fill_halo(block, f, full_lon=.false., full_lat=.false.)

  end subroutine latlon_damp_vtx

  subroutine div_damp(block, old_state, new_state, dt)

    type(block_type), intent(in) :: block
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state
    real(8), intent(in) :: dt

    type(mesh_type), pointer :: mesh
    integer i, j, k, cyc

    if (baroclinic) then
      mesh => old_state%mesh

      select case (div_damp_order)
      case (2)
        new_state%div = old_state%div
        do cyc = 1, div_damp_cycles
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
              do i = mesh%half_lon_ibeg, mesh%half_lon_iend
                new_state%u(i,j,k) = new_state%u(i,j,k) + dt / div_damp_cycles * cdiv_full_lat(j,k) * ( &
                  new_state%div(i+1,j,k) - new_state%div(i,j,k)) / mesh%de_lon(j)
              end do
            end do
          end do
          call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
              do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
                new_state%v(i,j,k) = new_state%v(i,j,k) + dt / div_damp_cycles * cdiv_half_lat(j,k) * ( &
                  new_state%div(i,j,k) - new_state%div(i,j-1,k)) / mesh%de_lat(j)
#else
                new_state%v(i,j,k) = new_state%v(i,j,k) + dt / div_damp_cycles * cdiv_half_lat(j,k) * ( &
                  new_state%div(i,j+1,k) - new_state%div(i,j,k)) / mesh%de_lat(j)
#endif
              end do
            end do
          end do
          call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)

          call calc_div(block, new_state)
        end do
      case (4)
        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              new_state%u(i,j,k) = new_state%u(i,j,k) - dt * cdiv_full_lat(j,k) * ( &
                old_state%div2(i+1,j,k) - old_state%div2(i,j,k)) / mesh%de_lon(j)
            end do
          end do
        end do
        call fill_halo(block, new_state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

        do k = mesh%full_lev_ibeg, mesh%full_lev_iend
          do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              new_state%v(i,j,k) = new_state%v(i,j,k) - dt * cdiv_half_lat(j,k) * ( &
                old_state%div2(i,j,k) - old_state%div2(i,j-1,k)) / mesh%de_lat(j)
#else
              new_state%v(i,j,k) = new_state%v(i,j,k) - dt * cdiv_half_lat(j,k) * ( &
                old_state%div2(i,j+1,k) - old_state%div2(i,j,k)) / mesh%de_lat(j)
#endif
            end do
          end do
        end do
#ifndef V_POLE
        if (mesh%has_south_pole()) then
          j = mesh%half_lat_ibeg
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%v(i,j,k) = new_state%v(i,j,k) + dt * 3.5d4 * ( &
                old_state%div(i,j+1,k) - old_state%div(i,j,k)) / mesh%de_lat(j)
            end do
          end do
        end if
        if (mesh%has_north_pole()) then
          j = mesh%half_lat_iend
          do k = mesh%full_lev_ibeg, mesh%full_lev_iend
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              new_state%v(i,j,k) = new_state%v(i,j,k) + dt * 3.5d4 * ( &
                old_state%div(i,j+1,k) - old_state%div(i,j,k)) / mesh%de_lat(j)
            end do
          end do
        end if
#endif
        call fill_halo(block, new_state%v, full_lon=.true., full_lat=.false., full_lev=.true.)
      end select
    end if

  end subroutine div_damp

end module damp_mod
