module damp_mod

  use const_mod
  use namelist_mod
  use parallel_mod
  use block_mod
  use operators_mod
  use reduce_mod

  implicit none

  private

  public damp_init
  public damp_final
  public zonal_damp
  public polar_damp
  public latlon_damp_lon
  public latlon_damp_lat
  public latlon_damp_cell
  public div_damp
  public shapiro_smooth

  integer, parameter :: diff_halo_width(2:8) = [1, 2, 2, 3, 3, 4, 4]

  real(r8), target :: diff_weights(9,2:8) = reshape([  &
     1, -2,   1,   0,   0,   0,   0,   0,  0,  & ! 2
    -1,  3, - 3,   1,   0,   0,   0,   0,  0,  & ! 3
     1, -4,   6, - 4,   1,   0,   0,   0,  0,  & ! 4
    -1,  5, -10,  10, - 5,   1,   0,   0,  0,  & ! 5
     1, -6,  15, -20,  15, - 6,   1,   0,  0,  & ! 6
    -1,  7, -21,  35, -35,  21, - 7,   1,  0,  & ! 7
     1, -8,  28, -56,  70, -56,  28, - 8,  1   & ! 8
  ], [9, 7])

  integer , allocatable :: ox_full_lat(:)
  integer , allocatable :: ox_half_lat(:)
  real(r8), allocatable :: cx_full_lat(:,:,:)
  real(r8), allocatable :: cy_full_lat(:,:,:)
  real(r8), allocatable :: cx_half_lat(:,:,:)
  real(r8), allocatable :: cy_half_lat(:,:,:)
  real(r8), allocatable :: sx_full_lat(:,:)   ! Shapiro filter strength coefficient
  real(r8), allocatable :: sy_full_lat(:,:)   !
  real(r8), allocatable :: sx_half_lat(:,:)   !
  real(r8), allocatable :: sy_half_lat(:,:)   !
  real(r8), allocatable :: cdiv_full_lat(:,:)
  real(r8), allocatable :: cdiv_half_lat(:,:)

contains

  subroutine damp_init()

    integer j, jr, jr0, k, order, r
    real(r8) c0, polar_damp_mul

    call damp_final()

    allocate(ox_full_lat(global_mesh%full_lat_ibeg_no_pole:global_mesh%full_lat_iend_no_pole))
    allocate(ox_half_lat(global_mesh%half_lat_ibeg_no_pole:global_mesh%half_lat_iend_no_pole))

    jr0 = 0
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
        if (reduce_factors(jr) > 1) jr0 = jr
      end if
    end do

    ! Set zonal damping orders.
    do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
      if (global_mesh%full_lat(j) <= 0) then
        jr = j - global_mesh%full_lat_ibeg_no_pole + 1
      else
        jr = global_mesh%full_lat_iend_no_pole - j + 1
      end if
      ox_full_lat(j) = merge(max(polar_damp_order, zonal_damp_orders(jr)), polar_damp_order, jr <= jr0)
    end do
    do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
      if (global_mesh%half_lat(j) <= 0) then
        jr = j - global_mesh%half_lat_ibeg_no_pole + 1
      else
        jr = global_mesh%half_lat_iend_no_pole - j + 1
      end if
#ifdef V_POLE
      ox_half_lat(j) = merge(max(polar_damp_order, min(zonal_damp_orders(jr), zonal_damp_orders(jr-1))), polar_damp_order, jr <= jr0)
#else
      ox_half_lat(j) = merge(max(polar_damp_order, min(zonal_damp_orders(jr), zonal_damp_orders(jr+1))), polar_damp_order, jr <= jr0)
#endif
    end do

    allocate(cx_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cy_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cx_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))
    allocate(cy_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend,2:4))

    do order = 2, 4, 2
      c0 = 0.5_r8**order / dt_in_seconds
      polar_damp_mul = 2
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
          ! Check if grids are reduced at this zonal circle.
          if (global_mesh%full_lat(j) <= 0) then
            jr = j - global_mesh%full_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%full_lat_iend_no_pole - j + 1
          end if
          cx_full_lat(j,k,order) = c0 * exp(-jr**2 / ((-jr0 / log(0.5)) * jr0))
          cy_full_lat(j,k,order) = c0 * exp(-jr**2 / ((-jr0 / log(0.5)) * jr0))
        end do
      end do
      do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
        do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
          if (global_mesh%half_lat(j) <= 0) then
            jr = j - global_mesh%half_lat_ibeg_no_pole + 1
          else
            jr = global_mesh%half_lat_iend_no_pole - j + 1
          end if
          cx_half_lat(j,k,order) = c0 * exp(-jr**2 / ((-jr0 / log(0.5)) * jr0))
          cy_half_lat(j,k,order) = c0 * exp(-jr**2 / ((-jr0 / log(0.5)) * jr0))
        end do
      end do
    end do

    allocate(sx_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend))
    allocate(sy_full_lat(global_mesh%full_lat_ibeg:global_mesh%full_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend))
    allocate(sx_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend))
    allocate(sy_half_lat(global_mesh%half_lat_ibeg:global_mesh%half_lat_iend,global_mesh%full_lev_ibeg:global_mesh%full_lev_iend))

    polar_damp_mul = 2
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = global_mesh%full_lat_ibeg_no_pole, global_mesh%full_lat_iend_no_pole
        ! Check if grids are reduced at this zonal circle.
        if (global_mesh%full_lat(j) <= 0) then
          jr = j - global_mesh%full_lat_ibeg_no_pole + 1
        else
          jr = global_mesh%full_lat_iend_no_pole - j + 1
        end if
        sx_full_lat(j,k) = exp(-jr**2 / ((-jr0 / log(0.5)) * jr0))
        sy_full_lat(j,k) = exp(-jr**2 / ((-jr0 / log(0.5)) * jr0))
      end do
    end do
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      do j = global_mesh%half_lat_ibeg_no_pole, global_mesh%half_lat_iend_no_pole
        if (global_mesh%half_lat(j) <= 0) then
          jr = j - global_mesh%half_lat_ibeg_no_pole + 1
        else
          jr = global_mesh%half_lat_iend_no_pole - j + 1
        end if
        sx_half_lat(j,k) = exp(-jr**2 / ((-jr0 / log(0.5)) * jr0))
        sy_half_lat(j,k) = exp(-jr**2 / ((-jr0 / log(0.5)) * jr0))
      end do
    end do

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

    if (allocated(ox_full_lat  )) deallocate(ox_full_lat  )
    if (allocated(ox_half_lat  )) deallocate(ox_half_lat  )
    if (allocated(cx_full_lat  )) deallocate(cx_full_lat  )
    if (allocated(cy_full_lat  )) deallocate(cy_full_lat  )
    if (allocated(cx_half_lat  )) deallocate(cx_half_lat  )
    if (allocated(cy_half_lat  )) deallocate(cy_half_lat  )
    if (allocated(sx_full_lat  )) deallocate(sx_full_lat  )
    if (allocated(sy_full_lat  )) deallocate(sy_full_lat  )
    if (allocated(sx_half_lat  )) deallocate(sx_half_lat  )
    if (allocated(sy_half_lat  )) deallocate(sy_half_lat  )
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

  subroutine polar_damp(block, dt, state)

    type(block_type), intent(in), target :: block
    real(r8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    type(reduced_mesh_type), pointer :: reduced_mesh
    type(reduced_state_type), pointer :: reduced_state
    type(reduced_tend_type), pointer :: reduced_tend
    real(r8), pointer, dimension(:,:,:) :: gx, gy
    integer i, io, j, jr, jb, k, move, ox, oy, nsx, nsy
    real(r8), pointer, dimension(:) :: wx(:), wy(:)

    mesh => block%mesh
    gx   => block%latlon_damp_lon_gx
    gy   => block%latlon_damp_lon_gy
    oy = polar_damp_order
    if (oy == 2) then
      nsy = diff_halo_width(oy)
      wy  => diff_weights(:,oy)
    else
      nsy = diff_halo_width(oy - 1)
      wy  => diff_weights(:,oy - 1)
    end if
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        ox = ox_full_lat(j)
        if (ox == 2) then
          nsx = diff_halo_width(ox)
          wx  => diff_weights(:,ox)
        else
          nsx = diff_halo_width(ox - 1)
          wx  => diff_weights(:,ox - 1)
        end if
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          reduced_mesh  => block%reduced_mesh (j)
          reduced_state => block%reduced_state(j)
          reduced_tend  => block%reduced_tend (j)
          gx(:,j,k) = 0.0_r8
          do move = 1, reduced_mesh%reduce_factor
            ! Calculate damping flux at interfaces.
            do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
              reduced_tend%gx(i,k) = sum(reduced_state%u(k,i-nsx:i+nsx-1,0,move) * wx(:2*nsx))
            end do
            ! Limit damping flux to avoid upgradient (Xue 2000).
            if (ox > 2) then
              reduced_tend%gx(:,k) = reduced_tend%gx(:,k) * (-1)**(ox / 2)
              do i = reduced_mesh%full_lon_ibeg, reduced_mesh%full_lon_iend
                reduced_tend%gx(i,k) = reduced_tend%gx(i,k) *        &
                  max(0.0_r8, sign(1.0_r8, -reduced_tend%gx(i,k) * ( &
                    reduced_state%u(k,i  ,0,move) -                  &
                    reduced_state%u(k,i-1,0,move)                    &
                  )))
              end do
            end if
            call reduce_append_array(move, reduced_mesh, reduced_tend%gx(:,k), mesh, gx(:,j,k))
          end do
          call overlay_inner_halo(block, gx(:,j,k), west_halo=.true.)
        else
          ! Calculate damping flux at interfaces.
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            !
            !                               i
            !   o      |      o      |      o      |      o      |      o
            !         i-2           i-1            i            i+1
            !                               ^
            !   o - cell   | - edge
            !
            gx(i,j,k) = sum(state%u(i-nsx:i+nsx-1,j,k) * wx(:2*nsx))
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (ox > 2) then
            gx(:,j,k) = gx(:,j,k) * (-1)**(ox / 2)
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * (state%u(i,j,k) - state%u(i-1,j,k))))
            end do
          end if
        end if
      end do
    end do
    call fill_halo(block, gx, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false.)
    ! - Merdional damping
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        ! Calculate damping flux at interfaces.
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          !
          !                               j
          !   |      o      |      o      |      o      |      o      |
          !         j-2           j-1            j            j+1
          !                               ^
          !   o - cell   | - edge
          !
          gy(i,j,k) = sum(state%u(i,j-nsy:j+nsy-1,k) * wy(:2*nsy))
#else
          !
          !                               j
          !   |      o      |      o      |      o      |      o      |
          !         j-1            j            j+1           j+2
          !                               ^
          !   o - cell   | - edge
          !
          gy(i,j,k) = sum(state%u(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
#endif
        end do
        ! Limit damping flux to avoid upgradient (Xue 2000).
        if (oy > 2) then
          gy(:,j,k) = gy(:,j,k) * (-1)**(oy / 2)
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (state%u(i,j,k) - state%u(i,j-1,k))))
#else
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (state%u(i,j+1,k) - state%u(i,j,k))))
#endif
          end do
        end if
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false.)
#else
    call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true., north_halo=.false.)
#endif
    ! Damp physical variable at last.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        ox = ox_full_lat(j)
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
          !
          !               - j+1
          !               |
          !               |
          !        o------|------o
          !        i      |     i+1
          !               |
          !               - j
          !
          state%u(i,j,k) = state%u(i,j,k) - dt / polar_damp_cycles * ( &
            (gx(i+1,j,k) - gx(i,j,k)) * cx_full_lat(j,k,ox) +          &
            (gy(i,j+1,k) - gy(i,j,k)) * cy_full_lat(j,k,oy)            &
          )
#else
          !
          !               - j
          !               |
          !               |
          !        o------|------o
          !        i      |     i+1
          !               |
          !               - j-1
          !
          state%u(i,j,k) = state%u(i,j,k) - dt / polar_damp_cycles * ( &
            (gx(i+1,j,k) - gx(i,j,k)) * cx_full_lat(j,k,ox) +          &
            (gy(i,j,k) - gy(i,j-1,k)) * cy_full_lat(j,k,oy)            &
          )
#endif
        end do
      end do
    end do
    call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

    gx => block%latlon_damp_lat_gx
    gy => block%latlon_damp_lat_gy
    oy = polar_damp_order
    if (oy == 2) then
      nsy = diff_halo_width(oy)
      wy  => diff_weights(:,oy)
    else
      nsy = diff_halo_width(oy - 1)
      wy  => diff_weights(:,oy - 1)
    end if
    ! - Zonal damping
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        ox = ox_half_lat(j)
        if (ox == 2) then
          nsx = diff_halo_width(ox)
          wx  => diff_weights(:,ox)
        else
          nsx = diff_halo_width(ox - 1)
          wx  => diff_weights(:,ox - 1)
        end if
#ifdef V_POLE
        if (mesh%half_lat(j) < 0) then
          jr = j - 1
          jb = 1
        else
          jr = j
          jb = 0
        end if
#else
        if (mesh%half_lat(j) < 0) then
          jr = j + 1
          jb = -1
        else
          jr = j
          jb = 0
        end if
#endif
        if (block%reduced_mesh(jr)%reduce_factor > 1) then
          reduced_mesh  => block%reduced_mesh (jr)
          reduced_state => block%reduced_state(jr)
          reduced_tend  => block%reduced_tend (jr)
          gx(:,j,k) = 0.0_r8
          do move = 1, reduced_mesh%reduce_factor
            ! Calculate damping flux at interfaces.
            do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
              reduced_tend%gx(i,k) = sum(reduced_state%v(k,i-nsx+1:i+nsx,jb,move) * wx(:2*nsx))
            end do
            ! Limit damping flux to avoid upgradient (Xue 2000).
            if (ox > 2) then
              reduced_tend%gx(:,k) = reduced_tend%gx(:,k) * (-1)**(ox / 2)
              do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
                reduced_tend%gx(i,k) = reduced_tend%gx(i,k) *        &
                  max(0.0_r8, sign(1.0_r8, -reduced_tend%gx(i,k) * ( &
                    reduced_state%v(k,i+1,jb,move) -                 &
                    reduced_state%v(k,i  ,jb,move)                   &
                  )))
              end do
            end if
            call reduce_append_array(move, reduced_mesh, reduced_tend%gx(:,k), mesh, gx(:,j,k))
          end do
          call overlay_inner_halo(block, gx(:,j,k), west_halo=.true.)
        else
          ! Calculate damping flux at interfaces.
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            !
            !                               i
            !   |      o      |      o      |      o      |      o      |
            !         i-1            i            i+1           i+2
            !                               ^
            !   o - cell   | - edge
            !
            gx(i,j,k) = sum(state%v(i-nsx+1:i+nsx,j,k) * wx(:2*nsx))
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (ox > 2) then
            gx(:,j,k) = gx(:,j,k) * (-1)**(ox / 2)
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * (state%v(i+1,j,k) - state%v(i,j,k))))
            end do
          end if
        end if
      end do
    end do
    call fill_halo(block, gx, full_lon=.false., full_lat=.false., full_lev=.true., east_halo=.false.)
    ! - Merdional damping
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        ! Calculate damping flux at interfaces.
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          !
          !                               j
          !   o      |      o      |      o      |      o      |      o
          !         j-1            j            j+1           j+2
          !                               ^
          !   o - cell   | - edge
          !
          gy(i,j,k) = sum(state%v(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
#else
          !
          !                               j
          !   o      |      o      |      o      |      o      |      o
          !         j-2           j-1            j            j+1
          !                               ^
          !   o - cell   | - edge
          !
          gy(i,j,k) = sum(state%v(i,j-nsy:j+nsy-1,k) * wy(:2*nsy))
#endif
        end do
        ! Limit damping flux to avoid upgradient (Xue 2000).
        if (oy > 2) then
          gy(:,j,k) = gy(:,j,k) * (-1)**(oy / 2)
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (state%v(i,j+1,k) - state%v(i,j,k))))
#else
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (state%v(i,j,k) - state%v(i,j-1,k))))
#endif
          end do
        end if
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, gy, full_lon=.true., full_lat=.true., full_lev=.true., north_halo=.false.)
#else
    call fill_halo(block, gy, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false.)
#endif
    ! Damp physical variable at last.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        ox = ox_half_lat(j)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          !
          !               o j
          !               |
          !               |
          !        |-------------|
          !       i-1     |      i
          !               |
          !               o j-1
          !
          state%v(i,j,k) = state%v(i,j,k) - dt / polar_damp_cycles * ( &
            (gx(i,j,k) - gx(i-1,j,k)) * cx_half_lat(j,k,ox) +          &
            (gy(i,j,k) - gy(i,j-1,k)) * cy_half_lat(j,k,oy)            &
          )
#else
          !
          !               o j+1
          !               |
          !               |
          !        |-------------|
          !       i-1     |      i
          !               |
          !               o j
          !
          state%v(i,j,k) = state%v(i,j,k) - dt / polar_damp_cycles * ( &
            (gx(i,j,k) - gx(i-1,j,k)) * cx_half_lat(j,k,ox) +          &
            (gy(i,j+1,k) - gy(i,j,k)) * cy_half_lat(j,k,oy)            &
          )
#endif
        end do
      end do
    end do
    call fill_halo(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)

    gx   => block%latlon_damp_cell_gx
    gy   => block%latlon_damp_cell_gy
    oy = polar_damp_order
    if (oy == 2) then
      nsy = diff_halo_width(oy)
      wy  => diff_weights(:,oy)
    else
      nsy = diff_halo_width(oy - 1)
      wy  => diff_weights(:,oy - 1)
    end if
    ! - Zonal damping
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        ox = ox_full_lat(j)
        if (ox == 2) then
          nsx = diff_halo_width(ox)
          wx  = diff_weights(:,ox)
        else
          nsx = diff_halo_width(ox - 1)
          wx  = diff_weights(:,ox - 1)
        end if
        if (block%reduced_mesh(j)%reduce_factor > 1) then
          reduced_mesh  => block%reduced_mesh (j)
          reduced_state => block%reduced_state(j)
          reduced_tend  => block%reduced_tend (j)
          gx(:,j,k) = 0.0_r8
          do move = 1, reduced_mesh%reduce_factor
            ! Calculate damping flux at interfaces.
            do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
              reduced_tend%gx(i,k) = sum(reduced_state%pt(k,i-nsx+1:i+nsx,0,move) * wx(:2*nsx))
            end do
            ! Limit damping flux to avoid upgradient (Xue 2000).
            if (ox > 2) then
              reduced_tend%gx(:,k) = reduced_tend%gx(:,k) * (-1)**(ox / 2)
              do i = reduced_mesh%half_lon_ibeg, reduced_mesh%half_lon_iend
                reduced_tend%gx(i,k) = reduced_tend%gx(i,k) *        &
                  max(0.0_r8, sign(1.0_r8, -reduced_tend%gx(i,k) * ( &
                    reduced_state%pt(k,i+1,0,move) -                 &
                    reduced_state%pt(k,i  ,0,move)                   &
                  )))
              end do
            end if
            call reduce_append_array(move, reduced_mesh, reduced_tend%gx(:,k), mesh, gx(:,j,k))
          end do
          call overlay_inner_halo(block, gx(:,j,k), west_halo=.true.)
        else
          ! Calculate damping flux at interfaces.
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            !
            !                               i
            !   |      o      |      o      |      o      |      o      |
            !         i-1            i            i+1           i+2
            !                               ^
            !   o - cell   | - edge
            !
            gx(i,j,k) = sum(state%pt(i-nsx+1:i+nsx,j,k) * wx(:2*nsx))
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (ox > 2) then
            gx(:,j,k) = gx(:,j,k) * (-1)**(ox / 2)
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * (state%pt(i+1,j,k) - state%pt(i,j,k))))
            end do
          end if
        end if
      end do
    end do
    call fill_halo(block, gx, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false.)
    ! - Merdional damping
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        ! Calculate damping flux at interfaces.
#ifdef V_POLE
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          !
          !                               j
          !   |      o      |      o      |      o      |      o      |
          !         j-2           j-1            j            j+1
          !                               ^
          !   o - cell   | - edge
          !
          gy(i,j,k) = sum(state%pt(i,j-nsy:j+nsy-1,k) * wy(:2*nsy))
        end do
#else
        if (mesh%is_outside_pole_full_lat(j-nsy+1)) then
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            io = i + mesh%num_full_lon / 2
            if (io > mesh%num_full_lon) io = io - mesh%num_full_lon
            state%pt(i,j-nsy+1,k) = state%pt(io,j,k)
            gy(i,j,k) = sum(state%pt(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
          end do
        else if (mesh%is_outside_pole_full_lat(j+nsy)) then
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            io = i + mesh%num_full_lon / 2
            if (io > mesh%num_full_lon) io = io - mesh%num_full_lon
            state%pt(i,j+nsy,k) = state%pt(io,j,k)
            gy(i,j,k) = sum(state%pt(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
          end do
        else
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            !
            !                               j
            !   |      o      |      o      |      o      |      o      |
            !         j-1            j            j+1           j+2
            !                               ^
            !   o - cell   | - edge
            !
            gy(i,j,k) = sum(state%pt(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
          end do
        end if
#endif
        ! Limit damping flux to avoid upgradient (Xue 2000).
        if (oy > 2) then
          gy(:,j,k) = gy(:,j,k) * (-1)**(oy / 2)
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (state%pt(i,j,k) - state%pt(i,j-1,k))))
#else
            gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (state%pt(i,j+1,k) - state%pt(i,j,k))))
#endif
          end do
        end if
      end do
    end do
#ifdef V_POLE
    call fill_halo(block, gy, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false.)
#else
    call fill_halo(block, gy, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)
#endif
    ! Damp physical variable at last.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        ox = ox_full_lat(j)
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
          !
          !               - j+1
          !               |
          !               |
          !        o------|------o
          !       i-1     |      i
          !               |
          !               - j
          !
          state%pt(i,j,k) = state%pt(i,j,k) - dt / polar_damp_cycles * ( &
            (gx(i,j,k) - gx(i-1,j,k)) * cx_full_lat(j,k,ox) +            &
            (gy(i,j+1,k) - gy(i,j,k)) * cy_full_lat(j,k,oy)              &
          )
#else
          !
          !               - j
          !               |
          !               |
          !        o------|------o
          !       i-1     |      i
          !               |
          !               - j-1
          !
          state%pt(i,j,k) = state%pt(i,j,k) - dt / polar_damp_cycles * ( &
            (gx(i,j,k) - gx(i-1,j,k)) * cx_full_lat(j,k,ox) +            &
            (gy(i,j,k) - gy(i,j-1,k)) * cy_full_lat(j,k,oy)              &
          )
#endif
        end do
      end do
    end do
    call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine polar_damp

  subroutine latlon_damp_lon(block, dt, f)

    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy
    integer i, j, k, cyc, ox, oy, nsx, nsy
    real(r8) wx(9), wy(9)

    mesh => block%mesh
    gx   => block%latlon_damp_lon_gx
    gy   => block%latlon_damp_lon_gy
    oy = polar_damp_order
    if (oy == 2) then
      nsy = diff_halo_width(oy)
      wy  = diff_weights(:,oy)
    else
      nsy = diff_halo_width(oy - 1)
      wy  = diff_weights(:,oy - 1)
    end if
    cycle_loop: do cyc = 1, polar_damp_cycles
      ! - Zonal damping
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          ox = ox_full_lat(j)
          ! Calculate damping flux at interfaces.
          if (ox == 2) then
            nsx = diff_halo_width(ox)
            wx  = diff_weights(:,ox)
          else
            nsx = diff_halo_width(ox - 1)
            wx  = diff_weights(:,ox - 1)
          end if
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            !
            !                               i
            !   o      |      o      |      o      |      o      |      o
            !         i-2           i-1            i            i+1
            !                               ^
            !   o - cell   | - edge
            !
            gx(i,j,k) = sum(f(i-nsx:i+nsx-1,j,k) * wx(:2*nsx))
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (ox > 2) then
            gx(:,j,k) = gx(:,j,k) * (-1)**(ox / 2)
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * (f(i,j,k) - f(i-1,j,k))))
            end do
          end if
        end do
      end do
      call fill_halo(block, gx, full_lon=.true., full_lat=.true., full_lev=.true., west_halo=.false.)
      ! - Merdional damping
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          ! Calculate damping flux at interfaces.
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            !
            !                               j
            !   |      o      |      o      |      o      |      o      |
            !         j-2           j-1            j            j+1
            !                               ^
            !   o - cell   | - edge
            !
            gy(i,j,k) = sum(f(i,j-nsy:j+nsy-1,k) * wy(:2*nsy))
#else
            !
            !                               j
            !   |      o      |      o      |      o      |      o      |
            !         j-1            j            j+1           j+2
            !                               ^
            !   o - cell   | - edge
            !
            gy(i,j,k) = sum(f(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
#endif
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (oy > 2) then
            gy(:,j,k) = gy(:,j,k) * (-1)**(oy / 2)
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j,k) - f(i,j-1,k))))
#else
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j+1,k) - f(i,j,k))))
#endif
            end do
          end if
        end do
      end do
#ifdef V_POLE
      call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false.)
#else
      call fill_halo(block, gy, full_lon=.false., full_lat=.false., full_lev=.true., north_halo=.false.)
#endif
      ! Damp physical variable at last.
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          ox = ox_full_lat(j)
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
#ifdef V_POLE
            !
            !               - j+1
            !               |
            !               |
            !        o------|------o
            !        i      |     i+1
            !               |
            !               - j
            !
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (    &
              (gx(i+1,j,k) - gx(i,j,k)) * cx_full_lat(j,k,ox) + &
              (gy(i,j+1,k) - gy(i,j,k)) * cy_full_lat(j,k,oy)   &
            )
#else
            !
            !               - j
            !               |
            !               |
            !        o------|------o
            !        i      |     i+1
            !               |
            !               - j-1
            !
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (    &
              (gx(i+1,j,k) - gx(i,j,k)) * cx_full_lat(j,k,ox) + &
              (gy(i,j,k) - gy(i,j-1,k)) * cy_full_lat(j,k,oy)   &
            )
#endif
          end do
        end do
      end do
      call fill_halo(block, f, full_lon=.false., full_lat=.true., full_lev=.true.)
    end do cycle_loop

  end subroutine latlon_damp_lon

  subroutine latlon_damp_lat(block, dt, f)

    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy
    integer i, j, k, cyc, ox, oy, nsx, nsy
    real(r8) wx(9), wy(9)

    mesh => block%mesh
    gx   => block%latlon_damp_lat_gx
    gy   => block%latlon_damp_lat_gy
    oy = polar_damp_order
    if (oy == 2) then
      nsy = diff_halo_width(oy)
      wy  = diff_weights(:,oy)
    else
      nsy = diff_halo_width(oy - 1)
      wy  = diff_weights(:,oy - 1)
    end if
    cycle_loop: do cyc = 1, polar_damp_cycles
      ! - Zonal damping
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          ox = ox_half_lat(j)
          ! Calculate damping flux at interfaces.
          if (ox == 2) then
            nsx = diff_halo_width(ox)
            wx  = diff_weights(:,ox)
          else
            nsx = diff_halo_width(ox - 1)
            wx  = diff_weights(:,ox - 1)
          end if
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            !
            !                               i
            !   |      o      |      o      |      o      |      o      |
            !         i-1            i            i+1           i+2
            !                               ^
            !   o - cell   | - edge
            !
            gx(i,j,k) = sum(f(i-nsx+1:i+nsx,j,k) * wx(:2*nsx))
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (ox > 2) then
            gx(:,j,k) = gx(:,j,k) * (-1)**(ox / 2)
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * (f(i+1,j,k) - f(i,j,k))))
            end do
          end if
        end do
      end do
      call fill_halo(block, gx, full_lon=.false., full_lat=.false., full_lev=.true., east_halo=.false.)
      ! - Merdional damping
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          ! Calculate damping flux at interfaces.
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            !
            !                               j
            !   o      |      o      |      o      |      o      |      o
            !         j-1            j            j+1           j+2
            !                               ^
            !   o - cell   | - edge
            !
            gy(i,j,k) = sum(f(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
#else
            !
            !                               j
            !   o      |      o      |      o      |      o      |      o
            !         j-2           j-1            j            j+1
            !                               ^
            !   o - cell   | - edge
            !
            gy(i,j,k) = sum(f(i,j-nsy:j+nsy-1,k) * wy(:2*nsy))
#endif
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (oy > 2) then
            gy(:,j,k) = gy(:,j,k) * (-1)**(oy / 2)
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j+1,k) - f(i,j,k))))
#else
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j,k) - f(i,j-1,k))))
#endif
            end do
          end if
        end do
      end do
#ifdef V_POLE
      call fill_halo(block, gy, full_lon=.true., full_lat=.true., full_lev=.true., north_halo=.false.)
#else
      call fill_halo(block, gy, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false.)
#endif
      ! Damp physical variable at last.
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          ox = ox_half_lat(j)
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            !
            !               o j
            !               |
            !               |
            !        |-------------|
            !       i-1     |      i
            !               |
            !               o j-1
            !
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (    &
              (gx(i,j,k) - gx(i-1,j,k)) * cx_half_lat(j,k,ox) + &
              (gy(i,j,k) - gy(i,j-1,k)) * cy_half_lat(j,k,oy)   &
            )
#else
            !
            !               o j+1
            !               |
            !               |
            !        |-------------|
            !       i-1     |      i
            !               |
            !               o j
            !
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (    &
              (gx(i,j,k) - gx(i-1,j,k)) * cx_half_lat(j,k,ox) + &
              (gy(i,j+1,k) - gy(i,j,k)) * cy_half_lat(j,k,oy)   &
            )
#endif
          end do
        end do
      end do
      call fill_halo(block, f, full_lon=.true., full_lat=.false., full_lev=.true.)
    end do cycle_loop

  end subroutine latlon_damp_lat

  subroutine latlon_damp_cell(block, dt, f)

    type(block_type), intent(in), target :: block
    real(8), intent(in) :: dt
    real(r8), intent(inout) :: f(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    type(mesh_type), pointer :: mesh
    real(r8), pointer, dimension(:,:,:) :: gx, gy
    integer i, j, k, cyc, ox, oy, nsx, nsy, io
    real(r8) wx(9), wy(9)

    mesh => block%mesh
    gx   => block%latlon_damp_cell_gx
    gy   => block%latlon_damp_cell_gy
    oy = polar_damp_order
    if (oy == 2) then
      nsy = diff_halo_width(oy)
      wy  = diff_weights(:,oy)
    else
      nsy = diff_halo_width(oy - 1)
      wy  = diff_weights(:,oy - 1)
    end if
    cycle_loop: do cyc = 1, polar_damp_cycles
      ! - Zonal damping
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          ox = ox_full_lat(j)
          ! Calculate damping flux at interfaces.
          if (ox == 2) then
            nsx = diff_halo_width(ox)
            wx  = diff_weights(:,ox)
          else
            nsx = diff_halo_width(ox - 1)
            wx  = diff_weights(:,ox - 1)
          end if
          do i = mesh%half_lon_ibeg, mesh%half_lon_iend
            !
            !                               i
            !   |      o      |      o      |      o      |      o      |
            !         i-1            i            i+1           i+2
            !                               ^
            !   o - cell   | - edge
            !
            gx(i,j,k) = sum(f(i-nsx+1:i+nsx,j,k) * wx(:2*nsx))
          end do
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (ox > 2) then
            gx(:,j,k) = gx(:,j,k) * (-1)**(ox / 2)
            do i = mesh%half_lon_ibeg, mesh%half_lon_iend
              gx(i,j,k) = gx(i,j,k) * max(0.0_r8, sign(1.0_r8, -gx(i,j,k) * (f(i+1,j,k) - f(i,j,k))))
            end do
          end if
        end do
      end do
      call fill_halo(block, gx, full_lon=.false., full_lat=.true., full_lev=.true., east_halo=.false.)
      ! - Merdional damping
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
          ! Calculate damping flux at interfaces.
#ifdef V_POLE
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
            !
            !                               j
            !   |      o      |      o      |      o      |      o      |
            !         j-2           j-1            j            j+1
            !                               ^
            !   o - cell   | - edge
            !
            gy(i,j,k) = sum(f(i,j-nsy:j+nsy-1,k) * wy(:2*nsy))
          end do
#else
          if (mesh%is_outside_pole_full_lat(j-nsy+1)) then
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              io = i + mesh%num_full_lon / 2
              if (io > mesh%num_full_lon) io = io - mesh%num_full_lon
              f(i,j-nsy+1,k) = f(io,j,k)
              gy(i,j,k) = sum(f(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
            end do
          else if (mesh%is_outside_pole_full_lat(j+nsy)) then
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              io = i + mesh%num_full_lon / 2
              if (io > mesh%num_full_lon) io = io - mesh%num_full_lon
              f(i,j+nsy,k) = f(io,j,k)
              gy(i,j,k) = sum(f(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
            end do
          else
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
              !
              !                               j
              !   |      o      |      o      |      o      |      o      |
              !         j-1            j            j+1           j+2
              !                               ^
              !   o - cell   | - edge
              !
              gy(i,j,k) = sum(f(i,j-nsy+1:j+nsy,k) * wy(:2*nsy))
            end do
          end if
#endif
          ! Limit damping flux to avoid upgradient (Xue 2000).
          if (oy > 2) then
            gy(:,j,k) = gy(:,j,k) * (-1)**(oy / 2)
            do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j,k) - f(i,j-1,k))))
#else
              gy(i,j,k) = gy(i,j,k) * max(0.0_r8, sign(1.0_r8, -gy(i,j,k) * (f(i,j+1,k) - f(i,j,k))))
#endif
            end do
          end if
        end do
      end do
#ifdef V_POLE
      call fill_halo(block, gy, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false.)
#else
      call fill_halo(block, gy, full_lon=.true., full_lat=.false., full_lev=.true., north_halo=.false.)
#endif
      ! Damp physical variable at last.
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
          ox = ox_full_lat(j)
          do i = mesh%full_lon_ibeg, mesh%full_lon_iend
#ifdef V_POLE
            !
            !               - j+1
            !               |
            !               |
            !        o------|------o
            !       i-1     |      i
            !               |
            !               - j
            !
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (    &
              (gx(i,j,k) - gx(i-1,j,k)) * cx_full_lat(j,k,ox) + &
              (gy(i,j+1,k) - gy(i,j,k)) * cy_full_lat(j,k,oy)   &
            )
#else
            !
            !               - j
            !               |
            !               |
            !        o------|------o
            !       i-1     |      i
            !               |
            !               - j-1
            !
            f(i,j,k) = f(i,j,k) - dt / polar_damp_cycles * (    &
              (gx(i,j,k) - gx(i-1,j,k)) * cx_full_lat(j,k,ox) + &
              (gy(i,j,k) - gy(i,j-1,k)) * cy_full_lat(j,k,oy)   &
            )
#endif
          end do
        end do
      end do
      call fill_halo(block, f, full_lon=.true., full_lat=.true., full_lev=.true.)
    end do cycle_loop

  end subroutine latlon_damp_cell

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
        cycle_loop: do cyc = 1, div_damp_cycles
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
        end do cycle_loop
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

  subroutine shapiro_smooth(block, state)

    type(block_type), intent(in) :: block
    type(state_type), intent(inout) :: state

    type(mesh_type), pointer :: mesh
    integer i, j, k

    mesh => state%mesh
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          state%u(i,j,k) = state%u(i,j,k) + sx_full_lat(j,k) * ( &
            -6 * state%u(i  ,j,k) + &
             4 * state%u(i-1,j,k) + &
             4 * state%u(i+1,j,k) + &
            -1 * state%u(i-2,j,k) + &
            -1 * state%u(i+2,j,k)   &
          ) / 16.0_r8
        end do
      end do
    end do
    call fill_halo(block, state%u, full_lon=.false., full_lat=.true., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg_no_pole, mesh%half_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%v(i,j,k) = state%v(i,j,k) + sx_half_lat(j,k) * ( &
            -6 * state%v(i  ,j,k) + &
             4 * state%v(i-1,j,k) + &
             4 * state%v(i+1,j,k) + &
            -1 * state%v(i-2,j,k) + &
            -1 * state%v(i+2,j,k)   &
          ) / 16.0_r8
        end do
      end do
    end do
    call fill_halo(block, state%v, full_lon=.true., full_lat=.false., full_lev=.true.)

    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          state%pt(i,j,k) = state%pt(i,j,k) + sx_full_lat(j,k) * ( &
            -6 * state%pt(i  ,j,k) + &
             4 * state%pt(i-1,j,k) + &
             4 * state%pt(i+1,j,k) + &
            -1 * state%pt(i-2,j,k) + &
            -1 * state%pt(i+2,j,k)   &
          ) / 16.0_r8
        end do
      end do
    end do
    call fill_halo(block, state%pt, full_lon=.true., full_lat=.true., full_lev=.true.)

  end subroutine shapiro_smooth

end module damp_mod
