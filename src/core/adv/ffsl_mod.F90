module ffsl_mod

  use const_mod
  use namelist_mod
  use block_mod
  use parallel_mod
  use adv_batch_mod

  implicit none

  private

  public ffsl_init
  public ffsl_calc_mass_hflx_cell
  public ffsl_calc_mass_hflx_vtx
  public ffsl_calc_tracer_hflx_cell
  public ffsl_calc_tracer_hflx_vtx
  public ffsl_calc_tracer_hval_vtx

  interface
    subroutine hflx_cell_interface(block, batch, u, v, mx, my, mfx, mfy)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: mx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: my (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    end subroutine hflx_cell_interface
    subroutine hflx_vtx_interface(block, batch, u, v, mx, my, mfx, mfy)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: u  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: v  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: mx (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: my (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mfy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    end subroutine hflx_vtx_interface
    subroutine hval_vtx_interface(block, batch, mx, my, mvx, mvy)
      import block_type, adv_batch_type, r8
      type(block_type    ), intent(in   ) :: block
      type(adv_batch_type), intent(inout) :: batch
      real(r8), intent(in ) :: mx (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(in ) :: my (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mvx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
      real(r8), intent(out) :: mvy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    end subroutine hval_vtx_interface
    pure real(r8) function slope_interface(fm1, f, fp1)
      import r8
      real(r8), intent(in) :: fm1, f, fp1
    end function slope_interface
  end interface

  procedure(hflx_cell_interface), pointer :: hflx_cell => null()
  procedure(hflx_vtx_interface ), pointer :: hflx_vtx  => null()
  procedure(hval_vtx_interface ), pointer :: hval_vtx  => null()
  procedure(slope_interface    ), pointer :: slope     => null()

contains

  subroutine ffsl_init()

    select case (ffsl_flux_type)
    case ('van_leer')
      hflx_cell => hflx_van_leer_cell
      hflx_vtx  => hflx_van_leer_vtx
      hval_vtx  => hval_van_leer_vtx
    case ('ppm')
      hflx_cell => hflx_ppm_cell
      hflx_vtx  => hflx_ppm_vtx
      hval_vtx  => hval_ppm_vtx
    end select

    select case (limiter_type)
    case ('none')
      slope => slope_simple
    case ('mono')
      slope => slope_mono
    case ('pd')
      slope => slope_pd
    end select

  end subroutine ffsl_init

  subroutine ffsl_calc_mass_hflx_cell(block, batch, m, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: m  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k
    real(r8) work(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend,block%mesh%num_full_lev)
    real(r8) pole(block%mesh%num_full_lev)

    associate (mesh => block%mesh, &
               dt   => batch%dt  , & ! in
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               mx   => batch%qx  , & ! work array
               my   => batch%qy)     ! work array
    ! Run inner advective operators.
    call hflx_cell(block, batch, u, v, m, m, mfx, mfy)
    call fill_halo(block, mfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ! Subtract divergence terms from flux to form advective operators.
          mx(i,j,k) = m(i,j,k) - 0.5_r8 * (          &
            (                                        &
              mfx(i,j,k) - mfx(i-1,j,k)              &
            ) * mesh%le_lon(j) / mesh%area_cell(j) - &
            divx(i,j,k) * m(i,j,k)                   &
          ) * dt
          my(i,j,k) = m(i,j,k) - 0.5_r8 * (     &
            (                                   &
              mfy(i,j  ,k) * mesh%le_lat(j  ) - &
              mfy(i,j-1,k) * mesh%le_lat(j-1)   &
            ) / mesh%area_cell(j) -             &
            divy(i,j,k) * m(i,j,k)              &
          ) * dt
        end do
      end do
    end do
    ! Handle the Pole boundary conditions.
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = mfy(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          mx(i,j,k) = m(i,j,k)
          my(i,j,k) = m(i,j,k) - 0.5_r8 * (pole(k) - divy(i,j,k) * m(i,j,k)) * dt
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = mfy(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          mx(i,j,k) = m(i,j,k)
          my(i,j,k) = m(i,j,k) + 0.5_r8 * (pole(k) - divy(i,j,k) * m(i,j,k)) * dt
        end do
      end do
    end if
    call fill_halo(block, mx, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, my, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx_cell(block, batch, u, v, my, mx, mfx, mfy)
    end associate

  end subroutine ffsl_calc_mass_hflx_cell

  subroutine ffsl_calc_mass_hflx_vtx(block, batch, m, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: m  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k

    associate (mesh => block%mesh, &
               dt   => batch%dt  , & ! in
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               mx   => batch%qx  , & ! work array
               my   => batch%qy)     ! work array
    ! Run inner advective operators.
    call hflx_vtx(block, batch, u, v, m, m, mfx, mfy)
    call fill_halo(block, mfx, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mfy, full_lon=.false., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ! Subtract divergence terms from flux to form advective operators.
          mx(i,j,k) = m(i,j,k) - 0.5_r8 * (          &
            (                                        &
              mfx(i+1,j,k) - mfx(i,j,k)              &
            ) * mesh%de_lat(j) / mesh%area_vtx(j) -  &
            divx(i,j,k) * m(i,j,k)                   &
          ) * dt
          my(i,j,k) = m(i,j,k) - 0.5_r8 * (     &
            (                                   &
              mfy(i,j+1,k) * mesh%de_lon(j+1) - &
              mfy(i,j  ,k) * mesh%de_lon(j  )   &
            ) / mesh%area_vtx(j) -              &
            divy(i,j,k) * m(i,j,k)              &
          ) * dt
        end do
      end do
    end do
    call fill_halo(block, mx, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, my, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx_vtx(block, batch, u, v, my, mx, mfx, mfy)
    end associate

  end subroutine ffsl_calc_mass_hflx_vtx

  subroutine ffsl_calc_tracer_hflx_cell(block, batch, q, qmfx, qmfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: q    (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qmfx (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qmfy (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k
    real(r8) work(block%mesh%full_lon_ibeg:block%mesh%full_lon_iend,block%mesh%num_full_lev)
    real(r8) pole(block%mesh%num_full_lev)

    associate (mesh => block%mesh, &
               dt   => batch%dt  , & ! in
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               mfx  => batch%mfx , & ! in
               mfy  => batch%mfy , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               qx   => batch%qx  , & ! work array
               qy   => batch%qy)     ! work array
    ! Run inner advective operators.
    call hflx_cell(block, batch, u, v, q, q, qmfx, qmfy)
    call fill_halo(block, qmfx, full_lon=.false., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, qmfy, full_lon=.true., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ! Subtract divergence terms from flux to form advective operators.
          qx(i,j,k) = q(i,j,k) - 0.5_r8 * (          &
            (                                        &
              qmfx(i,j,k) - qmfx(i-1,j,k)            &
            ) * mesh%le_lon(j) / mesh%area_cell(j) - &
            divx(i,j,k) * q(i,j,k)                   &
          ) * dt
          qy(i,j,k) = q(i,j,k) - 0.5_r8 * (      &
            (                                    &
              qmfy(i,j  ,k) * mesh%le_lat(j  ) - &
              qmfy(i,j-1,k) * mesh%le_lat(j-1)   &
            ) / mesh%area_cell(j) -              &
            divy(i,j,k) * q(i,j,k)               &
          ) * dt
        end do
      end do
    end do
    ! Handle the Pole boundary conditions.
    if (mesh%has_south_pole()) then
      j = mesh%full_lat_ibeg
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = qmfy(i,j,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          qx(i,j,k) = q(i,j,k)
          qy(i,j,k) = q(i,j,k) - 0.5_r8 * (pole(k) - divy(i,j,k) * q(i,j,k)) * dt
        end do
      end do
    end if
    if (mesh%has_north_pole()) then
      j = mesh%full_lat_iend
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          work(i,k) = qmfy(i,j-1,k)
        end do
      end do
      call zonal_sum(proc%zonal_circle, work, pole)
      pole = pole * mesh%le_lat(j-1) / global_mesh%num_full_lon / mesh%area_cell(j)
      do k = mesh%full_lev_ibeg, mesh%full_lev_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          qx(i,j,k) = q(i,j,k)
          qy(i,j,k) = q(i,j,k) + 0.5_r8 * (pole(k) - divy(i,j,k) * q(i,j,k)) * dt
        end do
      end do
    end if
    call fill_halo(block, qx, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, qy, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx_cell(block, batch, mfx, mfy, qy, qx, qmfx, qmfy)
    end associate

  end subroutine ffsl_calc_tracer_hflx_cell

  subroutine ffsl_calc_tracer_hflx_vtx(block, batch, q, qmfx, qmfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: q    (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qmfx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                   block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qmfy (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                   block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                   block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k

    associate (mesh => block%mesh, &
               dt   => batch%dt  , & ! in
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               mfx  => batch%mfx , & ! in
               mfy  => batch%mfy , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               qx   => batch%qx  , & ! work array
               qy   => batch%qy)     ! work array
    ! Run inner advective operators.
    call hflx_vtx(block, batch, u, v, q, q, qmfx, qmfy)
    call fill_halo(block, mfx, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mfy, full_lon=.false., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ! Subtract divergence terms from flux to form advective operators.
          qx(i,j,k) = q(i,j,k) - 0.5_r8 * (         &
            (                                       &
              qmfx(i+1,j,k) - qmfx(i,j,k)           &
            ) * mesh%de_lat(j) / mesh%area_vtx(j) - &
            divx(i,j,k) * q(i,j,k)                  &
          ) * dt
          qy(i,j,k) = q(i,j,k) - 0.5_r8 * (      &
            (                                    &
              qmfy(i,j+1,k) * mesh%de_lon(j+1) - &
              qmfy(i,j  ,k) * mesh%de_lon(j  )   &
            ) / mesh%area_vtx(j) -               &
            divy(i,j,k) * q(i,j,k)               &
          ) * dt
        end do
      end do
    end do
    call fill_halo(block, qx, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, qy, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hflx_vtx(block, batch, mfx, mfy, qy, qx, qmfx, qmfy)
    end associate

  end subroutine ffsl_calc_tracer_hflx_vtx

  subroutine ffsl_calc_tracer_hval_vtx(block, batch, q, qvx, qvy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: q  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qvx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: qvy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k

    associate (mesh => block%mesh, &
               dt   => batch%dt  , & ! in
               u    => batch%u   , & ! in
               v    => batch%v   , & ! in
               divx => batch%divx, & ! in
               divy => batch%divy, & ! in
               qx   => batch%qx  , & ! work array
               qy   => batch%qy)     ! work array
    ! Run inner advective operators.
    call hflx_vtx(block, batch, u, v, q, q, qvx, qvy)
    call fill_halo(block, qvx, full_lon=.true., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, qvy, full_lon=.false., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    ! Calculate intermediate tracer density due to advective operators.
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ! Subtract divergence terms from flux to form advective operators.
          qx(i,j,k) = q(i,j,k) - 0.5_r8 * (         &
            (                                       &
              qvx(i+1,j,k) - qvx(i,j,k)             &
            ) * mesh%de_lat(j) / mesh%area_vtx(j) - &
            divx(i,j,k) * q(i,j,k)                  &
          ) * dt
          qy(i,j,k) = q(i,j,k) - 0.5_r8 * (     &
            (                                   &
              qvy(i,j+1,k) * mesh%de_lon(j+1) - &
              qvy(i,j  ,k) * mesh%de_lon(j  )   &
            ) / mesh%area_vtx(j) -              &
            divy(i,j,k) * q(i,j,k)              &
          ) * dt
        end do
      end do
    end do
    call fill_halo(block, qx, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, qy, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    ! Run outer flux form operators.
    call hval_vtx(block, batch, qy, qx, qvx, qvy)
    end associate

  end subroutine ffsl_calc_tracer_hval_vtx

  subroutine hflx_van_leer_cell(block, batch, u, v, mx, my, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: mx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: my (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k, iu, ju, ci
    real(r8) cf, dm

    associate (mesh => block %mesh, &
               cflx => batch %cflx, & ! in
               cfly => batch %cfly)   ! in
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          if (cflx(i,j,k) > 0) then
            iu = i - ci
            dm = slope(mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k))
            mfx(i,j,k) = u(i,j,k) * (cf * (mx(iu,j,k) + dm * 0.5_r8 * (1 - cf)) + sum(mx(i+1-ci:i,j,k))) / cflx(i,j,k)
          else if (cflx(i,j,k) < 0) then
            iu = i - ci + 1
            dm = slope(mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k))
            mfx(i,j,k) = u(i,j,k) * (cf * (mx(iu,j,k) - dm * 0.5_r8 * (1 + cf)) - sum(mx(i+1:i-ci,j,k))) / cflx(i,j,k)
          else
            mfx(i,j,k) = 0
          end if
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          cf = cfly(i,j,k)
          ju = merge(j, j + 1, cf > 0)
          dm = slope(my(i,ju-1,k), my(i,ju,k), my(i,ju+1,k))
          mfy(i,j,k) = v(i,j,k) * (my(i,ju,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf))
        end do
      end do
    end do
    end associate

  end subroutine hflx_van_leer_cell

  subroutine hflx_van_leer_vtx(block, batch, u, v, mx, my, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: u  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: v  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: mx (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: my (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k, iu, ju, ci
    real(r8) cf, dm

    associate (mesh => block%mesh, &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly)   ! in
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          if (cflx(i,j,k) > 0) then
            iu = i - ci - 1
            dm = slope(mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k))
            mfx(i,j,k) = u(i,j,k) * (cf * (mx(iu,j,k) + dm * 0.5_r8 * (1 - cf)) + sum(mx(i-ci:i-1,j,k))) / cflx(i,j,k)
          else if (cflx(i,j,k) < 0) then
            iu = i - ci
            dm = slope(mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k))
            mfx(i,j,k) = u(i,j,k) * (cf * (mx(iu,j,k) - dm * 0.5_r8 * (1 + cf)) - sum(mx(i:i-ci-1,j,k))) / cflx(i,j,k)
          else
            mfx(i,j,k) = 0
          end if
        end do
      end do
      ! Along y-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          cf = cfly(i,j,k)
          ju = merge(j - 1, j, cf > 0)
          dm = slope(my(i,ju-1,k), my(i,ju,k), my(i,ju+1,k))
          mfy(i,j,k) = v(i,j,k) * (my(i,ju,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf))
        end do
      end do
    end do
    end associate

  end subroutine hflx_van_leer_vtx

  subroutine hval_van_leer_vtx(block, batch, mx, my, mvx, mvy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: mx (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: my (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mvx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mvy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k, iu, ju, ci
    real(r8) cf, dm

    associate (mesh => block%mesh, &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly)   ! in
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          if (cflx(i,j,k) > 0) then
            iu = i - ci - 1
            dm = slope(mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k))
            mvx(i,j,k) = (cf * (mx(iu,j,k) + dm * 0.5_r8 * (1 - cf)) + sum(mx(i-ci:i-1,j,k))) / cflx(i,j,k)
          else if (cflx(i,j,k) < 0) then
            iu = i - ci
            dm = slope(mx(iu-1,j,k), mx(iu,j,k), mx(iu+1,j,k))
            mvx(i,j,k) = (cf * (mx(iu,j,k) - dm * 0.5_r8 * (1 + cf)) - sum(mx(i:i-ci-1,j,k))) / cflx(i,j,k)
          else
            mvx(i,j,k) = (mx(i-1,j,k) + mx(i,j,k)) * 0.5_r8
          end if
        end do
      end do
      ! Along y-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          cf = cfly(i,j,k)
          ju = merge(j - 1, j, cf > 0)
          dm = slope(my(i,ju-1,k), my(i,ju,k), my(i,ju+1,k))
          mvy(i,j,k) = my(i,ju,k) + dm * 0.5_r8 * (sign(1.0_r8, cf) - cf)
        end do
      end do
    end do
    end associate

  end subroutine hval_van_leer_vtx

  subroutine hflx_ppm_cell(block, batch, u, v, mx, my, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: u  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: v  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: mx (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: my (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k, iu, ju, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3

    associate (mesh => block%mesh, &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly, & ! in
               mlx  => batch%qlx , & ! work array
               mly  => batch%qly , & ! work array
               dmx  => batch%dqx , & ! work array
               dmy  => batch%dqy , & ! work array
               m6x  => batch%q6x , & ! work array
               m6y  => batch%q6y )   ! work array
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          call ppm(mx(i-2,j,k), mx(i-1,j,k), mx(i,j,k), mx(i+1,j,k), mx(i+2,j,k), mlx(i,j,k), dmx(i,j,k), m6x(i,j,k))
          call ppm(my(i,j-2,k), my(i,j-1,k), my(i,j,k), my(i,j+1,k), my(i,j+2,k), mly(i,j,k), dmy(i,j,k), m6y(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, mlx, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, dmx, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, m6x, full_lon=.true., full_lat=.true., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mly, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, dmy, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, m6y, full_lon=.true., full_lat=.true., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          if (cflx(i,j,k) > 0) then
            iu = i - ci
            s1 = 1 - cf
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfx(i,j,k) =  u(i,j,k) * (sum(mx(i+1-ci:i,j,k)) + mlx(iu,j,k) * ds1 + 0.5_r8 * dmx(iu,j,k) * ds2 + m6x(iu,j,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx(i,j,k)
          else if (cflx(i,j,k) < 0) then
            iu = i - ci + 1
            s1 = 0
            s2 = -cf
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfx(i,j,k) = -u(i,j,k) * (sum(mx(i+1:i-ci,j,k)) + mlx(iu,j,k) * ds1 + 0.5_r8 * dmx(iu,j,k) * ds2 + m6x(iu,j,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx(i,j,k)
          else
            mfx(i,j,k) = 0
          end if
        end do
      end do
      ! Along y-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          if (cfly(i,j,k) > 0) then
            ju = j
            s1 = 1 - cfly(i,j,k)
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfy(i,j,k) =  v(i,j,k) * (mly(i,ju,k) * ds1 + 0.5_r8 * dmy(i,ju,k) * ds2 + m6y(i,ju,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly(i,j,k)
          else if (cfly(i,j,k) < 0) then
            ju = j + 1
            s1 = 0
            s2 = -cfly(i,j,k)
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfy(i,j,k) = -v(i,j,k) * (mly(i,ju,k) * ds1 + 0.5_r8 * dmy(i,ju,k) * ds2 + m6y(i,ju,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly(i,j,k)
          else
            mfy(i,j,k) = 0
          end if
        end do
      end do
    end do
    end associate

  end subroutine hflx_ppm_cell

  subroutine hflx_ppm_vtx(block, batch, u, v, mx, my, mfx, mfy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: u  (block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: v  (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: mx (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: my (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mfy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k, iu, ju, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3

    associate (mesh => block%mesh, &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly, & ! in
               mlx  => batch%qlx , & ! work array
               mly  => batch%qly , & ! work array
               dmx  => batch%dqx , & ! work array
               dmy  => batch%dqy , & ! work array
               m6x  => batch%q6x , & ! work array
               m6y  => batch%q6y )   ! work array
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          call ppm(mx(i-2,j,k), mx(i-1,j,k), mx(i,j,k), mx(i+1,j,k), mx(i+2,j,k), mlx(i,j,k), dmx(i,j,k), m6x(i,j,k))
          call ppm(my(i,j-2,k), my(i,j-1,k), my(i,j,k), my(i,j+1,k), my(i,j+2,k), mly(i,j,k), dmy(i,j,k), m6y(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, mlx, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, dmx, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, m6x, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mly, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, dmy, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, m6y, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          if (cflx(i,j,k) > 0) then
            iu = i - ci - 1
            s1 = 1 - cf
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfx(i,j,k) =  u(i,j,k) * (sum(mx(i-ci:i-1,j,k)) + mlx(iu,j,k) * ds1 + 0.5_r8 * dmx(iu,j,k) * ds2 + m6x(iu,j,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx(i,j,k)
          else if (cflx(i,j,k) < 0) then
            iu = i - ci
            s1 = 0
            s2 = -cf
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfx(i,j,k) = -u(i,j,k) * (sum(mx(i:i-ci-1,j,k)) + mlx(iu,j,k) * ds1 + 0.5_r8 * dmx(iu,j,k) * ds2 + m6x(iu,j,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx(i,j,k)
          else
            mfx(i,j,k) = 0
          end if
        end do
      end do
      ! Along y-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          if (cfly(i,j,k) > 0) then
            ju = j - 1
            s1 = 1 - cfly(i,j,k)
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfy(i,j,k) =  v(i,j,k) * (mly(i,ju,k) * ds1 + 0.5_r8 * dmy(i,ju,k) * ds2 + m6y(i,ju,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly(i,j,k)
          else if (cfly(i,j,k) < 0) then
            ju = j
            s1 = 0
            s2 = -cfly(i,j,k)
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mfy(i,j,k) = -v(i,j,k) * (mly(i,ju,k) * ds1 + 0.5_r8 * dmy(i,ju,k) * ds2 + m6y(i,ju,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly(i,j,k)
          else
            mfy(i,j,k) = 0
          end if
        end do
      end do
    end do
    end associate

  end subroutine hflx_ppm_vtx

  subroutine hval_ppm_vtx(block, batch, mx, my, mvx, mvy)

    type(block_type    ), intent(in   ) :: block
    type(adv_batch_type), intent(inout) :: batch
    real(r8), intent(in ) :: mx (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(in ) :: my (block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mvx(block%mesh%full_lon_lb:block%mesh%full_lon_ub, &
                                 block%mesh%half_lat_lb:block%mesh%half_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)
    real(r8), intent(out) :: mvy(block%mesh%half_lon_lb:block%mesh%half_lon_ub, &
                                 block%mesh%full_lat_lb:block%mesh%full_lat_ub, &
                                 block%mesh%full_lev_lb:block%mesh%full_lev_ub)

    integer i, j, k, iu, ju, ci
    real(r8) cf, s1, s2, ds1, ds2, ds3

    associate (mesh => block%mesh, &
               cflx => batch%cflx, & ! in
               cfly => batch%cfly, & ! in
               mlx  => batch%qlx , & ! work array
               mly  => batch%qly , & ! work array
               dmx  => batch%dqx , & ! work array
               dmy  => batch%dqy , & ! work array
               m6x  => batch%q6x , & ! work array
               m6y  => batch%q6y )   ! work array
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          call ppm(mx(i-2,j,k), mx(i-1,j,k), mx(i,j,k), mx(i+1,j,k), mx(i+2,j,k), mlx(i,j,k), dmx(i,j,k), m6x(i,j,k))
          call ppm(my(i,j-2,k), my(i,j-1,k), my(i,j,k), my(i,j+1,k), my(i,j+2,k), mly(i,j,k), dmy(i,j,k), m6y(i,j,k))
        end do
      end do
    end do
    call fill_halo(block, mlx, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, dmx, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, m6x, full_lon=.false., full_lat=.false., full_lev=.true., south_halo=.false., north_halo=.false.)
    call fill_halo(block, mly, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, dmy, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    call fill_halo(block, m6y, full_lon=.false., full_lat=.false., full_lev=.true.,  west_halo=.false.,  east_halo=.false.)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      ! Along x-axis
      do j = mesh%half_lat_ibeg, mesh%half_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          ci = int(cflx(i,j,k))
          cf = cflx(i,j,k) - ci
          if (cflx(i,j,k) > 0) then
            iu = i - ci - 1
            s1 = 1 - cf
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mvx(i,j,k) =  (sum(mx(i-ci:i-1,j,k)) + mlx(iu,j,k) * ds1 + 0.5_r8 * dmx(iu,j,k) * ds2 + m6x(iu,j,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx(i,j,k)
          else if (cflx(i,j,k) < 0) then
            iu = i - ci
            s1 = 0
            s2 = -cf
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mvx(i,j,k) = -(sum(mx(i:i-ci-1,j,k)) + mlx(iu,j,k) * ds1 + 0.5_r8 * dmx(iu,j,k) * ds2 + m6x(iu,j,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cflx(i,j,k)
          else
            mvx(i,j,k) = (mx(i-1,j,k) + mx(i,j,k)) * 0.5_r8
          end if
        end do
      end do
      ! Along y-axis
      do j = mesh%full_lat_ibeg_no_pole, mesh%full_lat_iend_no_pole
        do i = mesh%half_lon_ibeg, mesh%half_lon_iend
          if (cfly(i,j,k) > 0) then
            ju = j - 1
            s1 = 1 - cfly(i,j,k)
            s2 = 1
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mvy(i,j,k) =  (mly(i,ju,k) * ds1 + 0.5_r8 * dmy(i,ju,k) * ds2 + m6y(i,ju,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly(i,j,k)
          else if (cfly(i,j,k) < 0) then
            ju = j
            s1 = 0
            s2 = -cfly(i,j,k)
            ds1 = s2    - s1
            ds2 = s2**2 - s1**2
            ds3 = s2**3 - s1**3
            mvy(i,j,k) = -(mly(i,ju,k) * ds1 + 0.5_r8 * dmy(i,ju,k) * ds2 + m6y(i,ju,k) * (ds2 / 2.0_r8 - ds3 / 3.0_r8)) / cfly(i,j,k)
          else
            mvy(i,j,k) = (my(i,j-1,k) + my(i,j,k)) * 0.5_r8
          end if
        end do
      end do
    end do
    end associate

  end subroutine hval_ppm_vtx  

  subroutine ppm(fm2, fm1, f, fp1, fp2, fl, df, f6)

    real(r8), intent(in ) :: fm2
    real(r8), intent(in ) :: fm1
    real(r8), intent(in ) :: f
    real(r8), intent(in ) :: fp1
    real(r8), intent(in ) :: fp2
    real(r8), intent(out) :: fl
    real(r8), intent(out) :: df
    real(r8), intent(out) :: f6

    real(r8) dfl, dfr, fr

    ! Calculate values at left and right cell interfaces.
    dfl = slope(fm2, fm1, f  )
    df  = slope(fm1, f  , fp1)
    dfr = slope(f  , fp1, fp2)
    ! Why (B2) in Lin (2004) divide (dfl - df) and (df - dfr) by 3?
    fl = 0.5_r8 * (fm1 + f) + (dfl - df) / 6.0_r8
    fr = 0.5_r8 * (fp1 + f) + (df - dfr) / 6.0_r8
    ! Why (B3) and (B4) in Lin (2004) multiply df by 2?
    fl = f - sign(min(abs(df), abs(fl - f)), df)
    fr = f + sign(min(abs(df), abs(fr - f)), df)
    f6 = 6 * f - 3 * (fl + fr)
    df = fr - fl

  end subroutine ppm

  pure real(r8) function slope_simple(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    res = (fp1 - fm1) * 0.5_r8

  end function slope_simple

  pure real(r8) function slope_mono(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    real(r8) df, df_min, df_max

    df = (fp1 - fm1) * 0.5_r8 ! Initial guess
    df_min = 2 * (f - min(fm1, f, fp1))
    df_max = 2 * (max(fm1, f, fp1) - f)
    res = sign(min(abs(df), df_min, df_max), df)

  end function slope_mono

  pure real(r8) function slope_pd(fm1, f, fp1) result(res)

    real(r8), intent(in) :: fm1, f, fp1

    real(r8) df

    df = (fp1 - fm1) * 0.5_r8 ! Initial guess
    res = sign(min(abs(df), 2 * f), df)

  end function slope_pd

end module ffsl_mod
