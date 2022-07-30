module moist_mod

  use namelist_mod
  use block_mod
  use adv_mod

  implicit none

  private

  public moist_init
  public calc_qm

contains

  subroutine moist_init()

    call adv_add_tracer('moist', dt_adv, 'qv', 'water vapor' )
    ! call adv_add_tracer('moist', dt_adv, 'ql', 'cloud liquid')
    ! call adv_add_tracer('moist', dt_adv, 'qi', 'cloud ice'   )
    ! call adv_add_tracer('moist', dt_adv, 'qr', 'rain'        )
    ! call adv_add_tracer('moist', dt_adv, 'qs', 'snow'        )
    ! call adv_add_tracer('moist', dt_adv, 'qg', 'graupels'    )
    ! call adv_add_tracer('moist', dt_adv, 'qh', 'hail'        )

  end subroutine moist_init

  subroutine calc_qm(block, itime)

    type(block_type), intent(inout), target :: block
    integer, intent(in) :: itime

    real(r8), pointer, dimension(:,:,:) :: qv
    integer i, j, k

    associate (mesh => block%mesh,            &
               qm   => block%state(itime)%qm)   ! out
    qv => block%adv_batches(1)%q(:,:,:,1,block%adv_batches(1)%old)
    do k = mesh%full_lev_ibeg, mesh%full_lev_iend
      do j = mesh%full_lat_ibeg, mesh%full_lat_iend
        do i = mesh%full_lon_ibeg, mesh%full_lon_iend
          qm(i,j,k) = qv(i,j,k)
        end do
      end do
    end do
    end associate

  end subroutine calc_qm

end module moist_mod
