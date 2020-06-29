module vert_coord_mod

  use const_mod
  use hybrid_coord_mod

  implicit none

  real(r8), allocatable, dimension(:) :: eta
  real(r8), allocatable, dimension(:) :: eta_lev

  private

  public vert_coord_init
  public vert_coord_final
  public vert_coord_calc_ph_lev

contains

  subroutine vert_coord_init(num_lev, namelist_path)

    integer, intent(in) :: num_lev
    character(*), intent(in) :: namelist_path

    if (allocated(eta    )) deallocate(eta    )
    if (allocated(eta_lev)) deallocate(eta_lev)

    allocate(eta    (num_lev  ))
    allocate(eta_lev(num_lev+1))

    call hybrid_coord_mod(namelist_path)

  end subroutine vert_coord_init

  subroutine vert_coord_final()

    if (allocated(eta    )) deallocate(eta    )
    if (allocated(eta_lev)) deallocate(eta_lev)

  end subroutine vert_coord_final

  real(r8) function vert_coord_calc_ph_lev(k, phs) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: phs

    res = hybrid_coord_calc_ph_lev(k, phs)

  end function vert_coord_calc_ph_lev

end module vert_coord_mod
