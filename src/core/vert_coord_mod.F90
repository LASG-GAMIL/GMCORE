module vert_coord_mod

  use const_mod
  use sigma_coord_mod
  use hybrid_coord_mod

  implicit none

  private

  public vert_coord_init
  public vert_coord_final
  public vert_coord_calc_ph_lev
  public vert_coord_calc_dphdt_lev

contains

  subroutine vert_coord_init(num_lev, namelist_path)

    integer, intent(in) :: num_lev
    character(*), intent(in) :: namelist_path

#ifdef SIGMA
    call sigma_coord_init(num_lev, namelist_path)
#else
    call hybrid_coord_init(num_lev, namelist_path)
#endif

  end subroutine vert_coord_init

  subroutine vert_coord_final()

#ifdef SIGMA
    call sigma_coord_final()
#else
    call hybrid_coord_final()
#endif

  end subroutine vert_coord_final

  pure real(r8) function vert_coord_calc_ph_lev(k, phs) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: phs

#ifdef SIGMA
    res = sigma_coord_calc_ph_lev(k, phs)
#else
    res = hybrid_coord_calc_ph_lev(k, phs)
#endif

  end function vert_coord_calc_ph_lev

  pure real(r8) function vert_coord_calc_dphdt_lev(k, dphsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dphsdt

#ifdef SIGMA
    res = sigma_coord_calc_dphdt_lev(k, dphsdt)
#else
    res = hybrid_coord_calc_dphdt_lev(k, dphsdt)
#endif

  end function vert_coord_calc_dphdt_lev

end module vert_coord_mod
