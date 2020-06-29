module hybrid_coord_mod

  use const_mod

  implicit none

  private

  public hybrid_coord_init
  public hybrid_coord_final
  public hybrid_coord_calc_ph_lev

  real(r8), allocatable, dimension(:) :: hyai
  real(r8), allocatable, dimension(:) :: hybi
  real(r8), allocatable, dimension(:) :: hyam
  real(r8), allocatable, dimension(:) :: hybm

  real(r8) :: p0 = 1d5 ! Reference pressure (Pa)

  namelist /hybrid_coord_params/ &
    hyai, hybi, hyam, hybm, p0

contains

  subroutine hybrid_coord_init(num_lev, namelist_path)

    integer, intent(in) :: num_lev
    character(*), intent(in) :: namelist_path

    if (allocated(hyai)) deallocate(hyai)
    if (allocated(hybi)) deallocate(hybi)
    if (allocated(hyam)) deallocate(hyam)
    if (allocated(hybm)) deallocate(hybm)

    allocate(hyai(num_lev+1))
    allocate(hybi(num_lev+1))
    allocate(hyam(num_lev  ))
    allocate(hybm(num_lev  ))

    open(10, file=namelist_path, status='old')
    read(10, nml=hybrid_coord_params)
    close(10)

  end subroutine hybrid_coord_init

  subroutine hybrid_coord_final()

    deallocate(hyai)
    deallocate(hybi)
    deallocate(hyam)
    deallocate(hybm)

  end subroutine hybrid_coord_final

  real(r8) function hybrid_coord_calc_ph_lev(k, phs) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: phs

    res = hyai(k) * p0 + hybi(k) * phs

  end function hybrid_coord_calc_ph_lev

end module hybrid_coord_mod
