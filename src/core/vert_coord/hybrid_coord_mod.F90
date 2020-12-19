module hybrid_coord_mod

  use flogger
  use namelist_mod
  use const_mod
  use hybrid_coord_test_mod
  use hybrid_coord_ecmwf_mod
  use mesh_mod
  use process_mod

  implicit none

  private

  public hybrid_coord_init
  public hybrid_coord_final
  public hybrid_coord_calc_ph
  public hybrid_coord_calc_ph_lev
  public hybrid_coord_calc_dphdt_lev

  real(r8), allocatable, dimension(:) :: hyai
  real(r8), allocatable, dimension(:) :: hybi
  real(r8), allocatable, dimension(:) :: hyam
  real(r8), allocatable, dimension(:) :: hybm

  character(10) :: template = '' ! Template:
                                 ! - ecmwf_l50

  real(r8) :: p0 = 1d5 ! Reference pressure (Pa)

  namelist /hybrid_coord/ &
    hyai, hybi, hyam, hybm, p0, template

contains

  subroutine hybrid_coord_init(num_lev, namelist_file, template_)

    integer, intent(in) :: num_lev
    character(*), intent(in), optional :: namelist_file
    character(*), intent(in), optional :: template_

    integer ierr, k

    if (allocated(hyai)) deallocate(hyai)
    if (allocated(hybi)) deallocate(hybi)
    if (allocated(hyam)) deallocate(hyam)
    if (allocated(hybm)) deallocate(hybm)

    allocate(hyai(num_lev+1)); hyai = 0
    allocate(hybi(num_lev+1)); hybi = 0
    allocate(hyam(num_lev  )); hyam = 0
    allocate(hybm(num_lev  )); hybm = 0

    if (present(namelist_file)) then
      open(10, file=namelist_file, status='old')
      read(10, nml=hybrid_coord, iostat=ierr)
      close(10)
    else if (present(template_)) then
      template = template_
    end if

    if (ierr /= 0) then
      if (.not. baroclinic) then
        if (is_root_proc()) call log_notice('Run shallow-water model.')
        return
      else
        call log_error('No hybrid_coord parameters in ' // trim(namelist_file) // '!')
      end if
    end if

    select case (template)
    case ('test_l15')
      call hybrid_coord_test_l15(p0, hyai, hybi)
    case ('test_l26')
      call hybrid_coord_test_l26(p0, hyai, hybi)
    case ('test_l30')
      call hybrid_coord_test_l30(p0, hyai, hybi)
    case ('ecmwf_l50')
      call hybrid_coord_ecmwf_l50(p0, hyai, hybi)
    case default
    end select

    if (all(hyam == 0)) then
      do k = 1, num_lev
        hyam(k) = 0.5d0 * (hyai(k) + hyai(k+1))
        hybm(k) = 0.5d0 * (hybi(k) + hybi(k+1))
      end do
    end if

    do k = 1, num_lev
      global_mesh%full_lev(k) = hyam(k) + hybm(k)
    end do
    do k = 1, num_lev + 1
      global_mesh%half_lev(k) = hyai(k) + hybi(k)
    end do

  end subroutine hybrid_coord_init

  subroutine hybrid_coord_final()

    deallocate(hyai)
    deallocate(hybi)
    deallocate(hyam)
    deallocate(hybm)

  end subroutine hybrid_coord_final

  pure real(r8) function hybrid_coord_calc_ph(k, phs) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: phs

    res = hyai(k) * p0 + hybm(k) * phs

  end function hybrid_coord_calc_ph

  pure real(r8) function hybrid_coord_calc_ph_lev(k, phs) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: phs

    res = hyai(k) * p0 + hybi(k) * phs

  end function hybrid_coord_calc_ph_lev

  pure real(r8) function hybrid_coord_calc_dphdt_lev(k, dphsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dphsdt

    res = hybi(k) * dphsdt

  end function hybrid_coord_calc_dphdt_lev

end module hybrid_coord_mod
