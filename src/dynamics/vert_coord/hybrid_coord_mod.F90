module hybrid_coord_mod

  use flogger
  use const_mod, only: r8
  use namelist_mod, p0 => hybrid_coord_p0
  use hybrid_coord_test_mod
  use hybrid_coord_ecmwf_mod
  use hybrid_coord_mars_mod
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

  real(r8) :: local_ptop = 0

  namelist /hybrid_coord/ &
    hyai, hybi, hyam, hybm

contains

  subroutine hybrid_coord_init(num_lev, namelist_file, template)

    integer, intent(in) :: num_lev
    character(*), intent(in), optional :: namelist_file
    character(*), intent(in), optional :: template

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
    end if

    if (ierr /= 0 .and. .not. (baroclinic .or. advection) .and. is_root_proc()) then
      call log_notice('Run shallow-water model.')
    end if

    if (present(template)) then
      select case (template)
      case ('test_l15')
        call hybrid_coord_test_l15(p0, ptop, hyai, hybi)
      case ('test_l26')
        call hybrid_coord_test_l26(p0, ptop, hyai, hybi)
      case ('test_l30')
        call hybrid_coord_test_l30(p0, ptop, hyai, hybi)
      case ('test_l95')
        call hybrid_coord_test_l95(p0, ptop, hyai, hybi)
      case ('ecmwf_l50')
        call hybrid_coord_ecmwf_l50(p0, ptop, hyai, hybi)
      case ('ecmwf_l90')
        call hybrid_coord_ecmwf_l90(p0, ptop, hyai, hybi)
      case ('wrf_l16')
        call hybrid_coord_wrf_l16(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case ('wrf_l26')
        call hybrid_coord_wrf_l26(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case ('wrf_l32')
        call hybrid_coord_wrf_l32(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case ('wrf_l60')
        call hybrid_coord_wrf_l60(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case ('wrf_l64')
        call hybrid_coord_wrf_l64(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case ('schar_l40')
        call hybrid_coord_schar_l40(p0, ptop, hyai, hybi)
      case ('dcmip21_l60')
        call hybrid_coord_dcmip21_l60(p0, ptop, hyai, hybi)
      case ('dcmip31_l10')
        call hybrid_coord_dcmip31_l10(p0, ptop, hyai, hybi)
      case ('waccm_l70')
        call hybrid_coord_waccm_l70(p0, ptop, hyai, hybi)
      case ('emars28')
        call hybrid_coord_mars_emars28(p0, ptop, hyai, hybi)
      case ('dcmip_l60')
        call hybrid_coord_dcmip_l60(p0, ptop, hyai, hybi)
      case ('equal_eta')
        call hybrid_coord_equal_eta(p0, ptop, hyai, hybi)
        local_ptop = ptop
      case default
        if (baroclinic .and. template == 'N/A' .and. is_root_proc()) then
          call log_error('Hybrid vertical coordinate template "' // trim(template) // '" is invalid!')
        end if
      end select
    end if

    if (is_root_proc()) then
      call log_notice('Model top pressure is ' // to_str(ptop, '(ES10.2)') // 'Pa.')
    end if

    if ((baroclinic .or. advection) .and. global_mesh%num_full_lev > 1 .and. all(hyai == 0)) then
      call log_error('Hybrid coordinate parameters are not set!', pid=proc%id)
    end if

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

    res = hyam(k) * (p0 - local_ptop) + hybm(k) * (phs - local_ptop) + local_ptop

  end function hybrid_coord_calc_ph

  pure real(r8) function hybrid_coord_calc_ph_lev(k, phs) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: phs

    res = hyai(k) * (p0 - local_ptop) + hybi(k) * (phs - local_ptop) + local_ptop

  end function hybrid_coord_calc_ph_lev

  pure real(r8) function hybrid_coord_calc_dphdt_lev(k, dphsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dphsdt

    res = hybi(k) * dphsdt

  end function hybrid_coord_calc_dphdt_lev

end module hybrid_coord_mod
