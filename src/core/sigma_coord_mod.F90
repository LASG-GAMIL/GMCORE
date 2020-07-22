module sigma_coord_mod

  use flogger
  use namelist_mod
  use const_mod
  use mesh_mod
  use process_mod

  implicit none

  private

  public sigma_coord_init
  public sigma_coord_final
  public sigma_coord_calc_ph_lev
  public sigma_coord_calc_dphdt_lev

  real(r8), allocatable, dimension(:) :: sigi
  real(r8), allocatable, dimension(:) :: sig

  real(r8) :: pt = 1d5 ! Model top pressure (Pa)

  namelist /sigma_coord/ &
    sig, sigi, pt

contains

  subroutine sigma_coord_init(num_lev, namelist_path)

    integer, intent(in) :: num_lev
    character(*), intent(in) :: namelist_path

    integer ierr, k

    if (allocated(sigi)) deallocate(sigi)
    if (allocated(sig )) deallocate(sig )

    allocate(sigi(num_lev+1)); sigi = 0
    allocate(sig (num_lev  )); sig  = 0

    open(10, file=namelist_path, status='old')
    read(10, nml=sigma_coord, iostat=ierr)
    close(10)

    if (ierr /= 0) then
      if (.not. baroclinic) then
        if (is_root_proc()) call log_notice('Run shallow-water model.')
        return
      !else
      !  call log_error('No sigma_coord parameters in ' // trim(namelist_path) // '!')
      end if
    end if

    sigi = [                                                             &
      0.0, 0.00270708138944839, 0.00770525729190707, 0.0158928127640089, &
      0.0277039568735249, 0.0425225717851423, 0.0595424438370873,        &
      0.0764861789945283, 0.09037000846462, 0.106703636692459,           &
      0.125919292600401, 0.148525517478954, 0.175120586853872,           &
      0.206408247111974, 0.243216666654805, 0.286519798881535,           &
      0.33746367336785, 0.397396487729323, 0.467904355709719,            &
      0.550853249115629, 0.648438367513145, 0.743820793246579,           &
      0.830649646952043, 0.903087677425118, 0.955900674543326,           &
      0.985079453568069, 1.0                                             &
    ]
    sig = [                                                                     &
      0.00135354069472419, 0.00520616934067773, 0.011799035027958,              &
      0.0217983848187669, 0.0351132643293336, 0.0510325078111148,               &
      0.0680143114158078, 0.0834280937295741, 0.0985368225785395,               &
      0.11631146464643, 0.137222405039678, 0.161823052166413,                   &
      0.190764416982923, 0.22481245688339, 0.26486823276817, 0.311991736124692, &
      0.367430080548587, 0.432650421719521, 0.509378802412674,                  &
      0.599645808314387, 0.696129580379862, 0.787235220099311,                  &
      0.86686866218858, 0.929494175984222, 0.970490064055697, 0.992539726784035 &
    ]
    pt = 219.406699761748

    do k = 1, num_lev
      global_mesh%full_lev(k) = sig(k)
    end do
    do k = 1, num_lev + 1
      global_mesh%half_lev(k) = sigi(k)
    end do

  end subroutine sigma_coord_init

  subroutine sigma_coord_final()

    deallocate(sigi)
    deallocate(sig )

  end subroutine sigma_coord_final

  pure real(r8) function sigma_coord_calc_ph_lev(k, phs) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: phs

    res = sigi(k) * (phs - pt) + pt

  end function sigma_coord_calc_ph_lev

  pure real(r8) function sigma_coord_calc_dphdt_lev(k, dphsdt) result(res)

    integer, intent(in) :: k
    real(r8), intent(in) :: dphsdt

    res = sigi(k) * dphsdt

  end function sigma_coord_calc_dphdt_lev

end module sigma_coord_mod
