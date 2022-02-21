module vert_coord_mod

  use flogger
  use const_mod
  use namelist_mod
  use sigma_coord_mod
  use hybrid_coord_mod
  use mesh_mod
  use process_mod

  implicit none

  private

  public vert_coord_init
  public vert_coord_final
  public vert_coord_calc_ph
  public vert_coord_calc_ph_lev
  public vert_coord_calc_dphdt_lev

  interface
    pure real(r8) function vert_coord_calc_ph_interface(k, phs)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: phs
    end function vert_coord_calc_ph_interface

    pure real(r8) function vert_coord_calc_ph_lev_interface(k, phs)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: phs
    end function vert_coord_calc_ph_lev_interface

    pure real(r8) function vert_coord_calc_dphdt_lev_interface(k, dphsdt)
      import r8
      integer, intent(in) :: k
      real(r8), intent(in) :: dphsdt
    end function vert_coord_calc_dphdt_lev_interface
  end interface

  procedure(vert_coord_calc_ph_interface), pointer :: vert_coord_calc_ph
  procedure(vert_coord_calc_ph_lev_interface), pointer :: vert_coord_calc_ph_lev
  procedure(vert_coord_calc_dphdt_lev_interface), pointer :: vert_coord_calc_dphdt_lev

contains

  subroutine vert_coord_init(num_lev, namelist_file, scheme, template)

    integer, intent(in) :: num_lev
    character(*), intent(in), optional :: namelist_file
    character(*), intent(in), optional :: scheme
    character(*), intent(in), optional :: template

    integer k

    if (present(scheme)) then
      if (vert_coord_scheme /= scheme .and. is_root_proc()) then
        call log_notice('Change vert_coord_scheme to ' // trim(scheme) // '.')
      end if
      vert_coord_scheme = scheme
    end if
    if (present(template)) then
      if (vert_coord_template /= template .and. is_root_proc()) then
        call log_notice('Change vert_coord_template to ' // trim(template) // '.')
      end if
      vert_coord_template = template
    else if (baroclinic .and. (vert_coord_scheme == 'hybrid' .and. vert_coord_template == 'N/A')) then
      call log_error('Parameter vert_coord_template is not set!', pid=proc%id)
    end if

    select case (vert_coord_scheme)
    case ('sigma')
      call sigma_coord_init(num_lev, namelist_file, vert_coord_template)
      vert_coord_calc_ph => sigma_coord_calc_ph
      vert_coord_calc_ph_lev => sigma_coord_calc_ph_lev
      vert_coord_calc_dphdt_lev => sigma_coord_calc_dphdt_lev
    case ('hybrid')
      call hybrid_coord_init(num_lev, namelist_file, vert_coord_template)
      vert_coord_calc_ph => hybrid_coord_calc_ph
      vert_coord_calc_ph_lev => hybrid_coord_calc_ph_lev
      vert_coord_calc_dphdt_lev => hybrid_coord_calc_dphdt_lev
    end select

    ! Set vertical level intervals.
    do k = global_mesh%full_lev_ibeg, global_mesh%full_lev_iend
      global_mesh%full_dlev(k) = global_mesh%half_lev(k+1) - global_mesh%half_lev(k)
    end do
    do k = global_mesh%half_lev_ibeg + 1, global_mesh%half_lev_iend - 1
      global_mesh%half_dlev(k) = global_mesh%full_lev(k) - global_mesh%full_lev(k-1)
      global_mesh%half_dlev_upper(k) = global_mesh%half_lev(k) - global_mesh%full_lev(k-1)
      global_mesh%half_dlev_lower(k) = global_mesh%full_lev(k) - global_mesh%half_lev(k)
    end do
    global_mesh%half_dlev(1) = global_mesh%full_lev(1) - global_mesh%half_lev(1)
    global_mesh%half_dlev_lower(1) = global_mesh%half_dlev(1)
    global_mesh%half_dlev(global_mesh%half_lev_iend) = global_mesh%half_lev(global_mesh%half_lev_iend) - global_mesh%full_lev(global_mesh%full_lev_iend)
    global_mesh%half_dlev_upper(global_mesh%half_lev_iend) = global_mesh%half_dlev(global_mesh%half_lev_iend)

  end subroutine vert_coord_init

  subroutine vert_coord_final()

    select case (vert_coord_scheme)
    case ('sigma')
      call sigma_coord_final()
    case ('hybrid')
      call hybrid_coord_final()
    end select

  end subroutine vert_coord_final

end module vert_coord_mod
