module pgf_mod

  use flogger
  use pgf_swm_mod
  use pgf_lin97_mod

  implicit none

  private

  public pgf_init
  public pgf_prepare
  public pgf_run

  interface
    subroutine pgf_prepare_interface(block, state)
      import block_type, state_type
      type(block_type), intent(in) :: block
      type(state_type), intent(inout) :: state
    end subroutine pgf_prepare_interface

    subroutine pgf_run_interface(block, state, tend)
      import block_type, state_type, tend_type
      type(block_type), intent(inout) :: block
      type(state_type), intent(in) :: state
      type(tend_type), intent(inout) :: tend
    end subroutine pgf_run_interface
  end interface

  procedure(pgf_prepare_interface), pointer :: pgf_prepare
  procedure(pgf_run_interface), pointer :: pgf_run

contains

  subroutine pgf_init()

    if (baroclinic) then
      select case (pgf_scheme)
      case ('lin97')
        pgf_prepare => pgf_lin97_prepare
        pgf_run => pgf_lin97_run
      case default
        if (is_root_proc()) call log_error('Unknown PGF scheme ' // trim(pgf_scheme) // '!')
      end select
    else
      pgf_prepare => pgf_swm_prepare
      pgf_run => pgf_swm_run
    end if

  end subroutine pgf_init

end module pgf_mod
