module test_forcing_mod

  use namelist_mod
  use held_suarez_test_mod
  use block_mod

  private

  public test_forcing_run

contains

  subroutine test_forcing_run(block, dt, state)

    type(block_type), intent(in) :: block
    real(8), intent(in) :: dt
    type(state_type), intent(inout) :: state

    select case (test_case)
    case ('held_suarez')
      call held_suarez_test_apply_forcing(block, dt, state)
    end select

  end subroutine test_forcing_run

end module test_forcing_mod
