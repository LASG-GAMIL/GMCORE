module tridiag_mod

  use flogger
  use const_mod
  use namelist_mod, only: zonal_tridiag_solver
  use process_mod
#ifdef HAS_MKL
  use tridiag_mkl_mod
#endif
  use tridiag_spk_mod

  implicit none

  private

  public tridiag_solver_type

  type tridiag_solver_type
    integer :: solver = 0
#ifdef HAS_MKL
    type(tridiag_mkl_solver_type), allocatable :: mkl_solver
#endif
    type(tridiag_spk_solver_type), allocatable :: spk_solver
  contains
    procedure :: init_sym_const => tridiag_solver_init_sym_const
    procedure :: clone => tridiag_solver_clone
    procedure :: solve => tridiag_solver_solve
    final :: tridiag_solver_final
  end type tridiag_solver_type

contains

  subroutine tridiag_solver_init_sym_const(this, n, a, b, solver)

    class(tridiag_solver_type), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: a
    real(r8), intent(in) :: b
    character(*), intent(in) :: solver

    select case (solver)
#ifdef HAS_MKL
    case ('mkl')
      this%solver = 1
      allocate(this%mkl_solver)
      call this%mkl_solver%init_sym_const(n, a, b)
#endif
    case ('spk')
      this%solver = 2
      allocate(this%spk_solver)
      call this%spk_solver%init_sym_const(n, a, b)
    case default
      if (is_root_proc()) call log_error('Unknown tridiag solver ' // trim(solver) // '!')
    end select

  end subroutine tridiag_solver_init_sym_const

  subroutine tridiag_solver_clone(this, other)

    class(tridiag_solver_type), intent(inout) :: this
    type(tridiag_solver_type), intent(in) :: other

    this%solver = other%solver
    select case (this%solver)
#ifdef HAS_MKL
    case (1)
      allocate(this%mkl_solver)
      call this%mkl_solver%clone(other%mkl_solver)
#endif
    case (2)
      allocate(this%spk_solver)
      call this%spk_solver%clone(other%spk_solver)
    case default
      if (is_root_proc()) call log_error('Uninitialized tridiag solver!', __FILE__, __LINE__)
    end select

  end subroutine tridiag_solver_clone

  subroutine tridiag_solver_solve(this, rhs, x)

    class(tridiag_solver_type), intent(inout) :: this
    real(r8), intent(in) :: rhs(:)
    real(r8), intent(out) :: x(:)

    select case (this%solver)
#ifdef HAS_MKL
    case (1)
      call this%mkl_solver%solve(rhs, x)
#endif
    case (2)
      call this%spk_solver%solve(rhs, x)
    case default
      if (is_root_proc()) call log_error('Uninitialized tridiag solver!', __FILE__, __LINE__)
    end select

  end subroutine tridiag_solver_solve

  subroutine tridiag_solver_final(this)

    type(tridiag_solver_type), intent(inout) :: this

#ifdef HAS_MKL
    if (allocated(this%mkl_solver)) deallocate(this%mkl_solver)
#endif
    if (allocated(this%spk_solver)) deallocate(this%spk_solver)

  end subroutine tridiag_solver_final

end module tridiag_mod