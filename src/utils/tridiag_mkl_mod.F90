#ifdef HAS_MKL

include 'mkl_dss.f90'

module tridiag_mkl_mod

  use mkl_dss
  use flogger
  use const_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public tridiag_mkl_solver_type

  type sparse_matrix_type
    integer :: nrow = 0
    integer :: ncol = 0
    integer :: nval = 0
    integer, allocatable :: row_idx(:)
    integer, allocatable :: cols(:)
    real(r8), allocatable :: vals(:)
    ! Interval working variables
    integer :: row_  = 0
    integer :: count_ = 0
  contains
    procedure :: init => sparse_matrix_init
    procedure :: set => sparse_matrix_set
    final :: sparse_matrix_final
  end type sparse_matrix_type

  type tridiag_mkl_solver_type
    type(MKL_DSS_HANDLE) handle
    type(sparse_matrix_type) A
  contains
    procedure :: init_sym_const => tridiag_mkl_solver_init_sym_const
    procedure :: clone => tridiag_mkl_solver_clone
    procedure :: solve => tridiag_mkl_solver_solve
  end type tridiag_mkl_solver_type

contains

  subroutine tridiag_mkl_solver_init_sym_const(this, n, a, b)

    class(tridiag_mkl_solver_type), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: a
    real(r8), intent(in) :: b

    integer i, j, ierr
    integer perm(1)

    ! Create coefficient sparse matrix for upper triangle.
    call this%A%init(n, n, 2 * n)
    do i = 1, n
      do j = 1, n
        if (i == j) then
          call this%A%set(i, j, a)
        else if (i == j - 1 .or. (i == 1 .and. j == n)) then
          call this%A%set(i, j, b)
        end if
      end do
    end do

    ierr = DSS_CREATE(this%handle, MKL_DSS_DEFAULTS)
    if (ierr /= MKL_DSS_SUCCESS) then
      if (is_root_proc()) call log_error('Failed to call DSS_CREATE!', __FILE__, __LINE__)
    end if

    ierr = DSS_DEFINE_STRUCTURE(this%handle, MKL_DSS_SYMMETRIC, this%A%row_idx, this%A%nrow, this%A%ncol, this%A%cols, this%A%nval)
    if (ierr /= MKL_DSS_SUCCESS) then
      if (is_root_proc()) call log_error('Failed to call DSS_DEFINE_STRUCTURE!', __FILE__, __LINE__)
    end if

    ierr = DSS_REORDER(this%handle, MKL_DSS_DEFAULTS, perm)
    if (ierr /= MKL_DSS_SUCCESS) then
      if (is_root_proc()) call log_error('Failed to call DSS_REORDER!', __FILE__, __LINE__)
    end if

    ierr = DSS_FACTOR_REAL(this%handle, MKL_DSS_DEFAULTS, this%A%vals)
    if (ierr /= MKL_DSS_SUCCESS) then
      if (is_root_proc()) call log_error('Failed to call DSS_FACTOR_REAL!', __FILE__, __LINE__)
    end if

  end subroutine tridiag_mkl_solver_init_sym_const

  subroutine sparse_matrix_init(this, nrow, ncol, nval)

    class(sparse_matrix_type), intent(out) :: this
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: nval

    this%nrow = nrow
    this%ncol = ncol
    this%nval = nval

    if (allocated(this%row_idx)) deallocate(this%row_idx)
    allocate(this%row_idx(nrow+1))
    if (allocated(this%cols)) deallocate(this%cols)
    allocate(this%cols(nval))
    if (allocated(this%vals)) deallocate(this%vals)
    allocate(this%vals(nval))

    this%row_idx(1) = 1
    this%row_idx(nrow+1) = nval + 1

  end subroutine sparse_matrix_init

  subroutine sparse_matrix_set(this, i, j, val)

    class(sparse_matrix_type), intent(inout) :: this
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(r8), intent(in) :: val

    this%count_ = this%count_ + 1
    this%vals(this%count_) = val
    this%cols(this%count_) = j
    if (this%row_ == 0) this%row_ = i ! Here we go.
    if (this%row_ /= i) then
      ! Jump to new row.
      this%row_idx(i) = this%count_
      this%row_ = i
    end if

  end subroutine sparse_matrix_set

  subroutine sparse_matrix_final(this)

    type(sparse_matrix_type), intent(inout) :: this

    if (allocated(this%vals)) deallocate(this%vals)
    if (allocated(this%cols)) deallocate(this%cols)
    if (allocated(this%row_idx)) deallocate(this%row_idx)

  end subroutine sparse_matrix_final

  subroutine tridiag_mkl_solver_clone(this, other)

    class(tridiag_mkl_solver_type), intent(inout) :: this
    type(tridiag_mkl_solver_type), intent(in) :: other

    this%handle = other%handle

  end subroutine tridiag_mkl_solver_clone

  subroutine tridiag_mkl_solver_solve(this, rhs, x)

    class(tridiag_mkl_solver_type), intent(inout) :: this
    real(r8), intent(in) :: rhs(:)
    real(r8), intent(out) :: x(:)

    integer ierr

    ierr = DSS_SOLVE_REAL(this%handle, MKL_DSS_DEFAULTS, rhs, 1, x)
    if (ierr /= MKL_DSS_SUCCESS) then
      if (is_root_proc()) call log_error('Failed to call DSS_SOLVE_REAL!', __FILE__, __LINE__)
    end if

  end subroutine tridiag_mkl_solver_solve

end module tridiag_mkl_mod

#endif