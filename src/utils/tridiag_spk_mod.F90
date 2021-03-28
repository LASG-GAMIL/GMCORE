module tridiag_spk_mod

  use mpi
  use flogger
  use const_mod
  use parallel_mod
  use block_mod

  implicit none

  private

  public tridiag_spk_solver_type

#if (REAL_KIND == 4)
  integer, parameter :: md_mpi_real = MPI_REAL4
#elif (REAL_KIND == 8)
  integer, parameter :: md_mpi_real = MPI_REAL8
#endif

  type md_matrix_type
    integer n_size
    integer nnz
    integer lnz
    real(r8), allocatable, dimension(:,:) :: matrix
    real(r8), allocatable, dimension(:,:) :: w, v
    real(r8), allocatable, dimension(:,:) :: w_next, v_prev
  contains
    procedure :: init => md_matrix_init
    procedure :: clear => md_matrix_clear
    final :: md_matrix_final
  end type md_matrix_type

  type tridiag_spk_solver_type
    logical :: cloned = .false.
    type(md_matrix_type), pointer :: A(:)
    integer :: comm    = MPI_COMM_NULL
    integer :: id      = MPI_COMM_NULL
    integer :: prev_id = MPI_PROC_NULL
    integer :: next_id = MPI_PROC_NULL
    ! Working variables
    integer :: nrhs = 0
    integer buf_count
    integer buf_size
    integer pos1
    integer pos2
    real(r8), allocatable, dimension(:) :: send_buf1, recv_buf1
    real(r8), allocatable, dimension(:) :: send_buf2, recv_buf2
    real(r8), allocatable, dimension(:,:) :: x_prev, x_next
  contains
    procedure :: init_sym_const => tridiag_spk_solver_init_sym_const
    procedure :: cyclic_factor => tridiag_spk_solver_cyclic_factor
    procedure :: cyclic_solve => tridiag_spk_solver_cyclic_solve
    procedure :: solve => tridiag_spk_solver_solve
    procedure :: clone => tridiag_spk_solver_clone
    final :: tridiag_spk_solver_final
  end type tridiag_spk_solver_type

contains

  subroutine md_matrix_init(this, n_size, nnz)

    class(md_matrix_type), intent(inout) :: this
    integer, intent(in) :: n_size
    integer, intent(in) :: nnz

    call this%clear()

    this%n_size = n_size
    this%nnz    = nnz
    this%lnz    = (nnz - 1) / 2
    allocate(this%matrix(this%nnz,n_size))
    allocate(this%w     (this%lnz,n_size))
    allocate(this%v     (this%lnz,n_size))
    allocate(this%w_next(this%lnz,this%lnz))
    allocate(this%v_prev(this%lnz,this%lnz))

  end subroutine md_matrix_init

  subroutine md_matrix_clear(this)

    class(md_matrix_type), intent(inout) :: this

    if (allocated(this%matrix)) deallocate(this%matrix)
    if (allocated(this%w     )) deallocate(this%w     )
    if (allocated(this%v     )) deallocate(this%v     )
    if (allocated(this%w_next)) deallocate(this%w_next)
    if (allocated(this%v_prev)) deallocate(this%v_prev)

  end subroutine md_matrix_clear

  subroutine md_matrix_final(this)

    type(md_matrix_type), intent(inout) :: this

    call this%clear()

  end subroutine md_matrix_final

  subroutine tridiag_spk_solver_init_sym_const(this, n, a, b)

    class(tridiag_spk_solver_type), intent(inout) :: this
    integer, intent(in) :: n
    real(r8), intent(in) :: a
    real(r8), intent(in) :: b

    allocate(this%A(1))

    call this%A(1)%init(n, 3)

    this%A(1)%matrix(1,:) = b
    this%A(1)%matrix(2,:) = a
    this%A(1)%matrix(3,:) = b

    this%comm    = proc%zonal_circle%comm
    this%id      = proc%zonal_circle%id
    this%prev_id = proc%zonal_circle%west_ngb_id
    this%next_id = proc%zonal_circle%east_ngb_id

    call this%cyclic_factor()

  end subroutine tridiag_spk_solver_init_sym_const

  subroutine tridiag_spk_solver_cyclic_factor(this)

    class(tridiag_spk_solver_type), intent(inout) :: this

    integer n, i, j, k
    real(r8) lij
    integer neqs, n_size, nnz, lnz
    integer status(MPI_STATUS_SIZE), ierr

    neqs   = size(this%A)
    n_size = this%A(1)%n_size
    nnz    = this%A(1)%nnz
    lnz    = this%A(1)%lnz

    this%buf_count = lnz * lnz * neqs
    this%buf_size  = this%buf_count * REAL_KIND

    allocate(this%send_buf1(this%buf_count))
    allocate(this%send_buf2(this%buf_count))
    allocate(this%recv_buf1(this%buf_count))
    allocate(this%recv_buf2(this%buf_count))

    ! LU factorization
    this%pos1 = 0
    this%pos2 = 0
    do n = 1, neqs
      do j = 1, n_size - 1
        do i = j + 1, min(j + lnz, n_size)
          lij = this%A(n)%matrix(lnz+1+j-i,i) / this%A(n)%matrix(lnz+1,j)
          this%A(n)%matrix(lnz+1+j-i,i) = lij
          do k = 1, nnz - 1 + j - i
            this%A(n)%matrix(lnz+1+j-i+k,i) = this%A(n)%matrix(lnz+1+j-i+k,i) - &
              lij * this%A(n)%matrix(lnz+1+k,j)
          end do
        end do
      end do
      this%A(n)%v = 0
      this%A(n)%w = 0
      do i = 1, lnz
        do j = 1, i
          this%A(n)%v(j,n_size-lnz+i) = this%A(n)%matrix(nnz-i+j,n_size-lnz+i)
        end do
        do j = 1, lnz - i + 1
          this%A(n)%w(j,i) = this%A(n)%matrix(j,i)
        end do
      end do

      call struct_sptrsv(n_size, nnz, lnz, this%A(n)%matrix, this%A(n)%v)
      call struct_sptrsv(n_size, nnz, lnz, this%A(n)%matrix, this%A(n)%w)

      call MPI_PACK(this%A(n)%v(1, n_size-lnz+1), lnz * lnz, md_mpi_real, this%send_buf1, this%buf_size, this%pos1, this%comm, ierr)
      call MPI_PACK(this%A(n)%w(1, 1), lnz * lnz, md_mpi_real, this%send_buf2, this%buf_size, this%pos2, this%comm, ierr)
    end do

    if (this%prev_id == this%id .and. this%next_id == this%id) then
      this%recv_buf1(:) = this%send_buf1(:)
      this%recv_buf2(:) = this%send_buf2(:)
    else
      call MPI_SENDRECV(this%send_buf1, this%buf_count, md_mpi_real, this%next_id, 0, &
                        this%recv_buf1, this%buf_count, md_mpi_real, this%prev_id, 0, &
                        this%comm, status, ierr)
      call MPI_SENDRECV(this%send_buf2, this%buf_count, md_mpi_real, this%prev_id, 0, &
                        this%recv_buf2, this%buf_count, md_mpi_real, this%next_id, 0, &
                        this%comm, status, ierr)
    end if

    this%pos1 = 0
    this%pos2 = 0
    do n = 1, neqs
      call MPI_UNPACK(this%recv_buf1, this%buf_size, this%pos1, this%A(n)%v_prev(1, 1), lnz * lnz, md_mpi_real, this%comm, ierr)
      call MPI_UNPACK(this%recv_buf2, this%buf_size, this%pos2, this%A(n)%w_next(1, 1), lnz * lnz, md_mpi_real, this%comm, ierr)
    end do

    deallocate(this%send_buf1)
    deallocate(this%send_buf2)
    deallocate(this%recv_buf1)
    deallocate(this%recv_buf2)

  end subroutine tridiag_spk_solver_cyclic_factor

  subroutine tridiag_spk_solver_cyclic_solve(this, rhs, x)

    class(tridiag_spk_solver_type), intent(inout) :: this
    real(r8), intent(in) :: rhs(:,:,:)
    real(r8), intent(out) :: x(:,:,:)

    integer n, i
    integer neqs, nrhs, n_size, lnz, nnz
    integer status(mpi_status_size), ierr
    real(r8) b1, b2, g1(size(x, 3)), g2(size(x, 3))

    x = rhs

    neqs   = size(this%A)
    nrhs   = size(x, 3)
    n_size = this%A(1)%n_size
    lnz    = this%A(1)%lnz
    nnz    = this%A(1)%nnz

    do n = 1, neqs
      call struct_sptrsv(n_size, nnz, nrhs, this%A(n)%matrix, x(:,:,n))
    end do

    if (this%nrhs /= nrhs) then
      if (allocated(this%send_buf1)) deallocate(this%send_buf1)
      if (allocated(this%recv_buf1)) deallocate(this%recv_buf1)
      if (allocated(this%send_buf2)) deallocate(this%send_buf2)
      if (allocated(this%recv_buf2)) deallocate(this%recv_buf2)
      if (allocated(this%x_next   )) deallocate(this%x_next   )
      if (allocated(this%x_prev   )) deallocate(this%x_prev   )
      this%nrhs = nrhs
      this%buf_count = nrhs * lnz * neqs
      this%buf_size = this%buf_count * REAL_KIND
      allocate(this%send_buf1(this%buf_count))
      allocate(this%recv_buf1(this%buf_count))
      allocate(this%send_buf2(this%buf_count))
      allocate(this%recv_buf2(this%buf_count))
      allocate(this%x_next(nrhs,lnz))
      allocate(this%x_prev(nrhs,lnz))
    end if

    this%pos1 = 0
    this%pos2 = 0
    do n = 1, neqs
      call MPI_PACK(x(1,1,n), nrhs * lnz, md_mpi_real, this%send_buf2, this%buf_size, this%pos1, this%comm, ierr)
      call MPI_PACK(x(1,n_size-lnz+1,n), nrhs * lnz, md_mpi_real, this%send_buf1, this%buf_size, this%pos2, this%comm, ierr)
    end do

    if (this%prev_id == this%id .and. this%next_id == this%id) then
      this%recv_buf1(:) = this%send_buf1(:)
      this%recv_buf2(:) = this%send_buf2(:)
    else
      call MPI_SENDRECV(this%send_buf1, this%buf_count, md_mpi_real, this%next_id, 0, &
                        this%recv_buf1, this%buf_count, md_mpi_real, this%prev_id, 0, &
                        this%comm, status, ierr)
      call MPI_SENDRECV(this%send_buf2, this%buf_count, md_mpi_real, this%prev_id, 0, &
                        this%recv_buf2, this%buf_count, md_mpi_real, this%next_id, 0, &
                        this%comm, status, ierr)
    end if

    if (lnz == 1) then
      do n = 1, neqs
        ! Bottom
        b1 = this%A(n)%v(1,n_size)
        b2 = this%A(n)%w_next(1,1)
        g1 = x(:,n_size,n)
        g2 = this%recv_buf2((n-1)*nrhs+1:n*nrhs)
        x(:,n_size,n) = (g1 - b1 * g2) / (1 - b1 * b2)
        this%x_next(:,1) = (g2 - b2 * g1) / (1 - b1 * b2)
        ! Top
        b1 = this%A(n)%v_prev(1,1)
        b2 = this%A(n)%w(1,1)
        g1 = this%recv_buf1((n-1)*nrhs+1:n*nrhs)
        g2 = x(:,1,n)
        this%x_prev(:,1) = (g1 - b1 * g2) / (1 - b1 * b2)
        x(:,1,n) = (g2 - b2 * g1) / (1 - b1 * b2)
        ! Calculate X'
        do i = 2, n_size - 1
          x(:,i,n) = x(:,i,n) - this%A(n)%v(1,i) * this%x_next(:,1) - this%A(n)%w(1,i) * this%x_prev(:,1)
        end do
      end do
    else
      call log_error('Not implemented!', __FILE__, __LINE__)
    end if

  end subroutine tridiag_spk_solver_cyclic_solve

  subroutine struct_sptrsv(n_size, nnz, nrhs, a, x)

    integer, intent(in) :: n_size
    integer, intent(in) :: nnz
    integer, intent(in) :: nrhs
    real(r8), intent(in) :: a(nnz,n_size)
    real(r8), intent(inout) :: x(nrhs,n_size)

    integer i, j, k
    integer lnz

    lnz = (nnz - 1) / 2

    ! forward  LY = RHS
    do i = 1, n_size
      do j = 1, lnz
        k = i - lnz - 1 + j
        if (k >= 1) x(:,i) = x(:,i) - a(j,i) * x(:,k)
      end do
    end do
    ! backward UX = Y
    do i = n_size, 1, -1
      do j = 1, lnz
        k = i + j
        if (k <= n_size) x(:,i) = x(:,i) - a(lnz+1+j,i) * x(:,k)
      end do
      x(:,i) = x(:,i) / a(lnz+1,i)
    end do

  end subroutine struct_sptrsv

  subroutine tridiag_spk_solver_clone(this, other)

    class(tridiag_spk_solver_type), intent(inout) :: this
    type (tridiag_spk_solver_type), intent(in   ) :: other

    this%cloned   = .true.
    this%A        => other%A
    this%comm     = other%comm
    this%id       = other%id
    this%prev_id  = other%prev_id
    this%next_id  = other%next_id

  end subroutine tridiag_spk_solver_clone

  subroutine tridiag_spk_solver_solve(this, rhs, x)

    class(tridiag_spk_solver_type), intent(inout) :: this
    real(r8), intent(in) :: rhs(:)
    real(r8), intent(out) :: x(:)

    real(r8) rhs_array(1,this%A(1)%n_size,1)
    real(r8) x_array  (1,this%A(1)%n_size,1)

    rhs_array(1,:,1) = rhs
    x_array  (1,:,1) = x

    call this%cyclic_solve(rhs_array, x_array)

    x = x_array(1,:,1)

  end subroutine tridiag_spk_solver_solve

  subroutine tridiag_spk_solver_final(this)

    type(tridiag_spk_solver_type), intent(inout) :: this

    if (.not. this%cloned) then
      call this%A(1)%clear()
      deallocate(this%A)
    end if

    if (allocated(this%send_buf1)) deallocate(this%send_buf1)
    if (allocated(this%recv_buf1)) deallocate(this%recv_buf1)
    if (allocated(this%send_buf2)) deallocate(this%send_buf2)
    if (allocated(this%recv_buf2)) deallocate(this%recv_buf2)
    if (allocated(this%x_next   )) deallocate(this%x_next   )
    if (allocated(this%x_prev   )) deallocate(this%x_prev   )

  end subroutine tridiag_spk_solver_final

end module tridiag_spk_mod