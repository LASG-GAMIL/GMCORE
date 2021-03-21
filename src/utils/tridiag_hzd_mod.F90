module tridiag_hzd_mod

  use mpi
  use flogger
  use const_mod
  use parallel_mod
  use block_mod

#define md_type real(8)
#define md_size 8
#if (md_size == 8)
#define MD_MPI_TYPE MPI_REAL8
#else
#define MD_MPI_TYPE MPI_REAL4
#endif

  implicit none

  private

  ! process id
  integer :: pid
  integer :: next_id
  integer :: prev_id

  ! whether cyclic solver is on a single process
  logical :: single_process

  ! communication group
  integer :: comm

  public hzd_tridiag_solver_type
	public zonal_comm_init

  ! 多对角矩阵
  type md_matrix
		integer :: n_size !矩阵规模
		integer :: nnz !总对角元个数
		integer :: lnz !副对角元个数
		md_type, allocatable, dimension(:, :) :: matrix 
		!矩阵数据(nnz, n_size)，1~lnz为左副对角元，第lnz+1个为主对角元，lnz+2~nnz为右副对角元
		md_type, allocatable, dimension(:, :) :: w, v
		!SPIKE分解用数据(lnz, n_size)，分别代表SPIKE分解A=DS后，S矩阵单位阵左/右的条状数据
		!单进程中，用v和w分别存储循环三对角LU分解后L矩阵右端的竖状数据和L矩阵下端的横状数据
		md_type, allocatable, dimension(:, :) :: w_next, v_prev
		!SPIKE分解用数据(lnz, lnz)，w_next为下个划分的wt，v_prev为上个划分的vb
  end type md_matrix

  ! 三对角求解器
  type hzd_tridiag_solver_type
		type(md_matrix) A(1) ! 特殊处理后的系数矩阵
  contains
		procedure :: init_sym_const => tridiag_solver_init_sym_const
		procedure :: solve => tridiag_solver_solve
		procedure :: clone => tridiag_solver_clone
		final     ::          tridiag_solver_final
  end type hzd_tridiag_solver_type

  !type(hzd_tridiag_solver_type) zonal_div_damp_solver

contains

  subroutine zonal_comm_init(comm_in, rank, next_rank, prev_rank)

		integer, intent(in) :: comm_in !纬圈三对角通信子
		integer, intent(in) :: rank, next_rank, prev_rank!当前进程id，纬圈下一个/上一个进程id

		comm = comm_in

		pid = rank
		next_id = next_rank
		prev_id = prev_rank

		if (pid == next_id .and. pid == prev_rank) then
			single_process = .true.
		else if (pid .ne. next_id .and. pid .ne. prev_rank) then
			single_process = .false.
		else
			print *, "Error: Invalid init args!"
		end if

		print *, "Pid:", proc%id, " finished zonal_comm_init: comm is", comm, "next_id and prev_id are ", next_id, prev_id, "single_proc:", single_process

  end subroutine zonal_comm_init

	subroutine md_cyclic_fact_kernel(neqs, ma)
		
		integer, intent(in) :: neqs !纬圈三对角方程个数
    type(md_matrix), intent(inout) :: ma(neqs) !预处理结构体，每个方程对应一个
    ! logical, intent(in) :: ifprint

    integer :: n, i, j, k
    md_type :: lij
    integer :: n_size, lnz, nnz

    integer :: buf_count, buf_size, pos_v, pos_w
    integer :: status(mpi_status_size), ierr
    md_type, allocatable, dimension(:) :: send_buf_v, recv_buf_v
    md_type, allocatable, dimension(:) :: send_buf_w, recv_buf_w

    ! LU factorization
    lnz = ma(1)%lnz
    nnz = ma(1)%nnz
    n_size = ma(1)%n_size

    do n = 2, neqs
      if (ma(n)%nnz .ne. nnz .or. ma(n)%n_size .ne. n_size) then
        print *, "Error: matrices of md_cyclic_fact must have the same size!"
        return
      end if
    end do

    buf_count = lnz * lnz * neqs
    buf_size = buf_count * md_size
    pos_v = 0; pos_w = 0

    allocate(send_buf_v(lnz*lnz*neqs))
    allocate(send_buf_w(lnz*lnz*neqs))
    allocate(recv_buf_v(lnz*lnz*neqs))
    allocate(recv_buf_w(lnz*lnz*neqs))

    do n = 1, neqs
      ! if (pid == 0) then
      !     print *, "matrix", ma(n)%matrix
      ! end if
      do j = 1, n_size - 1
        do i = j + 1, min(j + lnz, n_size)
          lij = ma(n)%matrix(lnz + 1 + j - i, i)/ ma(n)%matrix(lnz + 1, j)
          ma(n)%matrix(lnz + 1 + j - i, i) = lij
          do k = 1, nnz - 1 + j - i
            ma(n)%matrix(lnz + 1 + j - i + k, i) = ma(n)%matrix(lnz + 1 + j - i + k, i) - &
              lij * ma(n)%matrix(lnz + 1 + k, j)
          end do
        end do
      end do
      ! if (pid == 0) then
      !     print *, "a", ma(n)%matrix
      !     print *, "a(2, 2)", ma(n)%matrix(1, n_size) * ma(n)%matrix(3, n_size-1) + ma(n)%matrix(2, n_size)
      ! end if
      ma(n)%v(:, :) = 0
      ma(n)%w(:, :) = 0
      do i = 1, lnz
        do j = 1, i
          ma(n)%v(j, n_size-lnz+i) = ma(n)%matrix(nnz-i+j, n_size-lnz+i)
        end do
        do j = 1, lnz - i + 1
          ma(n)%w(j, i) = ma(n)%matrix(j, i)
        end do
      end do

      call struct_sptrsv(n_size, nnz, lnz, ma(n)%matrix, ma(n)%v)
      call struct_sptrsv(n_size, nnz, lnz, ma(n)%matrix, ma(n)%w)

      ! if (ifprint) then
      !   print *, "before packing"
      !   print *, "send_buf_v"
      !   print *, send_buf_v, lnz, neqs
      !   print *, "send_buf_w"
      !   print *, send_buf_w
      ! end if

      ! send_buf_v(n) = ma(n)%v(1, n_size-lnz+1)
      ! send_buf_w(n) = ma(n)%w(1, 1)
      call mpi_pack(ma(n)%v(1, n_size-lnz+1), lnz*lnz, MD_MPI_TYPE, send_buf_v(1), buf_size, pos_v, comm, ierr)
      call mpi_pack(ma(n)%w(1, 1), lnz*lnz, MD_MPI_TYPE, send_buf_w(1), buf_size, pos_w, comm, ierr)
        
      ! if (ifprint) then
      !   print *, "after packing before communicating"
      !   print *, "send_buf_v"
      !   print *, send_buf_v
      !   print *, "send_buf_w"
      !   print *, send_buf_w
      ! end if
    end do

    ! MPI communication
    ! Send vb & wt
    if (single_process == .true.) then
        ! recv_buf_v(:) = send_buf_w(:)
        ! recv_buf_w(:) = send_buf_v(:)
        recv_buf_v(:) = send_buf_v(:)
        recv_buf_w(:) = send_buf_w(:)
    else
      call mpi_sendrecv(send_buf_v, buf_count, MD_MPI_TYPE, next_id, 0, &
      recv_buf_v, buf_count, MD_MPI_TYPE, prev_id, 0, comm, status, ierr)

      call mpi_sendrecv(send_buf_w, buf_count, MD_MPI_TYPE, prev_id, 0, &
      recv_buf_w, buf_count, MD_MPI_TYPE, next_id, 0, comm, status, ierr)
    end if

    ! if (ifprint) then
    !   print *, "after communicating"
    !   print *, "recv_buf_v"
    !   print *, recv_buf_v
    !   print *, "recv_buf_w"
    !   print *, recv_buf_w
    ! end if

    pos_v = 0
    pos_w = 0
    do n = 1, neqs
        ! ma(n) = ma(n)
        call mpi_unpack(recv_buf_v, buf_size, pos_v, ma(n)%v_prev(1, 1), lnz*lnz, MD_MPI_TYPE, comm, ierr)
        call mpi_unpack(recv_buf_w, buf_size, pos_w, ma(n)%w_next(1, 1), lnz*lnz, MD_MPI_TYPE, comm, ierr)
    end do

    ! do n = 1, neqs
    !   ma(n)%v_prev(1, 1) = recv_buf_v(n)
    !   ma(n)%w_next(1, 1) = recv_buf_w(n)
    ! end do

    deallocate(send_buf_v)
    deallocate(send_buf_w)
    deallocate(recv_buf_v)
    deallocate(recv_buf_w)

    ! if (pid <= 3) then
    !   ! print *, pid, "wt", ma(1)%w_next
    !   ! print *, pid, "vb", ma(1)%v_prev
    !   ! print *, pid, "w", ma(1)%w
    !   ! print *, pid, "v", ma(1)%v
    ! end if

	end subroutine md_cyclic_fact_kernel

	subroutine md_cyclic_solver_kernel(neqs, n_size, nrhs, ma, rhs, x)

		integer, intent(in) :: neqs !纬圈三对角方程个数
    integer, intent(in) :: n_size !当前进程对应三对角方程的行数
    integer, intent(in) :: nrhs !每个三对角方程的右端向量个数
    type(md_matrix), intent(in) :: ma(neqs) !预处理结构体，每个方程对应一个
    md_type, intent(in) :: rhs(nrhs, n_size, neqs) !求解右端向量
    md_type, intent(out) :: x(nrhs, n_size, neqs) !求解结果

    integer :: buf_count, buf_size, pos_b, pos_t
    integer :: status(mpi_status_size), ierr
    md_type, allocatable, dimension(:) :: send_buf_b, recv_buf_b
    md_type, allocatable, dimension(:) :: send_buf_t, recv_buf_t

    md_type, allocatable, dimension(:, :) :: x_next, x_prev

    !!! For trd solver       !!!
    !!! | 1 b1 | |x1| = |g1| !!!
    !!! | b2 1 | |x2| = |g2| !!!
    md_type :: b1, b2, g1(nrhs), g2(nrhs)

    integer :: n, i
    integer :: lnz, nnz

    x(:, :, :) = rhs(:, :, :)

    lnz = ma(1) % lnz
    nnz = ma(1) % nnz

    do n = 2, neqs
      if (ma(n)%nnz .ne. nnz .or. ma(n)%n_size .ne. n_size) then
        print *, "Error: matrices of md_cyclic_solver must have the same size!"
        return
      end if
    end do

    ! DY = B
    do n = 1, neqs
      call struct_sptrsv(n_size, nnz, nrhs, ma(n)%matrix, x(1, 1, n))
    end do

    allocate(send_buf_b(nrhs*lnz*neqs))
    allocate(recv_buf_b(nrhs*lnz*neqs))
    allocate(send_buf_t(nrhs*lnz*neqs))
    allocate(recv_buf_t(nrhs*lnz*neqs))
    allocate(x_next(nrhs, lnz))
    allocate(x_prev(nrhs, lnz))

    buf_count = nrhs * lnz * neqs
    buf_size = buf_count * md_size
    pos_b = 0; pos_t = 0

    ! MPI communication
    ! Send rhs_b & rhs_t
    do n = 1, neqs
      call mpi_pack(x(1, 1, n), nrhs*lnz, MD_MPI_TYPE, send_buf_t, buf_size, pos_t, comm, ierr)
      call mpi_pack(x(1, n_size-lnz+1, n), nrhs*lnz, MD_MPI_TYPE, send_buf_b, buf_size, pos_b, comm, ierr)
    end do
    ! do n = 1, neqs
    !   send_buf_t(nrhs*(n-1)+1:nrhs*lnz*n) = x(1:nrhs, 1, n)
    !   send_buf_b(nrhs*(n-1)+1:nrhs*lnz*n) = x(1:nrhs, n_size-lnz+1, n)
    ! end do

    if (single_process == .true.) then
      ! recv_buf_b(:) = send_buf_t(:)
      ! recv_buf_t(:) = send_buf_b(:)
      recv_buf_b(:) = send_buf_b(:)
      recv_buf_t(:) = send_buf_t(:)
    else
      call mpi_sendrecv(send_buf_b, buf_count, MD_MPI_TYPE, next_id, 0, &
      recv_buf_b, buf_count, MD_MPI_TYPE, prev_id, 0, comm, status, ierr)

      call mpi_sendrecv(send_buf_t, buf_count, MD_MPI_TYPE, prev_id, 0, &
      recv_buf_t, buf_count, MD_MPI_TYPE, next_id, 0, comm, status, ierr)
    end if

    if (lnz == 1) then
      do n = 1, neqs
        ! bottom
        b1 = ma(n)%v(1, n_size)
        b2 = ma(n)%w_next(1, 1)
        g1(:) = x(:, n_size, n)
        g2(:) = recv_buf_t((n-1)*nrhs+1:n*nrhs)
        x(:, n_size, n) = (g1(:) - b1 * g2(:)) / (1 - b1 * b2)
        x_next(:, 1) = (g2(:) - b2 * g1(:)) / (1 - b1 * b2)

        ! top
        b1 = ma(n)%v_prev(1, 1)
        b2 = ma(n)%w(1, 1)
        g1(:) = recv_buf_b((n-1)*nrhs+1:n*nrhs)
        g2(:) = x(:, 1, n)
        x_prev(:, 1) = (g1(:) - b1 * g2(:)) / (1 - b1 * b2)
        x(:, 1, n) = (g2(:) - b2 * g1(:)) / (1 - b1 * b2)
            
        ! Calculate X'
        do i = 2, n_size - 1
          x(:, i, n) = x(:, i, n) - &
          ma(n)%v(1, i) * x_next(:, 1) - ma(n)%w(1, i) * x_prev(:, 1)
        end do
      end do
    else
      print *, "Not implemented!"
    end if

    ! print *, pid, "b1", rhs(1, 1, 1)
    ! print *, pid, "test1", x_prev(1, 1)*test_a(1, 1, 1) + x(1, 1, 1)*test_a(2, 1, 1) + x(1, 2, 1)*test_a(3, 1, 1)
    ! print *, pid, "bn", rhs(1, n_size, 1)
    ! print *, pid, "test2", x(1, n_size-1, 1)*test_a(1, n_size, 1) + x(1, n_size, 1)*test_a(2, n_size, 1) + &
    ! x_next(1, 1)*test_a(3, n_size, 1)

    deallocate(send_buf_b)
    deallocate(recv_buf_b)
    deallocate(send_buf_t)
    deallocate(recv_buf_t)
    deallocate(x_next)
    deallocate(x_prev)

	end subroutine md_cyclic_solver_kernel

	subroutine struct_sptrsv(n_size, nnz, nrhs, a, x)

		integer, intent(in) :: n_size
		integer, intent(in) :: nnz
		integer, intent(in) :: nrhs
		md_type, intent(in) :: a(nnz, n_size)
		md_type, intent(inout) :: x(nrhs, n_size)

		integer :: i, j, k
		integer :: lnz

		lnz = (nnz - 1) / 2
		!md_type :: res(nrhs, n_size)

		! forward  LY = RHS
		do i = 1, n_size
			do j = 1, lnz
				k = i - lnz - 1 + j
				if (k >= 1) x(:, i) = x(:, i) - a(j, i) * x(:, k)
			end do
		end do
		! backward UX = Y
		do i = n_size, 1, -1
			do j = 1, lnz
				k = i + j
				if (k <= n_size) x(:, i) = x(:, i) - a(lnz + 1 + j, i) * x(:, k)
			end do
			x(:, i) = x(:, i) / a(lnz + 1, i)
		end do

	end subroutine struct_sptrsv

	! structure management
	subroutine allocate_md_matrix(n_size, nnz, ma)

		integer, intent(in) :: n_size
		integer, intent(in) :: nnz
		type(md_matrix), intent(inout) :: ma

		ma%n_size = n_size
		ma%nnz = nnz
		ma%lnz = (nnz - 1) / 2
		allocate(ma%matrix(ma%nnz, n_size))
		allocate(ma%w(ma%lnz, n_size))
		allocate(ma%v(ma%lnz, n_size))
		allocate(ma%w_next(ma%lnz, ma%lnz))
		allocate(ma%v_prev(ma%lnz, ma%lnz))
		
	end subroutine allocate_md_matrix

	subroutine deallocate_md_matrix(ma)

		type(md_matrix), intent(inout) :: ma
		if(allocated(ma%matrix)) deallocate(ma%matrix)
		if(allocated(ma%w))      deallocate(ma%w)
		if(allocated(ma%v))      deallocate(ma%v)
		if(allocated(ma%w_next)) deallocate(ma%w_next)
		if(allocated(ma%v_prev)) deallocate(ma%v_prev)

	end subroutine deallocate_md_matrix

  subroutine tridiag_solver_init_sym_const(this, n, a, b) ! 系数是对称的

		class(hzd_tridiag_solver_type), intent(inout) :: this
		integer, intent(in) :: n
		real(r8), intent(in) :: a
		real(r8), intent(in) :: b
    ! logical, intent(in) :: ifprint

		!call zonal_comm_init(proc%zonal_circle%comm, proc%id, proc%ngb(east)%id, proc%ngb(west)%id)

		call allocate_md_matrix(n, 3, this%A(1))

		! 稀疏矩阵赋值
		this%A(1)%matrix(1,:) = b
		this%A(1)%matrix(2,:) = a
		this%A(1)%matrix(3,:) = b

		call md_cyclic_fact_kernel(1, this%A(1)) ! neqs=1: 单member

  end subroutine tridiag_solver_init_sym_const

	subroutine tridiag_solver_clone(this, cpy_solver) ! 系数是对称的

		class(hzd_tridiag_solver_type), intent(inout) :: this
		type (hzd_tridiag_solver_type), intent(in   ) :: cpy_solver

		!call zonal_comm_init(proc%zonal_circle%comm, proc%id, proc%ngb(east)%id, proc%ngb(west)%id)

		call allocate_md_matrix(cpy_solver%A(1)%n_size, cpy_solver%A(1)%nnz, this%A(1))

		! 全部都赋值
		this%A(1)%matrix = cpy_solver%A(1)%matrix
    this%A(1)%w      = cpy_solver%A(1)%w
    this%A(1)%w_next = cpy_solver%A(1)%w_next
    this%A(1)%v      = cpy_solver%A(1)%v
    this%A(1)%v_prev = cpy_solver%A(1)%v_prev

  end subroutine tridiag_solver_clone

  subroutine tridiag_solver_solve(this, rhs, x)

		class(hzd_tridiag_solver_type), intent(inout) :: this
		real(r8), intent(in) :: rhs(:)
		real(r8), intent(out) :: x(:)

		real(r8) rhs_array(1, this%A(1)%n_size, 1), x_array(1, this%A(1)%n_size, 1)
		rhs_array(1,:,1) = rhs

		call md_cyclic_solver_kernel(1, this%A(1)%n_size, 1, this%A(1), rhs_array, x_array)

		x = x_array(1,:,1)

  end subroutine tridiag_solver_solve

  subroutine tridiag_solver_final(this)

		type(hzd_tridiag_solver_type), intent(inout) :: this

		call deallocate_md_matrix(this%A(1))

  end subroutine tridiag_solver_final

end module tridiag_hzd_mod
