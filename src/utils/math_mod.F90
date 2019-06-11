module math_mod

  use flogger
  use const_mod

  private
  
  public cross_product
  public math_inv_matrix

contains

  function cross_product(x, y) result(res)

    real(real_kind), intent(in) :: x(3)
    real(real_kind), intent(in) :: y(3)
    real(real_kind) res(3)

    res(1) = x(2) * y(3) - x(3) * y(2)
    res(2) = x(3) * y(1) - x(1) * y(3)
    res(3) = x(1) * y(2) - x(2) * y(1)

  end function cross_product
    
  subroutine math_inv_matrix(n, A, B)

    integer, intent(in) :: n
    real(real_kind), intent(inout) :: A(n,n)
    real(real_kind), intent(out)   :: B(n,n)

    real(real_kind) C(n,n)
    integer ipiv(n)
    integer i, j, k

    C = 0.0d0
    do i = 1, n
      C(i,i) = 1.0d0
    end do

    call gaussian_elimination(n, A, ipiv)

    do i = 1, n - 1
      do j = i + 1, n
        do k = 1, n
          C(ipiv(j),k) = C(ipiv(j),k) - A(ipiv(j),i) * C(ipiv(i),k)
        end do
      end do
    end do

    do i = 1, n
      B(n,i) = C(ipiv(n),i) / A(ipiv(n),n)
      do j = n - 1, 1, -1
        B(j,i) = C(ipiv(j),i)
        do k = j + 1, n
          B(j,i) = B(j,i) - A(ipiv(j),k) * B(k,i)
        end do
        B(j,i) = B(j,i) / A(ipiv(j),j)
      end do
    end do

  end subroutine math_inv_matrix

  subroutine gaussian_elimination(n, A, ipiv)

    integer, intent(in) :: n
    real(real_kind), intent(inout) :: A(n,n)
    integer, intent(out) :: ipiv(n)

    real(real_kind) c1, c(n), pi, pi1, pj
    integer i, j, k, itmp

    do i = 1, n
      ipiv(i) = i
    end do

    do i = 1, n
      c1 = 0.0d0
      do j = 1, n
        c1 = max(c1, abs(A(i,j)))
      end do
      c(i) = c1
    end do

    do j = 1, n - 1
      pi1 = 0.0
      do i = j, n
        pi = abs(A(ipiv(i),j)) / c(ipiv(i))
        if (pi > pi1) then
          pi1 = pi
          k   = i
        end if
      end do

      itmp    = ipiv(j)
      ipiv(j) = ipiv(k)
      ipiv(k) = itmp
      do i = j + 1, n
        pj = A(ipiv(i),j) / A(ipiv(j),j)
        A(ipiv(i),j) = pj
        do k = j + 1, n
          A(ipiv(i),k) = A(ipiv(i),k) - pj * A(ipiv(j),k)
        end do
      end do
    end do

  end subroutine gaussian_elimination

end module math_mod
