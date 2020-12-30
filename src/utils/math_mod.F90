module math_mod

  use flogger
  use const_mod

  private
  
  public cross_product
  public math_inv_matrix
  public tridiag_thomas

contains

  function cross_product(x, y) result(res)

    real(16), intent(in) :: x(3)
    real(16), intent(in) :: y(3)
    real(16) res(3)

    res(1) = x(2) * y(3) - x(3) * y(2)
    res(2) = x(3) * y(1) - x(1) * y(3)
    res(3) = x(1) * y(2) - x(2) * y(1)

  end function cross_product
    
  subroutine math_inv_matrix(n, A, B)

    integer , intent(in   ) :: n
    real(r8), intent(inout) :: A(n,n)
    real(r8), intent(  out) :: B(n,n)

    real(r8) C(n,n)
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

    integer , intent(in   ) :: n
    real(r8), intent(inout) :: A(n,n)
    integer , intent(  out) :: ipiv(n)

    real(r8) c1, c(n), pi, pi1, pj
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

  subroutine tridiag_thomas(a, b, c, d, x)

    real(r8), intent(inout) :: a(:)
    real(r8), intent(inout) :: b(:)
    real(r8), intent(inout) :: c(:)
    real(r8), intent(inout) :: d(:)
    real(r8), intent(out) :: x(:)

    real(r8) denominator
    integer n, i
    !  _                                                _   _      _     _      _
    ! |  b(1)  c(1)                                      | | x(1  ) |   | d(1  ) |
    ! |  a(2)  b(2)  c(2)                                | | x(2  ) |   | d(2  ) |
    ! |        a(3)  b(3)  c(3)                          | | x(3  ) |   | d(3  ) |
    ! |           ...  ...  ...                          | | ...    | = | ...    |
    ! |                ...  ...  ...                     | | ...    |   | ...    |
    ! |                          a(n-1)  b(n-1)  c(n-1)  | | x(n-1) |   | d(n-1) |
    ! |_                                 a(n  )  b(n  ) _| |_x(n  )_|   |_d(n  )_|

    n = size(x)

    c(1) = c(1) / b(1)
    d(1) = d(1) / b(1)
    do i = 2, n
      denominator = b(i) - a(i) * c(i-1)
      c(i) = c(i) / denominator
      d(i) = (d(i) - a(i) * d(i - 1)) / denominator
    end do

    x(n) = d(n)
    do i = n - 1, 1, -1
      x(i) = d(i) - c(i) * x(i+1)
    end do

  end subroutine tridiag_thomas

end module math_mod
