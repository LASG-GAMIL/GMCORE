module damp_mod

  use const_mod
  use parallel_mod

  implicit none

  integer, parameter :: diff_halo_width(2:8) = [1, 2, 2, 3, 3, 4, 4]

  real(r8) :: diff_weights(9,2:8) = reshape([ &
     1, -2,   1,   0,   0,   0,   0,   0,  0,  & ! 2
    -1,  3, - 3,   1,   0,   0,   0,   0,  0,  & ! 3
     1, -4,   6, - 4,   1,   0,   0,   0,  0,  & ! 4
    -1,  5, -10,  10, - 5,   1,   0,   0,  0,  & ! 5
     1, -6,  15, -20,  15, - 6,   1,   0,  0,  & ! 6
    -1,  7, -21,  35, -35,  21, - 7,   1,  0,  & ! 7
     1, -8,  28, -56,  70, -56,  28, - 8,  1   & ! 8
  ], [9, 7])

contains

  subroutine damp_run(order, dt, dx, lb, ub, n, f)

    integer , intent(in   ) :: order
    real(r8), intent(in   ) :: dt
    real(r8), intent(in   ) :: dx
    integer , intent(in   ) :: lb
    integer , intent(in   ) :: ub
    integer , intent(in   ) :: n
    real(r8), intent(inout) :: f(lb:ub)

    integer i, ns
    real(r8) g(lb:ub)
    real(r8) df(lb:ub)
    real(r8) a, w(9)

    a  = (dx / 2.0_r8)**order / dt

    call parallel_fill_halo(1 - lb, f)
    if (order == 2) then
      ns = diff_halo_width(order)
      w  = diff_weights(:,order)
      do i = 1, n
        g(i) = sum(f(i-ns:i+ns) * w(:2*ns+1))
      end do
      g = g * (-1)**(order / 2 + 1) * a / dx**order
      do i = 1, n
        f(i) = f(i) + dt * g(i)
      end do
    else
      ns = diff_halo_width(order - 1)
      w  = diff_weights(:,order - 1)
      do i = 1, n
        g (i) = sum(f(i-ns+1:i+ns) * w(:2*ns))
        df(i) = f(i+1) - f(i)
      end do
      g = g * (-1)**(order / 2) * a / dx**order
      do i = 1, n
        g(i) = g(i) * max(0.0_r8, sign(1.0_r8, -g(i) * df(i)))
      end do
      call parallel_fill_halo(1 - lb, g)
      do i = 1, n
        f(i) = f(i) - dt * (g(i) - g(i-1))
      end do
    end if

  end subroutine damp_run

end module damp_mod