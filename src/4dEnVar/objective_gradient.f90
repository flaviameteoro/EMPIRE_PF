subroutine objective_gradient(n,x,g)
implicit none
integer, parameter :: rk = kind(1.0d0)
integer, intent(in) :: n
real(kind=rk), dimension(n), intent(in) :: x
real(kind=rk), dimension(n), intent(out) :: g
! Rosenbrock function gradient
g(1) = 200*(x(2) - x(1)**2)*(-2*x(1)) - 2*(1 - x(1))
g(2) = 200*(x(2) - x(1)**2)

end subroutine objective_gradient
