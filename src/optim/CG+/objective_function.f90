subroutine objective_function(n,x,f)
implicit none
integer, parameter :: rk = kind(1.0d0)
integer, intent(in) :: n
real(kind=rk), dimension(n), intent(in) :: x
real(kind=rk), intent(out) :: f
! Rosenbrock function
f = 100.*((x(2) - x(1)**2)**2) + (1. - x(1))**2
end subroutine objective_function
