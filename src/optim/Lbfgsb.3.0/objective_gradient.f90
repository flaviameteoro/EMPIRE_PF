subroutine objective_gradient(n,x,g)
implicit none
integer, parameter :: rk = kind(1.0d0)
integer, intent(in) :: n
real(kind=rk), dimension(n), intent(in) :: x
real(kind=rk), dimension(n), intent(out) :: g

!     Declare a few additional variables for this sample problem
integer                :: i
real(rk)               :: t1, t2


! test function gradient

t1   = x(2) - x(1)**2
g(1) = 2.d0*(x(1) - 1.d0) - 1.6d1*x(1)*t1
do 22 i=2, n-1
   t2   = t1
   t1   = x(i+1) - x(i)**2
   g(i) = 8.d0*t2 - 1.6d1*x(i)*t1
22 continue
   g(n) = 8.d0*t1

end subroutine objective_gradient
