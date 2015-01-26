subroutine objective_function(n,x,f)
implicit none
integer, parameter :: rk = kind(1.0d0)
integer, intent(in) :: n
real(kind=rk), dimension(n), intent(in) :: x
real(kind=rk), intent(out) :: f

integer :: i
! test function

f=.25d0*( x(1)-1.d0 )**2
do 20 i=2, n
   f = f + ( x(i)-x(i-1 )**2 )**2
20 continue
   f = 4.d0*f


end subroutine objective_function
