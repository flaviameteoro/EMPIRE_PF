Subroutine UniformRandomNumbers1D(minv, maxv, n,phi)
!use random
implicit none
integer, parameter :: rk = kind(1.0D0)
integer, intent(in) :: n
real(kind=rk), intent(in) :: minv,maxv
real(kind=rk), dimension(n), intent(out) :: phi 

call random_number(phi)

phi = minv + (maxv-minv)*phi
end Subroutine UniformRandomNumbers1D


Subroutine NormalRandomNumbers1D(mean,stdev,n,phi)
use random
IMPLICIT NONE
integer, parameter :: rk = kind(1.0D0)
integer, intent(in) :: n
real(kind=rk), INTENT(IN) :: mean, stdev
real(kind=rk), dimension(n), INTENT(OUT) :: phi
integer :: i

do i = 1,n
   phi(i) = mean+stdev*random_normal()
end do

End Subroutine NormalRandomNumbers1D
