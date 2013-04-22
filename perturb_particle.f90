subroutine perturb_particle(x)
use sizes
integer, parameter :: rk=kind(1.0D0)
real(kind=rk), dimension(state_dim), intent(inout) :: x
real(kind=rk), dimension(state_dim) :: rdom

call NormalRandomNumbers1D(0.0D0,sqrt(2.0D0),state_dim,rdom)
x = x + rdom

end subroutine perturb_particle
