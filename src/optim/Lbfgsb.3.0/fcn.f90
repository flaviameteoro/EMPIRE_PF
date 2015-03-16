subroutine fcn( n, x, f, g )
implicit none
integer :: n
real(kind=kind(1.0d0)), dimension(n), intent(in) :: x
real(kind=kind(1.0d0)), intent(out) :: f
real(kind=kind(1.0d0)), dimension(n), intent(out) :: g

! Rosenbrock 

call objective_function(n,x,f) !get the objective function evaluated at x

call objective_gradient(n,x,g) !get gradient of the objective function at x

end subroutine fcn

