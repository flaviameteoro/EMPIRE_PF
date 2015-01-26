program call
implicit none
integer :: i

integer :: n !size of optimization state vector
integer :: method ! the type of nonlinear cg

real(kind=kind(1.0d0)), allocatable, dimension(:) :: x0

n = 2 !rosenbrock test function




!method = 1 ! FLETCHER-REEVES 
method = 2 ! POLAK-RIBIERE
!method = 3 ! POSITIVE POLAK-RIBIERE ( BETA=MAX{BETA,0} )

allocate(x0(n))
do i=1,n
   x0(i) = -2.d+00
end do




call subroutine_cg(method,n,x0)


deallocate(x0)
end program call
