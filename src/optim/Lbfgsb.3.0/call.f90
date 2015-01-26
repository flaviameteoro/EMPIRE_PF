program call
implicit none
integer :: i

integer :: n !size of optimization state vector
integer :: method ! the type of nonlinear cg

real(kind=kind(1.0d0)), allocatable, dimension(:) :: x0
integer, allocatable, dimension(:) :: nbd
real(kind=kind(1.0d0)), allocatable, dimension(:) :: l
real(kind=kind(1.0d0)), allocatable, dimension(:) :: u

n = 25 ! test function




!method = 1 ! FLETCHER-REEVES 
method = 2 ! POLAK-RIBIERE
!method = 3 ! POSITIVE POLAK-RIBIERE ( BETA=MAX{BETA,0} )

allocate(x0(n))
do i=1,n
   x0(i) = 3.0d0
end do




call lbfgs_sub(n,x0)
print*,'result of lbfgs = ',x0


allocate(nbd(n),l(n),u(n))

do i=1,n
   x0(i) = 3.0d0
end do


nbd = 0
do 10 i=1, n, 2
   nbd(i) = 2
   l(i)   = 1.0d0
   u(i)   = 1.0d2
10 continue

!     Next set bounds on the even-numbered variables.
   
do 12 i=2, n, 2
   nbd(i) =  2
   l(i)   = -1.0d2
   u(i)   =  1.0d2
12 continue

call lbfgsb_sub(n,x0,nbd,l,u)
print*,'result of lbfgsb = ',x0

deallocate(x0)
end program call
