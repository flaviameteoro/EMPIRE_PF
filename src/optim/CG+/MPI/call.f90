program call
implicit none
integer :: i

integer :: n !size of optimization state vector
integer :: method ! the type of nonlinear cg

real(kind=kind(1.0d0)), allocatable, dimension(:) :: x0
real(kind=kind(1.0d0)) :: epsin=0.1d0

include 'mpif.h'
integer :: mpi_err,mpi_size


call mpi_init(mpi_err)
call mpi_comm_size(MPI_COMM_WORLD,mpi_size,mpi_err)


n = 1 !rosenbrock test function, one variable per process




!method = 1 ! FLETCHER-REEVES 
method = 2 ! POLAK-RIBIERE
!method = 3 ! POSITIVE POLAK-RIBIERE ( BETA=MAX{BETA,0} )

allocate(x0(n))
do i=1,n
   x0(i) = -2.d+00
end do




call subroutine_cg(method,n,epsin,x0,MPI_COMM_WORLD,mpi_size)


deallocate(x0)

call mpi_finalize(mpi_err)

end program call
