subroutine objective_gradient(n,x,g)
implicit none
integer, parameter :: rk = kind(1.0d0)
integer, intent(in) :: n
real(kind=rk), dimension(n), intent(in) :: x
real(kind=rk), dimension(n), intent(out) :: g

real(kind=rk), dimension(2) :: fullx
integer :: mpi_rank,mpi_err
include 'mpif.h'
integer :: status(MPI_STATUS_SIZE)

! Rosenbrock function gradient
call mpi_comm_rank(MPI_COMM_WORLD,mpi_rank,mpi_err)

select case(mpi_rank)
case(0)
   fullx(1) = x(1)
   call mpi_send(fullx(1),1,MPI_DOUBLE_PRECISION,1,1,MPI_COMM_WORLD&
        &,mpi_err)
   call mpi_recv(fullx(2),1,MPI_DOUBLE_PRECISION,1,1,MPI_COMM_WORLD&
        &,status,mpi_err)
   g(1) = 200*(fullx(2) - fullx(1)**2)*(-2*fullx(1)) - 2*(1 - fullx(1))
case(1)
   fullx(2) = x(1)
   call mpi_recv(fullx(1),1,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD&
        &,status,mpi_err)
   call mpi_send(fullx(2),1,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD&
        &,mpi_err)
   g(1) = 200*(fullx(2) - fullx(1)**2)
case default
   print*,'too many processes launched for the model. Error'
   stop
end select



end subroutine objective_gradient
