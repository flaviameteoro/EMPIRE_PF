!> Subroutine to initialise mpi in a special way if the model is
!! weird like HadCM3 for example
subroutine user_initialise_mpi
  use comms
  use pf_control
  implicit none

end subroutine user_initialise_mpi


subroutine user_mpi_send(stateDim,nrhs,x,tag)
  implicit none
  integer, intent(in) :: stateDim
  integer, intent(in) :: nrhs
  real(kind=kind(1.0d0)), intent(in), dimension(stateDim,nrhs) :: x
  integer, intent(in) :: tag
end subroutine user_mpi_send

subroutine user_mpi_recv(stateDim,nrhs,x)
  implicit none
  integer, intent(in) :: stateDim
  integer, intent(in) :: nrhs
  real(kind=kind(1.0d0)), intent(out), dimension(stateDim,nrhs) :: x
end subroutine user_mpi_recv

subroutine user_mpi_irecv(stateDim,nrhs,x,requests)
  implicit none
  integer, intent(in) :: stateDim
  integer, intent(in) :: nrhs
  real(kind=kind(1.0d0)), intent(out), dimension(stateDim,nrhs) :: x
  integer, dimension(nrhs), intent(inout) :: requests
end subroutine user_mpi_irecv

