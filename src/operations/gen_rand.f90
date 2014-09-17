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

Subroutine NormalRandomNumbers2D(mean,stdev,n,k,phi)
use random
IMPLICIT NONE
integer, parameter :: rk = kind(1.0D0)
integer, intent(in) :: n,k
real(kind=rk), INTENT(IN) :: mean, stdev
real(kind=rk), dimension(n,k), INTENT(OUT) :: phi
integer :: i,j

do j = 1,k
   do i = 1,n
      phi(i,j) = mean+stdev*random_normal()
   end do
end do
End Subroutine NormalRandomNumbers2D

subroutine MixtureRandomNumbers1D(mean,stdev,ufac,epsi,n,phi,uniform)
use random
implicit none
real(kind=kind(1.0D0)), intent(in) :: mean,stdev,ufac,epsi
integer, intent(in) :: n
real(kind=kind(1.0D0)), dimension(n), intent(out) :: phi
logical, intent(out) :: uniform
real(kind=kind(1.0D0)) :: draw

call random_number(draw)

if(draw .gt. epsi) then
   call UniformRandomNumbers1D(mean-ufac,mean+ufac, n,phi)
   uniform = .true.
else
   call NormalRandomNumbers1D(mean,stdev,n,phi)
   uniform = .false.
end if

end subroutine MixtureRandomNumbers1D

subroutine MixtureRandomNumbers2D(mean,stdev,ufac,epsi,n,k,phi,uniform)
use random
implicit none
real(kind=kind(1.0D0)), intent(in) :: mean,stdev,ufac,epsi
integer, intent(in) :: n,k
real(kind=kind(1.0D0)), dimension(n,k), intent(out) :: phi
logical, dimension(k), intent(out) :: uniform
real(kind=kind(1.0D0)) :: draw
integer :: i

do  i = 1,k
   call random_number(draw)

   if(draw .gt. epsi) then
      call UniformRandomNumbers1D(mean-ufac,mean+ufac,n,phi(:,i))
      uniform(i) = .true.
   else
      call NormalRandomNumbers1D(mean,stdev,n,phi(:,i))
      uniform(i) = .false.
   end if
end do
end subroutine MixtureRandomNumbers2D

subroutine random_seed_mpi(pfid)
  integer, intent(in) :: pfid

  integer :: n
  integer, allocatable, dimension(:) :: seed

  call random_seed(SIZE=n)
  allocate(seed(n))
  call random_seed(GET=seed)
  !add the particle filter id to the seed to make it
  !independent on each process
  seed = seed+pfid
  call random_seed(PUT=seed)
  deallocate(seed)

end subroutine random_seed_mpi
