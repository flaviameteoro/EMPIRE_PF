subroutine nudge_data
  use comms
  
  implicit none
  integer, parameter :: rk= kind(1.0D0)
  real(kind=rk), allocatable, dimension(:) :: randomVec  !vector to store uncorrelated random error
  real(kind=rk), allocatable, dimension(:) :: corrRandom !vector to store sqrtQ correlated random error
  integer :: i,j,k,n
  integer :: alloc_error

  allocate(randomVec(nxn*nyn), STAT=alloc_error)
  allocate(corrRandom(3*nxn*nyn), STAT=alloc_error)
  

!  call NormalRandomNumbers1D(0.0,1.0,randomVec)
!  call MultiplySqrtQ(randomVec,corrRandom)  

  !This is an example of nudging the data
  do n=1,nens
     do j=1,levels
        do i=1,nxn
           !           psiGrand(i,j,n)=1./qrel*corrRandom((i-1)*nyn+j)
        enddo
     enddo
  enddo
  
end subroutine nudge_data
