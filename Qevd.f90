module Qevd
implicit none
integer :: Qnev
real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: QU
real(kind=kind(1.0D0)), allocatable, dimension(:) :: QD
contains
  subroutine loadQevd
    use sizes
    open(2,file='Qevd.dat',action='read',form='unformatted')
    
    read(2) Qnev
    allocate(QU(state_dim,Qnev),QD(Qnev))
    print*,'allocation of QU and QD done'
    
    read(2) QD
    read(2) QU
    close(2)
    print*,'loaded Qev'
!    QD = QD/1.0D2


  end subroutine loadQevd

  subroutine killQevd
    if(allocated(QU)) deallocate(QU)
    if(allocated(QD)) deallocate(QD)
  end subroutine killQevd

end module QEVD
