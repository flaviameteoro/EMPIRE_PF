module arse
implicit none
integer :: Qn,Qne
integer, allocatable, dimension(:) :: Qrow,Qcol
real(kind=kind(1.0D0)), allocatable, dimension(:) :: Qval,Qdiag
contains
  subroutine loadQ
!    use pf_control
    integer :: i,state_dim
!    real(kind=kind(1.0d0)) :: eig
!!$    Qn = state_dim
!!$    Qne = state_dim
!!$    allocate(Qrow(Qne),Qcol(Qne),Qval(Qne))
!!$    Qrow = (/ (i, i = 1,Qne) /)
!!$    Qcol = (/ (i, i = 1,Qne) /)
!!$    Qval = 1.0D-8

    print*,'fuck sake'
    open(2,file='Qdata.dat',action='read',form='unformatted')
    
    read(2) state_dim
    print*,'state_dim =',state_dim
    Qn = state_dim
    read(2) Qne
    print*,'Qne = ',Qne
    allocate(Qrow(Qne),Qcol(Qne),Qval(Qne))
    print*,'allocation of Qrow,Qcol and Qval done'
    
    do i = 1,Qne
!       print*,i
       read(2) Qval(i)
    end do
    do i = 1,Qne
       read(2) Qrow(i)
    end do
    do i = 1,Qne
       read(2) Qcol(i)
    end do
    close(2)


  end subroutine loadQ
end module Arse


program matrixchange
  use arse
integer :: i,c,a
integer, allocatable, dimension(:) :: Srow,Scol
real(kind=kind(1.0d0)), allocatable, dimension(:) :: Sval
print*,'pissing fuck'
call loadQ

allocate(Qdiag(Qn))
allocate(Srow(Qne-Qn),Scol(Qne-Qn),Sval(Qne-Qn))
c = 0
a = 0
do i = 1,Qne
   if(Qrow(i) .eq. Qcol(i)) then
      a = a + 1
      Qdiag(a) = Qval(i)
   else
      c = c + 1
      Sval(c) = Qval(i)
      Srow(c) = Qrow(i)
      Scol(c) = Qcol(i)
   end if
end do
print*,'c = ',c


open(2,file='NEWQdata.dat',action='write',form='unformatted')
write(2) Qn
write(2) c
do i = 1,Qn
   write(2) Qdiag(i)
end do
do i = 1,c
   write(2) Sval(i)
end do
do i = 1,c
   write(2) Srow(i)
end do
do i = 1,c
   write(2) Scol(i)
end do
close(2)
end program matrixchange

