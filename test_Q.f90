program test_Q
use sizes
use pf_control
use Qdata
implicit none
integer, parameter :: rk = kind(1.0d0)

!integer :: ne
!integer, allocatable, dimension(:) :: row,col
real(kind=rk), allocatable, dimension(:) :: normaln,b
integer :: i,j
logical :: err
call set_pf_controls


call loadQ
print*,'Q is now loaded'

print*,'check the rows'
do i = 1,Qne
   if(Qrow(i) .lt. 1) print*,i,Qrow(i),Qcol(i),Qval(i)
end do

print*,'check the columns'
do i  = 1,Qne
   if(Qcol(i) .lt. 1) print*,i,Qrow(i),Qcol(i),Qval(i)
end do
stop




allocate(normaln(state_dim),b(state_dim))
print*,'Allocation of normaln and b done'

call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,normaln)
print*,'normaln generated'

err = .true.
j = 1
do 
   print*,'iteration = ',j
   call Qhalf(normaln,b,err)
   
   if(.not. err) then
      print*,'Qhalf succeeded'
      exit
   else
      do i = 1,Qne
         if(Qcol(i) .eq. Qrow(i)) Qval(i) = 2.5d0*Qval(i)
      end do
   end if
   j = j+1
end do
call killQ
end program test_Q
