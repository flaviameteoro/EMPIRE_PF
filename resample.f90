subroutine resample
!this takes pf%phi with corresponding weights pf%weight
!which are stored as -log(w_i)
!and returns pf%phi which has been resampled with completely equal weights



use pf_control
use random
implicit none
integer, parameter :: rk = kind(1.0D0)
integer :: i,old,j,dupepos,dupe,particle,k
real(kind=rk), dimension(pf%ngrand) :: cumweights
integer, dimension(pf%ngrand) :: new,tobereplaced,brandnew
logical :: in
real(kind=rk) :: point,draw


pf%weight = exp(-pf%weight)
pf%weight = pf%weight/sum(pf%weight)

cumweights = 0.0_rk

cumweights(1) = pf%weight(1)
do i = 2,pf%ngrand
   cumweights(i) = cumweights(i-1) + pf%weight(i)
end do

print*,cumweights

call random_number(draw)
draw = draw/pf%ngrand

new = 0
old = 1

do particle = 1,pf%ngrand
   point = real(particle-1,rk)/pf%ngrand + draw
   print*,particle,point
   do
      if(point .lt. cumweights(old)) then
         new(particle) = old
         exit
      else
         old = old + 1
      end if
   end do
end do

tobereplaced = -1


do i = 1,pf%ngrand
   in = .false.
   do j= 1,pf%ngrand
      if(i .eq. new(j)) then
         in = .true.
         exit
      end if
   end do
   if(.not. in) tobereplaced(i) = i
end do

dupepos = 1
brandnew = (/ (i,i=1,pf%ngrand) /)
do i = 1,pf%ngrand
   if(tobereplaced(i) .ne. -1) then
      do k = dupepos,pf%ngrand
         if(new(k) .eq. new(k+1)) then
            dupe = new(k+1)
            dupepos = k+1
            exit
         end if
      end do
      brandnew(tobereplaced(i)) = dupe
   end if
end do
print*,'##################################'
print*,brandnew
print*,'##################################'
do particle = 1,pf%ngrand
   if(brandnew(particle) .ne. particle) then
      print*,'replacing particle ',particle,'with particle ',brandnew(particle)
      pf%psi(:,particle) = pf%psi(:,brandnew(particle))
   end if
end do

pf%weight = 1.0_rk/real(pf%ngrand,rk)
pf%weight = -log(pf%weight)
print*,'finished resampling'
print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
end subroutine resample
