!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-04-03 10:48:56 pbrowne>
!!!
!!!    Collection of routines to perturb and update states
!!!    Copyright (C) 2014  Philip A. Browne
!!!
!!!    This program is free software: you can redistribute it and/or modify
!!!    it under the terms of the GNU General Public License as published by
!!!    the Free Software Foundation, either version 3 of the License, or
!!!    (at your option) any later version.
!!!
!!!    This program is distributed in the hope that it will be useful,
!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!!    GNU General Public License for more details.
!!!
!!!    You should have received a copy of the GNU General Public License
!!!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!!
!!!    Email: p.browne @ reading.ac.uk
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> Subroutine to perturb state vector 
!! governed by the \link pf_control::pf_control_type::init init\endlink option
subroutine perturb_particle(x)
use sizes
use comms
use pf_control
integer, parameter :: rk=kind(1.0D0)
real(kind=rk), dimension(state_dim), intent(inout) :: x
real(kind=rk), dimension(state_dim) :: rdom,y,kgain
character(14) :: filename

if(pf%init .eq. 'P') then
   call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,rdom)
   call Qhalf(1,rdom,y)
   rdom = 0.0_rk
   kgain = 0.0_rk
   call update_state(rdom,x,kgain,y)
   x = rdom
elseif(pf%init .eq. 'N') then

!!$   x = (/4.82160531136,-1.77565265252,2.10461993106,3.46854171351&
!!$        &,7.19135143229,7.12754206422,-0.796082655036,2.45172556459&
!!$        &,-0.531641795086,-0.19366799976,-4.17678736568,6.94259660694&
!!$        &,4.67019291167,-1.63118486372,0.830259694341,1.74204790657&
!!$        &,5.63332497432,1.38456126234,-1.67858158117,2.15882451035&
!!$        &,8.24797269894,-1.18246777768,-3.58939011025,3.26473996747&
!!$        &,8.83474281013,3.34455161902,2.31141994205,-2.84763733923&
!!$        &,1.87634007581,3.3191730343,7.36150967249,3.99624410686 &
!!$        &,-7.63041979612,1.43836147006,-0.955124293143,5.38338216121&
!!$        &,2.92386476709,1.14393995662,3.17942949736,6.00540694151/)

   call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,rdom)
   kgain = 0.0_rk
   call update_state(y,x,kgain,rdom)
   x = y

elseif(pf%init .eq. 'R') then
   !get ensemble member from the restart folder                  
   write(filename,'(A,i2.2,A)') 'rstrt/',pfrank,'.state'
   call get_state(x,filename)
elseif(pf%init .eq. 'S') then
   !get ensemble member from the start folder                  
   if(pf%gen_data) then
      write(filename,'(A,i2.2,A)') 'start/',32,'.state'
   else
      write(filename,'(A,i2.2,A)') 'start/',pfrank,'.state'
      print*,'pf #',pfrank,' starting from ',filename
   end if
   call get_state(x,filename)
elseif(pf%init .eq. 'B') then
   call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,rdom)
   call Bhalf(1,rdom,y)
   rdom = 0.0_rk
   kgain = 0.0_rk
   call update_state(rdom,x,kgain,y)
   x = rdom
else
   print*,'ERROR: incorrect pf%init selected in perturb_particle'
   stop
end if


end subroutine perturb_particle

