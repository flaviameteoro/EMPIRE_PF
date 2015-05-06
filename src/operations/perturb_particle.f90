!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-05-06 10:12:57 pbrowne>
!!!
!!!    Collection of routines to perturb states
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

  select case(pf%init)
  case('P')
     call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,rdom)
     call Qhalf(1,rdom,y)
     rdom = 0.0_rk
     kgain = 0.0_rk
     call update_state(rdom,x,kgain,y)
     x = rdom
  case('N')
     call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,rdom)
     kgain = 0.0_rk
     call update_state(y,x,kgain,rdom)
     x = y
  case('R')
     !get ensemble member from the restart folder                  
     write(filename,'(A,i2.2,A)') 'rstrt/',pfrank,'.state'
     call get_state(x,filename)
  case('S')
     !get ensemble member from the start folder                  
     if(pf%gen_data) then
        write(filename,'(A,i2.2,A)') 'start/',32,'.state'
     else
        write(filename,'(A,i2.2,A)') 'start/',pfrank,'.state'
        print*,'pf #',pfrank,' starting from ',filename
     end if
     call get_state(x,filename)
  case('B')
     call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,rdom)
     call Bhalf(1,rdom,y)
     rdom = 0.0_rk
     kgain = 0.0_rk
     call update_state(rdom,x,kgain,y)
     x = rdom
  case('U')
     call user_perturb_particle(state_dim,x)
  case('Z')
     !no perturbation. x remains as is
  case default
     print*,'ERROR: incorrect pf%init selected in perturb_particle'
     stop
  end select


end subroutine perturb_particle

