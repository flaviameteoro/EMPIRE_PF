!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:39:18 pbrowne>
!!!
!!!    subroutine be called by optimization codes
!!!    Copyright (C) 2015  Philip A. Browne
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

!> This is the subroutine which the optimization routines call
!! to get the objective function value and its gradient
subroutine fcn( n, x, f, g )
  use output_empire, only : emp_e
  use pf_control
  implicit none
  integer,intent(in) :: n !< the dimension of the optimzation problem
  real(kind=kind(1.0d0)), dimension(n), intent(in) :: x !< the
  !!current optimization state
  real(kind=kind(1.0d0)), intent(out) :: f !< the objective function value
  real(kind=kind(1.0d0)), dimension(n), intent(out) :: g !< the
  !!gradient of the objective function

  select case(pf%filter)
  case('3D')
     call threedvar_fcn(n,x,f,g)
  case default
     write(emp_e,*) 'wrong case in fcn'
     stop
  end select
!  print*,'function = ',f
!  print*,'gradient = ',g
end subroutine fcn
