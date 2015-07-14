!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-07-14 10:56:54 pbrowne>
!!!
!!!    subroutine to compute a localisation weighting based on distance
!!!    Copyright (C) 2015 Philip A. Browne
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


!> subroutine to compute a localisation weighting based on a distance
!! @todo include multiple localisation functions such as Gaspari-Cohn ones
subroutine loc_function(loctype,dis,scal,inc)
  use pf_control
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: loctype !< the choice of localisation
  !! function
  real(kind=rk), intent(in) :: dis !< the input distance
  real(kind=rk), intent(out) :: scal !< the localisation weighting to
  !! use
  logical, intent(out) :: inc !< logical signifying whether to
  !! compute with this sized weight or not.
  real(kind=rk), parameter :: minscal_gaussian = exp(-8.0d0)
  

  select case(loctype)
  case(1) ! exponential function
     scal = exp(-(dis**2)/(2.0_rk*pf%len**2))
     if(scal .gt. minscal_gaussian) then
        inc = .true.
     else
        inc = .false.
     end if
  case default
     print*,'EMP: ERROR: loctype not supported in subroutine loc_funct&
          &ion'
     stop '-1'
  end select


end subroutine loc_function
