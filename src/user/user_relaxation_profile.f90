!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-11-23 11:36:57 pbrowne>
!!!
!!!    Subroutine to compute the relaxation strength via something
!!!    the user has coded
!!!    Copyright (C) 2016 Philip A. Browne
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
!> subroutine to compute the relaxation strength that the user codes
subroutine user_relaxation_profile(tau,p,zero)
  use output_empire, only : emp_e
  use pf_control
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), intent(in) :: tau !< the pseudotime between
  !! observations.  \$f 0 \le \tau \le 1 \$f
  real(kind=rk), intent(out) :: p !< the relaxation strength at time tau
  logical, intent(out) :: zero !< switch if false is no relaxation
  write(emp_e,*) 'EMPIRE ERROR: YOU HAVE CALLED A USER DEFINED RELAXAT&
       &ION PROFILE'
  write(emp_e,*) 'EMPIRE ERROR: BUT YOU HAVE NOT CODED THIS YOURSELF.'
  write(emp_e,*) 'EMPIRE ERROR: STOPPING. PLEASE CODE src/user/user_re&
       &laxation_profile.f90'
  stop 'EMPIRE ERROR IN USER_RELAXATION_PROFILE'
end subroutine user_relaxation_profile
