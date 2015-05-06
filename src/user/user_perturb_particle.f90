!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-05-06 11:14:07 pbrowne>
!!!
!!!    Subroutine to implement user specific initial perturbations
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

!> Subroutine to perturb state vector as defined by the user
!! governed by the \link pf_control::pf_control_type::init init\endlink option
!!
!! This should be considered an example routine. Here I shall
!! implement a perturbation with a uniform variable on the interval
!! \f$[-10,15]\f$
subroutine user_perturb_particle(n,x)
  integer, parameter :: rk=kind(1.0D0)
  integer, intent(in) :: n !< the dimension of the state vector x
  real(kind=rk), dimension(n), intent(inout) :: x !< the state to be perturbed


  !specific variables needed for uniform perturbation example
  real(kind=rk), dimension(n) :: lower,upper,unif

  !set lower bound
  lower = -10.0d0
  
  !ser upper bound
  upper = 15.0d0

  !get random numbers in [0,1]
  call UniformRandomNumbers1D(0.0d0, 1.0d0, n,unif)
  
  !scale and shift unif in [0,1] to [lower,upper]
  unif = unif*(upper-lower) + lower

  !now add this perturbation to the state
  x = x + unif
 
  
end subroutine user_perturb_particle
