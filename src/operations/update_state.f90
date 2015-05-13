!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-05-13 12:01:23 pbrowne>
!!!
!!!    Routines to update states
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


!> Subroutine to update the state
!!
!! This subroutine is here because, mathematically, in a particle
!! filter \f$x^{k+1} = f(x^k) + A^k + \xi^k\f$
!!
!! However sometimes the result needs to be bounded, some variables
!! need to be exactly related or maybe even something else.
!!
!! This can be changed for the specific model
!! if it needs to be, in order to bound variables etc.
!!
!! NOTE this the theory not mathematically correct.
!!
!! @param[in] fpsi deterministic model update \f$f(x^{n-1})\f$
!! @param[in] kgain nudging term
!! @param[inout] betan Stochastic term
!! @param[out] state The updated state vector 

subroutine update_state(state,fpsi,kgain,betan)
use sizes
implicit none
integer, parameter :: rk=kind(1.0D0)
real(kind=rk), dimension(state_dim), intent(out) :: state
real(kind=rk), dimension(state_dim), intent(in) :: fpsi,kgain
real(kind=rk), dimension(state_dim), intent(inout) :: betan
real(kind=rk) :: dnrm2
logical, parameter :: norms=.false.



!do the addition
state = fpsi+kgain+betan

if(norms) print*,' |kgain|= ',dnrm2(state_dim,kgain,1),' |betan| = ',dnrm2(state_dim,betan,1)


end subroutine update_state
