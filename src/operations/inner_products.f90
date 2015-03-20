!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-03-16 10:19:28 pbrowne>
!!!
!!!    Collection of inner product wrappers
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


!>subroutine to compute the inner product with \f$R^{-1}\f$
!! @param[in] y multiple vectors in observation space (pf\%count of them)
!! @param[in] n length of each vector in y
!! @param[in] c number of vectors in y
!! @param[in] t current timestep
!! @param[out] w multiple scalars (pf\%count) where w(i) has the value
!! \f$y(:,i)^TR^{-1}y(:,i)\f$
subroutine innerR_1(y,n,c,w,t)
  !subroutine to take an observation vector y and return w = y^T R^(-1) y
!  use sizes
!  use pf_control
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: n
  integer, intent(in) :: c
  integer, intent(in) :: t
  real(kind=rk), dimension(n,c), intent(in) :: y
  real(kind=rk), dimension(n,c) :: v
  real(kind=rk), dimension(c), intent(out) :: w
  real(kind=rk) :: ddot
  integer :: i

  call solve_r(n,c,y,v,t)

  !this can defo be done better using BLAS PAB...
  do i = 1,c
!!$     w(i) = sum(y(:,i)*v(:,i))
     w(i) = ddot(n,y(:,i),1,v(:,i),1)
  end do

end subroutine innerR_1


!>subroutine to compute the inner product with \f$(HQH^T+R)^{-1}\f$      
!! @param[in] y vector in observation space
!! @param[out] w scalar with value \f$y^TR^{-1}y\f$
!! @param[in] t current timestep
subroutine innerHQHt_plus_R_1(y,w,t)
  !subroutine to take an observation vector y and return w = y^T (HQH^T+R)^(-1) y
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim) :: v
  real(kind=rk), intent(out) :: w
  integer, intent(in) :: t

  call solve_hqht_plus_r(obs_dim,y,v,t)

  !this can defo be done better using BLAS PAB...
  w = sum(y*v)

end subroutine innerHQHt_plus_R_1
