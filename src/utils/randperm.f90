!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-05-21 13:32:24 pbrowne>
!!!
!!!    Subroutine to create an array for random permutations
!!!    Copyright (C) 2015 Mengbin Zhu 
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
!!!    Email: zhumengbin @ gmail.com  
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> subroutine to create an array of a random permutations of the
!! natural numbers from 1 to N
subroutine randperm(N, p)
  implicit none
  integer, intent(in)                :: N !< length of array P
  integer, dimension(N), intent(out) :: p !< on output, this is a
  !!random permutation of the integers from 1 to N

  integer :: i, j, k
  integer :: temp
  real(kind=kind(1.0d0))     :: u

  p = (/ (i, i=1,N) /)

  do j=N,2,-1

     call random_number(u)

     k = floor(j*u) + 1

     ! exchange p(k) and p(j)
     temp = p(k)
     p(k) = p(j)
     p(j) = temp

  end do

end subroutine randperm
