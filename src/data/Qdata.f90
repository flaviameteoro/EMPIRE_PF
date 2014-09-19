!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-18 14:15:37 pbrowne>
!!!
!!!    {one line to give the program's name and a brief idea of what it does.}
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
module Qdata
implicit none
integer :: Qn,Qne
integer, allocatable, dimension(:) :: Qrow,Qcol
real(kind=kind(1.0D0)), allocatable, dimension(:) :: Qval,Qdiag
real(kind=kind(1.0d0)) :: Qscale
contains
  subroutine loadQ
    use sizes
    use pf_control

  end subroutine loadQ

  subroutine killQ
    if(allocated(Qrow)) deallocate(Qrow)
    if(allocated(Qcol)) deallocate(Qcol)
    if(allocated(Qval)) deallocate(Qval)
    if(allocated(Qdiag)) deallocate(Qdiag)
  end subroutine killQ

end module Qdata
