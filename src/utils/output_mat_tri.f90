!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-05-08 14:55:02 pbrowne>
!!!
!!!    Subroutine to output triangular matrix
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

!> subroutine to output triangluar matrix in packed storage form
subroutine output_mat_tri(n,A,filename)
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: n !< number of columns of matrix A
  real(kind=rk), dimension(n*(n+1)/2) :: A !< matrix to be output in
  !!                                          packed storage
  character(40), intent(in) :: filename !<   the name of the file to
  !!                                          be output

  integer :: outunit=99
  logical :: opend
  

  !ensure unit is not opened then open file=filename
  do
     inquire(unit=outunit,opened=opend)
     if(.not. opend) then
        open(outunit,file=filename,action='write',form='unformatted')
        exit
     else
        outunit=outunit-1
     end if
  end do

  write(outunit) A

  close(outunit)


end subroutine output_mat_tri
