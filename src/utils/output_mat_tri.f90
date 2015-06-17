!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-06-17 14:38:39 pbrowne>
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

!> subroutine to output triangluar matrix various formats
subroutine output_mat_tri(n,A,filename,output_type)
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: n !< number of columns of matrix A
  real(kind=rk), intent(in), dimension(n*(n+1)/2) :: A !< matrix to be output in
  !!                                          rectangular full packed
  !!                                          format (TF)
  character(40), intent(in) :: filename !<   the name of the file to
  !!                                          be output
  integer, intent(in) :: output_type !< output file type. 
                              !! -  0 - undefined
                              !! -  1 - standard packed format (TP)
                              !! -  2 - rectangular full packed
                              !!                        format (TF)
                              !! \n
                              !! Negative values will be formatted.
                              !! \n
                              !! Positive values will be unformatted.


  real(kind=rk), dimension(n*(n+1)/2) :: Aout
  integer :: outunit=99
  logical :: opend
  character(11) :: fm
  integer :: err

  if(output_type .eq. 0) then
     write(*,*) 'Error in output_mat_tri. output_type = 0&
          & unsupported. Stopping.'
     stop '-4'
  elseif(output_type .gt. 0) then
     fm='unformatted'
  else
     fm='formatted'
  end if

  


  !ensure unit is not opened then open file=filename
  do
     inquire(unit=outunit,opened=opend)
     if(.not. opend) then
        open(outunit,file=filename,action='write',form=fm)
        exit
     else
        outunit=outunit-1
     end if
  end do

  select case(abs(output_type))
  case(1) ! standard packed format (TP)
     call dtfttp('N','U',n,A,Aout,err)
  case(2) ! rectangular full packed format (TF)
     Aout = A
  case default
     write(*,*) 'Error in output_mat_tri, unsupported output_type'
     write(*,*) 'output_type=',output_type,'. Stopping'
     stop '-4'
  end select


  if(output_type .gt. 0) then
     write(outunit) Aout
  else
     write(outunit,*) Aout
  end if

  close(outunit)


end subroutine output_mat_tri
