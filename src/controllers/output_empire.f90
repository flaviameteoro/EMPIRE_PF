!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-07-14 14:22:25 pbrowne>
!!!
!!!    Module that stores information about outputting from empire
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

!> @brief Module that stores the information about the outputting
!! from empire
Module output_empire
  implicit none
  integer,parameter :: emp_o=6 !< the output stream number

contains

  !> subroutine to open the file for outputting
  subroutine open_emp_o(id_num)
    implicit none
    integer, intent(in) :: id_num
    character(14) :: filename
    logical :: opend
    
    !define filename appropriately
    write(filename,'(A,i0)') 'emp.out.',id_num
    

    inquire(unit=emp_o,opened=opend)
    if(opend) then
       close(emp_o)
    end if

    open(emp_o,file=filename,action='write',status='replace'&
         &,form='formatted')
        
  end subroutine open_emp_o


  !> subroutine to close the output file
  subroutine close_emp_o()
    close(emp_o)
  end subroutine close_emp_o
End Module Output_empire



