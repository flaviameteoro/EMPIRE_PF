!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-09-18 15:59:18 pbrowne>
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
  !>
  !> in order to redirect the STDOUT used by EMPIRE, this subroutine
  !> will read from the file 'empire.nml'. If it exists, it looks for
  !> the namelist &empire_output, which consists of a single string
  !> up to 10 characters called 'basename' which will be read, and the
  !> STDOUT redirected
  !> to that string appended with the MPI rank of the EMPIRE process.
  !>
  !> In order to suppress most of the STDOUT from EMPIRE, this path
  !> can be set to a platform specific Null device:
  !> - Unix: /dev/null
  !> - MS: nul
  !>
  !> If you are running on any other system, please let me know what
  !> Null Device you would like to use, and we can add a check for it
  subroutine open_emp_o(id_num)
    implicit none
    integer, intent(in) :: id_num
    character(14) :: filename
    character(10) :: basename
    logical :: opend
    character(3), parameter :: msnul='nul'
    character(9), parameter :: unixnul='/dev/null'
    character(10), parameter :: emp_out_name='empire.nml'
    logical :: file_exists
    integer :: ios
    
    namelist/empire_output/basename

    inquire(file=emp_out_name,exist=file_exists)
    if(file_exists) then
       open(32,file=emp_out_name,iostat=ios,action='read'&
            &,status='old',form='formatted')
       if(ios .ne. 0) then
          print*,'Cannot open ',emp_out_name
          stop 'open_emp_o ERROR' 
       end if
       read(32,nml=empire_output,iostat=ios)
       close(32)
       if(ios .ne. 0) basename ='emp.out'
    else
       basename = 'emp.out'
    end if
    
    select case(basename)
    case(msnul)
       filename = msnul
    case(unixnul)
       filename = unixnul
    case default
       !define filename appropriately
       write(filename,'(A,A,i0)') trim(basename),'.',id_num
    end select

    print*,"this process' output will henceforth be directed to ",filename

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



