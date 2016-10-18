!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:40:19 pbrowne>
!!!
!!!    Subroutine to output trajectories of state variables
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
!> module to hold data for trajectories
module traj_data
  integer :: trajn
  integer, allocatable, dimension(:) :: trajvar
  character(28), parameter :: traj_list='traj_list.dat'
  contains
    !> subroutine to read in which trajectories are required
    !!
    !! this requires that the directory traj/ exists before runtime.
    !!
    !! Then this reads the file @ref traj_list .
    !!
    !! The format for @ref traj_list is a list of K+1 integers,
    !! 
    !! where the first integer is K
    !!
    !! and the following K integers are the index in the state
    !! dimension for which the trajectories are required.    
    subroutine setup_traj
      use output_empire, only : unit_traj_read,emp_e
      use comms
      use sizes, only : state_dim_g
      logical :: dir_e,file_e
      integer :: i
      integer :: counter,lowerbound,upperbound
      integer, allocatable, dimension(:) :: tempvar
      !first check that the traj directory exists:
      ! a trick to be sure traj is a dir
      inquire( file="./traj/.", exist=dir_e )
      if ( .not. dir_e ) then
         write(emp_e,*) 'EMPIRE ERROR -559: ./traj/ directory does not exi&
              &st'
         write(emp_e,*) 'Please create the directory traj/ in the folder t&
              &hat you wish to run'
         stop '-559'
      end if

      !now we check to see that traj_list.dat exists:
      inquire( file=traj_list, exist=file_e )
      if ( .not. file_e ) then
         write(emp_e,*) 'EMPIRE ERROR -560: file ',traj_list,' does not exi&
              &st'
         write(emp_e,*) 'This file should contain a list of state variable&
              &s'
         write(emp_e,*) 'for which you want to output trajectories'
         stop '-560'
      end if

      open(unit_traj_read,file=traj_list,action='read')
      read(unit_traj_read,*) trajn
      
      allocate(trajvar(trajn))

      do i = 1,trajn
         read(unit_traj_read,*) trajvar(i)
         if(trajvar(i) .le. 0) then
            write(emp_e,*) 'EMPIRE ERROR -561: trajectory variable ',i,' less &
                 &than 0'
            write(emp_e,*) '                 : variable read as',trajvar(i),' &
                 &STOP.'
            stop '-561'
         elseif(trajvar(i) .gt. state_dim_g) then
            write(emp_e,*) 'EMPIRE ERROR -562: trajectory variable ',i,' larger &
                 &than the state dimension ',state_dim_g
            write(emp_e,*) '                 : variable read as',trajvar(i),' &
                 &STOP.'
            stop '-562'
         end if
      end do

      close(unit_traj_read)


      if(comm_version .eq. 3) then

         lowerbound = state_displacements(pf_member_rank+1)
         if(pf_member_rank .eq. pf_member_size-1) then
            upperbound = state_dim_g
         else
            upperbound = state_displacements(pf_member_rank+2)
         end if

         counter = 0
         do i = 1,trajn
            if(trajvar(i) .gt. lowerbound .and. trajvar(i) .le.&
                 & upperbound) then
               counter = counter + 1
            end if
         end do

         allocate(tempvar(trajn))
         tempvar = trajvar
         deallocate(trajvar)
         allocate(trajvar(counter))
         
         counter = 0
         do i = 1,trajn
            if(tempvar(i) .gt. lowerbound .and. tempvar(i) .le.&
                 & upperbound) then
               counter = counter + 1
               trajvar(counter) = tempvar(i)
            end if
         end do
         deallocate(tempvar)
         trajvar = trajvar - lowerbound

      end if
      
    end subroutine setup_traj
    
    subroutine deallocate_traj
      if(allocated(trajvar)) deallocate(trajvar)
    end subroutine deallocate_traj
end module traj_data


!> subroutine to output trajectories
!!
subroutine trajectories
  use output_empire, only : unit_traj_write
  use traj_data
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  integer :: particle,i,j
  character(28) :: filename

  if(pf%timestep .eq. 0) call setup_traj

  do i = 1,pf%count
     particle = pf%particles(i)
     
     do j = 1,trajn
        if(pf%gen_data) then
           write(filename,'(A,i7.7)') 'traj/truth_var',trajvar(j)
        else
           write(filename,'(A,i7.7,A,i5.5)') 'traj/var',trajvar(j),'particle',particle
        end if
        if(pf%timestep .eq. 0) then
           open(unit_traj_write,file=filename,action='write',status='replace')
        else
           open(unit_traj_write,file=filename,action='write',status='old',position='append')
        end if
        write(unit_traj_write,'(es22.15)') pf%psi(trajvar(j),i)
        close(unit_traj_write)
        
     end do
  end do
end subroutine trajectories
