!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-08-16 15:22:48 pbrowne>
!!!
!!!    Module to control what variables are used to generate rank histograms
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

!> Module to control what variables are used to generate rank histograms
module histogram_data
  integer, allocatable, dimension(:) :: rank_hist_list
  integer, allocatable, dimension(:) :: rank_hist_nums
  integer :: rhl_n,rhn_n

contains
  !> subroutine to read from variables_hist.dat which 
  !! holds the variables to be used to make the rank histograms
  !!
  !! In order for histograms to be output, the file
  !! "variables_hist.dat" must contain the following infomation:
  !!
  !! - rhn_n -- the number of different rank histograms to be output
  !!
  !! - the numbers of variables to be included in each rank histogram
  !!
  !! - the index of the state vector for each different variable in
  !! each different rank histogram, grouped by the different histograms
  !!
  !! So as an example, suppose we wanted to produce 3 rank histograms,
  !! the first relating to the 10th, and 16th variables in the
  !! state vector, the second containing the 1st, 2nd, 56th and 98th
  !! variables of the state vector and the final rank histogram
  !! relating to the 6th, 11th, 19th, 45th and 32nd variables. Then
  !! variables_hist.dat would look as follows:
  !!
  !! \verbatim
  !! 3
  !! 2
  !! 4
  !! 5
  !! 10
  !! 16
  !! 1
  !! 2
  !! 56
  !! 98
  !! 6
  !! 11
  !! 19
  !! 45
  !! 32
  !! \endverbatim

  subroutine load_histogram_data
    use output_empire, only : unit_hist_read
    use comms
    use sizes
    implicit none
    integer :: i
    integer :: counter,lowerbound,upperbound
    integer, allocatable, dimension(:) :: tempvar

!    rhn_n = 9
    open(unit_hist_read,file='variables_hist.dat',action='read',status='old')
    read(unit_hist_read,'(i7.7)') rhn_n
    allocate(rank_hist_nums(rhn_n))
    do i = 1,rhn_n
       read(unit_hist_read,'(i7.7)') rank_hist_nums(i)
    end do
    rhl_n = sum(rank_hist_nums)
    allocate(rank_hist_list(rhl_n))
    do i = 1,rhl_n
       read(unit_hist_read,'(i7.7)') rank_hist_list(i)
    end do
    close(unit_hist_read)


    if(comm_version .eq. 3) then

       lowerbound = state_displacements(pf_member_rank+1)
       if(pf_member_rank .eq. pf_member_size-1) then
          upperbound = state_dim_g
       else
          upperbound = state_displacements(pf_member_rank+2)
       end if
       
       counter = 0
       do i = 1,rhl_n
          if(rank_hist_list(i) .gt. lowerbound .and. rank_hist_list(i) .le.&
            & upperbound) then
             counter = counter + 1
          end if
       end do
       
       allocate(tempvar(rhl_n))
       tempvar = rank_hist_list
       deallocate(rank_hist_list)
       allocate(rank_hist_list(counter))
       
       counter = 0
       do i = 1,rhl_n
          if(tempvar(i) .gt. lowerbound .and. tempvar(i) .le.&
               & upperbound) then
             counter = counter + 1
             rank_hist_list(counter) = tempvar(i)
          end if
       end do
       
       deallocate(tempvar)
       rank_hist_list = rank_hist_list - lowerbound
    end if
    
  end subroutine load_histogram_data

  !>subroutine to clean up arrays used in rank histograms
  subroutine kill_histogram_data
    if(allocated(rank_hist_list)) deallocate(rank_hist_list)
    if(allocated(rank_hist_nums)) deallocate(rank_hist_nums)
  end subroutine kill_histogram_data
end module histogram_data
