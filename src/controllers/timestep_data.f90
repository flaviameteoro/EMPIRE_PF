!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-07-22 11:47:22 pbrowne>
!!!
!!!    Module that stores information about timestepping
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

!> @brief Module that stores the information about the timestepping
!!process
Module timestep_data
  implicit none
  type, public :: timestep_data_type
     integer :: total_timesteps     !< total number of timesteps that the
                                    !! model will run
     integer :: current_timestep    !< the current timestep that
                                    !! empire is running 
     integer :: completed_timesteps !< the number of timesteps that
                                    !! empire has so far finished
     logical :: is_analysis         !< if true, then the current ensemble is
                                    !! an analysis. If false then the
                                    !! current ensemble is not an analysis
     logical :: do_analysis         !< if true then on this timestep
                                    !! we are required to do an
                                    !! analysis. If false we do not
                                    !! have an observation at this timestep
     integer, allocatable, dimension(:) :: obs_times !< an integer
                                    !! array that will hold a mapping
                                    !! from observation number in
                                    !! time to model timesteps. I.e.
                                    !! obs_times(i) is the timestep
                                    !! of observation i in time.     
     integer :: tau                 !< the pseudotimestep between observations
  end type timestep_data_type
  type(timestep_data_type), save :: TSData !< the derived data type
  !!holding all timestep data

  contains
    !> subroutine to allocate space for obs_times array
    subroutine timestep_data_allocate_obs_times(n)
      implicit none
      integer, intent(in) :: n
      integer :: st
      allocate(TSData%obs_times(n),stat=st)
      if(st .ne. 0) then
         print*,'ERROR: allocation of obs_times in timestep data'
         print*,'ERROR: STOPPING.'
         stop '-7'
      end if
      TSData%obs_times = -1
    end subroutine timestep_data_allocate_obs_times

    !> subroutine to deallocate obs_times array
    subroutine timestep_data_deallocate_obs_times
      deallocate(TSData%obs_times)
    end subroutine timestep_data_deallocate_obs_times

    !> subroutine to set the timestep corresponding to the
    !! observation number in time
    subroutine timestep_data_set_obs_times(obs_num_in_time,timestep)
      implicit none
      integer, intent(in) :: obs_num_in_time
      integer, intent(in) :: timestep

      TSData%obs_times(obs_num_in_time) = timestep
    end subroutine timestep_data_set_obs_times

    !> subroutine to extract the timestep corresponding to the
    !! observation number in time
    subroutine timestep_data_get_obs_times(obs_num_in_time,timestep)
      implicit none
      integer, intent(in) :: obs_num_in_time
      integer, intent(out) :: timestep
      timestep = TSData%obs_times(obs_num_in_time)
    end subroutine timestep_data_get_obs_times
    
    !> subroutine to define if the current timestep should perform an analysis
    subroutine timestep_data_set_do_analysis
      implicit none
      TSData%do_analysis = .true.
    end subroutine timestep_data_set_do_analysis

    !> subroutine to define if the current timestep should not
    !! perform an analysis
    subroutine timestep_data_set_do_no_analysis
      implicit none
      TSData%do_analysis = .false.
    end subroutine timestep_data_set_do_no_analysis
    
    !> subroutine to define if the current ensemble is an analysis
    subroutine timestep_data_set_is_analysis
      implicit none
      TSData%is_analysis = .true.
    end subroutine timestep_data_set_is_analysis
    
    !> subroutine to define if the current ensemble is not an analysis
    subroutine timestep_data_set_no_analysis
      implicit none
      TSData%is_analysis = .false.
    end subroutine timestep_data_set_no_analysis

    !> subroutine to define the number of completed timesteps
    subroutine timestep_data_set_completed(t)
      implicit none
      integer, intent(in) :: t
      TSData%completed_timesteps = t
    end subroutine timestep_data_set_completed

    !> subroutine to define the current timestep
    subroutine timestep_data_set_current(t)
      implicit none
      integer, intent(in) :: t
      TSData%current_timestep = t
    end subroutine timestep_data_set_current

    !> subroutine to define the total number of timesteps that the
    !! model will run for
    subroutine timestep_data_set_total(t)
      implicit none
      integer, intent(in) :: t
      TSData%total_timesteps = t
    end subroutine timestep_data_set_total

    !> subroutine to define the current number of timesteps between observations
    subroutine timestep_data_set_tau(pseudotimestep)
      implicit none
      integer, intent(in) :: pseudotimestep
      TSData%tau = pseudotimestep
    end subroutine timestep_data_set_tau

    
End Module Timestep_data



