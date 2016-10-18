!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:36:02 pbrowne>
!!!
!!!    Subroutine to output RMSE
!!!    Copyright (C) 2015 Philip A. Browne
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
!> subroutine to output RMSEs
!>
subroutine output_spatial_rmse(mean)
  use output_empire, only : unit_spatial_rmse,emp_e
  use pf_control
  use timestep_data
  use sizes
  use comms, only : comm_version,pf_ens_rank
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(state_dim), intent(in) :: mean
  real(kind=rk), dimension(state_dim) :: se
  real(kind=rk), dimension(state_dim) :: truth
  real(kind=rk) :: mse,rmse
  integer :: ios
  
  character(256) :: filename
  
  call get_truth(truth)
  
  se = (mean-truth)**2
  mse = sum(se)/real(state_dim)
  rmse = sqrt(mse)
  
  if(pf%timestep .eq. 0) then

     select case(comm_version)
     case(1)
        write(filename,'(A)') trim(pf%rmse_filename)
     case(2)
        write(filename,'(A)') trim(pf%rmse_filename)
     case(3)
        write(filename,'(A,A,i0)') trim(pf%rmse_filename),'.',pf_ens_rank
     case default    
        write(emp_e,*) 'ERROR in output_spatial_rmse: comm_version ',comm_version,' &
             &is not supported'
        write(emp_e,*) 'Stopping'
        stop '-23'
     end select

     
     open(unit_spatial_rmse,file=trim(filename),iostat=ios,action='write',status='replace')
     if(ios .ne. 0)  then
        write(emp_e,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open &
             &file ',filename
        write(emp_e,*) 'Very strange that I couldnt open it. Im goin&
             &g to stop now.'
        stop
     end if
  end if
  
  write(unit_spatial_rmse,'(es24.16)') rmse
  
  if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) close(unit_spatial_rmse)
end subroutine output_spatial_rmse
