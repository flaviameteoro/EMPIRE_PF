!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-08-16 15:08:45 pbrowne>
!!!
!!!    Collection of subroutines to deal with i/o
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

!> Subroutine to read observation from a file
!! \n
!! Uses pf%timestep to determine which observation to read
!> @param[out] y The observation
!! @param[in] t the current timestep
subroutine default_get_observation_data(y,t)
  use output_empire, only : unit_obs
  use timestep_data
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  integer, intent(in) :: t
  real(kind=rk), dimension(obs_dim), intent(out) :: y
  integer :: ios
  character(14) :: filename

  write(filename,'(A,i7.7)') 'obs_ts_',TSData%next_ob_timestep

  open(unit_obs,file=filename,iostat=ios,action='read',status='old',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Check it exists. I need a lie down.'
     stop
  end if
  read(unit_obs) y
  close(unit_obs)
end subroutine default_get_observation_data

!> Subroutine to save observation to a file              
!! \n              
!! Uses pf%timestep to determine which observation to save   
!> @param[in] y The observation
subroutine save_observation_data(y)
  use output_empire, only : unit_obs
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  integer :: ios
  character(14) :: filename

  write(filename,'(A,i7.7)') 'obs_ts_',pf%timestep

  open(unit_obs,file=filename,iostat=ios,action='write',status='replace',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
     stop
  end if
  write(unit_obs) y
  close(unit_obs)

end subroutine save_observation_data

!> Subroutine to read truth from the file written by @ref save_truth              
!! \n                 
!> @param[out] x The state vector
subroutine get_truth(x)
  use output_empire, only : unit_truth
  use timestep_data
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(state_dim), intent(out) :: x
  integer :: ios
  if(pf%timestep .eq. 0) then
     print*,'opening pf_truth'
     open(unit_truth,file='pf_truth',iostat=ios,action='read')
     if(ios .ne. 0)  then
        write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_truth'
        write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
        stop
     end if
  end if
  read(unit_truth,*) x
  call flush(unit_truth)
  !if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
  if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) then
     close(unit_truth)
     print*,'closing pf_truth'
  end if
end subroutine get_truth

!> Subroutine to save truth to a file              
!! \n                 
!> @param[in] x The state vector
subroutine save_truth(x)
  use output_empire, only : unit_truth,unit_mean
  use timestep_data
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer :: ios
  if(pf%timestep .eq. 0) then
     print*,'opening pf_truth'
     open(unit_truth,file='pf_truth',iostat=ios,action='write',status='replace')
     if(ios .ne. 0)  then
        write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_truth'
        write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
        stop
     end if
  end if
  write(unit_truth,*) x
  call flush(unit_truth)
  !if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
  if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) then
     close(unit_truth)
     print*,'closing pf_truth'
  end if
end subroutine save_truth



!>subroutine to output data from the filter
subroutine output_from_pf
  use output_empire, only : unit_weight,unit_mean
  use timestep_data
  use matrix_pf
  use pf_control
  use sizes
  use comms
  implicit none
  include 'mpif.h'
  real(kind=kind(1.0D0)), dimension(state_dim) :: mean,mtemp
  integer :: ios,particle,mpi_err
  character(20) :: filename


  if(pf%timestep .eq. 0) then

     if(pf%output_weights) then
        write(filename,'(A,i2.2)') 'ensemble_weights_',pfrank
        open(unit_weight,file=filename,iostat=ios,action='write',status='replace')
        if(ios .ne. 0)  then
           write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_out'
           write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
           stop
        end if
     end if !end if output_weights
     call read_matrix_pf_information
  end if

  if(pf%output_weights) then
     write(unit_weight,'(i6.6,A)',advance='no') pf%timestep,' '
     do ios = 1,pf%count-1
        write(unit_weight,'(i6.6,A,e21.15,A)',advance='no') pf%particles(ios),' ',pf&
             &%weight(pf%particles(ios)),' '
     end do
     write(unit_weight,'(i6.6,A,e21.15)',advance='yes') pf%particles(pf%count),' ',pf&
          &%weight(pf%particles(pf%count))
     call flush(unit_weight)

     if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) close(unit_weight)  
  end if !end if output_weights
     
  if(pf%use_mean .and. pf_ens_rank .eq. 0) then
     if(pf%timestep .eq. 0) then
        open(unit_mean,file='pf_mean',iostat=ios,action='write',status='replace')
        if(ios .ne. 0)  then
           write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_mean'
           write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
           stop
        end if
     end if
  end if

  if(pf%use_mean .or. (pf%use_spatial_rmse .and. .not. pf%gen_data)) then
     mtemp = 0.0D0
     do particle = 1,pf%count
        !mean(:) = mean(:) + pf%psi(:,particle)*exp(-pf&
        !     &%weight(particles(particle)))
        mtemp(:) = mtemp(:) + pf%psi(:,particle)/real(pf%nens)
     end do

     call mpi_allreduce(mtemp,mean,state_dim,MPI_DOUBLE_PRECISION,MPI_SUM&
             &,pf_ens_comm,mpi_err)
    
  end if

  if(pf%use_mean .and. pf_ens_rank .eq. 0) then
     write(unit_mean,*) mean(:)
     call flush(unit_mean)
     if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) close(unit_mean)
  end if

  if(pf_ens_rank .eq. 0 .and. pf%use_spatial_rmse .and. .not. pf%gen_data) call output_spatial_rmse(mean)

  if(pf%use_mean .and. pf%use_variance) call output_variance(mean)

  if(comm_version .ne. 3) then !empire_v3 models will be too large!!!
     call matrix_pf_output(npfs-1,pf_mpi_comm,state_dim,cnt,pf%psi&
          &,pf%timestep,TSData%is_analysis)
  end if
  
end subroutine output_from_pf

!> subroutine to save the state vector to a named file
!! as an unformatted fortran file
subroutine save_state(state,filename)
  use output_empire, only : unit_state
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(state_dim), intent(in) :: state !< the
  !!state vector
  character(256), intent(in) :: filename !< the name of the file to
  !!save the state vector in
  integer :: ios

  open(unit_state,file=trim(filename),iostat=ios,action='write',status='replace'&
       &,form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file '&
          &,filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop&
          & now.'
     stop
  end if
  write(unit_state) state
  close(unit_state)
end subroutine save_state


!> subroutine to read the state vector from a named file
!! as an unformatted fortran file
subroutine get_state(state,filename)
  use output_empire, only : unit_state
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(state_dim), intent(out) :: state !< the
  !!state vector
  character(256), intent(in) :: filename!< the name of the file to 
  !!write the state vector in
  integer :: ios

  open(unit_state,file=trim(filename),iostat=ios,action='read',status='old',form='un&
       &formatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file '&
          &,filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop&
          & now.'
     stop
  end if
  read(unit_state) state
  close(unit_state)
end subroutine get_state
