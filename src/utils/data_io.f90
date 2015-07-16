!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-07-16 15:50:46 pbrowne>
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

  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  integer, intent(in) :: t
  real(kind=rk), dimension(obs_dim), intent(out) :: y
  integer :: ios
  character(14) :: filename

  write(filename,'(A,i7.7)') 'obs_ts_',t

  open(67,file=filename,iostat=ios,action='read',status='old',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Check it exists. I need a lie down.'
     stop
  end if
  read(67) y
  close(67)
end subroutine default_get_observation_data

!> Subroutine to save observation to a file              
!! \n              
!! Uses pf%timestep to determine which observation to save   
!> @param[in] y The observation
subroutine save_observation_data(y)

  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  integer :: ios
  character(14) :: filename

  write(filename,'(A,i6.6)') 'obs_ts_',pf%timestep

  open(67,file=filename,iostat=ios,action='write',status='replace',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
     stop
  end if
  write(67) y
  close(67)

end subroutine save_observation_data

!> Subroutine to save truth to a file              
!! \n                 
!> @param[in] x The state vector
subroutine save_truth(x)
  use timestep_data
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer :: ios
  if(pf%timestep .eq. 0) then
     print*,'opening pf_truth'
     open(62,file='pf_truth',iostat=ios,action='write',status='replace')
     if(ios .ne. 0)  then
        write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_truth'
        write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
        stop
     end if
  end if
  write(62,*) x
  call flush(62)
  !if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
  if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) then
     close(62)
     print*,'closing pf_truth'
  end if
end subroutine save_truth

!>subroutine to ouput data from the filter
subroutine output_from_pf
  use timestep_data
  use matrix_pf
  use pf_control
  use sizes
  use comms
  implicit none
  real(kind=kind(1.0D0)), dimension(state_dim) :: mean
  integer :: ios,particle
  character(9) :: filename
  if(pf%timestep .eq. 0) then
     write(filename,'(A,i2.2)') 'pf_out_',pfrank
     open(68,file=filename,iostat=ios,action='write',status='replace')
     if(ios .ne. 0)  then
        write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_out'
        write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
        stop
     end if

     call read_matrix_pf_information
  end if
!  print*,'output: ',pf%timestep,pf%weight
!  write(68,*) pf%timestep,pf%particles,pf%weight(:)
  write(68,'(i6.6,A)',advance='no') pf%timestep,' '
  do ios = 1,pf%count-1
     write(68,'(i6.6,A,e21.15,A)',advance='no') pf%particles(ios),' ',pf&
          &%weight(pf%particles(ios)),' '
  end do
  write(68,'(i6.6,A,e21.15)',advance='yes') pf%particles(pf%count),' ',pf&
          &%weight(pf%particles(pf%count))
  call flush(68)
!  if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) close(68)
  if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) close(68)  

  if(pf%use_mean) then
     if(pf%timestep .eq. 0) then
        open(61,file='pf_mean',iostat=ios,action='write',status='replace')
        if(ios .ne. 0)  then
           write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_mean'
           write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
           stop
        end if
     end if
     mean = 0.0D0
     do particle = 1,pf%nens
        mean(:) = mean(:) + pf%psi(:,particle)*exp(-pf%weight(particle))
     end do
     write(61,*) mean(:)
     call flush(61)
     !if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) close(61)
     if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) close(61)
  end if

end subroutine output_from_pf

!> subroutine to save the state vector to a named file
!! as an unformatted fortran file
subroutine save_state(state,filename)
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(state_dim), intent(in) :: state !< the
  !!state vector
  character(14), intent(in) :: filename !< the name of the file to
  !!save the state vector in
  integer :: ios

  open(16,file=filename,iostat=ios,action='write',status='replace'&
       &,form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file '&
          &,filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop&
          & now.'
     stop
  end if
  write(16) state
  close(16)
end subroutine save_state


!> subroutine to read the state vector from a named file
!! as an unformatted fortran file
subroutine get_state(state,filename)
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(state_dim), intent(out) :: state !< the
  !!state vector
  character(14), intent(in) :: filename!< the name of the file to 
  !!write the state vector in
  integer :: ios

  open(16,file=filename,iostat=ios,action='read',status='old',form='un&
       &formatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file '&
          &,filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop&
          & now.'
     stop
  end if
  read(16) state
  close(16)
end subroutine get_state
