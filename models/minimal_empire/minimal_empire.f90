!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-09-09 10:25:55 pbrowne>
!!!
!!!    minimal code to setup and test empire comms
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

!> the main program
!!
program minimal_empire
  use comms
  use pf_control
  use sizes
  implicit none
  include 'mpif.h'
  integer :: mpi_err,total_timesteps,i

  real(kind=kind(1.0d0)), allocatable, dimension(:,:) :: x


  print*,'RUNNING MINIMAL_EMPIRE'
  print*,'EMPIRE COUPLING VERSION ',comm_version


  print*,'Reading state_dim from file state_dim: '
  open(11,file='state_dim',action='read',status='old')
  read(11,*) state_dim
  close(11)


  !> set up EMPIRE coupling
  call initialise_mpi


  print*,'allocating space for state vector (',state_dim,',',cnt,'):'
  allocate(x(state_dim,cnt))
  

  print*,'Reading total_timesteps from file timesteps: '
  open(11,file='timesteps',action='read',status='old')
  read(11,*) total_timesteps
  close(11)
  print*,'Total timesteps read as ',total_timesteps

  print*,'receiving initial states'
  call recv_all_models(state_dim,cnt,x)
  
  do i = 1,total_timesteps
     print*,'starting timestep ',i
     call send_all_models(state_dim,cnt,x,1)
     call recv_all_models(state_dim,cnt,x)
  end do

  print*,'sending final states'
  call send_all_models(state_dim,cnt,x,3)



  call mpi_finalize(mpi_err)
  print*,'MINIMAL_EMPIRE got to the end without breaking!'

  


end program minimal_empire


