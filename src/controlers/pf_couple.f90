!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-11-26 15:21:46 pbrowne>
!!!
!!!    {one line to give the program's name and a brief idea of what it does.}
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
!!! Time-stamp: <2014-08-13 15:10:45 pbrowne>

!program to run the particle filter on the model HadCM3.
!this shall hopefully have minimal changes specific to the model.
!Simon Wilson and Philip Browne 2013
!----------------------------------------------------------------


!> the main program
!!
program empire
  use comms
  use pf_control
  use sizes
  implicit none
  include 'mpif.h'
  integer :: i,j,k
  integer :: mpi_err,particle,tag
  INTEGER, DIMENSION(:), ALLOCATABLE  :: requests
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: mpi_statuses
  logical :: mpi_flag
  logical, dimension(:), ALLOCATABLE :: received
  real(kind=kind(1.0d0)) :: start_t,end_t


  write(6,'(A)') 'PF: Starting PF code'
  call flush(6)
  !> set up EMPIRE coupling
  call initialise_mpi
  print*,'PF: setting controls'
  !> read in controlling data
  call set_pf_controls
  print*,'PF: configuring model'
  !> call user specific routine for initialisation
  call configure_model
  print*,'allocating pf'
  !> allocate space for the filter
  call allocate_pf

  !> ensure random seed is set across mpi processes
  call random_seed_mpi(pfrank)
  write(6,*) 'PF: starting to receive from model'


  tag = 1        
  allocate(requests(pf%count))
  allocate(mpi_statuses(mpi_status_size,pf%count))
  allocate(received(pf%count))

  !HADCM3 MODEL SPECIFIC...
  !let us spin the model up for one day, or 72 timesteps
  start_t = mpi_wtime() 
  pf%time=mpi_wtime()

  if(.not. pf%gen_Q) then

  DO k = 1,pf%count
     particle = pf%particles(k)
     print*,'receiving  from ',particle
     CALL MPI_IRECV(pf%psi(:,k), state_dim, MPI_DOUBLE_PRECISION, &
          particle-1, tag, CPL_MPI_COMM,requests(k), mpi_err)
  end DO
  print*,'receives launched'
  k = 0
  received = .false.
  do
     k = mod(k,pf%count)+1
!     print*,k
     if(.not. received(k)) then
        particle = pf%particles(k)
!        print*,particle ,'not received so testing it'
        call MPI_TEST(requests(k), mpi_flag, mpi_statuses(:,k), mpi_err)
        
        if(mpi_flag) then
!           PRINT*,'Particle filter ',pfrank,'has received initial state_v&
!                &ector over mpi from ensemble member ',particle
           received(k) = .true.
!           if(.not. pf%gen_data) call perturb_particle(pf%psi(:,k))
!           print*,pf%psi(:,k)
           call perturb_particle(pf%psi(:,k))
        end if
     end if
     if(all(received)) exit
  end do
  write(6,*) 'PF: All models received in pf couple' 
  call flush(6)

!  if(pf%gen_data) call save_truth(pf%psi(:,1))
  call output_from_pf
  if(pf%gen_data) call save_truth(pf%psi(:,1))
  if(pf%use_traj) call trajectories
  start_t = mpi_wtime()

  do j=1,pf%time_obs
     write(6,*) 'PF: observation counter = ',j
     do i = 1,pf%time_bwn_obs-1
        pf%timestep = pf%timestep + 1
        if(pf%type .eq. 'EW') then
           call proposal_filter
        elseif(pf%type .eq. 'SI') then
           call stochastic_model
        elseif(pf%type .eq. 'SE') then
           call stochastic_model
        elseif(pf%type .eq. 'ET') then
           !this may not need to be stochastic...
           call deterministic_model
        else
           print*,'Error -555: Incorrect pf%type'
        end if
!        write(6,*) 'PF: timestep = ',pf%timestep, 'after proposal filter'
        call flush(6)
        if(pf%use_traj) call trajectories
        call output_from_pf
     end do
           
     pf%timestep = pf%timestep + 1
     write(6,*) 'starting the equivalent weights filter step'
     call flush(6)


     if(pf%type .eq. 'EW') then
           call equivalent_weights_filter
        elseif(pf%type .eq. 'SI') then
           call sir_filter
        elseif(pf%type .eq. 'SE') then
           call stochastic_model
           call diagnostics
        elseif(pf%type .eq. 'ET') then
           print*,'starting the letkf'
           call deterministic_model
           call letkf_analysis
           print*,'finished the letkf'
        else
           print*,'Error -556: Incorrect pf%type'
        end if
     write(6,*) 'PF: timestep = ',pf%timestep, 'after equal weight filter'
     call flush(6)

     if(pf%gen_data) call save_truth(pf%psi(:,1))
     if(pf%use_traj) call trajectories
     call output_from_pf

     call reconfigure_model

  end do
  call diagnostics
  write(6,*) 'PF: finished the loop - now to tidy up'
  call flush(6)


  tag = 1        
  DO k = 1,pf%count
     particle = pf%particles(k)
     CALL MPI_ISEND(pf%psi(:,k), state_dim , MPI_DOUBLE_PRECISION, &
          particle-1, tag, CPL_MPI_COMM, requests(k), mpi_err)
     PRINT*,'Particle filter ',pfrank,'has sent final state_vector over mpi &
          &to ensemble member ',particle
  END DO
  CALL MPI_WAITALL(pf%count,requests,mpi_statuses, mpi_err)


  else

     call genQ

  end if
  end_t = mpi_wtime()


  
  call deallocate_data

  write(6,*) 'PF: finished deallocate_data - off to mpi_finalize'
  call flush(6)


  call MPI_Finalize(mpi_err)
  deallocate(requests)
  deallocate(mpi_statuses)
  deallocate(received)
  write(*,*) 'Program couple_pf terminated successfully.'
  write(*,*) 'Time taken in running the model = ',end_t-start_t


end program empire


