!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-18 10:09:54 pbrowne>
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
program couple_pf
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
  integer :: mpi_status(MPI_STATUS_SIZE)
  real(kind=kind(1.0d0)) :: start_t,end_t

  write(6,'(A)') 'PF: Starting PF code'
  call flush(6)
  call initialise_mpi
  print*,'PF: setting controls'
  call set_pf_controls
  print*,'PF: configuring model'
  call configure_model
  print*,'allocating pf'
  call allocate_pf

  call random_seed_mpi(pfrank)
  write(6,*) 'PF: starting to receive from model'
! 1st call to model to get psi
!  do k =1,pf%count
!     call receive_from_model(pf%psi(:,k),pf%particles(k))

     !lets add some random noise to the initial conditions
!     call perturb_particle(pf%psi(:,k))

!  enddo
  print*,'launching the receives'
 ! print*,pf%count
 ! print*,pf%particles
 ! print*,allocated(requests)
  tag = 1        
  allocate(requests(pf%count))
  allocate(mpi_statuses(mpi_status_size,pf%count))
  allocate(received(pf%count))

!!$  !HADCM3 MODEL SPECIFIC...
!!$  !let us spin the model up for one day, or 72 timesteps
!!$  start_t = mpi_wtime() 
!!$  do j = 1,72
!!$     do k = 1,pf%count
!!$        particle = pf%particles(k)
!!$        call mpi_recv(pf%psi(:,k), state_dim, MPI_DOUBLE_PRECISION, &
!!$             particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
!!$     end do
!!$     do k = 1,pf%count
!!$        particle = pf%particles(k)
!!$        call mpi_send(pf%psi(:,k), state_dim , MPI_DOUBLE_PRECISION, &
!!$             particle-1, tag, CPL_MPI_COMM, mpi_err)
!!$     end do
!!$  end do
!!$  end_t = mpi_wtime()
!!$  write(6,*) 'Time for initial day run = ',end_t-start_t,' seconds.'

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
           call perturb_particle(pf%psi(:,k))
        end if
     end if
     if(all(received)) exit
  end do
  write(6,*) 'PF: All models received in pf couple' 
  call flush(6)


!  if(pf%gen_data) call save_truth(pf%psi(:,1))
  call output_from_pf

  start_t = mpi_wtime()

  do j=1,pf%time_obs
     write(6,*) 'PF: observation counter = ',j
     do i = 1,pf%time_bwn_obs-1
        pf%timestep = pf%timestep + 1
        call proposal_filter
!        write(6,*) 'PF: timestep = ',pf%timestep, 'after proposal filter'
        call flush(6)
        
        call output_from_pf
     end do
           
     pf%timestep = pf%timestep + 1
     write(6,*) 'starting the equal weight filter step'
     call flush(6)
     call equal_weight_filter
     write(6,*) 'PF: timestep = ',pf%timestep, 'after equal weight filter'
     call flush(6)

     call output_from_pf

  enddo
  write(6,*) 'PF: finished the loop - now to tidy up'
  end_t = mpi_wtime()
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


  
  call deallocate_data

  write(6,*) 'PF: finished deallocate_data - off to mpi_finalize'
  call flush(6)


  call MPI_Finalize(mpi_err)
  deallocate(requests)
  deallocate(mpi_statuses)
  deallocate(received)
  write(*,*) 'Program couple_pf terminated successfully.'
  write(*,*) 'Time taken in running the model = ',end_t-start_t
  write(*,*) 'Which works out at ',(end_t-start_t)/real(pf%time_obs),'seconds per day (72 timesteps)'


end program couple_pf


