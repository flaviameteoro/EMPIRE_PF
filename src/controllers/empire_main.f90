!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-07-14 14:05:07 pbrowne>
!!!
!!!    The main program to run EMPIRE
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

!program to run the particle filter on the model HadCM3.
!this shall hopefully have minimal changes specific to the model.
!Simon Wilson and Philip Browne 2013
!----------------------------------------------------------------


!> the main program
!!
program empire
  use output_empire
  use comms
  use pf_control
  use sizes
  implicit none
  include 'mpif.h'
  integer :: i,j,k
  integer :: mpi_err
  real(kind=kind(1.0d0)) :: start_t,end_t


  write(6,'(A)') 'PF: Starting PF code'
  call flush(6)


  !> set up EMPIRE coupling
  call initialise_mpi

  !> define output files
  call open_emp_o(pfrank)

  write(emp_o,*) 'PF: setting controls'
  !> read in controlling data
  call set_pf_controls


  write(emp_o,*) 'PF: configuring model'
  !> call user specific routine for initialisation
  call configure_model


  write(emp_o,*) 'allocating pf'
  !> allocate space for the filter
  call allocate_pf


  !> ensure random seed is set across mpi processes
  call random_seed_mpi(pfrank)
  write(6,*) 'PF: starting to receive from model'



  start_t = mpi_wtime() 
  pf%time=mpi_wtime()

  if(.not. pf%gen_Q) then

  call recv_all_models(state_dim,pf%count,pf%psi)

  do k = 1,pf%count
      call perturb_particle(pf%psi(:,k))
  end do

  write(6,*) 'PF: All models received in pf couple' 
  call flush(6)

  call output_from_pf
  if(pf%gen_data) call save_truth(pf%psi(:,1))
  if(pf%use_traj) call trajectories
  start_t = mpi_wtime()

  
  !start the timestep loop
  do j=1,pf%time_obs
     write(6,*) 'PF: observation counter = ',j
     do i = 1,pf%time_bwn_obs-1
        pf%timestep = pf%timestep + 1
        
        select case(pf%filter)
        case('EW')
           call proposal_filter
        case('EZ')
           call proposal_filter
        case('SI')
           call stochastic_model
        case('SE')
           call stochastic_model
        case('LE')
           call stochastic_model
        case('LD')
           call deterministic_model
        case('DE')
           call deterministic_model
        case default
           write(emp_o,*) 'Error -555: Incorrect pf%filter'
           stop -555
        end select
        call flush(6)
        if(pf%use_traj) call trajectories
        call output_from_pf
     end do
           
     pf%timestep = pf%timestep + 1
     write(6,*) 'starting the observation timestep'
     call flush(6)

     select case(pf%filter)
     case('EW')
        call equivalent_weights_filter
     case('EZ')
        call equivalent_weights_filter_zhu
     case('SI')
        call sir_filter
     case('SE')
        call stochastic_model
        call diagnostics
     case('LE')
        call stochastic_model
        call letkf_analysis
     case('LD')
        call deterministic_model
        call letkf_analysis
     case('DE')
        call deterministic_model
     case default
        write(emp_o,*) 'Error -556: Incorrect pf%filter'
        stop -556
     end select
     

     write(6,*) 'PF: timestep = ',pf%timestep, 'after observation analysis'
     call flush(6)

     if(pf%gen_data) call save_truth(pf%psi(:,1))
     if(pf%use_traj) call trajectories
     call output_from_pf

     call reconfigure_model

  end do
  call diagnostics
  write(6,*) 'PF: finished the loop - now to tidy up'
  call flush(6)



  !send the final state to the model to allow it to finish cleanly
  call send_all_models(state_dim,pf%count,pf%psi,3)

  else

     call genQ

  end if
  end_t = mpi_wtime()


  
  call deallocate_data

  write(6,*) 'PF: finished deallocate_data - off to mpi_finalize'
  call flush(6)

  call MPI_Finalize(mpi_err)
  write(*,*) 'Program couple_pf terminated successfully.'
  write(*,*) 'Time taken in running the model = ',end_t-start_t
  
  call close_emp_o

end program empire


