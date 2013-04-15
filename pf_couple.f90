!program to run the particle filter on the model HadCM3.
!this shall hopefully have minimal changes specific to the model.
!Simon Wilson and Philip Browne 2013
!----------------------------------------------------------------

program couple_pf
  use comms
  use pf_control
  implicit none
  
  integer :: i,j,particle
  integer :: mpi_err
  real(kind=kind(1.0D0)),dimension(3) :: rdom
  write(6,'(A)') 'PF: Starting PF code'
  call flush(6)
  call initialise_mpi
  
!  call allocate_data
 
  call configure_model
 
  call set_pf_controls

  call allocate_pf
  write(6,*) 'PF: starting to recieve from model'
! 1st call to model to get psi
  do particle =1,pf%ngrand
     call recieve_from_model(pf%psi(:,particle),particle)

     !lets add some random noise to the initial conditions
     call NormalRandomNumbers1D(0.0D0,sqrt(2.0D0),3,rdom)
     pf%psi(:,particle) = pf%psi(:,particle) + rdom

  enddo
  write(6,*) 'PF: All models recieved in pf couple' 
  call flush(6)
  if(pf%gen_data) call save_truth(pf%psi(:,1))
  call output_from_pf

  do j=1,pf%time_obs
     write(6,*) 'PF: j = ',j
     do i = 1,pf%time_bwn_obs-1
        pf%timestep = pf%timestep + 1
        call proposal_filter
        write(6,*) 'PF: i = ',i, 'after proposal filter'
        call flush(6)
        
        call output_from_pf
     end do
           
     pf%timestep = pf%timestep + 1
     call equal_weight_filter
     write(6,*) 'PF: timestep = ',pf%timestep, 'after equal weight filter'
     call flush(6)

     call output_from_pf

  enddo
  write(6,*) 'PF: finished the loop - now to tidy up'
  call flush(6)


  do particle =1,pf%ngrand
     call send_to_model(pf%psi(:,particle),particle)
  enddo

  
  call deallocate_data

  write(6,*) 'PF: finished deallocate_data - off to mpi_finalize'
  call flush(6)


  call MPI_Finalize(mpi_err)
  
  write(*,*) 'Program couple_pf terminated successfully.'

end program couple_pf


