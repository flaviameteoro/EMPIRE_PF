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
  enddo
  write(6,*) 'PF: All models recieved in pf couple' 
  call flush(6)
  do j=1,pf%time_obs
     write(6,*) 'PF: j = ',j
     do i = 1,pf%time_bwn_obs-1
        pf%timestep = pf%timestep + 1
        call proposal_filter
        write(6,*) 'PF: i = ',i, 'after proposal filter'
        call flush(6)
     end do
           
     pf%timestep = pf%timestep + 1
     call equal_weight_filter
     write(6,*) 'PF: timestep = ',pf%timestep, 'after equal weight filter'
     call flush(6)

  enddo
  write(6,*) 'PF: finished the loop - now to tidy up'
  call flush(6)

  
  call deallocate_data

  write(6,*) 'PF: finished deallocate_data - off to mpi_finalize'
  call flush(6)


!  call MPI_Finalize(mpi_err)
  
  write(*,*) 'Program couple_pf terminated successfully.'

end program couple_pf


