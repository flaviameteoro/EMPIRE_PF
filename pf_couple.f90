!program to run the particle filter on the model HadCM3.
!this shall hopefully have minimal changes specific to the model.
!Simon Wilson and Philip Browne 2013
!----------------------------------------------------------------

program couple_pf
  use comms
  use pf_control
  implicit none
  
  integer :: i,j
  integer :: mpi_err
  
  call initialise_mpi
  
!  call allocate_data
 
  call configure_model
 
  call set_pf_controls

  call allocate_pf
  
  do j=1,pf%time_obs

     do i = 1,pf%time_bwn_obs-1
        pf%timestep = pf%timestep + 1
        call proposal_filter
     end do
           
     pf%timestep = pf%timestep + 1
     call equal_weight_filter

  enddo

  
  call deallocate_data
  
  call MPI_Finalize(mpi_err)
  
  write(*,*) 'Program couple_pf terminated successfully.'

end program couple_pf


