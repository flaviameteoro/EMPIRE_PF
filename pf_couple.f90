program couple_pf
  use comms

  implicit none

  integer i
  integer mpi_err

  call initialise_mpi

  call allocate_data

  do i=1,ntimesteps

     call gather_data

!    call Statistics(i)
!    call Analyse(i)
!    call Statistics(i) 

     call nudge_data

     call scatter_data

  enddo

  
  call deallocate_data

  call MPI_Finalize(mpi_err)
  
end program couple_pf


