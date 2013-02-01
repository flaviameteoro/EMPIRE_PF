program couple_pf
  use comms

  implicit none

  integer i
  integer mpi_err

  call initialise_mpi

  call allocate_data
  
  open(20,file='out_pf',action='write')
  write(20,*) 'Starting couple_pf. MyRank = ',myRank,'nProc = ',nProc

  do i=1,ntimesteps

     call gather_data

!    call Statistics(i)
!    call Analyse(i)
!    call Statistics(i) 

     call nudge_data

     call scatter_data
     write(20,*) 'finished timestep ',i,' of ',ntimesteps
  enddo

  
  call deallocate_data
  write(20,*) 'Finishing couple_pf. MyRank = ',myRank,'nProc = ',nProc
  close(20)
  call MPI_Finalize(mpi_err)
  
end program couple_pf


