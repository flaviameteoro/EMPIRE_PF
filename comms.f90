module comms
  use hadcm3_config
  use hadcm3_data
  integer, dimension(:), Allocatable :: requests
  integer :: nens
  integer, parameter :: num_fields=(2*levels)
  integer :: COUPLE_MPI_COMMUNICATOR,mype_id,myRank,nProc
contains
  
  subroutine allocate_data
    ! Enough room for U and V 
    use GrandField
    implicit none
    allocate( psiGrand(nxn*nyn,num_fields ,nEns))
    
    ! U and V have one fewer point NS
    allocate( u(nxn*(nyn-1),levels ,nEns))
    allocate( v(nxn*(nyn-1),levels ,nEns))
    allocate( thetal(nxn*nyn,levels ,nEns))
    allocate( qt(nxn*nyn,levels ,nEns))
    allocate( pstar(nxn*nyn,nEns))
    
    allocate( requests(nens))  
  end subroutine allocate_data
  
  subroutine deallocate_data
    use GrandField
    implicit none
    deallocate(psiGrand)
    deallocate(requests)
    deallocate(pstar)
    deallocate(u)
    deallocate(v)
    deallocate(thetal)
    deallocate(qt)
  end subroutine deallocate_data
  
  subroutine initialise_mpi
    implicit none
    include 'mpif.h'
    
    integer :: dummy_colour,mpi_err
    integer :: DUMMY_MPI_COMMUNICATOR,couple_colour
    integer :: couple_mype_id,couple_root
    integer*8 rtmp,ctmp
    
    couple_colour=9999
    dummy_colour=10000
    
    
    call MPI_Init (mpi_err)
    CALL MPi_COMM_RANK(MPI_COMM_WORLD,mype_id,mpi_err)
    CALL MPi_COMM_SPLIT(MPI_COMM_WORLD,dummy_colour, &
         mype_id,DUMMY_MPI_COMMUNICATOR,mpi_err)
    CALL MPi_COMM_SPLIT(MPI_COMM_WORLD,couple_colour, &
         mype_id,COUPLE_MPI_COMMUNICATOR,mpi_err)
    
    call MPI_Comm_Rank (COUPLE_MPI_COMMUNICATOR, myRank, mpi_err)
    call MPI_Comm_Size (COUPLE_MPI_COMMUNICATOR, nProc, mpi_err)
    
    rtmp=myrank
    call MPi_allreduce(rtmp,ctmp,1,MPi_INTeger8,MPi_MAX,&
         COUPLE_MPI_COMMUNICATOR, mpi_err)
    
    nens=nproc-1

    couple_root=ctmp
    write(6,*)'chello',mype_id,myrank,nproc,couple_root

  end subroutine initialise_mpi

  subroutine gather_data

    use hadcm3_config
    use GrandField

    implicit none
    include 'mpif.h'
    
    
    integer :: i,j,k,n
    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: mpi_err
    
    
    !Gather
    n=1
    do j=0,nens
       if (j.ne.myrank) then
          write(6,*)'ohello1',j,n,nens,nxn*nyn*levels
          call MPI_IRecv(psiGrand(:,:,n), nxn*nyn*num_fields, MPI_DOUBLE, &
               j, MPI_ANY_TAG, COUPLE_MPI_COMMUNICATOR,&
               requests(n), mpi_err)
          write(6,*)'ohello2',i,j,n,psiGrand(100,1,n), &
               psiGrand(100,1+levels,n)
          call flush(6)
          n=n+1
       endif
    enddo
    call mpi_waitall(nens,requests,mpi_status, mpi_err)
    
    do n=1,nens
       do j=1,levels
          do i=1,nxn*(nyn-1)
             u(i,j,n)=psiGrand(i,j,n)
             v(i,j,n)=psiGrand(i,j+levels,n)
          enddo
       enddo
    enddo
    
    
  end subroutine gather_data
  
  
  subroutine scatter_data
    use hadcm3_config
    use GrandField
    implicit none
    include 'mpif.h'
    
    
    integer :: i,j,k,n
    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: mpi_err
    
    
    !Scatter
    n=1
    do j=0,nens
       if (j.ne.myrank) then
          call MPI_ISend(psiGrand(:,:,n), nxn*nyn*num_fields, MPI_DOUBLE, &
               j, 1, COUPLE_MPI_COMMUNICATOR,&
               requests(n), mpi_err)
          n=n+1
       endif
    enddo
    call mpi_waitall(nens,requests,mpi_status, mpi_err)
    
    
  end subroutine scatter_data
end module comms
