module comms
  use hadcm3_config
  use hadcm3_data

  integer :: COUPLE_MPI_COMMUNICATOR,mype_id,myRank,nProc
contains
  
  subroutine allocate_data
    ! Enough room for U and V 
    implicit none
    
!    ! U and V have one fewer point NS
!    allocate( u(a_nxn*(a_nyn-1),a_levels ,nEns))
!    allocate( v(a_nxn*(a_nyn-1),a_levels ,nEns))

!    allocate( u(a_nxn,a_nyn,a_levels ,nEns))
!    allocate( v(a_nxn,a_nyn,a_levels ,nEns))
!    allocate( thetal(a_nxn,a_nyn,a_levels ,nEns))
!    allocate( qt(a_nxn,a_nyn,a_levels ,nEns))
!    allocate( pstar(a_nxn,a_nyn,nEns))
    
!Ocean
!    allocate( b_u(o_nxn,o_nyn,o_levels ,nEns))
!    allocate( b_v(o_nxn,o_nyn,o_levels ,nEns))
!    allocate( t_o(o_nxn,o_nyn,o_levels ,nEns))
!    allocate( sal(o_nxn,o_nyn,o_levels ,nEns))
!    allocate( mld(o_nxn,o_nyn,nEns))


  end subroutine allocate_data
  
  subroutine deallocate_data
    implicit none
!    deallocate(pstar)
!    deallocate(u)
!    deallocate(v)
!    deallocate(thetal)
!    deallocate(qt)
  end subroutine deallocate_data
  
  subroutine initialise_mpi
    use pf_control
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
    
!    nens=nproc-1
    pf%ngrand=nproc-1


    couple_root=ctmp
    write(6,*)'chello',mype_id,myrank,nproc,couple_root

  end subroutine initialise_mpi

  subroutine recieve_from_model(mdata,particle)

    use hadcm3_config
    use sizes

    implicit none
    include 'mpif.h'
    
    real(kind=kind(1.0D+0)), INTENT(OUT)     ::mdata(state_dim)
    integer, INTENT(IN) :: particle

    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: mpi_err
    
    
          call MPI_Recv(mdata, state_dim, MPI_DOUBLE, &
               particle-1, 1, COUPLE_MPI_COMMUNICATOR,&
               mpi_status, mpi_err)
     write(6,*)'phellor1',particle,state_dim,mdata(100)
    
!    do n=1,nens
!       do j=1,a_levels
!          do i=1,a_nxn*(a_nyn-1)
!             u(i,j,n)=psiGrand(i,j,n)
!             v(i,j,n)=psiGrand(i,j+levels,n)
!          enddo
!       enddo
!    enddo
    
    
  end subroutine recieve_from_model
  
  
  subroutine send_to_model(mdata,particle)
    use hadcm3_config
    use sizes
    implicit none
    include 'mpif.h'
    
    integer, INTENT(IN) :: particle
    real(kind=kind(1.0D+0)), INTENT(IN)     ::mdata(state_dim)

    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: mpi_err
    
     call MPI_Send(mdata, state_dim , MPI_DOUBLE, &
         particle-1, 1, COUPLE_MPI_COMMUNICATOR, mpi_err)
     write(6,*)'phellos1',particle,state_dim,mdata(100)
    
    
  end subroutine send_to_model
end module comms
