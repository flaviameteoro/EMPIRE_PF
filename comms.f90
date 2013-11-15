module comms
  use hadcm3_config
  use hadcm3_data

  integer :: CPL_MPI_COMM,mype_id,myRank,nProc
  integer :: pf_mpi_comm,pfrank
  integer*8 :: npfs
  integer, allocatable, dimension(:) :: gblcount,gbldisp
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
    
    integer :: mpi_err!dummy_colour,mpi_err
    integer :: couple_colour !DUMMY_MPI_COMMUNICATOR,couple_colour
    !integer :: couple_mype_id,couple_root
    integer :: rtmp!,ctmp

 !   integer :: tag!,state_dim!,iter
 !   integer :: num_iters
    integer :: particle,mype_id
    integer :: myrank !nproc,myrank
!    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: nens,i
    integer*8 :: da
    integer :: count,pf_colour,pf_id!,pf_mpi_comm

    
    pf_colour = 10000
    couple_colour=9999
    CALL MPI_INIT (mpi_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype_id,mpi_err)
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,pf_colour,pf_id,pf_mpi_comm&
         &,mpi_err)
    
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,couple_colour,mype_id&
         &,CPL_MPI_COMM,mpi_err)
    CALL MPI_COMM_RANK (CPL_MPI_COMM, myRank, mpi_err)
    CALL MPI_COMM_SIZE (CPL_MPI_COMM, nens, mpi_err)
    
    da = 1
    CALL MPI_ALLREDUCE(da,npfs,1,mpi_integer8,mpi_sum,cpl_mpi_comm&
         &,mpi_err)
    nens = nens-npfs
    
    pfrank = myrank-nens
    
    !lets find the particles:
    count = 0
    do particle = 1,nens
       if( real(particle-1) .ge. real(nens*(pfrank))/real(npfs) .and.&
            & real(particle-1) .lt. real(nens*(pfrank+1))/real(npfs)) then
          count = count + 1
       end if
    end do
    
    allocate(pf%particles(count))
    rtmp = 0
    do particle = 1,nens
       if(real(particle-1) .ge. real(nens*(pfrank))/real(npfs) .and.&
            & real(particle-1) .lt. real(nens*(pfrank+1))/real(npfs))&
            & then
          rtmp = rtmp + 1
          pf%particles(rtmp) = particle
       end if
    end do
    
    allocate(gblcount(npfs))
    allocate(gbldisp(npfs))
    print*,'woohoo allgather'
    print*,count
    print*,gblcount
    print*,pf_mpi_comm
    call mpi_allgather(count,1,mpi_integer,gblcount,1,mpi_integer&
         &,pf_mpi_comm,mpi_err)
    print*,'allgather did not break'
    gbldisp = 0
    if(npfs .gt. 1) then
       do i = 2,npfs
          gbldisp(i) = gbldisp(i-1) + gblcount(i-1)
       end do
    end if
    pf%count = count

    pf%nens = nens
    PRINT*,'PF_rank = ',pfrank,' and I own particles ',pf%particles

    
  end subroutine initialise_mpi

  subroutine receive_from_model(mdata,particle,request)

    use hadcm3_config
    use sizes

    implicit none
    include 'mpif.h'
    
    real(kind=kind(1.0D+0)), INTENT(OUT)     ::mdata(state_dim)
    integer, INTENT(IN) :: particle
    integer, intent(inout) :: request
    !integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: mpi_err
    

    call MPI_IRECV(mdata, state_dim, MPI_DOUBLE_PRECISION, &
         particle-1, 1, CPL_MPI_COMM,&
         request, mpi_err)
    
    
    
    
  end subroutine receive_from_model
  
  
  subroutine send_to_model(mdata,particle,request)
    use hadcm3_config
    use sizes
    implicit none
    include 'mpif.h'
    
    integer, INTENT(IN) :: particle
    integer, intent(inout) :: request
    real(kind=kind(1.0D+0)), INTENT(IN)     ::mdata(state_dim)

    integer :: mpi_err
    
    call MPI_ISend(mdata, state_dim , MPI_DOUBLE_PRECISION, &
         particle-1, 1, CPL_MPI_COMM,request,mpi_err)
    
    
  end subroutine send_to_model

  subroutine wait_mpi(requests)
    use pf_control
    implicit none
    include 'mpif.h'
    
    integer, intent(inout), dimension(pf%count) :: requests
    integer :: mpi_statuses(MPI_STATUS_SIZE,pf%count)
    integer :: mpi_err

    call MPI_WAITALL(pf%count,requests,mpi_statuses,mpi_err)
  end subroutine wait_mpi
    
end module comms
