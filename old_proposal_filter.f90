!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine proposal_filter based on subroutine IntegrateModel
!given from the Particle filter code of Mel and Peter Jan
!PAB 04-02-2013

subroutine proposal_filter
  use pf_control
  use Sizes
  use comms

  IMPLICIT NONE
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)

  !Mel-14|11|11-added to allow fullQ eror in weights
  real(kind=rk) :: pWeight!, qWeight

  !Mel-20|12|11-added to allow sqrtQ correlated random error
  real(kind=rk), dimension(state_dim) :: normaln     !vector to store uncorrelated random error
  real(kind=rk), dimension(state_dim) :: betan       !vector to store sqrtQ correlated random error

  real(kind=rk), dimension(obs_dim) :: y             !y, the observations
  real(kind=rk), dimension(obs_dim) :: Hpsi          !H(psi^(n-1))
  real(kind=rk), dimension(obs_dim,pf%count) :: y_Hpsin1      !y-H(psi^(n-1))
  real(kind=rk), dimension(state_dim,pf%count) :: fpsi        !f(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: kgain       !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  real(kind=rk), dimension(state_dim) :: Qkgain
!  real(kind=rk), dimension(state_dim) :: prop_diff   !psi^n -f(psi^(n-1))
  real(kind=rk) :: dnrm2
  integer :: particle,k,tag,mpi_err
  INTEGER, DIMENSION(:), ALLOCATABLE  :: requests
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: mpi_statuses
  logical :: mpi_flag
  logical, dimension(:), ALLOCATABLE :: received



  if(.not. pf%gen_data) call get_observation_data(y)

  ALLOCATE( requests(pf%count),received(pf%count) )
  ALLOCATE( mpi_statuses(MPI_STATUS_SIZE,pf%count) )
  
  do k =1,pf%count

     if(.not. pf%gen_data) then
        call H(pf%psi(:,k),Hpsi)
        
        y_Hpsin1(:,k) = y - Hpsi
        
     else
        y_Hpsin1(:,k) = 0.0_rk
     end if

!     call send_to_model(pf%psi(:,k),pf%particles(k))

     particle = pf%particles(k)
     tag = 1
     call MPI_ISEND(pf%psi(:,k),state_dim,MPI_DOUBLE_PRECISION&
          &,particle-1,tag,CPL_MPI_COMM,requests(k),mpi_err)


  enddo
  call MPI_WAITALL(pf%count,requests,mpi_statuses,mpi_err)
  !all particles sent to the model now.

  DO k = 1,pf%count
     particle = pf%particles(k)
     tag = 1
     CALL MPI_IRECV(fpsi(:,k), state_dim, MPI_DOUBLE_PRECISION, &
          particle-1, tag, CPL_MPI_COMM,&
          requests(k), mpi_err)
  END DO

  k = 0
  received = .false.
  do
     k = mod(k,pf%count)+1
     if(.not. received(k)) then
        particle = pf%particles(k)
        call MPI_TEST(requests(k),mpi_flag,mpi_statuses(:,k),mpi_err)
        if(mpi_flag) then
           received(k) = .true.
           !DO SOMETHING WITH THE RECEIVED MODEL STATE
           call Bprime(y_Hpsin1(:,k),kgain,Qkgain)

!           call Q(kgain,Qkgain)
           
           call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,normaln)
           call Qhalf(normaln,betan)

!           print*,particle,'min Q fpsi = ',minval(fpsi(406465:539616,k)),' before'
!           print*,particle,'min Q Qkgn = ',minval(Qkgain(406465:539616)),' before'
!           print*,particle,'min Q beta = ',minval(betan(406465:539616)),' before'
!           pf%psi(:,k) = fpsi(:,k) + Qkgain + betan
           print*,'|fpsi-psi|_2 = ',dnrm2(state_dim,(fpsi(:,k)-pf%psi(:,k)),1)
           call update_state(pf%psi(:,k),fpsi(:,k),Qkgain,betan)
!           print*,particle,'min Q = ',minval(pf%psi(406465:539616,k))&
!                &,maxval(pf%psi(406465:539616,k))
!           print*,particle,'min S = ',minval(pf%psi(997128:1454638,k)),maxval(pf%psi(997128:1454638,k))
!           prop_diff = kgain! + betan

           !now we assume that Q and hat{Q} are the same, so
           ! contribution of betan in pWeight
           ! and qWeight will cancel each other out.
           
           !now we can replace the following solve using Q to a solve
           !using R:
           !call innerQ_1(prop_diff,pWeight)
           !call innerR_1(y_Hpsin1(:,k),pWeight)

!           print*,'particle ',particle,'sum(QKgain*kgain) = ',sum(Qkgain&
!                &*kgain),'sum(betan*kgain) = ',sum(betan*kgain)
           pweight = sum(Qkgain*kgain)+2.0D0*sum(betan*kgain)!&
!                &+sum(normaln*normaln)

!           call innerQ_1(betan,qWeight)

           pf%weight(particle) = pf%weight(particle) + 0.5*pWeight! - 0.5*qWeight
           
        else
!           print*,particle,' was not ready'
        end if
     end if
     if(all(received)) exit
  end do
!  print*,'all models received by proposal filter'
  
  

  
!  print*,'pf%weight:',pf%weight


!  pf%weight = exp(-pf%weight)
!  pf%weight = pf%weight/sum(pf%weight)
!  pf%weight = -log(pf%weight)
  !############ only normalise before the equal weights step to
  ! reduce the communications

!  print*,'pf%weight:',pf%weight
!  if(pf%gen_data) call save_truth(pf%psi(:,1))

end subroutine proposal_filter
