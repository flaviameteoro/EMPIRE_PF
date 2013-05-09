!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine proposal_filter based on subroutine IntegrateModel
!given from the Particle filter code of Mel and Peter Jan
!PAB 04-02-2013

subroutine proposal_filter
  use pf_control
  use Sizes
  use comms

  IMPLICIT NONE

  integer, parameter :: rk = kind(1.0D0)

  !Mel-14|11|11-added to allow fullQ eror in weights
  real(kind=rk) :: pWeight, qWeight

  !Mel-20|12|11-added to allow sqrtQ correlated random error
  real(kind=rk), dimension(state_dim) :: normaln     !vector to store uncorrelated random error
  real(kind=rk), dimension(state_dim) :: betan       !vector to store sqrtQ correlated random error

  real(kind=rk), dimension(obs_dim) :: y             !y, the observations
  real(kind=rk), dimension(obs_dim) :: Hpsi          !H(psi^(n-1))
  real(kind=rk), dimension(obs_dim,pf%ngrand) :: y_Hpsin1      !y-H(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: fpsi        !f(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: kgain       !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  real(kind=rk), dimension(state_dim) :: prop_diff   !psi^n -f(psi^(n-1))
!  real(kind=rk) :: dnrm2
  integer :: particle

  if(.not. pf%gen_data) call get_observation_data(y)

  !put stuff from psi vector into a nice form to work with it
  
  do particle =1,pf%ngrand

     if(.not. pf%gen_data) then
        call H(pf%psi(:,particle),Hpsi)
        
        y_Hpsin1(:,particle) = y - Hpsi
        
     else
        y_Hpsin1(:,particle) = 0.0_rk
     end if

     call send_to_model(pf%psi(:,particle),particle)

  enddo

   
  do particle =1,pf%ngrand
     call receive_from_model(fpsi,particle)

     call B(y_Hpsin1(:,particle),kgain)

!     print*,'kgain',kgain

     call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,normaln)
     call Q(normaln,betan)

     betan = 5.0_rk*betan

!     print*,dnrm2(state_dim,fpsi-pf%psi(:,particle),1),dnrm2(state_dim,kgain,1)&
!          &,dnrm2(state_dim,betan,1)

     pf%psi(:,particle) = fpsi + kgain + betan

     prop_diff = kgain + betan

     call innerQ_1(prop_diff,pWeight)

     call innerQ_1(betan,qWeight)
!     print*,'pW =',pWeight,'qw = ',qWeight
     pf%weight(particle) = pf%weight(particle) + 0.5*pWeight - 0.5*qWeight
    
  end do
  
!  print*,'pf%weight:',pf%weight
  pf%weight = exp(-pf%weight)
  pf%weight = pf%weight/sum(pf%weight)
  pf%weight = -log(pf%weight)
!  print*,'pf%weight:',pf%weight
!  if(pf%gen_data) call save_truth(pf%psi(:,1))

end subroutine proposal_filter
