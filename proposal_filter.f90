!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine proposal_filter based on subroutine IntegrateModel
!given from the Particle filter code of Mel and Peter Jan
!PAB 04-02-2013

subroutine proposal_filter
  !--------------------------------------------------------------------------
  !  Integrate the system from time-step nbegin to time-step nend 
  !--------------------------------------------------------------------------
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
  real(kind=rk), dimension(obs_dim) :: y_Hpsin1      !y-H(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: fpsi        !f(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: kgain       !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  real(kind=rk), dimension(state_dim) :: prop_diff   !psi^n -f(psi^(n-1))
  
  integer :: particle

!!  call get_observation_data(y)




  !put stuff from psi vector into a nice form to work with it
  !$omp parallel do
  do particle =1,pf%ngrand

!!     call H(pf%psi(:,particle),Hpsi)

!!     y_Hpsin1 = y - Hpsi


     !call the model now to make one timestep.....
     !............................................
     !HELP PLEASE SIMON! HOW DO I CALL THE MODEL HERE?!
     !I would like to give the model the current particle
     !which is at psi(:,i) and return the variable fpsi
     !============================================
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !############################################
     !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !||||||||||||||||||||||||||||||||||||||||||||
     call send_to_model(pf%psi(:,particle),particle)
  enddo
  do particle =1,pf%ngrand
     call recieve_from_model(fpsi,particle)

     pf%psi(:,particle) = fpsi

!!     call K(y_Hpsin1,kgain)

     !Mel-20|12|11-changed to allow random error correlated by sqrt Q
!!     call NormalRandomNumbers1D(0.0,1.0,state_dim,normaln)
     call Qhalf(normaln,betan)

!!     pf%psi(:,particle) = fpsi + kgain + betan

!!     prop_diff = kgain + betan

!!     call innerQ_1(prop_diff,pWeight)

!!     call innerQ_1(betan,qWeight)

!!     pf%weight(particle) = pf%weight(particle) + 0.5*pWeight - 0.5*qWeight

  end do
  !$omp end parallel do


end subroutine proposal_filter
