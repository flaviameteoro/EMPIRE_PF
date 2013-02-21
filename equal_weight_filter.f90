subroutine equal_weight_filter
  use pf_control
  use sizes
  use random
  use comms
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(pf%ngrand) :: a,b,c,alpha
  real(kind=rk) :: cmax
  integer :: particle 
  real(kind=rk), dimension(obs_dim) :: y               !y, the observations
  real(kind=rk), dimension(obs_dim) :: Hpsi            !H(psi^(n-1))
  real(kind=rk), dimension(obs_dim) :: y_Hpsin1        !y-H(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: fpsi          !f(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: kgain         !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  real(kind=rk), dimension(state_dim) :: prop_diff     !psi^n -f(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: betan         !the mixture random variable
  real(kind=rk), dimension(state_dim) :: statev        !temporary state space vector 
  real(kind=rk), dimension(obs_dim) :: obsv,obsvv      !temporary  obs  space vector

  real(kind=rk) :: w,uu

  if(.not. pf%gen_data) then
     call get_observation_data(y)
     

     !$omp parallel do
     do particle =1,pf%ngrand


        call H(pf%psi(:,particle),Hpsi)
     
        y_Hpsin1 = y - Hpsi
        
        call send_to_model(pf%psi(:,particle),particle)

     enddo
     !$omp end parallel do 
     
     !$omp parallel do 
     do particle =1,pf%ngrand
        call recieve_from_model(fpsi,particle)
        
        !c(particle) = pf%weight(particle) + 0.5*(y-Hf(x_i^n-1))^T (HQH^T + R)^(-1) (y-Hf(x_i^n-1))
        call innerHQHt_plus_R_1(y_Hpsin1,w)
        c(particle) = pf%weight(particle) + 0.5*w
        
     end do
     !$omp end parallel do
     
     !here we can pick somehow the 80% level etc...
     cmax = maxval(c)
  end if


  !$omp parallel do
  do particle = 1,pf%ngrand
     if(.not. pf%gen_data) then
        !a(particle) = 0.5* (y-Hf(x_i^n-1))^T R^(-1) (HQH^T) (HQH^T + R)^(-1) (y-Hf(x_i^n-1))
        call K(y_Hpsin1,kgain)
        call H(kgain,obsv)
        call solve_r(obsv,obsvv)
        w = sum(obsv*y_Hpsin1)
        a(particle) = 0.5*w
        
        !b(particle) = 0.5* (y-Hf(x_i^n-1))^T R^(-1) (y-Hf(x_i^n-1)) - cmax - pf%weight(i)
        call innerR_1(y_Hpsin1,w)
        b(particle) = 0.5*w - cmax - pf%weight(particle)



        !note the plus sign in the below equation. See Ades & van Leeuwen 2012.
        alpha(particle) = 1.0 + sqrt(1.0 - b(particle)/a(particle) + 1.0D-10) 
     else
        kgain = 0.0_rk
        alpha(particle) = 0.0_rk
        y_Hpsin1 = 0.0_rk
     end if

     !generate beta from a mixture density
     call MixtureRandomNumbers1D(0.0D0,pf%nfac,pf%ufac,pf%efac,state_dim,statev)
     call Qhalf(statev,betan)

     !now do the following
     !x^n = f(x^(n-1)) + alpha(i) K (y-Hf(x_i^n-1)) + beta
     pf%psi(:,particle) = fpsi + a(particle)*kgain + betan

     prop_diff = a(particle)*kgain + betan

     !now calculate the new weights
     !pf%weight(particle) = 0.5*(y-Hf(x_i^n-1))^T R^(-1) (y-Hf(x_i^n-1)) + 0.5*(x_i^n-f(x_i^(n-1)))^T Q^(-1) (x_i^n-f(x_i^(n-1))) + pf^weight(i)
     call innerR_1(y_Hpsin1,w)
     call innerQ_1(prop_diff,uu)
     
     pf%weight(particle) = pf%weight(particle)  + 0.5*(w+uu)


  end do
  !$omp end parallel do

  !normalise the weights:
  pf%weight = exp(-pf%weight)
  pf%weight = pf%weight/sum(pf%weight)
  pf%weight = -log(pf%weight)

  if(pf%gen_data) then
     call H(pf%psi(:,1),y)
     call save_observation_data(y)
  end if
  

end subroutine equal_weight_filter
