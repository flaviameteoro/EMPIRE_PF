subroutine equal_weight_filter
  use pf_control
  use sizes
  use random
  use comms
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(pf%ngrand) :: a,b,c,alpha,csorted
  real(kind=rk) :: cmax
  integer :: particle 
  real(kind=rk), dimension(obs_dim) :: y     !y, the observations
  real(kind=rk), dimension(obs_dim) :: Hfpsi            !H(f(psi^(n-1)))
  real(kind=rk), dimension(obs_dim,pf%ngrand) :: y_Hfpsin1  !y-H(f(psi^(n-1)))
  real(kind=rk), dimension(state_dim,pf%ngrand) :: fpsi          !f(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: psimean
  real(kind=rk), dimension(state_dim) :: kgain         !QH^T(HQH^T+R)^(-1)(y-H(f(psi^(n-1))))
  real(kind=rk), dimension(state_dim) :: betan         !the mixture random variable
  real(kind=rk), dimension(state_dim) :: statev        !temporary state space vector 
  real(kind=rk), dimension(obs_dim) :: obsv,obsvv      !temporary  obs  space vector
  real(kind=rk) :: w,e                                 !e = d_i^t R^(-1) d_i
  real(kind=rk), parameter :: pi = 4.0D0*atan(1.0D0)
  logical :: uniform
  
  if(.not. pf%gen_data) then
     call get_observation_data(y)
 
     do particle =1,pf%ngrand
        call send_to_model(pf%psi(:,particle),particle)
     enddo
      
     do particle =1,pf%ngrand
        call recieve_from_model(fpsi(:,particle),particle)
        
        !c(particle) = pf%weight(particle) + 0.5*(y-Hf(x_i^n-1))^T (HQH^T + R)^(-1) (y-Hf(x_i^n-1))
        call H(fpsi(:,particle),Hfpsi)
        y_Hfpsin1(:,particle) = y - Hfpsi
        
        call innerHQHt_plus_R_1(y_Hfpsin1(:,particle),w)
        
        c(particle) = pf%weight(particle) + 0.5*w
        
     end do
     
     !here we can pick somehow the 80% level etc...
     csorted = c
     call kb05ad(csorted,pf%ngrand)
     cmax = csorted(nint(pf%keep*pf%ngrand))

  else
      
     do particle =1,pf%ngrand
        call send_to_model(pf%psi(:,particle),particle)
     enddo
     
     do particle =1,pf%ngrand
        call recieve_from_model(fpsi(:,particle),particle)
     end do
     
  end if


  psimean = 0.0_rk

  
  do particle = 1,pf%ngrand
     if(c(particle) .le. cmax) then
     if(.not. pf%gen_data) then

        !a(particle) = 0.5* (y-Hf(x_i^n-1))^T R^(-1) (HQH^T) (HQH^T + R)^(-1) (y-Hf(x_i^n-1))
        
        call K(y_Hfpsin1(:,particle),kgain)
        
        call H(kgain,obsv)
        
        call solve_r(obsv,obsvv)

        e = sum(obsvv*y_Hfpsin1(:,particle))
        
        a(particle) = 0.5*e
        
           !b(particle) = 0.5* (y-Hf(x_i^n-1))^T R^(-1) (y-Hf(x_i^n-1))- cmax + pf%weight(i)
        
        call innerR_1(y_Hfpsin1(:,particle),e)
        
        b(particle) = 0.5*e - cmax + pf%weight(particle)
        
        !note the plus sign in the below equation. See Ades & van Leeuwen 2012.
        alpha(particle) = 1.0 + sqrt(1.0 - b(particle)/a(particle) + 1.0D-6) 
     else !if(.not. pf%gen_data) 
        kgain = 0.0_rk
        alpha(particle) = 0.0_rk
        y_Hfpsin1(:,particle) = 0.0_rk
     end if !if(.not. pf%gen_data)

     !generate beta from a mixture density
     call MixtureRandomNumbers1D(0.0D0,pf%nfac,pf%ufac,pf%efac,state_dim,statev,uniform)
     call Qhalf(statev,betan)

     if(uniform) then
        pf%weight(particle) = pf%weight(particle) +&
             (alpha(particle)**2.0_rk - 2.0_rk*alpha(particle))*a(particle) + & 
             0.5_rk*e
        
     else
      pf%weight(particle) = pf%weight(particle) +&
             (alpha(particle)**2.0_rk - 2.0_rk*alpha(particle))*a(particle) + &
             0.5_rk*e &
             + 2**(-real(state_dim,rk)/2.0_rk)*pi**(real(state_dim,rk)&
             &/2.0_rk)*pf%nfac*pf%ufac**(-real(state_dim,rk))*((1.0_rk&
             &-pf%efac)/pf%efac)*exp(0.5_rk*(sum(betan*betan)))
              
     end if !if(uniform)

     !now do the following
     !x^n = f(x^(n-1)) + alpha(i) K (y-Hf(x_i^n-1)) + beta

     pf%psi(:,particle) = fpsi(:,particle) + alpha(particle)*kgain + betan

     psimean = psimean + pf%psi(:,particle)

     !now calculate the new weights

!     pf%weight(particle) = pf%weight(particle) + (alpha(particle)&
!          &**2.0_rk - 2.0_rk*alpha(particle))*a(particle) + 0.5_rk*e
  else
     pf%weight(particle) = huge(1.0D0)
  end if !if(c(particle) .le. cmax)
  end do

  if(pf%use_talagrand) call diagnostics
  
  call resample
  
  if(pf%gen_data) then
     call H(pf%psi(:,1),y)
     call NormalRandomNumbers1D(0.0D0,sqrt(2.0D0),obs_dim,obsv)
     y = y + obsv
     call save_observation_data(y)
     call save_truth(pf%psi(:,1))
  end if !if(pf%gen_data)
  
  
end subroutine equal_weight_filter
