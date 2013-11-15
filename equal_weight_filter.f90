subroutine equal_weight_filter
  use pf_control
  use sizes
  use random
  use comms
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(pf%count) :: a,b,alpha,c
  real(kind=rk), dimension(pf%nens) :: csorted
  real(kind=rk) :: cmax
  integer :: particle,i,tag,mpi_err
  real(kind=rk), dimension(obs_dim) :: y     !y, the observations
  real(kind=rk), dimension(obs_dim,pf%count) :: Hfpsi            !H(f(psi^(n-1)))
  real(kind=rk), dimension(obs_dim,pf%count) :: y_Hfpsin1  !y-H(f(psi^(n-1)))
  real(kind=rk), dimension(state_dim,pf%count) :: fpsi          !f(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: psimean
  real(kind=rk), dimension(state_dim,pf%count) :: kgain         !QH^T(HQH^T+R)^(-1)(y-H(f(psi^(n-1))))
  real(kind=rk), dimension(state_dim,pf%count) :: betan         !the mixture random variable
  real(kind=rk), dimension(state_dim,pf%count) :: statev        !temporary state space vector 
  real(kind=rk), dimension(obs_dim,pf%count) :: obsv,obsvv      !temporary  obs  space vector
  real(kind=rk) :: w
  real(kind=rk), dimension(pf%count) :: e                     !e = d_i^t R^(-1) d_i
  real(kind=rk), parameter :: pi = 4.0D0*atan(1.0D0)
  logical, dimension(pf%count) :: uniform
  INTEGER, DIMENSION(MPI_STATUS_SIZE) :: mpi_status
  real(kind=rk), dimension(pf%count) :: weight_temp
  real(kind=rk) :: ddot
!  print*,'in equal weight filter the weights are:'
!  print*,pf%weight
  
  weight_temp = -huge(1.0d0)
  do i = 1,pf%count
     weight_temp(i) = pf%weight(pf%particles(i))
  end do
!  print*,'temporary weight = :'
!  print*,weight_temp

  call mpi_allgatherv(weight_temp,pf%count,mpi_double_precision,pf%weight,gblcount&
       &,gbldisp,mpi_double_precision,pf_mpi_comm,mpi_err)

!  print*,'after allgather, pf%weight = '
!  print*,pf%weight

  pf%weight = exp(-pf%weight+maxval(pf%weight))
  pf%weight = pf%weight/sum(pf%weight)
  pf%weight = -log(pf%weight)
!  print*,'weights should be normalised:'
!  print*,pf%weight
  if(.not. pf%gen_data) then
     call get_observation_data(y)
 
!     do particle =1,pf%count
!        call send_to_model(pf%psi(:,particle),particle)
!     enddo

     DO i = 1,pf%count
        particle = pf%particles(i)
        tag = 1
        CALL MPI_SEND(pf%psi(:,i), state_dim , MPI_DOUBLE_PRECISION, &
             particle-1, tag, CPL_MPI_COMM, mpi_err)
     END DO

     DO i = 1,pf%count
        particle = pf%particles(i)
        tag = 1
        CALL MPI_RECV(fpsi(:,i), state_dim, MPI_DOUBLE_PRECISION, &
             particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
     END DO
     

     call H(fpsi,Hfpsi)


     do i = 1,pf%count
        particle = pf%particles(i)
        y_Hfpsin1(:,i) = y - Hfpsi(:,i)
        
        call innerHQHt_plus_R_1(y_Hfpsin1(:,i),w)
        
        c(i) = pf%weight(particle) + 0.5*w
     end do
     
     !print*,'all models received by particle filter'
     !  end do
!     print*,'allgatherv in eq'
!     print*,c
!     print*,csorted
     !print*,gblcount
     !print*,gbldisp
     !print*,pf_mpi_comm
     !print*,pf%count
     !here we can pick somehow the 80% level etc...
     !print*,'launching mpi_allgatherv pfrank =',pfrank
     call mpi_allgatherv(c,pf%count,mpi_double_precision,csorted,gblcount&
          &,gbldisp,mpi_double_precision,pf_mpi_comm,mpi_err)
!     print*,'after allgatherv',mpi_err
!     print*,csorted
     call quicksort_d(csorted,pf%nens)
     cmax = csorted(nint(pf%keep*pf%nens))
!     print*,'cmax = ',cmax
  else
     
     DO i = 1,pf%count
        particle = pf%particles(i)
        tag = 1
        CALL MPI_SEND(pf%psi(:,i), state_dim , MPI_DOUBLE_PRECISION, &
             particle-1, tag, CPL_MPI_COMM, mpi_err)
        !PRINT*,'Particle filter ',pfrank,'has sent state_vector over mpi at iteratio&
        !    &n',iter,' to ensemble member ',particle
     END DO
     DO i = 1,pf%count
        particle = pf%particles(i)
        tag = 1
        CALL MPI_RECV(fpsi(:,i), state_dim, MPI_DOUBLE_PRECISION, &
             particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
     END DO

  end if

!  call mpi_barrier(pf_mpi_comm,mpi_err)
  psimean = 0.0_rk

  if(.not. pf%gen_data) then
     call K(y_Hfpsin1,kgain)
     call H(kgain,obsv)
     call solve_r(obsv,obsvv)
     
     do i = 1,pf%count
        a(i) = 0.5*ddot(obs_dim,obsvv(:,i),1,y_Hfpsin1(:,i),1)
     end do
     
     call innerR_1(y_Hfpsin1,e)
     
     do i = 1,pf%count
        particle = pf%particles(i)
        b(i) = 0.5*e(i) - cmax + pf%weight(particle)
        !note the plus sign in the below equation. See Ades & van Leeuwen 2012.
        alpha(i) = 1.0 + sqrt(1.0 - b(i)/a(i) + 1.0D-6)
        
!        print*,'i a  b alpha'
!        print*,i,a(i),b(i),alpha(i)
     end do
     
  else !if(.not. pf%gen_data) 
     kgain = 0.0_rk
     alpha = 0.0_rk
     y_Hfpsin1 = 0.0_rk
  end if !if(.not. pf%gen_data)

  call MixtureRandomNumbers2D(0.0D0,pf%nfac,pf%ufac,pf%efac,state_dim,pf%count,statev,uniform)
  call Qhalf(pf%count,statev,betan)

  do i = 1,pf%count
     if(c(i) .le. cmax) then
        particle = pf%particles(i)
        if(uniform(i)) then
           pf%weight(particle) = pf%weight(particle) +&
                (alpha(i)**2.0_rk - 2.0_rk*alpha(i))*a(i) + & 
                0.5_rk*e(i)
        else
           pf%weight(particle) = pf%weight(particle) +&
                (alpha(i)**2.0_rk - 2.0_rk*alpha(i))*a(i) + &
                0.5_rk*e(i) &
                + 2**(-real(state_dim,rk)/2.0_rk)*pi**(real(state_dim,rk)&
                &/2.0_rk)*pf%nfac*pf%ufac**(-real(state_dim,rk))*((1.0_rk&
                &-pf%efac)/pf%efac)*exp(0.5_rk*(sum(betan(:,i)*betan(:,i))))        
        end if !if(uniform)
        
        !now do the following
        !x^n = f(x^(n-1)) + alpha(i) K (y-Hf(x_i^n-1)) + beta
        call update_state(pf%psi(:,i),fpsi(:,i),alpha(i)*kgain,betan)
        psimean = psimean + pf%psi(:,i)
     else
        pf%weight(pf%particles(i)) = huge(1.0D0)
     end if
  end do
  

!=========================================================================

  if(pf%gen_data) then
     write(6,*) 'generating the data'
     call flush(6)
     call H(pf%psi,y)
     !     write(6,*) 'after H'
     !     call flush(6)
     call NormalRandomNumbers1D(0.0D0,1.0D0,obs_dim,obsv)
     call rhalf(obsv,obsvv)
     y = y + obsvv(:,1)
     call save_observation_data(y)
  else 
     if(pf%use_talagrand) call diagnostics
     !print*,'entering resample step'
     print*,'time until resample = ',mpi_wtime()-pf%time
     call flush(6)
     call resample
  end if !if(pf%gen_data)
  
  
end subroutine equal_weight_filter
