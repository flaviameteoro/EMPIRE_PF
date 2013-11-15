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
  real(kind=rk), dimension(obs_dim) :: Hfpsi            !H(f(psi^(n-1)))
  real(kind=rk), dimension(obs_dim,pf%count) :: y_Hfpsin1  !y-H(f(psi^(n-1)))
  real(kind=rk), dimension(state_dim,pf%count) :: fpsi          !f(psi^(n-1))
  real(kind=rk), dimension(state_dim) :: psimean
  real(kind=rk), dimension(state_dim) :: kgain         !QH^T(HQH^T+R)^(-1)(y-H(f(psi^(n-1))))
  real(kind=rk), dimension(state_dim) :: betan         !the mixture random variable
  real(kind=rk), dimension(state_dim) :: statev        !temporary state space vector 
  real(kind=rk), dimension(obs_dim) :: obsv,obsvv      !temporary  obs  space vector
  real(kind=rk) :: w,e                                 !e = d_i^t R^(-1) d_i
  real(kind=rk), parameter :: pi = 4.0D0*atan(1.0D0)
  logical :: uniform
  INTEGER, DIMENSION(pf%count) :: requests
  INTEGER, DIMENSION(MPI_STATUS_SIZE,pf%count) :: mpi_statuses
  logical :: mpi_flag
  logical, dimension(pf%count) :: received
  real(kind=rk), dimension(pf%count) :: weight_temp
  real(kind=rk) :: dnrm2
  print*,'in equal weight filter the weights are:'
  print*,pf%weight
  
  weight_temp = -huge(1.0d0)
  do i = 1,pf%count
     weight_temp(i) = pf%weight(pf%particles(i))
  end do
  print*,'temporary weight = :'
  print*,weight_temp

  call mpi_allgatherv(weight_temp,pf%count,mpi_double_precision,pf%weight,gblcount&
       &,gbldisp,mpi_double_precision,pf_mpi_comm,mpi_err)

  print*,'after allgather, pf%weight = '
  print*,pf%weight

  pf%weight = exp(-pf%weight+maxval(pf%weight))
  pf%weight = pf%weight/sum(pf%weight)
  pf%weight = -log(pf%weight)
!  print*,'now they should be normalised:'
!  print*,pf%weight
  if(.not. pf%gen_data) then
     call get_observation_data(y)
 
!     do particle =1,pf%count
!        call send_to_model(pf%psi(:,particle),particle)
!     enddo

     DO i = 1,pf%count
        particle = pf%particles(i)
        tag = 1
        CALL MPI_ISEND(pf%psi(:,i), state_dim , MPI_DOUBLE_PRECISION, &
             particle-1, tag, CPL_MPI_COMM, requests(i), mpi_err)
        !PRINT*,'Particle filter ',pfrank,'has sent state_vector over mpi at iteratio&
        !    &n',iter,' to ensemble member ',particle
     END DO
     CALL MPI_WAITALL(pf%count,requests,mpi_statuses, mpi_err)
     

     DO i = 1,pf%count
        particle = pf%particles(i)
        tag = 1
        CALL MPI_IRECV(fpsi(:,i), state_dim, MPI_DOUBLE_PRECISION, &
             particle-1, tag, CPL_MPI_COMM,&
             requests(i), mpi_err)
     END DO
     
     received = .false.
     i = 0
     do
        i = mod(i,pf%count)+1
        if(.not. received(i)) then
           particle = pf%particles(i)
           call MPI_TEST(requests(i), mpi_flag, mpi_statuses(:,i),&
                & mpi_err)
           
           if(mpi_flag) then
              received(i) = .true.
              !     do particle =1,pf%count
              !       call receive_from_model(fpsi(:,particle),particle)
              
              !c(particle) = pf%weight(particle) + 0.5*(y-Hf(x_i^n
              !-1))^T (HQH^T + R)^(-1) (y-Hf(x_i^n-1))
              call H(fpsi(:,i),Hfpsi)
              y_Hfpsin1(:,i) = y - Hfpsi
              
              call innerHQHt_plus_R_1(y_Hfpsin1(:,i),w)
              
              c(i) = pf%weight(particle) + 0.5*w
           else
              !print*,particle,' was not ready'
           end if
        end if
        if(all(received)) exit
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
     print*,csorted
     call quicksort_d(csorted,pf%nens)
     cmax = csorted(nint(pf%keep*pf%nens))
     print*,'cmax = ',cmax
  else
     
     !     do particle =1,pf%count
     !        call send_to_model(pf%psi(:,particle),particle)
     !     enddo
     
     DO i = 1,pf%count
        particle = pf%particles(i)
        tag = 1
        CALL MPI_ISEND(pf%psi(:,i), state_dim , MPI_DOUBLE_PRECISION, &
             particle-1, tag, CPL_MPI_COMM, requests(i), mpi_err)
        !PRINT*,'Particle filter ',pfrank,'has sent state_vector over mpi at iteratio&
        !    &n',iter,' to ensemble member ',particle
     END DO
     CALL MPI_WAITALL(pf%count,requests,mpi_statuses, mpi_err)

   
!     do particle =1,pf%count
!        call receive_from_model(fpsi(:,particle),particle)
!     end do
     DO i = 1,pf%count
        particle = pf%particles(i)
        tag = 1
        CALL MPI_IRECV(fpsi(:,i), state_dim, MPI_DOUBLE_PRECISION, &
             particle-1, tag, CPL_MPI_COMM,&
             requests(i), mpi_err)
     END DO
     CALL MPI_WAITALL(pf%count,requests,mpi_statuses, mpi_err)


  end if

  call mpi_barrier(pf_mpi_comm,mpi_err)
  psimean = 0.0_rk

  
  do i = 1,pf%count
     particle = pf%particles(i)
     if(c(i) .le. cmax) then
        if(.not. pf%gen_data) then

        !a(particle) = 0.5* (y-Hf(x_i^n-1))^T R^(-1) (HQH^T) (HQH^T + R)^(-1) (y-Hf(x_i^n-1))
        
        call K(y_Hfpsin1(:,i),kgain)
        print*,'i = ',i,' ||kgain||_2 = ',dnrm2(state_dim,kgain,1)
        call flush(6)
        call H(kgain,obsv)
        print*,'i = ',i,' ||obsv||_2 = ',dnrm2(obs_dim,obsv,1)
        call flush(6)
        call solve_r(obsv,obsvv)
        print*,'i = ',i,' ||obsvv||_2 = ',dnrm2(obs_dim,obsvv,1)
        call flush(6)
        e = sum(obsvv*y_Hfpsin1(:,i))
        
        a(i) = 0.5*e
        
           !b(particle) = 0.5* (y-Hf(x_i^n-1))^T R^(-1) (y-Hf(x_i^n-1))- cmax + pf%weight(i)
        
        call innerR_1(y_Hfpsin1(:,i),e)
        
        b(i) = 0.5*e - cmax + pf%weight(particle)
        
        !note the plus sign in the below equation. See Ades & van Leeuwen 2012.
        alpha(i) = 1.0 + sqrt(1.0 - b(i)/a(i) + 1.0D-6)

        print*,'i a b alpha'
        print*,i,a(i),b(i),alpha(i)
     else !if(.not. pf%gen_data) 
        kgain = 0.0_rk
        alpha(i) = 0.0_rk
        y_Hfpsin1(:,i) = 0.0_rk
     end if !if(.not. pf%gen_data)

     !generate beta from a mixture density
     call MixtureRandomNumbers1D(0.0D0,pf%nfac,pf%ufac,pf%efac,state_dim,statev,uniform)
     call Qhalf(statev,betan)

     if(uniform) then
        pf%weight(particle) = pf%weight(particle) +&
             (alpha(i)**2.0_rk - 2.0_rk*alpha(i))*a(i) + & 
             0.5_rk*e
        
     else
      pf%weight(particle) = pf%weight(particle) +&
             (alpha(i)**2.0_rk - 2.0_rk*alpha(i))*a(i) + &
             0.5_rk*e &
             + 2**(-real(state_dim,rk)/2.0_rk)*pi**(real(state_dim,rk)&
             &/2.0_rk)*pf%nfac*pf%ufac**(-real(state_dim,rk))*((1.0_rk&
             &-pf%efac)/pf%efac)*exp(0.5_rk*(sum(betan*betan)))
              
     end if !if(uniform)

     !now do the following
     !x^n = f(x^(n-1)) + alpha(i) K (y-Hf(x_i^n-1)) + beta

!     pf%psi(:,i) = fpsi(:,i) + alpha(i)*kgain + betan
     call update_state(pf%psi(:,i),fpsi(:,i),alpha(i)*kgain,betan)

     psimean = psimean + pf%psi(:,i)

     !now calculate the new weights

!     pf%weight(particle) = pf%weight(particle) + (alpha(particle)&
!          &**2.0_rk - 2.0_rk*alpha(particle))*a(particle) + 0.5_rk*e
  else
     pf%weight(particle) = huge(1.0D0)
  end if !if(c(particle) .le. cmax)
  end do
  !print*,'after da looooop'
!  print*,'at end of equal weight filter pf%weight = '
!  print*,pf%weight
  
!  call resample
  
  if(pf%gen_data) then
!     write(6,*) 'generating the data'
!     call flush(6)
     call H(pf%psi(:,1),y)
!     write(6,*) 'after H'
!     call flush(6)
     call NormalRandomNumbers1D(0.0D0,1.0D0,obs_dim,obsv)
     call rhalf(obsv,obsvv)
!     write(6,*) 'after rhalf'
!     call flush(6)
     y = y + obsvv
     call save_observation_data(y)
!     call save_truth(pf%psi(:,1))
  else 
     if(pf%use_talagrand) call diagnostics
     !print*,'entering resample step'
     print*,'time until resample = ',mpi_wtime()-pf%time
     call flush(6)
     open(4,file='done',action='write',status='replace')
     close(4)
     call resample
  end if !if(pf%gen_data)
  
  
end subroutine equal_weight_filter
