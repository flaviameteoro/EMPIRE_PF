!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-09-09 14:17:39 pbrowne>
!!!
!!!    Computes the equivalent weights step in the EWPF
!!!    Copyright (C) 2014  Philip A. Browne
!!!
!!!    This program is free software: you can redistribute it and/or modify
!!!    it under the terms of the GNU General Public License as published by
!!!    the Free Software Foundation, either version 3 of the License, or
!!!    (at your option) any later version.
!!!
!!!    This program is distributed in the hope that it will be useful,
!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!!    GNU General Public License for more details.
!!!
!!!    You should have received a copy of the GNU General Public License
!!!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!!
!!!    Email: p.browne @ reading.ac.uk
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> subroutine to do the equivalent weights step
!!
subroutine equivalent_weights_filter
  use timestep_data
  use pf_control
  use sizes
  use random
  use comms
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0) !<specify double precision
  real(kind=rk), dimension(pf%count) :: a,b,alpha,c
  real(kind=rk), dimension(pf%nens) :: csorted !<sorted vector of c
  real(kind=rk) :: cmax
  integer :: particle,i,mpi_err!,tag
  real(kind=rk), dimension(obs_dim) :: y     !y, !//!<the observations
  real(kind=rk), dimension(obs_dim,pf%count) :: Hfpsi           !H(f(psi^(n-1))) !< \f$H(f(x^{n-1}))\f$
  real(kind=rk), dimension(obs_dim,pf%count) :: y_Hfpsin1  !y-H(f(psi^(n-1))) !< \f$y-H(f(x^{n-1}))\f$
  real(kind=rk), dimension(state_dim,pf%count) :: fpsi     !f(psi^(n-1)) !< \f$f(x^{n-1})\f$
  real(kind=rk), dimension(state_dim) :: psimean !< the mean of the state vectors
  real(kind=rk), dimension(state_dim,pf%count) :: kgain !QH^T(HQH^T+R)^(-1)(y-H(f(psi^(n-1)))) !< \f$QH^T(HQH^T+R)^{-1}(y-H(f(x^{n-1})))\f$
  real(kind=rk), dimension(state_dim,pf%count) :: betan         !the mixture random variable
  real(kind=rk), dimension(state_dim,pf%count) :: statev        !<temporary state space vector 
  real(kind=rk), dimension(obs_dim,pf%count) :: obsv,obsvv      !<temporary  obs  space vector
  real(kind=rk) :: w
  real(kind=rk), dimension(pf%count) :: e                     !e = d_i^t R^(-1) d_i
  real(kind=rk), parameter :: pi = 4.0D0*atan(1.0D0)
  logical, dimension(pf%count) :: uniform

  real(kind=rk) :: ddot,wtemp,betanTbetan
  integer :: ensemble_comm

  if(comm_version .eq. 1 .or. comm_version .eq. 2) then
     ensemble_comm = pf_mpi_comm
  elseif(comm_version .eq. 3) then
     ensemble_comm = pf_ens_comm
  else
     print*,'EMPIRE VERSION ',comm_version,' NOT SUPPORTED IN proposal_filter'
     print*,'THIS IS AN ERROR. STOPPING'
     stop '-24'
  end if


  call send_all_models(state_dim,pf%count,pf%psi,1)
  
  !get the next observation and store it in vector y
  call get_observation_data(y,pf%timestep)
  
  call recv_all_models(state_dim,pf%count,fpsi)

  call H(obs_dim,pf%count,fpsi,Hfpsi,pf%timestep)
  
  !compute c for each particle on this mpi thread
  do i = 1,pf%count
     particle = pf%particles(i)
     y_Hfpsin1(:,i) = y - Hfpsi(:,i)
     
     call innerHQHt_plus_R_1(y_Hfpsin1(:,i),w,pf%timestep)

     if(comm_version .eq. 3) then
        !need to perform the sum across all parts of the state vector
        wtemp=w
        call mpi_allreduce(wtemp,w,1,MPI_DOUBLE_PRECISION,MPI_SUM&
             &,pf_member_comm,mpi_err)
     end if
     
     c(i) = pf%weight(particle) + 0.5d0*w
  end do
     
  
  !communicate c to all the mpi threads
  call mpi_allgatherv(c,pf%count,MPI_DOUBLE_PRECISION,csorted,gblcount&
       &,gbldisp,MPI_DOUBLE_PRECISION,ensemble_comm,mpi_err)
  
  !calculate cmax
  call quicksort_d(csorted,pf%nens)
  cmax = csorted(nint(pf%keep*pf%nens))

  psimean = 0.0_rk
  
  !compute the kalman gain
  call K(y_Hfpsin1,kgain)
  call H(obs_dim,pf%count,kgain,obsv,pf%timestep)
  call solve_r(obs_dim,pf%count,obsv,obsvv,pf%timestep)
  
  !compute a for each particle on this mpi thread
  do i = 1,pf%count
     a(i) = 0.5d0*ddot(obs_dim,obsvv(:,i),1,y_Hfpsin1(:,i),1)

     if(comm_version .eq. 3) then
        !need to perform the sum across all parts of the observation vector
        wtemp=a(i)
        call mpi_allreduce(wtemp,a(i),1,MPI_DOUBLE_PRECISION,MPI_SUM&
             &,pf_member_comm,mpi_err)
     end if
  end do
  
  call innerR_1(obs_dim,pf%count,y_Hfpsin1,e,pf%timestep)
  
  !compute alpha for each particle on this mpi thread
  do i = 1,pf%count
     particle = pf%particles(i)
     b(i) = 0.5d0*e(i) - cmax + pf%weight(particle)
     !note the plus sign in the below equation. See Ades & van Leeuwen 2012.
     alpha(i) = 1.0d0 + sqrt(1.0d0 - b(i)/a(i) + 1.0D-6)

     if(alpha(i) .ne. alpha(i)) alpha(i) = 1.0d0 !ensure that
     !alpha(i) is not a NaN because of roundoff errors.
  end do
  
  
  !draw from a mixture density for the random noise then correlate it
  call MixtureRandomNumbers2D(0.0D0,pf%nfac,pf%ufac,pf%efac,state_dim&
       &,pf%count,statev,uniform)
  call Qhalf(pf%count,statev,betan)
  
  !update the weights and the new state
  do i = 1,pf%count
     if(c(i) .le. cmax) then
        particle = pf%particles(i)
        if(uniform(i)) then
           pf%weight(particle) = pf%weight(particle) +&
                (alpha(i)**2.0_rk - 2.0_rk*alpha(i))*a(i) + & 
                0.5_rk*e(i)
        else
           betanTbetan = sum(betan(:,i)*betan(:,i))
           if(comm_version .eq. 3) then
              !need to perform the sum across all parts of the state vector
              wtemp=betanTbetan
              call mpi_allreduce(wtemp,betanTbetan,1,MPI_DOUBLE_PRECISION,MPI_SUM&
                   &,pf_member_comm,mpi_err)
           end if
           pf%weight(particle) = pf%weight(particle) +&
                (alpha(i)**2.0_rk - 2.0_rk*alpha(i))*a(i) + &
                0.5_rk*e(i) &
                + 2**(-real(state_dim_g,rk)/2.0_rk)*pi**(real(state_dim_g,rk)&
                &/2.0_rk)*pf%nfac*pf%ufac**(-real(state_dim_g,rk))*((1.0_rk&
                &-pf%efac)/pf%efac)*exp(0.5_rk*(betanTbetan))        
        end if !if(uniform)
        
        !now do the following
        !x^n = f(x^(n-1)) + alpha(i) K (y-Hf(x_i^n-1)) + beta
        call update_state(pf%psi(:,i),fpsi(:,i),alpha(i)*kgain(:,i)&
             &,betan(:,i))
        psimean = psimean + pf%psi(:,i)
     else
        pf%weight(pf%particles(i)) = huge(1.0D0)
     end if
  end do
  

!=========================================================================



  if(pf%use_talagrand) call diagnostics
  !print*,'entering resample step'
  print*,'time until resample = ',mpi_wtime()-pf%time
  call flush(6)
  call resample
  
  call timestep_data_set_is_analysis
end subroutine equivalent_weights_filter
