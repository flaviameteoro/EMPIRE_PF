!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-07-17 11:11:24 pbrowne>
!!!
!!!    Computes the equal weights step in the New Schem Equal Weights Particle Filter
!!!    Copyright (C) 2015  Mengbin Zhu
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
!!!    Email: p.browne @ reading.ac.uk  EMPIRE
!!!           zhumengbin@gmail.com      new scheme equal weights particle filter
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> subroutine to do the new scheme equal weights last time-step
!! 
!! structure of code loosely based on original equivalent weights
!! scheme @ref equivalent_weights_filter
subroutine equivalent_weights_filter_zhu
  use timestep_data
  use pf_control
  use sizes
  use random
  use comms
  implicit none
  include 'mpif.h'
  integer, parameter                             :: rk = kind(1.0D0) !<specify double precision
  real(kind=rk), dimension(pf%count)             :: alpha,c
  real(kind=rk), dimension(pf%count)             :: epsiloni,epsilon_n,epsilon_p,gammai,ce,win,wsolve
  real(kind=rk), dimension(pf%nens)              :: csorted !<sorted vector of c
  real(kind=rk), dimension(pf%nens)              :: weight_an !<weight for the analysis steps
  real(kind=rk)                                  :: cmax,xs1,xs2
  integer                                        :: particle,i,mpi_err,rc!,tag
  real(kind=rk)                                  :: percentage,randc
  real(kind=rk), dimension(obs_dim)              :: y     !y, !//!<the observations
  real(kind=rk), dimension(obs_dim,pf%count)     :: Hfpsi           !H(f(psi^(n-1))) !< \f$H(f(x^{n-1}))\f$
  real(kind=rk), dimension(obs_dim,pf%count)     :: y_Hfpsin1  !y-H(f(psi^(n-1))) !< \f$y-H(f(x^{n-1}))\f$
  real(kind=rk), dimension(state_dim,pf%count)   :: fpsi     !f(psi^(n-1)) !< \f$f(x^{n-1})\f$
  real(kind=rk), dimension(state_dim)            :: psimean !< the mean of the state vectors
  real(kind=rk), dimension(state_dim,pf%count)   :: kgain !QH^T(HQH^T+R)^(-1)(y-H(f(psi^(n-1)))) !< \f$QH^T(HQH^T+R)^{-1}(y-H(f(x^{n-1})))\f$
  real(kind=rk), dimension(state_dim,pf%count)   :: betan         !the random variable
  real(kind=rk), dimension(state_dim,pf%count)   :: statev        !<temporary state space vector 
  real(kind=rk)                                  :: w
  ! real(kind=rk), dimension(pf%count)             :: e          !e = d_i^t R^(-1) d_i
  real(kind=rk), parameter                       :: pi = 4.0D0*atan(1.0D0)
  ! INTEGER, DIMENSION(MPI_STATUS_SIZE)           :: mpi_status
  real(kind=rk), dimension(pf%count)             :: weight_temp2
  real(kind=rk)                                  :: ddot

  percentage = 0.5
  rc=1

  !start the model integration
  call send_all_models(state_dim,pf%count,pf%psi,1)

  !get the next observation and store it in vector y
  call get_observation_data(y,pf%timestep)

  !recv the models after integration
  call recv_all_models(state_dim,pf%count,fpsi)

  call H(obs_dim,pf%count,fpsi,Hfpsi,pf%timestep)

  !compute c for each particle on this mpi thread
  do i = 1,pf%count
     particle = pf%particles(i)
     y_Hfpsin1(:,i) = y - Hfpsi(:,i)

     call innerHQHt_plus_R_1(y_Hfpsin1(:,i),w,pf%timestep)

     c(i) = 2*pf%weight(particle) + w
  end do


  !communicate c to all the mpi threads
  call mpi_allgatherv(c,pf%count,mpi_double_precision,csorted,gblcount&
       &,gbldisp,mpi_double_precision,pf_mpi_comm,mpi_err)

  !calculate cmax
  call quicksort_d(csorted,pf%nens)
  cmax = csorted(pf%nens)


  psimean = 0.0_rk


  !compute the kalman gain
  call K(y_Hfpsin1,kgain)

  !generate the unit normal random numbers for all particles ~ xi
  call normalrandomnumbers2d(0.0D0,1.0D0,state_dim,pf%count,statev)

  !calculate epsiloni for the particles
  xs1=-0.90D0
  xs2=2.1D0
  do i = 1,pf%count
     gammai(i) = ddot(state_dim,statev(:,i),1,statev(:,i),1)
     ce(i) = cmax - c(i)

     win(i) = (-gammai(i)/state_dim)*exp(-gammai(i)/state_dim)*exp(-ce(i)/state_dim)
     call LambertW(0,win(i),wsolve(i))
     if(ce(i) .eq. 0.0) then
        wsolve(i) = -(gammai(i)/state_dim)
     end if
     epsilon_n(i) = (-state_dim/gammai(i))*wsolve(i) - 1

     !call newton(xs1,gammai(i),state_dim,c(i),cmax,epsilon_n(i))

     call LambertW(-1,win(i),wsolve(i))
     epsilon_p(i) = (-state_dim/gammai(i))*wsolve(i) - 1

     !call newton(xs2,gammai(i),state_dim,c(i),cmax,epsilon_p(i))
  end do

  do i = 1,pf%count
     call random_number(randc)
     if(randc .ge. 0.5) then
        epsiloni(i) = epsilon_p(i)
     else
        epsiloni(i) = epsilon_n(i)
     end if
  end do


  !compute alpha for each particle on this mpi thread
  do i = 1,pf%count
     alpha(i) = 1.0 + epsiloni(i)
  end do

  !compute the random term of the analysis 
  call Phalf(pf%count,statev,betan)

  !update the weights and the new state
  do i = 1,pf%count
     particle = pf%particles(i)
     weight_an(particle) = 0.0
     if(alpha(i) .eq. 0.0) then
        alpha(i) = alpha(i) + 1.0E-16
     end if
     weight_an(particle) = gammai(i)*epsiloni(i)-state_dim*log(alpha(i))+c(i)
     weight_an(particle) = 0.5*weight_an(particle)
     if(alpha(i) .eq. 1.0E-16) then
        weight_an(particle) = 0.5*cmax
     end if
     pf%Weight(particle) = weight_an(particle)

     !now do the following
     !x^n = f(x^(n-1)) + K (y-Hf(x_i^n-1)) + alpha(i)^(1/2) P^(1/2) xi
     call update_state(pf%psi(:,i),fpsi(:,i),kgain(:,i),sqrt(alpha(i))*betan(:,i))
     psimean = psimean + pf%psi(:,i)
  end do

  !store in weight_temp2 only those weights on this mpi thread
  weight_temp2 = -huge(1.0d0)
  do i = 1,pf%count
     weight_temp2(i) = pf%weight(pf%particles(i))
  end do

  !communicate the weights of all particles to each mpi thread
  call mpi_allgatherv(weight_temp2,pf%count,mpi_double_precision,pf%weight,gblcount&
       &,gbldisp,mpi_double_precision,pf_mpi_comm,mpi_err)

  !normalise the weights
  pf%weight = exp(-pf%weight+maxval(pf%weight))
  pf%weight = pf%weight/sum(pf%weight)
  pf%weight = -log(pf%weight)
  !=========================================================================


  call timestep_data_set_is_analysis
  if(pf%use_talagrand) call diagnostics

end subroutine equivalent_weights_filter_zhu
