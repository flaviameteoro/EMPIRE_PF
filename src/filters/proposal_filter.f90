!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-09-09 13:37:50 pbrowne>
!!!
!!!    Subroutine to perform nudging in the proposal step of EWPF
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine proposal_filter based on subroutine IntegrateModel
!given from the Particle filter code of Mel and Peter Jan
!PAB 04-02-2013

!> Subroutine to perform nudging in the proposal step of EWPF
subroutine proposal_filter
  use pf_control
  use Sizes
  use comms
  
  IMPLICIT NONE
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)

  real(kind=rk) :: pWeight
  real(kind=rk), dimension(state_dim,pf%count) :: normaln     !vector to store uncorrelated random error
  real(kind=rk), dimension(state_dim,pf%count) :: betan       !vector to store sqrtQ correlated random error
  real(kind=rk), dimension(obs_dim) :: y             !y, the observations
  real(kind=rk), dimension(obs_dim,pf%count) :: Hpsi          !H(psi^(n-1))
  real(kind=rk), dimension(obs_dim,pf%count) :: y_Hpsin1      !y-H(psi^(n-1))
  real(kind=rk), dimension(state_dim,pf%count) :: fpsi        !f(psi^(n-1))
  real(kind=rk), dimension(state_dim,pf%count) :: kgain       !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  real(kind=rk), dimension(state_dim,pf%count) :: Qkgain
  integer :: particle,k
  integer :: mpi_err
  real(kind=rk) :: pweighttemp,ddot

  !get the model to provide f(x)
  call send_all_models(state_dim,pf%count,pf%psi,1)


  !get the next observations and store it in vector y
  call get_observation_data(y,pf%timestep)


  !compute y - H(x)
  call H(obs_dim,pf%count,pf%psi,Hpsi,pf%timestep)

  !$omp parallel do
  do k = 1,pf%count
     y_Hpsin1(:,k) = y - Hpsi(:,k)
  end do
  !$omp end parallel do


  !draw from a Gaussian for the random noise
  call NormalRandomNumbers2D(0.0D0,1.0D0,state_dim,pf%count,normaln)

  call recv_all_models(state_dim,pf%count,fpsi)

  
  !compute the relaxation term Qkgain, the intermediate
  !term kgain and apply correlation to noise
  call Bprime(y_Hpsin1,kgain,Qkgain,normaln,betan)


  !update the new state and weights based on these terms
  !$omp parallel do private(particle,pweight)
  DO k = 1,pf%count
     particle = pf%particles(k)
     !     print*,'|fpsi-psi|_2 = ',dnrm2(state_dim,(fpsi(:,k)-pf%psi(:,k)),1)
     call update_state(pf%psi(:,k),fpsi(:,k),Qkgain(:,k),betan(:,k))

!     pweight = sum(Qkgain(:,k)*kgain(:,k))+2.0D0*sum(betan(:,k)*kgain(:,k))
     pweight = ddot(state_dim,Qkgain(:,k),1,kgain(:,k),1) + 2.0D0&
             &*ddot(state_dim, betan(:,k),1,kgain(:,k),1)

     if(empire_version .eq. 3) then
        !need to perform the sum across all parts of the state vector
        pweighttemp=pweight
        call mpi_allreduce(pweighttemp,pweight,1,MPI_DOUBLE_PRECISION,MPI_SUM&
             &,pf_member_comm,mpi_err)
     end if

     pf%weight(particle) = pf%weight(particle) + 0.5*pWeight
  end DO
  !$omp end parallel do

end subroutine proposal_filter
