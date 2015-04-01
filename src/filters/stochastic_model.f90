!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-04-01 22:19:49 pbrowne>
!!!
!!!    subroutine to simply move the model forward in time one timestep
!!!    then add model error
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
!> @brief
!> subroutine to simply move the model forward in time one 
!> timestep
!> PAB 21-05-2013

subroutine stochastic_model
  use pf_control
  use Sizes
  use comms

  IMPLICIT NONE
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)

  real(kind=rk), dimension(state_dim,pf%count) :: normaln     !vector to store uncorrelated random error
  real(kind=rk), dimension(state_dim,pf%count) :: betan       !vector to store sqrtQ correlated random error
  real(kind=rk), dimension(state_dim,pf%count) :: fpsi        !f(psi^(n-1))
!  real(kind=rk), dimension(state_dim,pf%count) :: kgain       !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  real(kind=rk), dimension(state_dim,pf%count) :: Qkgain
!  real(kind=rk) :: t,dnrm2
  integer :: particle,k!,tag,mpi_err
!  integer :: mpi_status( MPI_STATUS_SIZE )
  logical, parameter :: checkscaling=.true.

  real(kind=rk), dimension(obs_dim) :: obsv,obsvv,y

  call send_all_models(state_dim,pf%count,pf%psi,1)
!  do k =1,pf%count
!     particle = pf%particles(k)
!     tag = 1
!     call mpi_send(pf%psi(:,k),state_dim,MPI_DOUBLE_PRECISION&
!          &,particle-1,tag,CPL_MPI_COMM,mpi_err)
!  end do

  call NormalRandomNumbers2D(0.0D0,1.0D0,state_dim,pf%count,normaln)

  call Qhalf(pf%count,normaln,betan)
  Qkgain = 0.0_rk

!  DO k = 1,pf%count
!     particle = pf%particles(k)
!     tag = 1
!     CALL MPI_RECV(fpsi(:,k), state_dim, MPI_DOUBLE_PRECISION, &
!          particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
!  END DO
  call recv_all_models(state_dim,pf%count,fpsi)

  !$omp parallel do private(particle)
  DO k = 1,pf%count
     particle = pf%particles(k)
     call update_state(pf%psi(:,k),fpsi(:,k),Qkgain(:,k),betan(:,k))
  end DO
  !$omp end parallel do


  if(pf%gen_data .and. mod(pf%timestep,pf%time_bwn_obs) .eq. 0) then
     if(pf%count .ne. 1 .and. pf%nens .ne. 1) then
        print*,'OBS GEN ERROR -558: PLEASE RUN WITH ONLY A SINGLE &
             &ENSEMBLE MEMBER'
        stop -558
     end if

     write(6,*) 'generating the data'
     call flush(6)

     
     !get model equivalent of observations
     call H(obs_dim,1,pf%psi,y,pf%timestep)

     !generate uncorrelated random vector obsv
     call NormalRandomNumbers1D(0.0D0,1.0D0,obs_dim,obsv)

     !turn this into correlated random observation noise
     call rhalf(obs_dim,1,obsv,obsvv,pf%timestep)

     !add the noise to the model equivalent obs
     y = y + obsvv

     !output the observation to a file
     call save_observation_data(y)

     !call diagnostics to save things like the truth etc
     call diagnostics
  end if

end subroutine stochastic_model


