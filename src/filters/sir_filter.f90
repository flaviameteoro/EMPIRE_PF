!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-01-27 10:19:44 pbrowne>
!!!
!!!    Subroutine to perform SIR filter
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
!> Subroutine to perform SIR filter (Sequential Importance Resampling)
subroutine sir_filter
  use comms
  use sizes
  use pf_control
  implicit none
  include 'mpif.h'
  integer :: k,tag,particle
  real(kind=kind(1.0d0)), dimension(state_dim) :: zeros
  real(kind=kind(1.0d0)), dimension(pf%count) :: w
  real(kind=kind(1.0d0)), dimension(state_dim,pf%count) :: fpsi,normaln,betan
  real(kind=kind(1.0d0)), dimension(obs_dim) :: y
  real(kind=kind(1.0d0)), dimension(obs_dim,pf%count) :: y_Hfpsi,Hfpsi
  integer, dimension(mpi_status_size) :: mpi_status
  integer :: mpi_err
  zeros = 0.0d0

  !get the model to provide f(x)
  do k =1,pf%count
     particle = pf%particles(k)
     tag = 1
     call mpi_send(pf%psi(:,k),state_dim,MPI_DOUBLE_PRECISION&
          &,particle-1,tag,CPL_MPI_COMM,mpi_err)
  end do

  DO k = 1,pf%count
     particle = pf%particles(k)
     tag = 1
     CALL MPI_RECV(fpsi(:,k), state_dim, MPI_DOUBLE_PRECISION, &
          particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
  END DO


  !draw from a Gaussian for the random noise
  call NormalRandomNumbers2D(0.0D0,1.0D0,state_dim,pf%count,normaln)
  
  !compute the relaxation term Qkgain, the intermediate
  !term kgain and apply correlation to noise
  call Qhalf(pf%count,normaln,betan)

  
  !update the new state and weights based on these terms
  !$omp parallel do
  DO k = 1,pf%count
!     print*,'|fpsi-psi|_2 = ',dnrm2(state_dim,(fpsi(:,k)-pf%psi(:,k)),1)
     call update_state(pf%psi(:,k),fpsi(:,k),zeros,betan(:,k))
  end DO
  !$omp end parallel do







if(mod(pf%timestep,pf%time_bwn_obs) .eq. 0) then
   call get_observation_data(y)
   !this is the analysis step.
   
   call H(obs_dim,pf%count,fpsi,Hfpsi,pf%timestep)

   !$omp parallel do
   do k = 1,pf%count
      y_Hfpsi(:,k) = y - Hfpsi(:,k)
   end do
   !$omp end parallel do
    
   call innerR_1(obs_dim,pf%count,y_Hfpsi,w,pf%timestep)

   do k = 1,pf%count
      particle = pf%particles(k)
      pf%weight(particle) = 0.5*w(k)
   end do

   call resample
end if


end subroutine sir_filter
