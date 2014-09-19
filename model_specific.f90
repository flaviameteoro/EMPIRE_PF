!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-19 10:33:19 pbrowne>
!!!
!!!    This file must be adapted to the specific model in use.
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
subroutine configure_model
  !called initially to get up details and data for model specific functions
  use pf_control
  use sizes
  use Qdata
  use Rdata
  implicit none
  include 'mpif.h'
  
  real(kind=kind(1.0d0)) :: t1
!this is for Lorenz 63
!  state_dim=3
!  obs_dim = 1

  !this is for hadcm3
  state_dim = 2314430
  obs_dim = 27370

  print*,'#################################'
  print*,'######### SANITY CHECK ##########'
  print*,'#################################'
  print*,'## STATE DIMENSION = ',state_dim
  print*,'##  OBS  DIMENSION = ',obs_dim
  print*,'#################################'
  if(.not. pf%gen_Q) then
     t1 = mpi_wtime()
     call loadQ
     print*,'load Q     took ',mpi_wtime()-t1,' seconds'
     t1 = mpi_wtime()
     call loadR
     print*,'load R     took ',mpi_wtime()-t1,' seconds'
     t1 = mpi_wtime()
!     call load_HQHTR
     print*,'load HQHTR took ',mpi_wtime()-t1,' seconds'
     
  end if
end subroutine configure_model



subroutine solve_r(y,v,t)
  !subroutine to take an observation vector y and return v
  !in observation space.
  use sizes
  use pf_control
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: t !the timestep
  real(kind=rk), dimension(obs_dim,pf%count), intent(in) :: y
  real(kind=rk), dimension(obs_dim,pf%count), intent(out) :: v

  v = y/(0.3d0**2)
  !stop 'Solve_r not yet implemented'
  
end subroutine solve_r


subroutine solve_hqht_plus_r(y,v,t)
  !subroutine to take an observation vector y and return v
  !in observation space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: t !the timestep
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim), intent(out) :: v

  stop 'solve_hqht_plus_r not yet implemented'


end subroutine solve_hqht_plus_r

subroutine Q(nrhs,x,Qx)
  !subroutine to take a full state vector x and return Qx
  !in state space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: Qx
  real(kind=rk), dimension(state_dim,nrhs) :: temp

  call Qhalf(nrhs,x,temp)

  call Qhalf(nrhs,temp,Qx)
  

end subroutine Q

subroutine Qhalf(nrhs,x,Qx)
  !subroutine to take a full state vector x and return Q^{1/2}x
  !in state space.
  use sizes
  use Qdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: qx

  qx = 5.3d3*x
!  stop 'Qhalf not yet implemented'
  
end subroutine Qhalf


subroutine R(nrhs,y,Ry,t)
  !subroutine to take an observation vector x and return Rx
  !in observation space.
  use sizes
  use Rdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs
  integer, intent(in) :: t !the timestep
  real(kind=rk), dimension(obs_dim,nrhs), intent(in) :: y
  real(kind=rk), dimension(obs_dim,nrhs), intent(out) :: Ry

!  stop 'R not yet implemented'
  Ry = 0.3d0**2*y

end subroutine R

subroutine Rhalf(nrhs,y,Ry,t)
  !subroutine to take an observation vector x and return Rx
  !in observation space.
  use sizes
  use Rdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs
  integer, intent(in) :: t !the timestep
  real(kind=rk), dimension(obs_dim,nrhs), intent(in) :: y
  real(kind=rk), dimension(obs_dim,nrhs), intent(out) :: Ry

!  stop 'Rhalf not yet implemented'
  Ry = 0.3d0*y

end subroutine RHALF

subroutine H(x,hx,t)
  !subroutine to take a full state vector x and return H(x)
  !in observation space.
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: t !the timestep
  real(kind=rk), dimension(state_dim,pf%count), intent(in) :: x
  real(kind=rk), dimension(obs_dim,pf%count), intent(out) :: hx

  stop 'H not yet implemented'
  !hx(:,:) = x(539617:566986,:)

end subroutine H

subroutine HT(y,x,t)
  !subroutine to take an observation vector y and return x = H^T(y)
  !in full state space.
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: t !the timestep
  real(kind=rk), dimension(state_dim,pf%count), intent(out) :: x
  real(kind=rk), dimension(obs_dim,pf%count), intent(in) :: y

  stop 'HT not yet implemented'
  !x = 0.0_rk
  !x(539617:566986,:) = y(:,:)

end subroutine HT

