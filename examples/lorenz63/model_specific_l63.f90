!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-10-13 13:20:18 pbrowne>
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

!> subroutine called initially to set up details and data
!> for model specific functions
subroutine configure_model
  use pf_control
  use timestep_data
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
  state_dim = 3
  obs_dim = 3

  call timestep_data_set_total(pf%time_bwn_obs*pf%time_obs)

  print*,'#################################'
  print*,'######### SANITY CHECK ##########'
  print*,'#################################'
  print*,'## STATE DIMENSION = ',state_dim
  print*,'##  OBS  DIMENSION = ',obs_dim
  print*,'#################################'
!  if(.not. pf%gen_Q) then
     t1 = mpi_wtime()
     call loadQ
     print*,'load Q     took ',mpi_wtime()-t1,' seconds'
     t1 = mpi_wtime()
     call loadR
     print*,'load R     took ',mpi_wtime()-t1,' seconds'
     t1 = mpi_wtime()
!     call load_HQHTR
     print*,'load HQHTR took ',mpi_wtime()-t1,' seconds'
     
!  end if

     pf%time_bwn_obs = 8 !set this so the obs are picked up corrently
end subroutine configure_model


!>subroutine to reset variables that may change when the observation
!!network changes
subroutine reconfigure_model
end subroutine reconfigure_model


!>subroutine to take an observation vector y and return v
!!  in observation space.
!!
!! Given \f$y\f$ find \f$v\f$ such that \f$Rv=y\f$
subroutine solve_r(obsDim,nrhs,y,v,t)
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: obsDim !< the dimension of the observations
  integer, intent(in) :: nrhs   !< the number of right hand sides
  real(kind=rk), dimension(obsdim,nrhs), intent(in) :: y !<
                                                              !!input vector
  real(kind=rk), dimension(obsdim,nrhs), intent(out) :: v!<
  !!result vector where \f$v=R^{-1}y\f$
  integer, intent(in) :: t !<the timestep

  !R==1
  v = y

  
end subroutine solve_r

!>subroutine to take an observation vector y and return v
!!  in observation space.
!!
!! Given \f$y\f$ find \f$v\f$ such that \f$R^{\frac{1}{2}}v=y\f$
subroutine solve_rhalf(obsdim,nrhs,y,v,t)
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: obsDim !< the dimension of the observations
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(obsdim,nrhs), intent(in) :: y !<
                                                              !!input vector
  real(kind=rk), dimension(obsdim,nrhs), intent(out) :: v!<
  !!result vector where \f$v=R^{-\frac{1}{2}}y\f$
  integer, intent(in) :: t !<the timestep

  !R==1
  v = y

  
end subroutine solve_rhalf



!>subroutine to take an observation vector y and return v
!>in observation space.
!!
!! Given \f$y\f$ find \f$v\f$ such that \f$(HQH^T+R)v=y\f$
subroutine solve_hqht_plus_r(obsdim,y,v,t)
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: obsdim !< the dimension of the observations
  real(kind=rk), dimension(obsdim), intent(in) :: y !<the input vector
  real(kind=rk), dimension(obsdim), intent(out) :: v !< the result
  !! where \f$v = (HQH^T+R)^{-1}y\f$
  integer, intent(in) :: t !<the timestep
  
  !Q==0, R==1, H=I
  v = y

end subroutine solve_hqht_plus_r

!> subroutine to take a full state vector x and return Qx
!> in state space.
!!
!! Given \f$x\f$ compute \f$Qx\f$
subroutine Q(nrhs,x,Qx)

  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x !< the
  !!input vector
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: Qx !< the
  !!resulting vector where Qx \f$= Qx\f$
  real(kind=rk), dimension(state_dim,nrhs) :: temp

  Qx = 0.0d0
  
end subroutine Q



!> subroutine to take a full state vector x and return \f$Q^{1/2}x\f$
!> in state space.
!!
!! Given \f$x\f$ compute \f$Q^{\frac{1}{2}}x\f$
subroutine Qhalf(nrhs,x,Qx)
  use sizes
  use Qdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x !< the
  !!input vector
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: qx !< the
  !!resulting vector where Qx \f$= Q^{\frac{1}{2}}x\f$

  qx = 0.0d0
  
end subroutine Qhalf

!> subroutine to take an observation vector x and return Rx
!> in observation space.
!!
!! Given \f$y\f$ compute \f$Ry\f$
subroutine R(obsDim,nrhs,y,Ry,t)
  use Rdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: obsDim !< the dimension of the observations
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(obsDim,nrhs), intent(in) :: y !< the
  !!input vector
  real(kind=rk), dimension(obsDim,nrhs), intent(out) :: Ry !< the
  !!resulting vectors where Ry \f$= Ry\f$
  integer, intent(in) :: t !< the timestep


  !R ==1
  Ry = y

end subroutine R

!> subroutine to take an observation vector x and return Rx
!> in observation space.
!!
!! Given \f$y\f$ compute \f$R^{\frac{1}{2}}y\f$
subroutine Rhalf(obsDim,nrhs,y,Ry,t)
  use Rdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: obsDim !< the dimension of the observations
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(obsDim,nrhs), intent(in) :: y !< the
  !!input vector
  real(kind=rk), dimension(obsDim,nrhs), intent(out) :: Ry !< the
  !!resulting vector where Ry \f$= R^{\frac{1}{2}}y\f$
  integer, intent(in) :: t !<the timestep


  Ry = y

end subroutine RHALF


!> subroutine to take a full state vector x and return H(x)
!> in observation space.
!!
!! Given \f$x\f$ compute \f$Hx\f$
subroutine H(obsDim,nrhs,x,hx,t)
!  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: obsDim !< the dimension of the observations
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x !< the
  !!input vectors in state space
  real(kind=rk), dimension(obsDim,nrhs), intent(out) :: hx !< the
  !!resulting vector in observation space  where hx \f$= Hx\f$
  integer, intent(in) :: t !< the timestep


  !H = I
  hx = x

end subroutine H

!> subroutine to take an observation vector y and return x \f$= H^T(y)\f$
!> in full state space.
!!
!! Given \f$y\f$ compute \f$x=H^T(y)\f$
subroutine HT(obsDim,nrhs,y,x,t)
!  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: obsDim !< the dimension of the observations
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(obsDim,nrhs), intent(in) :: y !< the
  !!input vectors in observation space
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: x !< the
  !!resulting vector in state space where x \f$= H^Ty\f$
  integer, intent(in) :: t !< the timestep

  !H = I
  x = y
  

end subroutine HT

!> subroutine to compute the distance between the variable
!> in the state vector and the variable in the observations
!>
!> Compute \f$\mathrm{dist}(x(xp),y(yp))\f$
subroutine dist_st_ob(xp,yp,dis,t)
  use sizes
  implicit none
  integer, intent(in) :: xp !<the index in the state vector
  integer, intent(in) :: yp !<the index in the observation vector
  real(kind=kind(1.0d0)), intent(out) :: dis !<the distance between
                                             !!x(xp) and y(yp)
  integer, intent(in) :: t  !<the current time index for observations
  
  dis = 0.0d0
  
end subroutine dist_st_ob

!> subroutine to take a full state vector x and return Qx
!> in state space.
!!
!! Given \f$x\f$ compute \f$b^{1/2}x\f$
subroutine Bhalf(nrhs,x,Qx)

  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x !< the
  !!input vector
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: Qx !< the
  !!resulting vector where Qx \f$= B^{1/2}x\f$
  real(kind=rk), dimension(state_dim,nrhs) :: temp
  
  real(kind=rk), dimension(3,3) :: bhalfdat

  bhalfdat(1:3,1) = (/ 2.11314166d-1, 1.33392833d-1,-5.82026139d-3/)
  bhalfdat(1:3,2) = (/ 1.33392833d-1, 3.79075017d-1,-2.30954705d-4/)
  bhalfdat(1:3,3) = (/-5.82026139d-3,-2.30954705d-4, 4.26705685d-1/)

  Qx = matmul(bhalfdat,x)

  
end subroutine Bhalf


!> Subroutine to read observation from a file              
!! \n              
!> @param[out] y The observation              
!! @param[in] t the current timestep              
subroutine get_observation_data(y,t)

  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  integer, intent(in) :: t
  real(kind=rk), dimension(obs_dim), intent(out) :: y


  !This is set up tp call the routine written which will              
  !work to do twin experiments. If you want to use your own  
  !observations you should implement your own method of reading      
  !in the observations              
  call default_get_observation_data(y,t)
end subroutine get_observation_data



!> subroutine to take a full state vector x and return B^{-1}x
!> in state space.
!!
!! Given \f$x\f$ compute \f$B^{-1}x\f$
subroutine solve_b(nrhs,x,Qx)

  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x !< the
  !!input vector
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: Qx !< the
  !!resulting vector where Qx \f$= B^{1/2}x\f$
  real(kind=rk), dimension(state_dim,nrhs) :: temp
  
  real(kind=rk), dimension(3,3) :: binvdat

  binvdat(1:3,1) = (/ 41.639551051239188,-20.302009390334526, 0.748818628880694/)
  binvdat(1:3,2) = (/-20.302009390334526, 16.091028535138122,-0.332370273316068/)
  binvdat(1:3,3) = (/  0.748818628880694, -0.332370273316068, 5.504775516290841/)

  Qx = matmul(binvdat,x)

  
end subroutine Solve_b
