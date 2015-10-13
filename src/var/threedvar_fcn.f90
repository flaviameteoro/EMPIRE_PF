!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-10-13 15:25:35 pbrowne>
!!!
!!!    subroutine to compute 3DVar objective function and gradient
!!!    Copyright (C) 2015  Philip A. Browne
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


!> subroutine to provide the objective function and gradient for 3dvar
!!
!!
!! Let \f$x\f$ be the state we wish to find using Var.
!!
!! The objective function considered is
!!
!! \f$J(x) = \frac{1}{2}(x-x_b)^TB^{-1}(x-x_b) +
!! \frac{1}{2}(y-H(x))^T R^{-1} (y - H(x)) )\f$
!!
!! where \f$x_b\f$ is a background guess, \f$B\f$ the background
!! error covariance matrix,
!!
!! \f$y\f$ are the observations, and \f$H\f$ the
!! corresponding observation operator with associated observation
!! error covariance matrix \f$R\f$.
!!
!! The gradient of the objective function can then be written
!!
!! \f$g = \nabla J(x) \approx B^{-1}(x-x_b) - H^TR^{-1}(y-H(x))\f$
!!
!! which is exact if \f$H\f$ is linear 
!>
!> NOTE: this will only currently work for EMPIRE VERSION 1 of 2.
!> @todo update 3dvar to work with EMPIRE VERSION 3!
subroutine threedvar_fcn(n, x, f, g )
  use sizes
  use timestep_data
  use comms
  use threedvar_data
  include 'mpif.h'
  integer, intent(in) :: n !< the dimension of the state
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(n), intent(in) :: x !< current guess
  real(kind=rk), dimension(n), intent(out) :: g !< gradient of
  !!objective function
  real(kind=rk), intent(out) :: f !< the objective function
  
  real(kind=rk), dimension(obs_dim) :: y,Hx
  real(kind=rk), dimension(n) :: gtemp
  real(kind=rk), dimension(n) :: xtemp
  real(kind=rk) :: ftemp,ftemp3
  integer :: mpi_err
  !start with observation component:
  ! Hx = H(x)
  call H(obs_dim,1,x,Hx,TSData%current_timestep)
  
  call get_observation_data(y,TSData%current_timestep)

  ! USE BLAS to perform
  ! y = y - Hx
  call daxpy(obs_dim,-1.0d0,hx,1,y,1)
  
  ! hx = R^{-1}(y-Hx)
  call solve_r(obs_dim,1,y,Hx,TSData%current_timestep)

  ! compute the obs component of the objective function
  ! f = 0.5 (y-Hx)^T R^{-1} (y-Hx)
  ftemp = 0.5_rk*sum(y*Hx)
  if(empire_version .eq. 3) then
     !need to perform the sum across all parts of the obs vector
     ftemp3 = ftemp
     call mpi_allreduce(ftemp3,ftemp,1,MPI_DOUBLE_PRECISION,MPI_SUM&
          &,pf_member_comm,mpi_err)
  end if
  
  ! gtemp = H^TR^{-1}(y-Hx)
  call HT(obs_dim,1,y,gtemp,TSData%current_timestep)

  

  !now lets do the background component
  
  !x - xb
  xtemp = x - xb
  
  !g = B^{1}(x-xb)
  call solve_b(1,xtemp,g)

  ! compute the background component of the objective function
  ! f = 0.5 (x-x_b)^T B^{-1} (x-x_b)
  f = 0.5_rk*sum(xtemp*g)
  if(empire_version .eq. 3) then
     !need to perform the sum across all parts of the obs vector
     ftemp3 = f
     call mpi_allreduce(ftemp3,f,1,MPI_DOUBLE_PRECISION,MPI_SUM&
          &,pf_member_comm,mpi_err)
  end if
  
  ! put the background and observation components of the objective
  ! function together
  f = f + ftemp

  ! put the background and observation components of the gradient
  ! together. Note the minus sign in the obs component.
  ! Use blas to perform g = g - gtemp
  call daxpy(state_dim,-1.0d0,gtemp,1,g,1)

end subroutine threedvar_fcn
