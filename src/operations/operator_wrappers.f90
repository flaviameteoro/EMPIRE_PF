!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-11-23 11:07:26 pbrowne>
!!!
!!!    Collection of combinations of other subroutines
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

!> Subroutine to apply \f$K\f$ to a vector y in observation space
!! where \f$K := QH^T(HQH^T+R)^{-1}\f$
!! @param[in] y vector in observation space
!! @param[out] x vector in state space
subroutine K(y,x)
  !subroutine to apply the operator K to a vector y in obs space and return
  !the vector x in full state space.
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim,pf%count), intent(out) :: x
  real(kind=rk), dimension(obs_dim,pf%count), intent(in) :: y

  real(kind=rk), dimension(obs_dim,pf%count) :: v
  real(kind=rk), dimension(state_dim,pf%count) :: vv
  !  real(kind=rk) :: dnrm2
  integer :: i
  !  print*,'||y||_2 = ',dnrm2(obs_dim,y,1)

  do i = 1,pf%count
     call solve_hqht_plus_r(obs_dim,y(:,i),v(:,i),pf%timestep)
  end do

  !  print*,'||(HQHT+R)^(-1)y||_2 = ',dnrm2(obs_dim,v,1)
  call HT(obs_dim,pf%count,v,vv,pf%timestep)
  !  print*,'||HTv||_2 = ',dnrm2(state_dim,vv,1)
  call Q(pf%count,vv,x)
  !  print*,'||Qvv||_2 = ',dnrm2(state_dim,x,1)
  call flush(6)
end subroutine K




!!$subroutine B(y,x)
!!$use pf_control
!!$use sizes
!!$implicit none
!!$integer, parameter :: rk = kind(1.0D0)
!!$real(kind=rk), dimension(obs_dim), intent(in) :: y
!!$real(kind=rk), dimension(state_dim), intent(out) :: x
!!$real(kind=rk), dimension(obs_dim) :: R_1y
!!$real(kind=rk), dimension(state_dim) :: HtR_1y,QHtR_1y
!!$real(kind=rk) :: freetime,p,tau
!!$
!!$freetime = 0.6_rk
!!$
!!$tau = real(modulo(pf%timestep,pf%time_bwn_obs),rk)/real(pf&
!!$     &%time_bwn_obs,rk)
!!$
!!$if(tau .lt. freetime) then
!!$   x = 0.0_rk
!!$else
!!$   
!!$   call solve_r(y,R_1y)
!!$   
!!$   call HT(R_1y,HtR_1y)
!!$   
!!$   call Q(HtR_1y,QHtR_1y)
!!$
!!$   p = pf%nudgefac*(tau-freetime)/(1.0_rk-freetime)
!!$
!!$   x = p*QHtR_1y
!!$end if
!!$
!!$end subroutine B


!> subroutine to calculate nudging term and correlated random errors
!!efficiently
!! @param[in] y (obs_dim,pf\%count) vectors of innovations
!! \f$y-H(x^{n-1})\f$
!! @param[out] x (state_dim,pf\%count) vectors of
!! \f$\rho H^TR^{-1}[y-H(x^{n-1})]\f$
!! @param[out] QHtR_1y (state_dim,pf\%count) vectors of
!! \f$\rho QH^TR^{-1}[y-H(x^{n-1})]\f$ 
!! @param[in] normaln (state_dim,pf\%count) uncorrelated random vectors such that
!! normaln(:,i) \f$\sim \mathcal{N}(0,I)\f$
!! @param[out] betan (state_dim,pf\%count) correlated random vectors such that        
!! betan(:,i) \f$\sim \mathcal{N}(0,Q)\f$
subroutine Bprime(y,x,QHtR_1y,normaln,betan)
  !this is B but with separate Q multiplication
  use timestep_data
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(obs_dim,pf%count), intent(in) :: y
  real(kind=rk), dimension(state_dim,pf%count), intent(out) :: x
  real(kind=rk), dimension(obs_dim,pf%count) :: R_1y
  real(kind=rk), dimension(state_dim,pf%count) :: HtR_1y
  real(kind=rk), dimension(state_dim,pf%count), intent(out) :: QHtR_1y
  real(kind=rk), dimension(state_dim,pf%count), intent(in) :: normaln
  real(kind=rk), dimension(state_dim,pf%count), intent(out) :: betan
  real(kind=rk), dimension(state_dim,2*pf%count) :: temp1,temp2
  real(kind=rk) :: freetime,p,tau
  real(kind=rk) :: t
  real(kind=rk), dimension(7) :: ti
  logical, parameter :: time = .false.
  logical :: zero
  include 'mpif.h'
  if(time) t = mpi_wtime()
  freetime = 0.6_rk

  tau = real(TSData%tau,rk)/real(pf%time_bwn_obs,rk)

  call relaxation_profile(tau,p,zero)
  
  !this is when the nudging term is zero
  if(zero) then
     x = 0.0_rk
     if(time) ti(1:3) = mpi_wtime()
     QHtR_1y = 0.0_rk
     if(time) ti(4) = mpi_wtime()
     call Qhalf(pf%count,normaln,betan)
     if(time) ti(5:7) = mpi_wtime()
  else !this is when the nudging term is nonzero

     call solve_r(obs_dim,pf%count,y,R_1y,pf%timestep)
     if(time) ti(1) = mpi_wtime()
     call HT(obs_dim,pf%count,R_1y,HtR_1y,pf%timestep)
     if(time) ti(2) = mpi_wtime()
     
     !p = pf%nudgefac*(tau-freetime)/(1.0d0-freetime)

     x = p*HtR_1y
     if(time) ti(3) = mpi_wtime()
     temp1(:,1:pf%count) = x
     temp1(:,pf%count+1:2*pf%count) = normaln
     if(time) ti(4) = mpi_wtime()
     call Qhalf(2*pf%count,temp1,temp2)
     if(time) ti(5) = mpi_wtime()
     betan = temp2(:,pf%count+1:2*pf%count)
     !x = temp2(:,1:pf%count)
     if(time) ti(6) = mpi_wtime()
     call Qhalf(pf%count,temp2(:,1:pf%count),QHtR_1y)
     if(time) ti(7) = mpi_wtime()

  end if
  if(time) then
     ti(7) = ti(7) - ti(6)
     ti(6) = ti(6) - ti(5)
     ti(5) = ti(5) - ti(4)
     ti(4) = ti(4) - ti(3)
     ti(3) = ti(3) - ti(2)
     ti(2) = ti(2) - ti(1)
     ti(1) = ti(1) - t

     print*,'Bprime times =',ti
     print*,'_______'
  end if


end subroutine Bprime

!!$subroutine Qhalf(x,y)
!!$
!!$  ! Simple code to illustrate row entry to hsl_ea20
!!$!  use pf_control
!!$  use HSL_EA20_double
!!$  use sizes
!!$  implicit none
!!$
!!$  ! Derived types
!!$  type (ea20_control) :: cntl
!!$  type (ea20_info)    :: info
!!$  type (ea20_reverse) :: rev
!!$
!!$  ! Parameters
!!$
!!$  integer, parameter :: wp = kind(0.0d0)    
!!$  integer :: ido
!!$  real(kind=kind(1.0D0)), dimension(state_dim), intent(out) :: y
!!$  real(kind=kind(1.0D0)), dimension(state_dim), intent(in) :: x
!!$  real(kind=kind(1.0D0)), allocatable  :: w(:,:)
!!$  real(kind=kind(1.0D0))               :: s
!!$
!!$  !! set u
!!$  y = x 
!!$
!!$  !! set data
!!$  s = 0.5d0
!!$
!!$  !! set cntl
!!$  cntl%d     = 3         !! delay
!!$  cntl%tol   = 1.d-2     !! convergece tolerance
!!$  cntl%maxit = 20        !! max number iteration
!!$
!!$  cntl%diagnostics_level = 1 !! full error check
!!$
!!$  ido = -1
!!$
!!$  do while ( ido .ne. 0 .and. info%flag == 0)
!!$
!!$     call EA20(state_dim,y,s,ido,w,cntl,info,rev)
!!$
!!$     select case ( ido )
!!$
!!$     case (0)
!!$        exit
!!$
!!$     case (1) !! Matrix-Vector product w_out = A w(:,1)
!!$        call Q(w(:,1),w(:,2))
!!$
!!$!        if(maxval(abs(w(:,1)-w(:,2))) .gt. 1.0D-10) then
!!$!           print*,'shaft! ',maxval(abs(w(:,1)-w(:,2))) 
!!$!        end if
!!$        
!!$     case (2) !! Matrix-Vector product w(:,2) = M w(:,1)
!!$        w(:,2) = w(:,1)
!!$
!!$     case (3) !! Matrix-Vector product w_out = inv(M) w(:,1)
!!$        w(:,2) = w(:,1)
!!$
!!$     end select
!!$
!!$  end do
!!$
!!$
!!$  if (info%flag .ge. 0) then
!!$!     write(*,'(a,i5)') 'error code  = ',info%flag
!!$     !! print the final solution
!!$!     write(*,'(a,i5)') 'number of iterations  = ',info%iter
!!$!     write(*,'(a,1x,1pe14.2)') 'estimated final error = ',info%error
!!$!     write(*,'(a)') '    i       X(i)'
!!$!     write(*,'(i5,1x,1pe14.3)') (i,y(i), i=1,state_dim)
!!$
!!$  else
!!$!     write(*,'(i5,1x,1pe14.3)') (i,x(i), i=1,state_dim) 
!!$     print*,'EA20 broke: max x = ',maxval(x),' min x = ',minval(x)
!!$!     write(filename,'(A,i0)') 'data',pf%timestep
!!$!     open(20,file=filename,action='write')
!!$!     do ido = 1,state_dim
!!$!        write(20,*) x(ido)
!!$!     end do
!!$!     close(20)
!!$
!!$!     stop
!!$  end if
!!$
!!$
!!$end subroutine Qhalf



