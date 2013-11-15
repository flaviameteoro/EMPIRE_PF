subroutine K(y,x)
  !subroutine to apply the operator K to a vector y in obs space and return
  !the vector x in full state space.
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim), intent(out) :: x
  real(kind=rk), dimension(obs_dim), intent(in) :: y

  real(kind=rk), dimension(obs_dim) :: v
  real(kind=rk), dimension(state_dim) :: vv
  real(kind=rk) :: dnrm2

  print*,'||y||_2 = ',dnrm2(obs_dim,y,1)
  call solve_hqht_plus_r(y,v)
  print*,'||(HQHT+R)^(-1)y||_2 = ',dnrm2(obs_dim,v,1)
  call HT(v,vv)
  print*,'||HTv||_2 = ',dnrm2(state_dim,vv,1)
  call Q(vv,x)
  print*,'||Qvv||_2 = ',dnrm2(state_dim,x,1)
  call flush(6)
end subroutine K

subroutine innerR_1(y,w)
  !subroutine to take an observation vector y and return w = y^T R^(-1) y
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim) :: v
  real(kind=rk), intent(out) :: w

  call solve_r(y,v)

  !this can defo be done better using BLAS PAB...
  w = sum(y*v)

end subroutine innerR_1

subroutine innerHQHt_plus_R_1(y,w)
  !subroutine to take an observation vector y and return w = y^T (HQH^T+R)^(-1) y
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim) :: v
  real(kind=rk), intent(out) :: w

  call solve_hqht_plus_r(y,v)

  !this can defo be done better using BLAS PAB...
  w = sum(y*v)

end subroutine innerHQHt_plus_R_1


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

subroutine Bprime(y,x,QHtR_1y)
!this is B but with separate Q multiplication
use pf_control
use sizes
implicit none
integer, parameter :: rk = kind(1.0D0)
real(kind=rk), dimension(obs_dim), intent(in) :: y
real(kind=rk), dimension(state_dim), intent(out) :: x
real(kind=rk), dimension(obs_dim) :: R_1y
real(kind=rk), dimension(state_dim) :: HtR_1y
real(kind=rk), dimension(state_dim), intent(out) :: QHtR_1y
real(kind=rk) :: freetime,p,tau,atmos,ocean

freetime = 0.6_rk
atmos = 2.0d0/3.0d0
ocean = 1.0d0/3.0d0

tau = real(modulo(pf%timestep,pf%time_bwn_obs),rk)/real(pf&
     &%time_bwn_obs,rk)

if(tau .le. atmos) then
   !this is the atmosphere section
   if(tau .le. freetime*atmos) then
      x = 0.0_rk
      QHtR_1y = 0.0_rk
   else
      
      call solve_r(y,R_1y)
      
      call HT(R_1y,HtR_1y)
      
      !comment out Qhalf to make this subroutine Bprime  
      !   call Qhalf(HtR_1y,QHtR_1y)
      
      p = pf%nudgefac*(tau-freetime*atmos)/(atmos-freetime*atmos)
      
      x = p*HtR_1y
      
      call Q(x,QHtR_1y)
   end if
else
   !this is the ocean section
   if(tau .le. ocean*freetime+atmos) then
      x = 0.0_rk
      QHtR_1y = 0.0_rk
   else
      
      call solve_r(y,R_1y)
      
      call HT(R_1y,HtR_1y)
      
      !comment out Qhalf to make this subroutine Bprime  
      !   call Qhalf(HtR_1y,QHtR_1y)
      
      p = pf%nudgefac*(tau-freetime*ocean-atmos)/(ocean-ocean*freetime)
      
      x = p*HtR_1y
      
      call Q(x,QHtR_1y)
   end if
end if


!print*,'Bprime = ',dnrm2(state_dim,x,1),dnrm2(state_dim,QHtR_1y,1)

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


subroutine rhalf(x,y)

  ! Simple code to illustrate row entry to hsl_ea20
  use pf_control
  use HSL_EA20_double
  use sizes
  implicit none

  ! Derived types
  type (ea20_control) :: cntl
  type (ea20_info)    :: info
  type (ea20_reverse) :: rev

  ! Parameters

  integer, parameter :: wp = kind(0.0d0)    
  integer :: ido
  real(kind=kind(1.0D0)), dimension(obs_dim), intent(out) :: y
  real(kind=kind(1.0D0)), dimension(obs_dim), intent(in) :: x
  real(kind=kind(1.0D0)), allocatable  :: w(:,:)
  real(kind=kind(1.0D0))               :: s
  character(20) :: filename

  !! set u
  y = x 

  !! set data
  s = 0.5d0

  !! set cntl
  cntl%d     = 3         !! delay
  cntl%tol   = 1.d-2     !! convergece tolerance
  cntl%maxit = 10        !! max number iteration

  cntl%diagnostics_level = 1 !! full error check

  ido = -1

  do while ( ido .ne. 0 .and. info%flag == 0)

     call EA20(obs_dim,y,s,ido,w,cntl,info,rev)

     select case ( ido )

     case (0)
        exit

     case (1) !! Matrix-Vector product w_out = A w(:,1)
        call R(w(:,1),w(:,2))

!        if(maxval(abs(w(:,1)-w(:,2))) .gt. 1.0D-10) then
!           print*,'shaft! ',maxval(abs(w(:,1)-w(:,2))) 
!        end if
        
     case (2) !! Matrix-Vector product w(:,2) = M w(:,1)
        w(:,2) = w(:,1)

     case (3) !! Matrix-Vector product w_out = inv(M) w(:,1)
        w(:,2) = w(:,1)

     end select

  end do


  if (info%flag .ge. 0) then
!     write(*,'(a,i5)') 'error code  = ',info%flag
     !! print the final solution
!     write(*,'(a,i5)') 'number of iterations  = ',info%iter
!     write(*,'(a,1x,1pe14.2)') 'estimated final error = ',info%error
!     write(*,'(a)') '    i       X(i)'
!     write(*,'(i5,1x,1pe14.3)') (i,y(i), i=1,state_dim)
  else
!     write(*,'(i5,1x,1pe14.3)') (i,x(i), i=1,state_dim) 
     print*,'EA20 broke: max x = ',maxval(x),' min x = ',minval(x)
     write(filename,'(A,i0)') 'data',pf%timestep
     open(20,file=filename,action='write')
     do ido = 1,obs_dim
        write(20,*) x(ido)
     end do
     close(20)
!     stop
  end if


end subroutine rhalf


