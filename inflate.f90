subroutine inflate

  ! Simple code to illustrate row entry to hsl_ea20
!  use pf_control
  use HSL_EA20_double
  use sizes
  use Qdata
  implicit none

  ! Derived types
  type (ea20_control) :: cntl
  type (ea20_info)    :: info
  type (ea20_reverse) :: rev

  ! Parameters

  integer, parameter :: wp = kind(0.0d0)    
  integer :: ido,i
  real(kind=kind(1.0D0)), dimension(state_dim) :: y
  real(kind=kind(1.0D0)), dimension(state_dim) :: x
  real(kind=kind(1.0D0)), allocatable  :: w(:,:)
  real(kind=kind(1.0D0))               :: s
  character(20) :: filename
  logical :: err
  
  call NormalRandomNumbers1d(0.0d0,1.0d0,state_dim,x)

  err = .true.

  do
     if(.not. err) exit

     !! set u
     y = x 
     
     !! set data
     s = 0.5d0
     
     !! set cntl
     cntl%d     = 3         !! delay
     cntl%tol   = 1.d-2     !! convergece tolerance
     cntl%maxit = 20        !! max number iteration
     
     cntl%diagnostics_level = 1 !! full error check
     
     ido = -1
     
     do while ( ido .ne. 0 .and. info%flag == 0)
        
        call EA20(state_dim,y,s,ido,w,cntl,info,rev)
        
        select case ( ido )
           
        case (0)
           exit
           
        case (1) !! Matrix-Vector product w_out = A w(:,1)
           call Q(w(:,1),w(:,2))
           
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
        write(*,'(a,i5)') 'error code  = ',info%flag
        !! print the final solution
        !     write(*,'(a,i5)') 'number of iterations  = ',info%iter
        !     write(*,'(a,1x,1pe14.2)') 'estimated final error = ',info%error
        !     write(*,'(a)') '    i       X(i)'
        !     write(*,'(i5,1x,1pe14.3)') (i,y(i), i=1,state_dim)
        err = .false.
        print*,'AHA! WE HAVE AN SPD Q MATRIX' 
     else
        !     write(*,'(i5,1x,1pe14.3)') (i,x(i), i=1,state_dim) 
        print*,'EA20 broke:',info%flag
!!$     write(filename,'(A,i0)') 'data',pf%timestep
!!$     open(20,file=filename,action='write')
!!$     do ido = 1,state_dim
!!$        write(20,*) x(ido)
!!$     end do
!!$     close(20)
        err = .true.
        !now do some inflation
        print*,'INFLATING diag(Q) BY 0.5'
        do i = 1,Qne
           if(Qrow(i) .eq. Qcol(i)) Qval(i) = 0.5d0+Qval(i)
        end do
        
        !     stop
     end if
     
  end do
  

end subroutine inflate
