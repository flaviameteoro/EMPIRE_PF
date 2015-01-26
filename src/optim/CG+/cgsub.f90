!
!       --------------------------------------------------------------------
!       Main program for running the conjugate gradient methods described in 
!       the paper:
!
!       Gilbert, J.C. and Nocedal, J. (1992). "Global Convergence Properties 
!       of Conjugate Gradient Methods", SIAM Journal on Optimization, Vol. 2,
!       pp. 21-42.
!
!       A web-based Server which solves unconstrained nonlinear optimization
!       problems using this Conjugate Gradient code can be found at:
!
!       http://www-neos.mcs.anl.gov/neos/solvers/UCO:CGPLUS/
!
!       Written by G. Liu, J. Nocedal and R. Waltz
!       October 1998
!
!--------------------------------------------------------------------
! modified to be a callable subroutine by Philip A Browne Jan 2015
subroutine subroutine_cg(method,n,x)
!
! Change the maximum size of the problem dimension here
!
implicit none
integer, parameter :: rk = kind(1.0d0)

integer, intent(in) :: method
integer, intent(in) :: n
real(kind=rk), dimension(n), intent(inout) :: x

real(kind=rk), dimension(n) :: g,d,gold,w
real(kind=rk) :: f,eps,tlev
real(kind=rk) :: time1,time2,tottime
real(kind=rk), parameter :: one=1.0d0
logical :: finish
integer :: iprint(2),iflag,icall,i!mp,lp,i
integer :: iter,nfun,irest




FINISH= .FALSE.
!
! Read problem input information
!
irest = 0 !no restarts
!irest = 1 !restarts every n steps  

!     IPRINT =  FREQUENCY AND TYPE OF PRINTING
!               IPRINT(1) < 0 : NO OUTPUT IS GENERATED
!               IPRINT(1) = 0 : OUTPUT ONLY AT FIRST AND LAST ITERATION
!               IPRINT(1) > 0 : OUTPUT EVERY IPRINT(1) ITERATIONS
!               IPRINT(2)     : SPECIFIES THE TYPE OF OUTPUT GENERATED;
!                               THE LARGER THE VALUE (BETWEEN 0 AND 3),
!                               THE MORE INFORMATION
!               IPRINT(2) = 0 : NO ADDITIONAL INFORMATION PRINTED
! 		IPRINT(2) = 1 : INITIAL X AND GRADIENT VECTORS PRINTED
!		IPRINT(2) = 2 : X VECTOR PRINTED EVERY ITERATION
!		IPRINT(2) = 3 : X VECTOR AND GRADIENT VECTOR PRINTED 
!				EVERY ITERATION 
iprint(1) = 1 
iprint(2) = 0  
!
! Check for correct dimension value n
!
if (n .lt. 0) then
   iflag = -3
   write(*,850)
   go to 50
end if

!
! Print parameters
!
if (iprint(1) .ge. 0) then
   write (*,820)
   write (*,840) n, method, irest
end if

ICALL=0
!
! This is the convergence constant 
!
EPS= 1.0D-5

! IFLAG=0 indicates an initial entry to program

IFLAG=0

!
! Begin counting CPU time. 
! (Note: This function may not work on all operating systems.)
!


20 CONTINUE
!
! Calculate the function and gradient values here
!
! Rosenbrock test function
call fcn(n,x,f,g)

30 CONTINUE
!
! Call the main optimization code
!
CALL CGFAM(N,X,F,G,D,GOLD,IPRINT,EPS,W,IFLAG,IREST,METHOD,FINISH )
!
!       IFLAG=
!              0 : successful termination
!              1 : return to evaluate F and G
!              2 : return with a new iterate, try termination test
!             -i : error
!
IF(IFLAG.LE.0.OR.ICALL.GT.10000) GO TO 50
IF(IFLAG.EQ.1) THEN
   ICALL=ICALL + 1	   
   GO TO 20
ENDIF
IF(IFLAG.EQ.2) THEN
!
! Termination Test.  The user may replace it by some other test. However, 
! the parameter 'FINISH' must be set to 'TRUE' when the test is satisfied.
!
   TLEV= EPS*(ONE + DABS(F))
   I=0
40 I=I+1
   IF(I.GT.N) THEN
      FINISH = .TRUE.
      GO TO 30
   ENDIF
   IF(DABS(G(I)).GT.TLEV) THEN
      GO TO 30
   ELSE
      GO TO 40
   ENDIF
   
ENDIF
        
50 continue
        
!
! Code has terminated; print final results
!
if (iprint(1).ge.0.and.iflag.ge.0) then
   write (*,890) f
end if

!
! Formatting
!
800 format (12x, i3)
820 format (//, ' Conjugate Gradient Minimization Routine', /)
840 format (/, ' n	  =', i6, /,' method   =', i6,/,' irest    =', i6,/)
850 format (/,'  Error: negative N value'/)
890 format (/,' f(x*) =', 1pd16.8)


end subroutine subroutine_cg





















