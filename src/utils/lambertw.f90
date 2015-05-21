!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-05-21 13:35:03 pbrowne>
!!!
!!!    Subroutine to implement the lambertw function
!!!    Copyright (C) 2015 Mengbin Zhu 
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
!!!    Email: zhumengbin @ gmail.com  
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> subroutine to implement the lambertw function
!! see @ref https://en.wikipedia.org/wiki/Lambert_W_function
SUBROUTINE LAMBERTW(K,X,W)

!######################################################################
!Lambert W Function Numerical Solver
!In mathematics, the Lambert W function, also called the omega function 
!or product logarithm, is a set of functions, namely the branches of 
!the inverse relation of the function z = f(W) = We^W 
!where e^W is the exponential function and W is any complex number.
!
!University of Reading, 
!Dept. of Meteorology, 
!Mengbin Zhu 14/04/2015
!######################################################################

IMPLICIT NONE

INTEGER, INTENT(IN)  :: K
REAL(KIND=KIND(1.0d0)),    INTENT(IN)  :: X
REAL(KIND=KIND(1.0d0)),    INTENT(OUT) :: W

INTEGER              :: I,J,Q,R,Q2,Q3

REAL(KIND=KIND(1.0d0)),PARAMETER       :: EPS=4.0E-16
REAL(KIND=KIND(1.0d0)),PARAMETER       :: EM1=0.3678794411714423215955237701614608

REAL(KIND=KIND(1.0d0))                 :: PP,EE,TT2,L1,L2

IF(K==0) THEN     !Lambert W-Function 0 Branch
    IF(X<-EM1) THEN
        PRINT*,'Lambert W-Function: Bad Argument! STOP!'
        RETURN
    END IF
    IF(X==0.0) THEN
        W=0.0
        RETURN
    END IF
    IF(X<(-EM1+1E-4)) THEN
        Q=X+EM1
        R=SQRT(REAL(Q))
        Q2=Q*Q
        Q3=Q2*Q
        W=-1.0 &
        & +2.331643981597124203363536062168*R    &
        & -1.812187885639363490240191647568*Q    &
        & +1.936631114492359755363277457668*R*Q  &
        & -2.353551201881614516821543561516*Q2   &
        & +3.066858901050631912893148922704*R*Q2 &
        & -4.175335600258177138854984177460*Q3   &
        & +5.858023729874774148815053846119*R*Q3 &
        & -8.401032217523977370984161688514*Q3*Q
        RETURN
    END IF
    IF(X<1.0) THEN
        PP=SQRT(2.0*(2.7182818284590452353602874713526625*X+1.0))
        W=-1.0+PP*(1.0+PP*(-0.333333333333333333333+PP*0.152777777777777777777777))
    ELSE
        W=LOG(X)
    END IF
    !IF(X>3.0) THEN
    !    W=W-LOG(W)
    !END IF

    TT2=1.0
    I=0
    DO WHILE (ABS(TT2)>(EPS*(1.0+ABS(W))))
        !PRINT*,"Cannot go out of the cycle 0 Branch."
        I=I+1
        EE=EXP(W)
        TT2=W*EE-X
        PP=W+1.0
        TT2=TT2/(EE*PP-0.5*(PP+1.0)*TT2/PP)
        W=W-TT2
        !PRINT*,W
        IF(I==20) THEN
            W = W + 1.0E-6
        END IF
        IF(I==40) EXIT
    END DO

ELSE IF (K==-1) THEN  !Lambert W-Function -1 Branch.
    !IF(X<-EM1 .OR. X>0.0) THEN
    !    PRINT*,'Lambert W-Function -1 Branch: Bad Argument! STOP!'
    !    !STOP
    !    RETURN
    !END IF
    !PRINT*,"THE -1 BRANCH"
    IF(X==0.0) THEN
        W=0.0
    ELSE 
        !PRINT*,X
        PP=1.0
        IF(X<(-1e-6)) THEN
            !PRINT*,"IS THIS REAL?"
            PP=-SQRT(2.0*(2.7182818284590452353602874713526625*X+1.0))
            W=-1.0+PP*(1.0+PP*(-0.333333333333333333333+PP*0.152777777777777777777777))
        ELSE
            L1=LOG(-X)
            L2=LOG(-L1)
            W=L1-L2+L2/L1
        END IF
        IF(ABS(PP)<1E-4) THEN
            RETURN
        END IF
        TT2=1.0
        I=0
        DO WHILE(ABS(TT2)>EPS*(1.0+ABS(W)))
            !PRINT*,"Cannot go out of the cycle"
            I=I+1
            EE=EXP(W)
            TT2=W*EE-X
            PP=W+1.0
            TT2=TT2/(EE*PP-0.5*(PP+1.0)*TT2/PP)
            W=W-TT2
            !PRINT*,W
            IF(I==20) THEN
                W = W + 1.0E-6
            END IF
            IF(I==40) EXIT
        END DO
    ENDIF

ELSE
    PRINT*,"Lambert W-Function Branch Should Be -1 OR 0!"
    !STOP
    RETURN
END IF

END SUBROUTINE LAMBERTW
