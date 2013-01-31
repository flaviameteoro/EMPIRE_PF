subroutine InitiateModel (psi)
!-------------------------------------------------------------------------- 
! This is the multi-layer primitive equation version of the double gyre problem
!
!        domain:     10e  to  35e
!                    45s  to  30s
!
!        resolution: 0.1 degree (= approx. 10 km)
!
!--------------------------------------------------------------------------

use Sizes
use ErrorsAndVariances


 IMPLICIT NONE

  real, dimension(:, :)     :: psi
  integer  STAT, ntime, alloc_error,i,j
  real  t

  !  Initialize psi
  write(*,*)'initialize psi'

  open (50, file='Data/PE.dat', form = 'unformatted')
  read (50) psi
  close (50)

  5 format(e16.10)
  
  write(*,*)'end initialise'
  
end subroutine InitiateModel




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine IntegrateModel (nbegin, nend, psi,time,weight)
!Mel-08|08|11-changed to check calculation of weightQ 
subroutine IntegrateModel (nbegin, nend, psi,ttime,weight,rank)
!--------------------------------------------------------------------------
!  Integrate the system from time-step nbegin to time-step nend 
!--------------------------------------------------------------------------

use Sizes
use ErrorsAndVariances
use EnsembleOfFields
use ObservationField
use Output
use FullQ
use TimeInfo


IMPLICIT NONE

Interface

  subroutine NormalRandomNumbers1D (mean, stdev, phi)
    implicit none
    real, intent(in) :: mean, stdev
    real, intent(out), dimension(:) :: phi
  end subroutine NormalRandomNumbers1D
!  subroutine model(u,v,e,k)
!    use Sizes
!    implicit none
!    integer k
!    real, intent(inout), dimension(:,:,:) :: u,v,e
!  end subroutine model
  subroutine testQ(inVector)
    use Sizes
    use ErrorsAndVariances
    use fullQ
    IMPLICIT NONE
    real, dimension(:), intent(in) :: inVector
  end subroutine testQ

End Interface

integer :: nbegin, nend, ttime
real, dimension(:, :) :: psi
real, dimension(3) :: emax,umax,vmax
real ::  t, temp, thulp, total
character filename*20, timeString*6, rankNo*3

! Variables:
real, allocatable, dimension(:) :: qobs
real, allocatable, dimension(:,:,:)   :: u,v,e

!Mel-08|08|11-added to check calculation of weightQ
real :: weightQ=0.

!Mel-14|11|11-added to allow fullQ eror in weights
real, allocatable, dimension(:) :: exVector !vector to store (x^n-f(x^{n-1})) and to calculate Q^{-1}(x^n-f(x^{n-1}))
real :: pWeight, qWeight

!Mel-14|11|11-added to allow mutiplcation of nudge by full Q
real, allocatable, dimension(:) :: qhelp    !vector to store temp*pseudoH*(qobs-q)
real, allocatable, dimension(:) :: qcorr    !vector to store Q*qhelp

!Mel-20|12|11-added to allow sqrtQ correlated random error
real, allocatable, dimension(:) :: randomVec  !vector to store uncorrelated random error
real, allocatable, dimension(:) :: corrRandom !vector to store sqrtQ correlated random error

!Mel-14|08|12-added to test control variable transform
real, allocatable, dimension(:) :: psiTest


! Model setup:
real    :: b                                   ! Extrapolation parameter
real  :: weight,factor

integer          nstep,i,j,r, status, maxl,alloc_error,nex,rank,k,z,index

write(*,*)'start integration'
write(rankNo,'(i3.3)') rank

b=1.0
nex=4

write(*,*)'time is ',nbegin,' to ',nend

allocate(u(nxn,nyn,3) , STAT=alloc_error)
allocate(v(nxn,nyn,3) , STAT=alloc_error)
allocate(e(nxn,nyn,3) , STAT=alloc_error)

allocate(qobs(3*nxn*nyn) , STAT=alloc_error)

!Mel-14|11|11-added to allow fullQ error in weights
allocate(exVector(3*nxn*nyn),STAT=alloc_error)

!Mel-14|11|11-added to allow mutiplcation of nudge by full Q
allocate(qhelp(3*nxn*nyn), STAT=alloc_error)
allocate(qcorr(3*nxn*nyn), STAT=alloc_error)

!Mel-20|12|11-added for sqrt Q correlated error
allocate(randomVec(nxn*nyn), STAT=alloc_error)
allocate(corrRandom(3*nxn*nyn), STAT=alloc_error)

!Mel-14|08|12-added to test control variable transform
allocate(psiTest(3*nxn*nyn), STAT=alloc_error)

IF (alloc_error /= 0)  PRINT *, "Couldn't allocate space for matrices in enssmm.f90"

! Read data from files
write(*,*)'reading data...'
write (timeString, '(i6.6)') ttime
filename = 'Data/datan' // timeString
open (50, file=filename, form = 'unformatted') 
read(50)qobs
close(50)

do i=1,nxn
  do j=1,nyn
      u(i,j,2) = psi((i-1)*nyn+j,1)
      v(i,j,2) = psi(nxn*nyn+(i-1)*nyn+j,1)
      e(i,j,2) = psi(2*nxn*nyn+(i-1)*nyn+j,1)
  enddo
enddo

weightQ = 0.

!Mel-20|12|11-changed to allow sqrt Q correlated error
if (nbegin .eq. 0) then
    call NormalRandomNumbers1D(0.0,1.0,randomVec)
    call MultiplySqrtQ(randomVec,corrRandom)
    do i=1,nxn
      do j=1,nyn
        u(i,j,2)=u(i,j,2) +  1./qrel*corrRandom((i-1)*nyn+j) !1/qrel required to negate scaling of Q
        v(i,j,2)=v(i,j,2) +  1./qrel*corrRandom(nxn*nyn+(i-1)*nyn+j)
        e(i,j,2)=e(i,j,2) +  1./qrel*corrRandom(2*nxn*nyn+(i-1)*nyn+j) 
      enddo
    enddo
endif


write(*,*)'Start integration now'

do nstep = nbegin+1, nend

   relaxStrength(:,index)=0;
   !Mel-23|09|11-for plotting 3D movement
   if (indTraj .eq. 1) then
       open(35,file='output/close'//rankNo,position="append")
       write(35,*) e(114,128,1),e(128,128,1),e(128,142,1)
       close(35)
    
       open(36,file='output/apart'//rankNo,position="append")
       write(36,*) e(1,1,1),e(86,86,1),e(173,173,1)
       close(36)
  
       open(37,file='output/closeObs'//rankNo,position="append")
       write(37,*) e(119,127,1),e(127,127,1),e(127,135,1)
       close(37)
    
       open(38,file='output/apartObs'//rankNo,position="append")
       write(38,*) e(0,0,1),e(87,87,1),e(171,171,1)
       close(38) 

       open(38,file='output/trajObs2'//rankNo,position="append")
       write(38,*) e(127,63,1)
       close(38)

       open(38,file='output/trajUnObs2'//rankNo,position="append")
       write(38,*) e(128,192,1)
       close(38)
   endif

   !Mel-27|09|11-for plotting full Movie
   if (fullMovie .eq. 1) then
       if (mod(nstep,6) .eq. 0) then
           open(39,file='output/full'//rankNo,position="append")
           do i=1,nxn-1
             do j=1,nyn-1
               write(39,*) e(i,j,2)
             enddo
           enddo
           close(39)
        endif
    endif

    t=nstep
    if(nstep .eq. nbegin+1)then
!SSW don't call model
!      call model(u,v,e,1)
    else
!      call model(u,v,e,2)
    endif
 
    temp = nstep-nbegin
    thulp=(nend-nbegin)
    temp = nudgeFac*temp/thulp    !qd*qrel*temp/thulp*0.125 !Added to return to 0.001 overall
    !write(*,*)'temp is ',temp,nstep
    do i=1,nxn
      do j=1,nyn
        qhelp((i-1)*nyn+j)=temp*qd((i-1)*nyn+j)*pseudoH((i-1)*nyn+j)*(qobs((i-1)*nyn+j)-u(i,j,2))
        qhelp(nxn*nyn+(i-1)*nyn+j)=temp*qd(nxn*nyn+(i-1)*nyn+j)*pseudoH(nxn*nyn+(i-1)*nyn+j)*(qobs(nxn*nyn+(i-1)*nyn+j)-v(i,j,2))
        qhelp(2*nxn*nyn+(i-1)*nyn+j)=temp*qd(2*nxn*nyn+(i-1)*nyn+j)*pseudoH(2*nxn*nyn+(i-1)*nyn+j)*(qobs(2*nxn*nyn+(i-1)*nyn+j)-e(i,j,2))
      enddo
    enddo

    !Mel-14|11|11-added to allow multiplication by fullQ
    qcorr(:)=0.0
    call MultiplyQ(qhelp,qcorr)

    !Mel-20|12|11-changed to allow random error correlated by sqrt Q
    call NormalRandomNumbers1D(0.0,1.0,randomVec)
    call MultiplySqrtQ(randomVec,corrRandom)

    do i=1,nxn
      do j=1,nyn
        !Mel-14|11|11-changed to allow multiplication by fullQ
         u(i,j,2) = u(i,j,2) + qcorr((i-1)*nyn+j) + corrRandom((i-1)*nyn+j)
         v(i,j,2) = v(i,j,2) + qcorr(nxn*nyn+(i-1)*nyn+j) + corrRandom(nxn*nyn+(i-1)*nyn+j)
         e(i,j,2) = e(i,j,2) + qcorr(2*nxn*nyn+(i-1)*nyn+j) + corrRandom(2*nxn*nyn+(i-1)*nyn+j) 
      enddo
    enddo

! PJ- Not sure what happens here...
!          do i=1,nxn
!              do j=1,nyn
!                  relaxStrength(1,index) = relaxStrength(1,index) + corrRandom((i-1)*nxn+j)*corrRandom((i-1)*nxn+j)
!                  relaxStrength(2,index) = relaxStrength(2,index) + qcorr((i-1)*nxn+j)*qcorr((i-1)*nxn)/(qd*qd)
!                  relaxStrength(3,index) = relaxStrength(3,index) + (q(i-1,j-1)-qstart(i-1,j-1))*(q(i-1,j-1)-qstart(i-1,j-1))
!              enddo
!          enddo
!          relaxStrength(1,index) = sqrt(relaxStrength(1,index))
!          relaxStrength(2,index) = sqrt(relaxStrength(2,index))
!          relaxStrength(3,index) = sqrt(relaxStrength(3,index))
!       qcorr(:)=0.0
!       corrRandom(:)=0.0

   write(*,*)'max e is ',maxval(e(1:,:,2)),nstep

   exVector = qcorr + corrRandom

   call solveQ(exVector,pWeight)
   call solveQ(corrRandom,qWeight)

   !write(*,*) 'pWeight: ',  pWeight
   !write(*,*) 'qWeight: ',  qWeight

   weight = weight + 0.5*pWeight - 0.5*qWeight

   !write(*,*) 'weight is:  ',weight,'t= ',nstep

enddo

do i=1,nxn
  do j=1,nyn
      psi((i-1)*nyn+j,1) = u(i,j,2)
      psi(nxn*nyn+(i-1)*nyn+j,1) = v(i,j,2)
      psi(2*nxn*nyn+(i-1)*nyn+j,1) = e(i,j,2)
    
      psiTest((i-1)*nyn+j) = u(i,j,2)
      psiTest(nxn*nyn+(i-1)*nyn+j) = v(i,j,2)
      psiTest(2*nxn*nyn+(i-1)*nyn+j) = e(i,j,2)
   enddo
enddo

!if (nend .eq. 1920) then
!   call testQ(psiTest)
!endif

5 format(e16.10)

write(*,*)'Final weight for ',rank,' is ',weight

if(allocated(u))  deallocate(u)
if(allocated(v))  deallocate(v)
if(allocated(e))  deallocate(e)
if(allocated(qobs))  deallocate(qobs)
if(allocated(exVector))  deallocate(exVector)
if(allocated(qhelp))        deallocate(qhelp)
if(allocated(qcorr))        deallocate(qcorr)
if(allocated(randomVec))  deallocate(randomVec)
if(allocated(corrRandom))  deallocate(corrRandom)

end subroutine IntegrateModel


Subroutine testQ(inVector)

  use Sizes
  use ErrorsAndVariances
  use fullQ

  IMPLICIT NONE

     real, dimension(:), intent(in) :: inVector

     real, dimension(3*nxn*nyn) :: exVector, psiTest
     real, allocatable, dimension(:) :: psiTest1, psiTest2, psiTransform, psiSin, psiSin2
     integer :: alloc_error, i
     real :: total

     allocate(psiTest1(3*nxn*nyn), STAT=alloc_error)
     allocate(psiTest2(3*nxn*nyn), STAT=alloc_error)
     allocate(psiTransform(nxn*nyn),STAT=alloc_error)
     allocate(psiSin(nxn*nyn),STAT=alloc_error)
     allocate(psiSin2(nxn*nyn),STAT=alloc_error)

     psiTest=inVector
     exVector=0.

     !-----------------------test of SolveQ----------------------------------
     !write(*,*)'SolveQ test'
     !open(73,file='output/test1.dat',form='formatted')
     !write(73,5) psiTest
     !close(73)
        
     !call SolveQ(psiTest,exVector)
     !call MultiplyQ(psiTest,exVector)
     !call uAdjoint(psiTest,psiSin)
     !call uTransform(psiSin,exVector)
        
     !write(*,*)'SolveQ test'
     !open(73,file='output/test2.dat',form='formatted')
     !write(73,5) exVector
     !close(73)


     !-----------------------test of SolveQR----------------------------------
     !write(*,*)'SolveQR test'
     !open(73,file='output/test1.dat',form='formatted')
     !write(73,5) psiTest
     !close(73)
        
     !call SolveQR(psiTest,exVector)
     !call MultiplyQ(exVector,psiTest)

     !psiTest(:)=psiTest(:)+qd(:)*qd(:)*exVector(:)
        
     !write(*,*)'SolveQR test'
     !open(73,file='output/test2.dat',form='formatted')
     !write(73,5) psiTest
     !close(73)


     !-----------------------inverse test--------------------------------------
     !call uInverse(psiTest,psiTransform)
     !call uHatInverse(psiTransform,psiSin)
     !call uHatTransform(psiSin,psiTransform)
     !call uTransform(psiTransform,psiTest)     

     !write(*,*)'Inverse test'
     !open(73,file='output/test1.dat',form='formatted')
     !write(73,5) psiTest
     !close(73)

     !call uInverse(psiTest,psiTransform)
     !call uHatInverse(psiTransform,psiSin)
     !call uAdjoint(psiTest,psiTransform)
     !call uHatAdjoint(psiTransform,psiSin)
     !call uHatTransform(psiSin,psiTransform)
     !call uTransform(psiTransform,exVector)

     !open(74,file='output/test2.dat',form='formatted')
     !write(74,5) exVector
     !close(74)


     !----------------------adjoint u test------------------------------------
     write(*,*)'U^TU adjoint test'
     call uInverse(psiTest,psiTransform)
     call uTransform(psiTransform,psiTest1)

     total=0.
     do i=1,3*nxn*nyn
        total=total+psiTest1(i)*psiTest1(i)
     enddo
     write(*,*) 'true uhat version: ',total
   
     call uAdjoint(psiTest1,psiSin)
     
     total=0.
     do i=1,nxn*nyn
        total=total+psiTransform(i)*psiSin(i)
     enddo
     write(*,*) 'adjoint uhat version: ',total

     !----------------------2nd adjoint u test------------------------------------
     write(*,*)'UU^T adjoint test'
     call uAdjoint(psiTest,psiSin)

     total=0.
     do i=1,nxn*nyn
        total=total+psiSin(i)*psiSin(i)
     enddo
     write(*,*) 'true uhat version: ',total
   
     call uTransform(psiSin,psiTest1)
     
     total=0.
     do i=1,3*nxn*nyn
        total=total+psiTest(i)*psiTest1(i)
     enddo
     write(*,*) 'adjoint uhat version: ',total

    
     !----------------------Q^{1/2} test-------------------------------------------
     
     !psiTest=0.
     !psiTest(9900)=1.
     
     !call uAdjoint(psiTest,psiTest1)
     !call uHatAdjoint(psiTest1,psiTest2)
     
     !do i=1,nxn*nyn
     !   psiTest2(i)=qv(i)*psiTest2(i)
     !enddo
     
     !call uHatTransform(psiTest2,psiTest1)
     !call uTransform(psiTest1,psiTest)
     
     !open(72,file='output/testQ.dat',form='formatted')
     !write(72,5) psiTest
     !close(72)

     5 format(e16.10)

End Subroutine testQ
