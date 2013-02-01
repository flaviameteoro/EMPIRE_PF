Subroutine Analyse (time)

! In this routine the analysis is calculated
! Requires - a file with data
!          - the ensemble of current states: psiGrand
! Result: an ensemble of new states psiGrand

  use Sizes               ! nxx, nyy,  nEns, nGrand       
  use Grandfield          ! psiGrand
  use ErrorsAndVariances
  use ObservationField
  use Output
  use FullQ

  IMPLICIT NONE

!!$  Interface
!!$    Subroutine NormalRandomNumbers1D (mean, stdev, phi)
!!$       implicit none
!!$       real, intent(in) :: mean, stdev
!!$       real, intent(out), dimension(:) :: phi
!!$    End Subroutine NormalRandomNumbers1D
!!$    Subroutine UniformRandomNumbers1D(minval, maxval, phi)
!!$       implicit none
!!$       real, intent(in)::minval,maxval
!!$       real, intent(out), dimension(:) :: phi
!!$    End Subroutine UniformRandomNumbers1D
!!$  End Interface

  integer :: time         ! time is measured in timesteps since the beginning
  integer :: i, j,  nData, iolen, status, k, z, talCount
  integer :: nmax, ii, maxq(1),n2,count,nmax1, n1, num, ifail,nnmax,nxrandom,mbest,n3
  integer :: nrobs(nxx,nyy)
  integer :: lEnKF(nGrand)
  integer :: ndim, nnmin, nn
  integer :: clevel,index
  integer, parameter :: rk = kind(1.0D0)
  
  character  filename*20, timeString*6
  
  real(kind=rk) :: dataMean, lpsiMeanMean
  real(kind=rk) :: alpha, gamma,bh
  real(kind=rk) :: r2, dist2,sb, d2, d22, erange, dmin
  real(kind=rk) :: mindist,var
  real(kind=rk) :: hulp1,hulp2,xnum(1) 
  real(kind=rk) :: xxmin, priorxmin
  real(kind=rk) :: MI
  real(kind=rk) :: ccmax,qinv,qninv,aaa,bbb,bmin,qq,dy2,dy1
  real(kind=rk) :: storeVal
  real(kind=rk) :: qWeight

  real(kind=rk), dimension(:),   Allocatable :: data  ! Size nData
  real(kind=rk), dimension(:),   Allocatable :: lpsiMean         ! Size nData
  real(kind=rk), dimension(:,:), Allocatable :: lpsi, dlpsi      ! Size nGrand * nObs
  real(kind=rk), dimension(:), Allocatable :: b
  real(kind=rk), dimension(:), Allocatable :: bnew
  real(kind=rk), dimension(:), Allocatable :: xmin
  integer, dimension(:), Allocatable ::  irank                
  real(kind=rk), dimension(:,:), Allocatable ::  psisuper             
  real(kind=rk), dimension(:,:,:), Allocatable ::  psinew                ! Size nGrand 
  real(kind=rk), dimension(:), Allocatable ::  cc 
  real(kind=rk), dimension(:), Allocatable ::  weights
  real(kind=rk), dimension(:,:), Allocatable ::  qobs
  real(kind=rk), dimension(:,:), Allocatable :: vortobs ! height of bottom.

  !Mel-25|08|11-added to allow calculation of weight distribution
  real(kind=rk), dimension(:,:), Allocatable :: barWeights

  !Mel-10|01|12-added to allow calculation of strength of aew movement
  real(kind=rk), dimension(:), allocatable :: aewStrength

  !Mel-07|06|12-added to allow comp C-minC to aew movement
  real(kind=rk),dimension(:,:),allocatable :: cDiff

  !Mel-05|09|11 - added for use with full Q vector
  real(kind=rk), dimension(:,:), Allocatable :: xVector !for storing nObs size vector (d-H(f(x^{n-1})))
  real(kind=rk), dimension(:), Allocatable :: exVector !for storing full vector (x^n-f(x^{n-1})) or H^T(HQH^T+R)^{-1}*xVector
  real(kind=rk), dimension(:,:), Allocatable :: QRxVector !for storing (HQH^T+R)^{-1}*xVector
  real(kind=rk), dimension(:), Allocatable :: hTemp  !for storing HQH^T(HQH^T+R)^{-1}*xVector
  real(kind=rk), dimension(:), Allocatable :: qTemp  !for storing QH^T(HQH^T+R)^{-1}*xVector
 
 !Mel-14|11|11-added for sampling from a mixture density
  real(kind=rk) :: eFac, uFac, nFac, help
  real(kind=rk), dimension(:), Allocatable :: randomVec
  real(kind=rk), dimension(:), Allocatable :: corrRandom
  real(kind=rk), dimension(:), Allocatable :: proposalWeight
  real(kind=rk), dimension(:,:), Allocatable :: perbWeights !for storing (H^T(HQH^T+R)^{-1}*exVector+corrRandom)
  real(kind=rk), dimension(:), allocatable :: psiVector
  real(kind=rk), dimension(:), allocatable :: psiSin

  write(*,*)'analyse', time

! Read data from files
  write(*,*)'reading data...'
  write (timeString, '(i6.6)') time
  filename = 'Data/datan' // timeString
  open (50, file=filename, form = 'unformatted') 

  nData=(nxx)*(nyy)

  !---------!
  clevel=10
  eFac=0.001/nGrand
  uFac=1.e-5
  nFac=1.e-5
  !---------!

  allocate( data(nObs), STAT = status)
  if (status /= 0) stop " Allocation of data failed"
  allocate( qobs(nxx,nyy), STAT = status)
  if (status /= 0) stop " Allocation of qobs failed"

  !Mel-25|08|11-added to allow calculation of weight distribution
  if (weightDiv .eq. 1) then
      allocate( barWeights(nGrand,3), STAT = status)
      if (status /= 0) stop " Allocation of barWeights failed"
  endif

  !Mel-10|01|12-added to allow calculation of strength of aew movement
  if (relaxDiv .eq.1) then
      allocate( aewStrength(nGrand), STAT = status)
      if (status /= 0) stop " Allocation of aewStrength failed"
      allocate( cDiff(4,(8*nGrand/10)), STAT = status)
      if (status /= 0) stop " Allocation of aew_c failed"
  endif
  
  !Mel-02|09|11 - added for use with full Q vector
  allocate( xVector(nGrand,nObs), STAT = status) !for storing nObs size vector (d-H(f(x^{n-1})))
  if (status /= 0) stop " Allocation of xVector failed"
  allocate( exVector(nData), STAT = status) !for storing full vector (x^n-f(x^{n-1})) or H^T(HQH^T+R)^{-1}*xVector
  if (status /= 0) stop " Allocation of exVector failed"
  allocate( QRxVector(nGrand,nObs), STAT = status) !for storing (HQH^T+R)^{-1}*xVector
  if (status /= 0) stop " Allocation of QRxVector failed"
  allocate( hTemp(nObs), STAT = status)  !for storing HQH^T(HQH^T+R)^{-1}*xVector
  if (status /= 0) stop " Allocation of hTemp failed"
  allocate( qTemp(nData), STAT = status) !for storing QH^T(HQH^T+R)^{-1}*xVector 
  if (status /= 0) stop " Allocation of qTemp failed"

  !Mel-20|12|11-added for use with mixture density
  allocate( randomVec(nxn*nyn), STAT = status) !for storing uncorrelated random variables (normal or uniform)
  if (status /= 0) stop " Allocation of randomVec failed"
  allocate( corrRandom(nData), STAT = status) !for storing correlated random variables (normal or pseudo uniform)
  if (status /= 0) stop " Allocation of corrRandom failed"
  allocate( proposalWeight(nGrand), STAT = status) !for storing proposal weights
  if (status /= 0) stop " Allocation of proposalWeight failed"
  allocate( perbWeights(nGrand,nxn*nyn),STAT = status) !for storing factor needed for calculating transition for perturbed weights
  if (status /= 0) stop " Allocation of perbWeight failed"
  allocate( psiVector(nxn*nyn),STAT = status) 
  if (status /= 0) stop " Allocation of psiVector failed"  
  allocate( psiSin(nxn*nyn),STAT = status) 
  if (status /= 0) stop " Allocation of psiVector failed"  

  read(50)qobs
  close(50)

  write(*,*)'data read. Fill arrays...'
  call Measure (qobs, data)
 
  print *, ' In analyse:', nData, ' datapoints are read in for time ', time

  write(*,*)'arrays filled. Start allocation of extra arrays'

! Allocate rest of matrices 

  allocate( lpsiMean(nObs), STAT = status)
  if (status /= 0) stop " Allocation of lpsiMean failed"
  allocate( lpsi(nGrand, nObs), STAT = status)
  if (status /= 0) stop " Allocation of lpsi failed"
  allocate( dlpsi(nGrand, nObs), STAT = status)
  if (status /= 0) stop " Allocation of dlpsi failed"
  allocate( b(nGrand), STAT = status)
  if (status /= 0) stop " Allocation of b failed"
  allocate( bnew(nGrand), STAT = status)
  if (status /= 0) stop " Allocation of b failed"
  allocate( irank(nGrand), STAT = status)
  if (status /= 0) stop " Allocation of irank failed"
  allocate( xmin(nGrand), STAT = status)
  if (status /= 0) stop " Allocation of xmin failed"
  allocate( psiSuper(nxx, nyy), STAT = status)
  if (status /= 0) stop " Allocation of psiSuper failed"
  allocate( psinew(nxx,nyy,nGrand), STAT = status)
  if (status /= 0) stop " Allocation of psinew failed"
  allocate( cc(nGrand), STAT = status)
  if (status /= 0) stop " Allocation of cc failed"
  allocate( weights(nGrand), STAT = status)
  if (status /= 0) stop " Allocation of weights failed"

! Measure the various ensemble members
  write(*,*)' start measurements '
  do n3 = 1, nGrand
     call Measure (psiGrand(1:nxx,:,n3), lpsi(n3,:))
  enddo

  write(*,*)' measurements complete. Start analysis'
! Subtract it from each ensemble member to get a deviation
  do n3 = 1, nGrand
    dlpsi(n3,:) = (data(:) - lpsi(n3,:)) * (data(:) - lpsi(n3,:))
  enddo

! Determine the mean of all ensemble members
  lpsiMean = sum(lpsi, dim=1) / nGrand

! Block with print statements for getting an idea of what is happening
! Mel-11|08|11-whole block changed to allow calculations with reduced observations
  dataMean = sum(data(:)) / nObs
  print *, " rms of data     = ", sqrt(sum ((data(:) - dataMean) * (data(:) - dataMean)) /nObs),  &
           ", mean of data     = ", dataMean 
  lpsiMeanMean = sum(lpsiMean) / nData
  print *, " rms of lpsiMean = ", sqrt(sum ((lpsiMean - lpsiMeanMean) * (lpsiMean - lpsiMeanMean)) /nData),&
           ", mean of lpsiMean = ", lpsiMeanMean   
  print *, " Quality of ensemble before analysis:"
  var=0.
  do n3 = 1, nGrand
    print *, " n3= ", n3,  &
             " rms(data-lpsi) = ", sqrt(sum((data(:)-lpsi(n3,:)) * (data(:)-lpsi(n3,:))) / nObs),  &
             " mean(data-lpsi) = ", sum(data-lpsi(n3,:)) / nObs      
    var = var + sum((data-lpsi(n3,:))*(data-lpsi(n3,:)))
  enddo
  write(*,*)' var prior is ',var

! Calculate weights of ensemble members

   !Mel-02|09|11-calculation of (HQH^T+R)^{-1}(d-H(f(x^{-1})))
   do n3=1,nGrand

       write(*,*) 'Solve (HQH^T+R)^{-1}(d-H(x_i^{n-1}))'
       xVector(n3,:)=(data(:)-lpsi(n3,:))
       call SolveQR(xVector(n3,:),QRxVector(n3,:))
       
       !Mel-07|09|11-test that inverse is correct
       write(*,*)'Checking result of HQH+R inverse calculation'
       hTemp(:)=0.0;
       call checkQRinverse(QRxVector(n3,:),hTemp)

       write(*,*) 'Check of SolveQR'
       write(*,*) xVector(n3,5), hTemp(5)
   enddo

   weights(:) = psiGrand(nxx+1,1,:)
   write(*,*)'weights before equal weights',weights

   do n3=1,nGrand
     b(n3)=weights(n3)
     do ii=1,nObs
       !Mel-02|09|11-changed to allow for calculation with full Q vector
       b(n3)=b(n3) + 0.5*xVector(n3,ii)*QRxVector(n3,ii)
     enddo
   enddo

   write(*,*)'weights b before equal weights',b
   call M01DAF(b,1,nGrand,'A',irank,ifail)
   cc=b
   call M01EAF(cc,1,nGrand,irank,ifail)
   write(*,*)'Ranked weights',cc
   ccmax=cc(clevel*nGrand/10)
   write(*,*)'ccmax is ',ccmax

   psinew(:,:,:)=psiGrand(1:nxx,:,:)    ! Note that psinew is now f(psi^{n-1})

   index=1
   do n3=1,nGrand

      !Mel-02|09|11-calculation of QH^T(HQH^T+R)^{-1}x and HQH^T(HQH^T+R)^{-1}x
      hTemp(:)=0.0
      call Expand(QRxVector(n3,:),exVector(:))
      call MultiplyQ(exVector,qTemp)
      call MeasureVec(qtemp,hTemp) 

      if (b(n3).lt. ccmax) then
        write(*,*)'high weights n is ',n3
        dy1=0
        dy2=0
        do ii=1,nObs
            dy1 = dy1 + xVector(n3,ii)*hTemp(ii)/(qdData(ii)*qdData(ii))
            dy2 = dy2 + dlpsi(n3,ii)/(qdData(ii)*qdData(ii))
        enddo
        write(*,*)'dy1, dy2 ', dy1, dy2
        write(*,*)'ccmax, -log_wi ', ccmax, weights(n3)
        aaa = 0.5 * dy1                          
        bbb = 0.5 * dy2 - ccmax + weights(n3)
        write(*,*)'aaa bbb',aaa,bbb

        !Mel-change to 1 + sqrt(1.-bbb/aaa + 0.000000001)
        alpha = 1 + sqrt(1.-bbb/aaa + 0.000000001)
        !alpha = 1 - sqrt(1.-bbb/aaa + 0.000000001)
        !call UniformRandomNumbers1D(0.,1.,xnum)
        !Positive alpha
        !if (xnum(1) .gt. 0.5) then
        !    alpha = 1 + sqrt(1.-bbb/aaa + 0.000000001)
        !Negative alpha
        !else
        !    alpha = 1 - sqrt(1.-bbb/aaa + 0.000000001)
        !endif
        write(*,*)'alpha', alpha

      else
        write(*,*)'low weights n is ',n3
        alpha=1.0

      endif

      do j=1,nyy
        do i=1,nxx
          psiGrand(i,j,n3) = psiGrand(i,j,n3) + alpha *qTemp(i+(j-1)*nxx)
        enddo
      enddo    

      !Value needed for calculating transition densities stored here
      call uAdjoint(exVector,psiVector)
      call uHatAdjoint(psiVector,psiSin)
      do i=1,nxn*nyn
         psiSin(i)=psiSin(i)*sqrtQv(i)
      enddo
      psiSin=qrel*qo*psiSin
      perbWeights(n3,:)=alpha*psiSin(:)

      !Mel-10|01|12-added to calculate the strength of the aew movement
     ! if (relaxDiv .eq. 1) then
     !    if (b(n3)<ccmax) then
     !       cDiff(4,index)=0.
     !       do i=1,nxx
     !          do j=1,nyy
     !             cDiff(4,index)=cDiff(4,index)+(psiGrand(i,j,n3)-psinew(i,j,n3))*(psiGrand(i,j,n3)-psinew(i,j,n3))
     !          enddo
     !       enddo
     !
     !       cDiff(4,index)=sqrt(cDiff(4,index))
     !       cDiff(1,index)=weights(n3)
     !       cDiff(2,index)=b(n3)
     !       cDiff(3,index)=ccmax
     !       index=index+1
     !    endif
     ! endif

  enddo

  !if (relaxDiv .eq. 1) then
  !   open(10, file='output/Nudge/ewInfo'//timeString,form='unformatted')
  !   write(10) cDiff
  !   close(10)
  !endif

  
  !BEGIN Test of scheme - are all weights the same before jitter?--------------------------------------------------------------------------------
  do n3 = 1, nGrand
    call Measure (psiGrand(1:nxx,:,n3), lpsi(n3,:))
    dlpsi(n3,:) = (data(:) - lpsi(n3,:)) * (data(:) - lpsi(n3,:))
  enddo

  do n3=1,nGrand

    !-log w_{rest}
    b(n3)=weights(n3)                              ! Weights from integration -log w_{rest}
    if (weightDiv .eq. 1) then
        barWeights(n3,1)=weights(n3)                   !for calculating distribution of weights
    endif
    write(*,*) '-log w_{rest}: ', b(n3)


    !Likelihood
    do ii=1,nObs
      b(n3)=b(n3)+0.5*dlpsi(n3,ii)/(qdData(ii)*qdData(ii))         ! Weights from observations d-H(psi^n)
    enddo
    if (weightDiv .eq. 1) then
        barWeights(n3,3)=b(n3) - barWeights(n3,1)      !for calculating distribution of weights
    endif
    write(*,*) 'likelihood: ', b(n3)


    !Transition
    do i=1,nxn*nyn
       b(n3)=b(n3) + 0.5*perbWeights(n3,i)*perbWeights(n3,i)
    enddo
    if (weightDiv .eq. 1) then
       barWeights(n3,2)=b(n3) - barWeights(n3,2) - barWeights(n3,1)      !for calculating distribution of weights
    endif
    write(*,*)'transition weights: ', b(n3)
  enddo

  write(*,*)'b test ',b
 ! END test of scheme-----------------------------------------------------------------------------------------------------------------------


  if (weightDiv .eq. 1) then
      open(35,file='output/barWeights-t'//timeString)
        do n3=1,nGrand
          write(35,*) barWeights(n3,1),barWeights(n3,2),barWeights(n3,3)
        enddo
      close(35)
  endif


 !add random error--------------------------------------------------------------------------------------------------------------------------
 !Mel-20|12|11-add jitter using Q correlated random variables
 do n3=1,nGrand

    !Mel-14|11|11-added for mixture density randomness
    call UniformRandomNumbers1D(0.,1.,size(xnum),xnum)
 
    !Uniform error
    if (xnum(1) .gt. eFac) then
        call UniformRandomNumbers1D(-uFac,uFac,nxn*nyn,randomVec)
        call MultiplySqrtQ(randomVec,corrRandom)
        do i=1,nxx
           do j=1,nyy
               psiGrand(i,j,n3)=psiGrand(i,j,n3) +  corrRandom(i+(j-1)*nxx)
           enddo
        enddo
        proposalWeight = 0.0
        !Value needed for final perturbed weights stored here
        perbWeights(n3,:)=perbWeights(n3,:) + randomVec(:)

    !Normal error
    else
       call NormalRandomNumbers1D(0.,1.,nxn*nyn,randomVec)
       call MultiplySqrtQ(randomVec,corrRandom)
       do i=1,nxx
           do j=1,nyy
               psiGrand(i,j,n3)=psiGrand(i,j,n3) +  nFac * corrRandom((i-1)*nxx+j)
           enddo
       enddo
       help = 2**(-nData/2)*pi**(nData/2)*nFac*uFac**(-nData)*((1-eFac)/eFac)
       do i=1,nxn*nyn
          proposalWeight(n3) = proposalWeight(n3) + help*exp(-0.5*randomVec(i)*randomVec(i))
       enddo
        !Value needed for final perturbed weights stored here
        perbWeights(n3,:)=perbWeights(n3,:) + nFac * randomVec(:)

    endif

    write(*,*)'Proposal weight: ',proposalWeight(n3)
 enddo

 !Calculate final weights-------------------------------------------------------------------------------------------------------------------
 do n3 = 1, nGrand
   call Measure (psiGrand(1:nxx,:,n3), lpsi(n3,:))
   dlpsi(n3,:) = (data(:) - lpsi(n3,:)) * (data(:) - lpsi(n3,:))
 enddo

  !Calculate new weights
  do n3=1,nGrand
    
    !-log w_{rest}
    b(n3)=weights(n3)
    write(*,*)'-log w_i^{rest}: ',b(n3) 

    !Likelihood weights
    do ii=1,nObs
      b(n3)=b(n3)+0.5*dlpsi(n3,ii)/(qdData(ii)*qdData(ii))   !Likelihood weights
    enddo 
    write(*,*)'likelihood: ',b(n3)

    !Transition weights
    do i=1,nxn*nyn
         b(n3)=b(n3) + 0.5 * perbWeights(n3,i) * perbWeights(n3,i)
    enddo
    write(*,*)'transition: ',b(n3)

    !Proposal weights
    b(n3)=b(n3)+proposalWeight(n3)
  enddo

  bmin=minval(b)
  write(*,*)'b after ',b,bmin
  do n3=1,nGrand
    if(bmin-b(n3) .lt. -100)then
      b(n3)=1.e-10                     ! These were abandoned anyway
    else
      b(n3)=exp(bmin-b(n3))
    endif
  enddo
!END calculation of final weights-----------------------------------------------------------------------------------------------------------

!Mel-10|01|12-added to calculate the strength of the aew movement
if (relaxDiv .eq. 1) then
   aewStrength(:)=0.
   do n3=1,nGrand
      do i=1,nxx
         do j=1,nyy
             aewStrength(n3)=aewStrength(n3)+(psiGrand(i,j,n3)-psinew(i,j,n3))*(psiGrand(i,j,n3)-psinew(i,j,n3))
         enddo
      enddo
      aewStrength(n3)=sqrt(aewStrength(n3))
   enddo

   open(10, file='output/Nudge/aewStrength'//timeString,form='unformatted')
   write(10)aewStrength
   close(10)
endif
       
! Resample global !!!

b=b/sum(b)
write(*,*)'b',b

bnew(1) = b(1)
do n3=2,nGrand
  bnew(n3)=bnew(n3-1)+b(n3)
enddo
call uniformrandomnumbers1D(0.,1./nGrand,xnum)
write(*,*)'x and bnew ',xnum(1),bnew
nn=1 
do n3=1,nGrand
  do while(xnum(1).gt.bnew(nn))
    nn=nn+1
  enddo 
  write(*,*)'newpsi ',n3,xnum(1),nn,bnew(nn)
  psinew(:,:,n3)=psiGrand(1:nxx,:,nn)
  xnum(1)=xnum(1)+1./nGrand
enddo


open(34,file='statusnew')
  write(34,*)nGrand,b(:),sum(b(:))
close(34)

open(34,file='prob-b'// timeString)
  write(34,*)' weights'
  do n3=1,nGrand
    write(34,*)b(n3)
  enddo
close(34)


! Entropy: we calculate the mutual information MI=-int post log (post/prior) d psi
 write(*,*)'Mutual Information'
 MI=0
do n3=1,nGrand
  MI=-b(n3)*log(nGrand*b(n3))
enddo
write(*,*)'Mutual Information ',MI 

write(*,*)'PsiSuper'
psiSuper = 0
do i=1,nxx
  do j=1,nyy
    do n3=1,nGrand
      psiSuper(i,j) = psiSuper(i,j) + b(n3) * psiGrand(i,j,n3)
    enddo
  enddo
enddo

open(78,file='psisuper',form='unformatted')
write(78)psinew(:,:,1)
close(78)

b = nGrand * b

write(*,*)'before time is ', time
do n3=1,nGrand
  write(*,*)psiGrand(9900,1,n3),b(n3)
enddo

psiGrand(1:nxx,:,:)=psinew(:,:,:) 

write(*,*)'after time is ', time
do n3=1,nGrand
  write(*,*)psiGrand(9900,1,n3)
enddo

! Measure the various ensemble members but now after the analysis
  do n3 = 1, nGrand
    call Measure (psiGrand(1:nxx,:,n3), lpsi(n3,:))
 enddo

!  print *, " Quality of ensemble after analysis: "
  var = 0.
  do n3 = 1, nGrand
    print *, " n3= ", n3,  &
             " New rms(data-lpsi) = ", sqrt(sum((data-lpsi(n3,:)) *          &
                                                (data-lpsi(n3,:))) / nObs), &
             " New mean(data-lpsi) = ",      sum((data-lpsi(n3,:))) / nObs      
   var = var + sum((data-lpsi(n3,:))*(data-lpsi(n3,:)))
  enddo
  write(*,*)' var posterior is ',var
 
write(*,*)' end of analyse'

! Deallocate space for remaining matrices that were allocated in this routine
deallocate (b, bnew, data)
deallocate (irank,psiSuper,psinew)
deallocate (lpsiMean)
deallocate (lpsi)
deallocate (dlpsi)
deallocate( xmin)
deallocate(cc,qobs,weights)
!Mel-25|08|11-added to calculate the distribution of the weights
if (weightDiv .eq. 1) then
    deallocate (barWeights)
endif
!Mel-02|09|11-added to allow full Q calculations
deallocate (xVector)
deallocate (exVector)
deallocate (QRxVector)
deallocate (hTemp)
deallocate (qTemp)
deallocate (randomVec)
deallocate (corrRandom)
deallocate (proposalWeight)
deallocate (perbWeights)
deallocate (psiVector)
deallocate (psiSin)

End Subroutine Analyse


