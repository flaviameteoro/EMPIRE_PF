Subroutine AllocateArrays

  use Sizes
  use EnsembleOfFields

  integer :: status

  ! Allocate space for the ensemble of fields 
  allocate( psiEns(nxx+1, nyy,  nEns), STAT = status)
  if (status /= 0) stop " Allocation of psiEns failed"

End Subroutine AllocateArrays


Subroutine Statistics (time)

  use Sizes
  use GrandField
  use ErrorsAndVariances
  use Output 
  use ObservationField

  IMPLICIT NONE

  Interface
     Subroutine Talagrand (ensemble, truth, res)
       real(kind=kind(1.0D0)), dimension(:), intent(IN) :: ensemble
       real(kind=kind(1.0D0)), intent(IN) :: truth
       real(kind=kind(1.0D0)), intent(OUT) :: res
     End Subroutine Talagrand
  End Interface
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk) :: totalvar, storeVal
  real(kind=rk), dimension(nxx,nyy) ::  psiMeane
  real(kind=rk), dimension(nxx,nyy) ::  psiVare
  integer time, nM, talCount, status, i, j, k, l
  character alll*5, filename*20, timeString*6, particleNo*3
  real(kind=rk), dimension(:), Allocatable :: talValues !for storing the array of independent talagrand results
  real(kind=rk), dimension(:), Allocatable :: perturbVals !for storing the array of observation error perturbed ensemble members
  real(kind=rk), dimension(:), Allocatable :: randomVec !for storing the array of normally distributed random error
  real(kind=rk), dimension(:,:), Allocatable :: pdfValues !for storing the array of values for calculating pdfs
  real(kind=rk), dimension(:,:), Allocatable ::  qTruth

  write(*,*)'statistics', time

  ! Calculate mean, variance, and totalvariance
  do j=1,nyy
     do i=1,nxx
        psiMeane(i,j)=sum(psiGrand(i,j,:))/nGrand
        psiVare(i,j)=sum(psiGrand(i,j,:)*psiGrand(i,j,:))/nGrand
     enddo
  enddo
  psiVare=psiVare-psiMeane*psiMeane
  totalvar = SUM(psiVare) / SIZE(psiVare)

  ! Write to files
  write(alll,'(i5.5)')time
  open(60,file='output/vare.'//alll,form='unformatted')
  write(60)psiVare
  close(60)
  open(70,file='output/meane.'//alll,form='unformatted')
  write(70)psiMeane
  close(70)

  open(80,file='status')
  write(80,*)time,totalvar
  close(80)
  open(80,file='status')
  write(80,*)time,totalvar
  close(80)

  write(55,*)time,totalvar

  !open(80, file='output/prob.'//alll)
  !do nM=1,nGrand
  !  write(80,*)psiGrand(5,5,nM)
  !enddo
  !close(80)

  !Mel-05|01|12-Prints out values for individual particles
  !if ((time .eq. 1150) .OR. (time .eq. 1151)) then
  !   do nM=1,nGrand
  !      write(particleNo,'(i3.3)') nM
  !      open(81, file='output/Particle/particle'//particleNo//'_'//alll)
  !      do i=1,nxx
  !         do j=1,nyy
  !            write(81,*)psiGrand(i,j,nM)
  !         enddo
  !      enddo
  !      close(81)
  !   enddo
  !endif

  !Mel-08|11|11-Calculations for Talagrand diagrams
  if (talaTruth .eq. 1) then
     write(*,*)'Outputting talaValues compared to Truth at ', time

     allocate( qTruth(nxx,nyy), STAT = status)
     if (status /= 0) stop " Allocation of qTruth failed"

     if (mod(time,2) .eq. 0) then
        write(timeString,'(i6.6)') time
     else
        write(timeString,'(i6.6)') time-1
     endif

     filename = 'Data/truth' // timeString
     open (51, file=filename, form = 'unformatted') 
     read (51) qTruth
     close(51)

     talCount = ceiling(nxx/(3.0*qSigma+1))
     talCount = talCount*talCount
     allocate( talValues(talCount), STAT = status)
     if (status /= 0) stop " Allocation of talValues failed"

     k=1
     do i=1,nxx,3*qSigma+1
        do j=1,nyy,3*qSigma+1
           call talagrand(psiGrand(i,j,:),qTruth(i,j),storeVal)
           talValues(k)=storeVal
           k=k+1
        enddo
     enddo

     open(43,file='output/talaTruth'//alll,form='unformatted')
     write(43) talValues
     close(43)

     deallocate (qTruth)
     deallocate (talValues)
  endif

  !Mel-08|11|11-Calculations for Talagrand diagrams
  if (talaObs .eq. 1) then
     write(*,*)'Outputting talaValues compared to Obs at ', time

     allocate( qTruth(nxx,nyy), STAT = status)
     if (status /= 0) stop " Allocation of qTruth failed"
     allocate( perturbVals(nGrand), STAT=status)
     if (status /= 0) stop " Allocation of peturbVals failed"
     allocate( randomVec(nGrand), STAT=status)
     if (status /= 0) stop " Allocation of randomVec failed"


     if (mod(time,2) .eq. 0) then
        write(timeString,'(i6.6)') time
     else
        write(timeString,'(i6.6)') time-1
     endif

     filename = 'Data/datan' // timeString
     open (51, file=filename, form = 'unformatted') 
     read (51) qTruth
     close(51)

     talCount=0
     do i=1,nxx,3*qSigma+1
        do j=1,nyy,3*qSigma+1
           if (pseudoH((i-1)*nxx+j) .eq. 1) then
              talCount=talCount+1
           endif
        enddo
     enddo

     allocate( talValues(talCount), STAT = status)
     if (status /= 0) stop " Allocation of talValues failed"

     k=1
     do i=1,nxx,3*qSigma+1
        do j=1,nyy,3*qSigma+1
           if (pseudoH((i-1)*nxx+j) .eq. 1) then
              call g05fdf(0.0, 1.0, nGrand, randomVec)
              perturbVals(:) = psiGrand(i,j,:) + qd(:)*randomVec(:)
              call talagrand(perturbVals(:),qTruth(i,j),storeVal)
              talValues(k)=storeVal
              k=k+1
           endif
        enddo
     enddo

     open(43,file='output/talaObs'//alll,form='unformatted')
     write(43) talValues
     close(43)

     deallocate (qTruth)
     deallocate (talValues)
  endif

  !Mel-08|11|11-Calculations for pdfs
  if (pdf .eq. 1) then
     write(*,*)'Outputting pdfValues at ', time

     allocate( pdfValues(nxx,nGrand), STAT = status)
     if (status /= 0) stop " Allocation of pdfValues failed"

     k=1
     do i=1,nxx,16
        do j=1,nyy,16
           pdfValues(k,:)=psiGrand(i,j,:)
           k=k+1
        enddo
     enddo

     open(42,file='output/pdfValues'//alll,form='unformatted') 
     write(42) pdfValues
     close(42)

     deallocate (pdfValues)
  endif

  write(*,*) 'Finish statistics ', time

End Subroutine Statistics



Subroutine talagrand(ensemble,truth,res)

  use Sizes


  IMPLICIT NONE
  integer, parameter :: rk=kind(1.0D0)
  real(kind=rk), dimension(:), intent(in) :: ensemble
  real(kind=rk), intent(in) :: truth
  real(kind=rk), intent(out) :: res

  real(kind=rk), dimension(:), allocatable :: sortedE
  integer, dimension(:), allocatable :: irank
  integer :: status, ifail, i, test

  allocate( sortedE(nGrand), STAT=status)
  if (status /= 0) stop " Allocation of sortedE failed"
  allocate( irank(nGrand), STAT=status)
  if (status /= 0) stop " Allocation of irank failed" 

  call M01DAF(ensemble,1,nGrand,'A',irank,ifail)
  sortedE=ensemble
  call M01EAF(sortedE,1,nGrand,irank,ifail)

  test=0
  i=1
  do while (test .eq. 0)
     if (i .gt. nGrand) then
        test = 1
        res = i-0.5
     else if (sortedE(i) .gt. truth) then
        test = 1
        res = i-0.5
     endif
     i = i+1
  enddo

  deallocate(sortedE,irank)  

End Subroutine talagrand
