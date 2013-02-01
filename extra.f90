
!Module RandomFields

!  real, dimension(:, :),  Allocatable  ::     psiRandom

!End Module RandomFields


Module Sizes
  
  integer        :: nxx, nyy, nEns, nGrand, nxn, nyn
  real(kind=kind(1.0D0)) :: h, pi, dt, dx, d, tau0, dy, aa, ld, gac, f0, beta, rd 
  real(kind=kind(1.0D0)), dimension(:), Allocatable :: f
  real(kind=kind(1.0D0)), dimension(:,:), Allocatable :: tau

Contains

  Subroutine GetSizes

    use hadcm3_config

    IMPLICIT NONE

    integer :: j,j0, ncorio, ntau, status

!    open  (99, file="ProblemSizes")
!    read  (99, *) nxn, nyn,  nEns
!    read  (99, *) dt, dx
!    read  (99, *) d
!    read  (99, *) ncorio
!    read  (99, *) tau0, ntau

    
!    print *, " Values given for the problemsize dependent parameters:"
!    print *, " nxn =  ", nxn
!    print *, " nyn =  ", nyn
!    print *, " nEns = ", nEns
 
    nGrand = nEns
! U and V
    nxx = nxn*nyn*2
    nyy = 1
!    print*, " nxx =  ", nxx
!    print*, " nyy =  ", nyy
    
    !pi = 3.14159265359
    pi = 3.1415927
    dy = dx
    aa = dt/dx
    ld = d * dt/(dx*dx)
    
    h = 500.
    gac = 0.002 * 9.81
    f0 = 1.e-4
    beta = 1.e-10
    rd = sqrt(gac * h) / abs(f0)

    allocate( f(nyn), STAT = status)
    if (status /= 0) stop " Allocation of f failed"
    if (ncorio .eq. 0) then
      f = f0
    else
      j0 = nint(nyn/2.)
      do j=1,nyn
        f(j) = f0 + beta * dy * (j-j0)
      enddo
    endif

    allocate( tau(nxn,nyn), STAT = status)
    if (status /= 0) stop " Allocation of tau failed"
    tau0 = tau0/1024.
    if(ntau .eq. 0) then
       tau(3:nxn-1, 2:nyn-1) = tau0
    else
      do j = 2, nyn-1
        !tau(3:nxn-1, j) = - tau0 * cos(pi * (j-2)/real(nyn-3))
        tau(3:nxn-1, j) = - tau0 * cos(2 * pi * (j-2)/real(nyn-3))
      enddo
    endif
    
  End Subroutine GetSizes

End Module Sizes

Module SizesOfRandomField

  use Sizes
  integer :: nxr, nyr

Contains

  Subroutine GetSizesOfRandomField

    IMPLICIT NONE

    nxr = nxx 
    nyr = nyy 

  End Subroutine GetSizesOfRandomField

End Module SizesOfRandomField


Module TimeInfo

  integer :: time, time_to_start, time_to_stop, time_increment, &
              time_to_analyse, time_between_analyses

Contains

  Subroutine GetTimeInfo

    IMPLICIT NONE

    read  (99, *) time_to_start,        &           
                  time_to_stop,         &           
                  time_increment,       &          
                  time_to_analyse,      &            
                  time_between_analyses            

    print *, " Values given  for the various timestepping parameters: "
    print *, " time_to_start         = ", time_to_start       
    print *, " time_to_stop          = ", time_to_stop         
    print *, " time_increment        = ", time_increment              
    print *, " time_to_analyse       = ", time_to_analyse            
    print *, " time_between_analyses = ", time_between_analyses         
    time = time_to_start

  End Subroutine GetTimeInfo

End Module TimeInfo


Module ErrorsAndVariances

  real(kind=kind(1.0D0)) :: qo,qrel, nudgeFac
  real(kind=kind(1.0D0)), dimension(:), Allocatable :: qv, sqrtQv, qd
  real(kind=kind(1.0D0)), dimension(:,:), Allocatable :: fftQ, sqrtfftQ
  integer :: qSigma

Contains

  Subroutine GetErrorsAndVariances

    use Sizes

    IMPLICIT NONE

    integer :: status, ifail=0, i, j, p, q
    real(kind=kind(1.0D0)) :: qq, totalWave, help, qdu, qdv, qde                          !scaling of covariance matrix Q to allow for timestepping
    integer :: SOAR
    real(kind=kind(1.0D0)), dimension(nxn*nyn) :: waveValues
    
    read  (99, *) qo, qrel, qdu, qdv,qde, qSigma, SOAR, nudgeFac

    print *, " Values for errors and variances, correlations: .."
    print *, " qo    = ", qo
    print *, " qrel  = ", qrel
    print *, " qdu  = ", qdu
    print *, " qdv  = ", qdv
    print *, " qde  = ", qde
    print *, " qSigma = ", qSigma
    print *, " SOAR = ",SOAR
    print *, "nudgeFac = ", nudgeFac

     !Mel-02|09|11-added to allow calculation with full q vector
     allocate( qv(nxn*nyn), STAT = status)
     if (status /= 0) stop " Allocation of qv failed"
     allocate( sqrtQv(nxn*nyn), STAT = status)
     if (status /= 0) stop " Allocation of sqrtQv failed"

     !Mel-15|10|12-added to allow calculation with full qd vector
     allocate( qd(3*nxn*nyn), STAT = status)
     if (status /= 0) stop "Allocation of qd failed"

     qv=0.
     do i=3,nxn-1
        do j=3,nyn-1
           p=i-3
           q=j-3
           totalWave=sqrt(real(p*p+q*q))
           waveValues((i-1)*nyn+j)=totalWave
           !help = exp(-totalWave*totalWave/(2*qSigma*qSigma)) !Gaussian correlation function
           help = (1+totalWave/qSigma)*exp(-totalWave/qSigma) !SOAR correlation function
           if (help .gt. 1.e-12) then
              qv((i-1)*nyn+j) = help
           endif
           sqrtQv((i-1)*nyn+j) =sqrt(qv((i-1)*nyn+j))
        enddo
     enddo

     qd(:)=0.
     do i=2,nxn-1
        do j=2,nyn-1
           qd((i-1)*nyn+j) = qdu
           qd(nxn*nyn+(i-1)*nyn+j) = qdv
           qd(2*nxn*nyn+(i-1)*nyn+j) = qde
        enddo
     enddo

  End Subroutine GetErrorsAndVariances

  Subroutine ClearErrorsAndVariances
     
     if(allocated(qv))        deallocate (qv)
     if(allocated(sqrtQv))    deallocate (sqrtQv)     

  End Subroutine ClearErrorsAndVariances

End Module ErrorsAndVariances


Module EnsembleOfFields

  real(kind=kind(1.0D0)), dimension(:,:,:), Allocatable  :: psiEns   ! An ensemble of 2d-fields

End Module EnsembleOfFields


Module GrandField

  real(kind=kind(1.0D0)), dimension(:,:,:), Allocatable  :: psiGrand ! A lot of 2d-fields

End Module GrandField


Module RandomFilter

  real(kind=kind(1.0D0)), dimension(:,:),   Allocatable  :: xr         ! For random fields

End Module RandomFilter


Module ObservationField

   integer, dimension(:), Allocatable :: pseudoH !Specifies which elements of the array are observed
   real(kind=kind(1.0D0)), dimension(:), Allocatable :: qdData !Contains only the variances from observed data
   integer :: nObs, redObs, HStatus

 Contains

   Subroutine GetObsField

      use Sizes
      use ErrorsAndVariances

      IMPLICIT NONE

      integer :: status, i, j, k, loopEnd
      integer, dimension(:,:), Allocatable :: help

      read (99, *)  redObs, HStatus

      print *, "redObs =", redObs

      ! Allocate space for the observation field
      allocate( pseudoH(3*nxn*nyn), STAT = status)
      if (status /= 0) stop " Allocation of pseudoH failed"
      allocate( help(nxx, nyy), STAT = status)
      if (status /= 0) stop " Allocation of help failed"

      help(:,:)=0

      !Mel-22|08|11-staggered lines of observations-H1
      if (HStatus .eq. 1) then
         do i = 1,nxx,redObs
             do k = 0,redObs-1
                 do j = 1,nyy,redObs
                   if ((redObs .le. nxx) .AND. (j+k .le. nyy) .AND. (i+k .le. nxx)) then
                       help(i+k,j+k)=1
                   end if
               end do
             end do
         end do
         do j=1,nyy
            do i=1,nxx
               pseudoH(i+(j-1)*(nxx)) = help(i,j)
            enddo
         enddo
      endif

      !Mel-22|08|11-one observation
      !pseudoH(128,128)=1

      !Mel-17|09|11-regular pattern observations-H2
      if (HStatus .eq. 2) then
         pseudoH(:)=0.
         do i=1,nxn,redObs
             do j=1,nyn,redObs
                 pseudoH((i-1)*nyn+j) = 1
                 pseudoH(nxn*nyn+(i-1)*nyn+j) = 1
                 pseudoH(2*nxn*nyn+(i-1)*nyn+j) = 1
             enddo
          enddo
      endif

      !Mel-10|01|12-half the state observed-H3
      if (HStatus .eq. 3) then
         pseudoH(:)=0.
         loopEnd=nyy/2
         do i=1,nxx,redObs
             do j=1,loopEnd,redObs
                 pseudoH((i-1)*nxx+j) = 1
             enddo
         enddo
      endif     

      !Mel-17|01|12-3/4 of the state observed-H4
      if (HStatus .eq. 4) then
         pseudoH(:)=0.
         loopEnd=nyy/2
         do i=1,loopEnd,redObs
             do j=1,loopEnd,redObs
                 pseudoH((i-1)*nxx+j) = 1
             enddo
         enddo
         do i=loopEnd,nxx,redObs
             do j=1,nyy,redObs
                 pseudoH((i-1)*nxx+j) = 1
             enddo
         enddo
      endif
      
      !Full state observed minus boundary conditions
      if (HStatus .eq. 5) then
         pseudoH(:)=0.
         do i=2,nxn-1
            do j=2,nyn-1
               pseudoH((i-1)*nyn+j)=1.
               pseudoH(nxn*nyn+(i-1)*nyn+j)=1.
               pseudoH(2*nxn*nyn+(i-1)*nyn+j)=1.
            enddo
         enddo
      endif

      !Only 'e' observed every 30km = 3 grid point
      if (HStatus .eq. 6) then
         pseudoH(:)=0.
         do i=2,nxn-2,3
            do j=3,nyn-2,3
               pseudoH(2*nxn*nyn+(i-1)*nyn+j)=1
            enddo
         enddo
      endif



      nObs = count(pseudoH .eq. 1)
      print *, "nObs = ", nObs

      allocate( qdData(nObs), STAT=status) !storage for observed variable variances
      if (status /= 0) stop "Allocation of qdData failed"

      k=1
      do i=1,3*nxn*nyn
         if (pseudoH(i) .eq. 1) then
            qdData(k)=qd(i)
            k=k+1
         endif
      enddo

      deallocate(help)

   End Subroutine GetObsField


  Subroutine ClearPseudoH

      if(allocated(pseudoH))   deallocate (pseudoH)
      if(allocated(qdData))    deallocate (qdData)

  End Subroutine ClearPseudoH  
 
End Module ObservationField

Module Output

  integer :: talaTruth, talaObs, pdf, fullMovie, indTraj, weightDiv, relaxDiv
  real(kind=kind(1.0D0)), dimension(:,:), allocatable :: relaxStrength

Contains

  Subroutine GetOutput

    use Sizes
    use TimeInfo

    IMPLICIT NONE

    integer :: status=0, obNum, i
    character :: rankNo*3

    read  (99, *) talaTruth, talaObs, pdf, fullMovie, indTraj, weightDiv, relaxDiv
    print *, " Output required:"
    print *, " talaTruth =  ", talaTruth
    print *, " talaObs = ", talaObs
    print *, " pdf =  ", pdf
    print *, " fullMovie =  ", fullMovie
    print *, " indTraj =  ", indTraj
    print *, " weightDiv =  ", weightDiv
    print *, 'relaxDiv = ', relaxDiv


    if (pdf .eq. 1) then
       open (35, file='output/pdfValues',status='replace')
       close (35)
    endif

    if (fullMovie .eq. 1) then
       do i=0,nGrand
          write(rankNo,'(i3.3)') i
          open (36, file='output/full'//rankNo,status='replace')
          close (36)
       enddo
    endif

    if (indTraj .eq. 1) then
       do i=0,nGrand
          open(37,file='output/close'//rankNo,status='replace')
          close (37)
          open(38,file='output/apart'//rankNo,status='replace')
          close(38)
          open(39,file='output/closeObs'//rankNo,status='replace')
          close(39)
          open(40,file='output/apartObs'//rankNo,status='replace')
          close(40)
          open(41,file='output/trajObs'//rankNo,status='replace')
          close(41)
          open(42,file='output/trajUnObs'//rankNo,status='replace')
          close(42)
       enddo
    endif

    if (relaxDiv .eq. 1) then
       obNum = time_to_stop/time_between_analyses

       allocate( relaxStrength(3,obNum), STAT=status)
       if (status /= 0) stop 'Allocation of relaxStrength failed'
    endif
 
  End Subroutine GetOutput

  Subroutine ClearOutput

  IMPLICIT NONE

  if (allocated(relaxStrength)) deallocate (relaxStrength)

  End Subroutine ClearOutput


End Module Output
