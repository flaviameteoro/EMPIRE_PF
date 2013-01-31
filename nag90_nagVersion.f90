Module NagStuff
  !integer LSEED
  !integer SEED(1)
  !integer GENID
  !integer SUBID
  !integer IFAIL
  integer STATE(17)
  !integer LSTATE

End Module NagStuff

Subroutine SeedRandom (TheSeed)
  use NagStuff
  IMPLICIT NONE

  integer, INTENT(IN) :: TheSeed
  integer :: LSEED,SEED(1),GENID,SUBID,IFAIL,LSTATE
  integer, dimension(17) :: temp  

  LSEED = 1
  SEED(1) = TheSeed
  GENID = 1
  SUBID = 1
  LSTATE = 17
  IFAIL=0
  temp=state
  
  call g05kff(GENID,SUBID,SEED,LSEED,temp,LSTATE,IFAIL)
  
  state=temp
End Subroutine SeedRandom


Subroutine UniformRandomNumbers1D (minval, maxval, phi)
  use NagStuff
  IMPLICIT NONE

  real, INTENT(IN)                  :: minval, maxval
  real, dimension(:), INTENT(OUT)   :: phi
  integer                           :: n,IFAIL
  real ::A,B

  n = Size (phi, 1)
  A = MIN(minval,maxval)
  B = MAX(minval,maxval)
  IFAIL = 0
  CALL G05SQF(n,A,B,STATE,phi,IFAIL)

End Subroutine UniformRandomNumbers1D


Subroutine UniformRandomNumbers2D (minval, maxval, phi)
  Use NagStuff
  IMPLICIT NONE

  real, INTENT(IN)                  :: minval, maxval
  real, dimension(:,:), INTENT(OUT) :: phi
  integer                           :: nx, ny, IFAIL
  real:: A,B

  nx = Size (phi, 1)
  ny = Size (phi, 2)
!  call g05faf (minval, maxval, nx*ny, phi)
  A = MIN(minval,maxval)
  B = MAX(minval,maxval)
  IFAIL = 0
  CALL G05SQF(nx*ny,A,B,STATE,phi,IFAIL)

End Subroutine UniformRandomNumbers2D


Subroutine NormalRandomNumbers1D (mean, stdev, phi)
  Use NagStuff
  IMPLICIT NONE

  real, INTENT(IN)                  :: mean, stdev
  real, dimension(:), INTENT(OUT)   :: phi
  integer                           :: n, IFAIL
  real:: VAR

  n = Size (phi, 1)

  VAR = stdev**2
  IFAIL = 0
  CALL G05SKF(n,mean,VAR,STATE,phi,IFAIL)

End Subroutine NormalRandomNumbers1D


Subroutine NormalRandomNumbers2D (mean, stdev, phi)
  Use NagStuff
  IMPLICIT NONE

  real, INTENT(IN)                  :: mean, stdev
  real, dimension(:,:), INTENT(OUT) :: phi
  integer                           :: nx, ny, IFAIL
  real :: VAR
  nx = Size (phi, 1)
  ny = Size (phi, 2)

  VAR = stdev**2
  IFAIL = 0
  CALL G05SKF(nx*ny,mean,VAR,STATE,phi,IFAIL)

End Subroutine NormalRandomNumbers2D


Subroutine ComplexFFT2D (xr, xi)

  IMPLICIT NONE

  real, dimension(:,:), INTENT(INOUT) :: xr, xi

  integer                             :: nx, ny

  real, allocatable, dimension(:)     :: trigm, trign, work
  integer                             :: ifail, status

  nx = Size (xr, 1)
  ny = Size (xr, 2)
 
  if (nx /= Size (xi,1) .or. ny /= Size (xi, 2)) &
    stop "ComplexFFT2D called with different sizes for xr and xi"

  allocate (work(2 * nx * ny), Stat = status)
  if (status /= 0) stop "Allocation of work in ComplexFFT2D failed"

  allocate (trigm(2 * nx),     Stat = status)
  if (status /= 0) stop "Allocation of trigm in ComplexFFT2D failed"

  allocate (trign(2 * ny),     Stat = status)
  if (status /= 0) stop "Allocation of trign in ComplexFFT2D failed"

  call c06fuf (nx, ny, xr, xi, 'i', trigm, trign, work, ifail)

  deallocate (trigm, trign, work)

End Subroutine ComplexFFT2D


Subroutine SolveWithSVD (B, A, x)

! Purpose: Solve x, such that B = A . x, using SVD of A
! where A is a matrix of size dim * dim, and B is a matrix (dim * ncolb)
 
  IMPLICIT NONE

  real, dimension(:,:), INTENT(INOUT)  :: B
  real, dimension(:,:), INTENT(INOUT)  :: A
  real, dimension(:,:), INTENT(OUT)    :: x
     
  integer                              :: dim, ncolb
      
  integer  i, j, nn, ifail, status

  real     q(1,1)
       
  real, dimension(:,:), allocatable :: pT
  real, dimension(:),   allocatable :: work, sv

  real     tt

  dim   = size (b,1)
  ncolb = size (b,2)

  allocate (pT(dim, dim), STAT = status)
  if (status /= 0) stop "Allocation of pt in SolveWithSVD failed"
  
  allocate (work(dim * dim + 5 * dim), STAT = status)
  if (status /= 0) stop " Allocation of work in SolveWithSVD failed"
  
  allocate (sv(dim), STAT = status)
  if (status /= 0) stop " Allocation of sv in SolveWithSVD failed"
              
! Perform a Singular Value Decomposition of a real matrix
  call f02wef (dim, dim, A, dim, ncolb, B, dim, .true.,      &
       q, 1, sv, .true., pT, dim, work, ifail)

! Set tolerance tt to be the smallest allowed singular value 
  tt = sv(1) * sv(1) * 1.e-6
  nn = count(sv*sv .gt. tt)
 
! Calculate solution, note that on exit of f02wef B contains qT . B: 
  x = 0.
  do i = 1, dim
    do j = 1, nn
      x(i,:) = x(i,:) + pT(j,i) * 1/sv(j) * B(j,:)
    enddo
  enddo

!  write (24, *) "nn = ", nn
!  write (24, *) "sv = ", sv
!  write (24, *) "x = ", x

  deallocate (pT, work, sv)
      
End Subroutine SolveWithSVD
