Subroutine SeedRandom (seed)
  
  IMPLICIT NONE

  integer seed
  
  call g05cbf (seed)
  
End Subroutine SeedRandom


Subroutine UniformRandomNumbers1D (minval, maxval, phi)

  IMPLICIT NONE

  real, INTENT(IN)                  :: minval, maxval
  real, dimension(:), INTENT(OUT)   :: phi
  integer                           :: n

  n = Size (phi, 1)
  call g05faf (minval, maxval, n, phi)

End Subroutine UniformRandomNumbers1D


Subroutine UniformRandomNumbers2D (minval, maxval, phi)

  IMPLICIT NONE

  real, INTENT(IN)                  :: minval, maxval
  real, dimension(:,:), INTENT(OUT) :: phi
  integer                           :: nx, ny

  nx = Size (phi, 1)
  ny = Size (phi, 2)
  call g05faf (minval, maxval, nx*ny, phi)

End Subroutine UniformRandomNumbers2D


Subroutine NormalRandomNumbers1D (mean, stdev, phi)

  IMPLICIT NONE

  real, INTENT(IN)                  :: mean, stdev
  real, dimension(:), INTENT(OUT)   :: phi
  integer                           :: n

  n = Size (phi, 1)

  call g05fdf (mean, stdev, n, phi)

End Subroutine NormalRandomNumbers1D


Subroutine NormalRandomNumbers2D (mean, stdev, phi)

  IMPLICIT NONE

  real, INTENT(IN)                  :: mean, stdev
  real, dimension(:,:), INTENT(OUT) :: phi

  integer                           :: nx, ny

  nx = Size (phi, 1)
  ny = Size (phi, 2)

  call g05fdf (mean, stdev, nx*ny, phi)

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
