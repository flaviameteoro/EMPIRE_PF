
!Module RandomFields

!  real, dimension(:, :),  Allocatable  ::     psiRandom

!End Module RandomFields


Module Sizes
  
  integer :: obs_dim,state_dim
  real(kind=kind(1.0D0)) :: dt, dx, d, tau0, dy, aa, ld, gac, f0, beta, rd !h
  real(kind=kind(1.0D0)), dimension(:), Allocatable :: f
!  real(kind=kind(1.0D0)), dimension(:,:), Allocatable :: tau

End Module Sizes



