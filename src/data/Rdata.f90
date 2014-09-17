module Rdata
implicit none
integer :: Rn,Rne
integer, allocatable, dimension(:) :: Rrow,Rcol
real(kind=kind(1.0D0)), allocatable, dimension(:) :: Rval,Rdiag
contains
  subroutine loadR
    use sizes
    integer :: i
    Rn = obs_dim
    Rne = 1
    allocate(Rrow(Rne),Rcol(Rne),Rval(Rne))
    allocate(Rdiag(Rn))
    Rrow = (/ (i, i = 1,Rne) /)
    Rcol = (/ (i, i = 1,Rne) /)
    Rval = 0.0D0
    Rdiag = 0.3D0
  end subroutine loadR

  subroutine killR
    if(allocated(Rrow)) deallocate(Rrow)
    if(allocated(Rcol)) deallocate(Rcol)
    if(allocated(Rval)) deallocate(Rval)
    if(allocated(Rdiag)) deallocate(Rdiag)
  end subroutine killR

end module Rdata

module hqht_plus_r

  implicit none

contains

  subroutine load_HQHTR
    call HQHTR_factor
  end subroutine load_HQHTR

  subroutine HQHTR_factor

    
  end subroutine HQHTR_factor
  
  subroutine kill_HQHTR

  end subroutine kill_HQHTR




end module hqht_plus_r
