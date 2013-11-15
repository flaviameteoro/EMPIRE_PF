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
  use hsl_ma87_double
  implicit none
  integer :: HQHTRn,HQHTRne
  integer, allocatable, dimension(:) :: HQHTRrow,HQHTRcol,HQHTR_order
  real(kind=kind(1.0D0)), allocatable, dimension(:) :: HQHTRval
 
  type (ma87_keep) :: HQHTR_keep
  type (ma87_control) :: HQHTR_control
  type (ma87_info) :: HQHTR_info

contains

  subroutine load_HQHTR
    use Rdata
    use pf_control
    integer :: i
    
    open(2,file='HQHTRdata.dat',action='read',status='old',form='unfor&
         &matted')
    read(2) HQHTRn
    read(2) HQHTRne
    allocate(HQHTRrow(HQHTRne),HQHTRcol(HQHTRne),HQHTRval(HQHTRne))
    allocate(HQHTR_order(HQHTRn))
    read(2) HQHTRrow
    read(2) HQHTRcol
    read(2) HQHTRval
    read(2) HQHTR_order
    close(2)

    !scale HQHTR here as we scale Q...
    HQHTRval = HQHTRval/pf%Qscale


    do i = 1,HQHTRne
       if(HQHTRcol(i) .eq. HQHTRrow(i)) HQHTRval(i) = HQHTRval(i) + Rdiag(1)
    end do

    call HQHTR_factor
  end subroutine load_HQHTR

  subroutine HQHTR_factor

    use hsl_mc69_double
    integer, parameter :: wp = kind(0d0)
    

    integer :: lmap, flag
    integer, dimension(:), allocatable  :: ptr, row, map
    real(wp), dimension(:), allocatable :: val
    
    ! Read the first matrix in coordinate format
    ! Convert to HSL standard format

    allocate(ptr(HQHTRn+1))
    call mc69_coord_convert(HSL_MATRIX_REAL_SYM_PSDEF, HQHTRn, HQHTRn, HQHTRne, HQHTRrow, HQHTRcol, &
         ptr, row, flag, val_in=HQHTRval, val_out=val, lmap=lmap, map=map)

    print*,'ma87_analyse'
    ! Analyse
    call ma87_analyse(HQHTRn, ptr, row, HQHTR_order, HQHTR_keep, HQHTR_control, HQHTR_info)

    print*,'ma87_factor'
    ! Factor
    call ma87_factor(HQHTRn, ptr, row, val, HQHTR_order, HQHTR_keep, HQHTR_control&
         &, HQHTR_info)
    print*,'ma87 complete'
    
  end subroutine HQHTR_factor
  
  subroutine kill_HQHTR
    call ma87_finalise(HQHTR_keep, HQHTR_control)
  end subroutine kill_HQHTR




end module hqht_plus_r
