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
!!$  subroutine calc_HQHTR
!!$    use Qdata
!!$    use Rdata
!!$    use sizes
!!$
!!$    use hsl_mc68_double
!!$    use hsl_mc69_double
!!$    implicit none
!!$    integer, parameter :: wp = kind(0.0d0)
!!$    
!!$    type(mc68_control) :: control68
!!$    type(mc68_info) :: info68
!!$    
!!$    integer :: lmap, flag
!!$    integer, dimension(:), allocatable  :: ptr, row, map
!!$    real(wp), dimension(:), allocatable :: val
!!$
!!$    integer :: i,j
!!$    
!!$    allocate(HQHTRrow(Qne),HQHTRcol(Qne),HQHTRval(Qne))
!!$    
!!$    j = 0
!!$    do i = 1,Qne
!!$       if(Qrow(i) .ge. 539617 .and. Qrow(i) .le. 566986 .and. &
!!$            & Qcol(i) .ge. 539617 .and. Qcol(i) .le. 566986) then
!!$          j = j + 1
!!$          HQHTRrow(j) = Qrow(i)-539616
!!$          HQHTRcol(j) = Qcol(i)-539616
!!$          if(HQHTRrow(j) .eq. HQHTRcol(j)) then
!!$             HQHTRval(j) = Qval(i)+Rval(1)
!!$          else
!!$             HQHTRval(j) = Qval(i)
!!$          end if
!!$       end if
!!$    end do
!!$    print*,'min/max HQHTRrow = ',minval(HQHTRrow(1:j)),maxval(HQHTRrow(1:j))
!!$    print*,'min/max HQHTRcol = ',minval(HQHTRcol(1:j)),maxval(HQHTRcol(1:j))
!!$    HQHTRn = obs_dim
!!$    HQHTRne = j
!!$    allocate(ptr(HQHTRn+1))
!!$    call mc69_coord_convert(HSL_MATRIX_REAL_SYM_PSDEF, HQHTRn, &
!!$         HQHTRn, HQHTRne, HQHTRrow(1:HQHTRne), HQHTRcol(1:HQHTRne), &
!!$         ptr, row, flag, val_in=HQHTRval(1:HQHTRne), &
!!$         val_out=val, lmap=lmap, map=map)
!!$    
!!$    print*,'mc69 finished with flag (should read 0): ',flag
!!$    print*,'HQHTRn = ',HQHTRn
!!$    allocate(HQHTR_order(HQHTRn))
!!$    call mc68_order(1, HQHTRn, ptr, row, HQHTR_order, control68, info68)
!!$    
!!$    call HQHTR_factor
!!$
!!$!    call save_HQHTR
!!$    
!!$  end subroutine calc_HQHTR
!!$  
!!$!  subroutine save_HQHTR
!!$!    open(2,file='HQHTRdata.dat',action='write',status='replace',form='unformat&
!!$!         &ted')
!!$!    write(2) HQHTRn
!!$!    write(2) HQHTRne
!!$!    write(2) HQHTRrow
!!$!    write(2) HQHTRcol    
!!$!    write(2) HQHTRval
!!$!    write(2) HQHTR_order
!!$!    close(2)
!!$!  end subroutine save_HQHTR
!!$
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



!!$implicit none
!!$integer :: HQHTRnev
!!$real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: HQHTRU
!!$real(kind=kind(1.0D0)), allocatable, dimension(:) :: HQHTRD
!!$contains
!!$  subroutine loadHQHTRevd
!!$    use sizes
!!$    open(2,file='QevdHQHTR.dat',action='read',form='unformatted')
!!$    
!!$    read(2) HQHTRnev
!!$    allocate(HQHTRU(obs_dim,HQHTRnev),HQHTRD(HQHTRnev))
!!$    print*,'allocation of HQHTRU and HQHTRD done'
!!$    
!!$    read(2) HQHTRD
!!$    read(2) HQHTRU
!!$    close(2)
!!$    print*,'loaded HQHTRev'
!!$!    QD = QD/1.0D2
!!$
!!$
!!$  end subroutine loadHQHTRevd
!!$
!!$  subroutine killHQHTRevd
!!$    if(allocated(HQHTRU)) deallocate(HQHTRU)
!!$    if(allocated(HQHTRD)) deallocate(HQHTRD)
!!$  end subroutine killHQHTRevd













end module hqht_plus_r
