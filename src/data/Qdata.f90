!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-18 10:09:54 pbrowne>
!!!
!!!    {one line to give the program's name and a brief idea of what it does.}
!!!    Copyright (C) 2014  Philip A. Browne
!!!
!!!    This program is free software: you can redistribute it and/or modify
!!!    it under the terms of the GNU General Public License as published by
!!!    the Free Software Foundation, either version 3 of the License, or
!!!    (at your option) any later version.
!!!
!!!    This program is distributed in the hope that it will be useful,
!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!!    GNU General Public License for more details.
!!!
!!!    You should have received a copy of the GNU General Public License
!!!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!!
!!!    Email: p.browne @ reading.ac.uk
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Qdata
implicit none
integer :: Qn,Qne
integer, allocatable, dimension(:) :: Qrow,Qcol
real(kind=kind(1.0D0)), allocatable, dimension(:) :: Qval,Qdiag
real(kind=kind(1.0d0)) :: Qscale
contains
  subroutine loadQ
    use sizes
    use pf_control
    integer :: i
    real(kind=kind(1.0d0)) :: eig
!!$    Qne = state_dim
!!$    allocate(Qrow(Qne),Qcol(Qne),Qval(Qne))
!!$    Qrow = (/ (i, i = 1,Qne) /)
!!$    Qcol = (/ (i, i = 1,Qne) /)
!!$    Qval = 1.0D-8

    
    open(2,file='Qdata.dat',action='read',form='unformatted')
    
    read(2) state_dim
    print*,'state_dim =',state_dim
    Qn = state_dim
    read(2) Qne
    print*,'Qne = ',Qne
    allocate(Qrow(Qne),Qcol(Qne),Qval(Qne))
    allocate(Qdiag(Qn))
    print*,'allocation of Qrow,Qcola dn Qval done boyee'

    do i = 1,Qn
       read(2) Qdiag(i)
    end do
    print*,'qdiag finished loading'
    do i = 1,Qne
       read(2) Qval(i)
    end do
    print*,'qval finished loading'
    do i = 1,Qne
       read(2) Qrow(i)
    end do
    print*,'qrow finished loading'
    do i = 1,Qne
       read(2) Qcol(i)
    end do
    print*,'qcol finished loading'
    close(2)
    
    print*,'loaded Q'

    open(2,file='Qeig.dat',action='read',status='old')
    read(2,*) eig
    close(2)
    Qval = Qval/sqrt(eig)
    Qdiag = Qdiag/sqrt(eig)
    Qscale = pf%Qscale
    Qval = Qval/sqrt(Qscale)
    Qdiag = Qdiag/sqrt(Qscale)
!    Qval = Qval/1.0D2

!    do i = 1,Qne
!       if(Qrow(i) .eq. Qcol(i)) Qval(i) = Qval(i) + 2.2D1
!    end do

!    Qval = Qval/1.0D7
!    print*,'inflated Q'

  end subroutine loadQ

  subroutine killQ
    if(allocated(Qrow)) deallocate(Qrow)
    if(allocated(Qcol)) deallocate(Qcol)
    if(allocated(Qval)) deallocate(Qval)
    if(allocated(Qdiag)) deallocate(Qdiag)
  end subroutine killQ

end module Qdata
