!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-05-08 14:43:16 pbrowne>
!!!
!!!    Subroutine to generate Pf matrix given ensemble members on a communicator
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


!> subroutine to generate Pf matrix given ensemble members on a
!> communicator
subroutine generate_pf(stateDim,cnt,comm,x,pf)
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: stateDim !< size of the state vectors
  integer, intent(in) :: cnt      !< number of ensemble members on
  !<                                 this process
  integer, intent(in) :: comm     !< mpi communicator to use
  real(kind=rk), dimension(stateDim,cnt), intent(in) :: x !< the
  !< ensemble members on this process
  real(kind=rk), dimension(stateDim*(stateDim+1)/2), intent(out) :: Pf !< the Pf matrix
  !< i.e. the ensemble covariance matrix stored in packed storage in
  !! upper triangular form (see http://www.netlib.org/lapack/lug/node123.html)
  
  
  real(kind=rk), dimension(stateDim) :: mean
  real(kind=rk), dimension(stateDim,cnt) :: xp !< ensemble
  !perturbation matrix of those members only on this process

  integer :: m !< the total number of ensemble members
  integer :: i ! counter
  integer :: mpi_err !error flag for mpi

  ! get the total number of ensemble members used:
  call mpi_allreduce(cnt,m,1,MPI_INTEGER,MPI_SUM,comm,mpi_err)
  if(mpi_err .ne. MPI_SUCCESS) then
     print*,'ERROR in generate_pf: mpi_allreduce 1 failed with flag ',mpi_err
     stop
  end if
  
  !get the ensemble mean
  xp(:,1) = sum(x,dim=2)
  call mpi_allreduce(xp(:,1),mean,stateDim,MPI_DOUBLE_PRECISION&
       &,MPI_SUM,comm,mpi_err)
  if(mpi_err .ne. MPI_SUCCESS) then
     print*,'ERROR in generate_pf: mpi_allreduce 2 failed with flag ',mpi_err
     stop
  end if
  ! use BLAS to perform mean = mean/real(m,rk)
  call dscal(stateDim,1.0d0/real(m,rk),mean,1)


  !compute the unscaled ensemble perturbation matrix on this process
  do i = 1,cnt
     xp(:,i) = x(:,i)-mean
  end do

  !form pf with only ensemble members on this process
  call dsyrk('U',&      !UPLO  upper part of matrix formed
             'N',&      !TRANS no transposed done
             stateDim,& !N     dim of output matrix [square]
             cnt,&      !K     number of columns in input
             1.0d0,&    !ALPHA alpha scalar
             xp,&       !A     input matrix
             stateDim,& !LDA   leading dimension of A
             0.0d0,&    !BETA  add to empty array
             Pf,&       !C     output matrix
             stateDim&  !LDC   dimension of output matrix
             )

  ! get pf on all processors
  call mpi_allreduce(MPI_IN_PLACE,Pf,stateDim*(stateDim+1)/2,MPI_DOUBLE_PRECISION&
       &,MPI_SUM,comm,mpi_err)
  if(mpi_err .ne. MPI_SUCCESS) then
     print*,'ERROR in generate_pf: mpi_allreduce 3 failed with flag ',mpi_err
     stop
  end if
  
  !scale pf appropriately
  ! use BLAS to perform pf = pf/(m-1)
  call dscal(stateDim**2,1.0d0/real(m-1,rk),pf,1)
  
end subroutine generate_pf
