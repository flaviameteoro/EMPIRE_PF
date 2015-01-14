!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-01-14 18:11:17 pbrowne>
!!!
!!!    subroutine to simply move the model forward in time one timestep
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
!> @brief
!> subroutine to simply move the model forward in time one 
!> timestep
!>
!> PAB 21-05-2013

subroutine deterministic_model
  use pf_control
  use Sizes
  use comms

  IMPLICIT NONE
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)

  integer :: particle,k,tag,mpi_err
  integer :: mpi_status( MPI_STATUS_SIZE )
  logical, parameter :: nan_check=.false.

  if(nan_check) then
     if(.not. all(pf%psi .eq. pf%psi)) then
        stop 'NaN detected in deterministic_model before mpi_send to m&
             &odel'
     end if
  end if

  do k =1,pf%count
     particle = pf%particles(k)
     tag = 1
     call mpi_send(pf%psi(:,k),state_dim,MPI_DOUBLE_PRECISION&
          &,particle-1,tag,CPL_MPI_COMM,mpi_err)
  end do
  DO k = 1,pf%count
     particle = pf%particles(k)
     tag = 1
     CALL MPI_RECV(pf%psi(:,k), state_dim, MPI_DOUBLE_PRECISION, &
          particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
  END DO

  if(nan_check) then
     if(.not. all(pf%psi .eq. pf%psi)) then
        stop 'NaN detected in deterministic_model after mpi_recv from &
             &model'
     end if
  end if


end subroutine deterministic_model
