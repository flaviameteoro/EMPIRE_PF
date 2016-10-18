!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:05:01 pbrowne>
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
  use output_empire, only : emp_e
  use pf_control
  use Sizes
  use comms

  IMPLICIT NONE
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)

!  integer :: particle,k,tag,mpi_err
!  integer :: mpi_status( MPI_STATUS_SIZE )
  logical, parameter :: nan_check=.false.

  if(nan_check) then
     if(.not. all(pf%psi .eq. pf%psi)) then
        write(emp_e,*)  'NaN detected in deterministic_model before mp&
             &i_send to model'
        stop
     end if
  end if

  call send_all_models(state_dim,pf%count,pf%psi,1)

  call recv_all_models(state_dim,pf%count,pf%psi)

  if(nan_check) then
     if(.not. all(pf%psi .eq. pf%psi)) then
        write(emp_e,*)  'NaN detected in deterministic_model after mpi_recv from &
             &model'
        stop
     end if
  end if


end subroutine deterministic_model
