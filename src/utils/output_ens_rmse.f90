!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:32:38 pbrowne>
!!!
!!!    Subroutine to output RMSE
!!!    Copyright (C) 2015 Philip A. Browne
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
!> subroutine to output RMSEs
!>
subroutine output_ens_rmse()
  use output_empire, only : unit_ens_rmse,emp_e
  use pf_control
  use timestep_data
  use sizes
  use comms
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(state_dim,pf%count) :: e
  real(kind=rk), dimension(state_dim) :: mse
  real(kind=rk), dimension(state_dim) :: truth
  integer :: ios,i
  integer :: ensemble_comm,comm_root,mpi_err,rank
  
  character(256) :: filename
  
  call get_truth(truth)


  !calculate the error
  do i = 1,pf%count
     e(:,i) = pf%psi(:,i) - truth
  end do
  !square the error
  e = e*e

  !now let us sum the fields
  !for memory usage we we store the summed across the ensemble
  !members in the TRUTH variable
  truth = sum(e,dim=2)

  !now the sum is only over those ensemble members on this mpi
  !process. let us use an mpi reduce to sum up all the rest
  !and store the result on the highest ranking empire process


  select case(comm_version)
  case(1,2,5)
     ensemble_comm = pf_mpi_comm
     comm_root = npfs-1
     rank = pfrank
     write(filename,'(A,A,A,i7.7)') 'ens_',trim(pf%rmse_filename),'_',pf%timestep
  case(3)
     ensemble_comm = pf_ens_comm
     comm_root = pf_ens_size-1
     rank = pf_ens_rank
     write(filename,'(A,A,A,i7.7,A,i0)') 'ens_',trim(pf%rmse_filename),'_'&
          &,pf%timestep,'.',pf_ens_rank     
  case default
     write(emp_e,*) 'EMPIRE ERROR: output_ens_rmse comm_version'
     write(emp_e,*) 'EMPIRE ERROR: comm_version',comm_version
     write(emp_e,*) 'EMPIRE ERROR: is currently not supported. STOPPING'
     stop '-9'
  end select
  
  call mpi_reduce(truth,mse,state_dim,MPI_DOUBLE_PRECISION,MPI_SUM&
       &,comm_root,ensemble_comm,mpi_err)

  ! now need to ensure the mean is correct and take the root:
  if( rank == comm_root ) then
     mse = mse/real(nens,rk)
     mse = sqrt(mse)
     
  
     !mse is now the rmse that we want. Time to output it

     open(unit_ens_rmse,file=trim(filename),iostat=ios,action='write',status='replace'&
          &,form='formatted')
     if(ios .ne. 0)  then
        write(emp_e,*) 'EMPIRE ERROR: CANNOT open file ',filename
        write(emp_e,*) 'Very strange that I couldnt open it. Im going to stop&
             & now.'
        stop '-12'
     end if
     write(unit_ens_rmse,*) mse
     close(unit_ens_rmse)
     

  end if




end subroutine output_ens_rmse
