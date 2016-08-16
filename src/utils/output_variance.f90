!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-08-16 15:20:28 pbrowne>
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
!> subroutine to output ensemble variance
!>
subroutine output_variance(mean)
  use output_empire, only : unit_variance
  use pf_control
  use timestep_data
  use sizes
  use comms, only : pf_ens_size,pf_ens_comm,pf_ens_rank,pf_member_rank
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(state_dim), intent(in) :: mean !< the
  !ensemble mean
  real(kind=rk), dimension(state_dim) :: sum_perts
  real(kind=rk), dimension(state_dim) :: variance
  integer :: ios
  
  character(256) :: filename
  integer :: i,mpi_err
  
  sum_perts = 0.0_rk
  do i = 1,pf%count
     sum_perts = sum_perts + (pf%psi(:,i) - mean)*(pf%psi(:,i) - mean)
  end do

  !continue the sum across all da processes
  call mpi_reduce(sum_perts,variance,state_dim,MPI_DOUBLE_PRECISION,MPI_SUM&
       &,pf_ens_size-1,pf_ens_comm,mpi_err)
  
  ! divide by Ne - 1 to get the sample variance
  variance = variance/(pf%nens-1)

  if(pf_ens_rank .eq. pf_ens_size -1) then
  
     if(pf%timestep .eq. 0) then
        write(filename,'(A,i0)') 'ensemble_variance_',pf_member_rank
        open(unit_variance,file=trim(filename),iostat=ios,action='write',status='replace')
        if(ios .ne. 0)  then
           write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open &
                &file ',filename
           write(*,*) 'Very strange that I couldnt open it. Im goin&
                &g to stop now.'
           stop
        end if
     end if

     write(unit_variance,*) variance

     if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) close(unit_variance)
     
  end if !end if master process
end subroutine output_variance
