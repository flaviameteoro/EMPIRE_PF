!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:09:55 pbrowne>
!!!
!!!    Collection of subroutines to make multidimensional random arrays
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
module random_number_controls
  character(10) :: normal_generator='random_d'
  contains
    subroutine set_random_number_controls
      use output_empire, only : emp_o,unit_nml,emp_e
      integer :: ios
      logical :: file_exists
      character(14) :: empire_namelist='empire.nml'
      
      namelist/random_number_controls/normal_generator
      
      inquire(file=empire_namelist,exist=file_exists)
      if(file_exists) then
         open(unit_nml,file=empire_namelist,iostat=ios,action='read'&
              &,status='old',form='formatted')
         if(ios .ne. 0) then
            write(emp_e,*) 'Cannot open ',empire_namelist
            write(emp_e,*) 'open_emp_o ERROR'
            stop
         end if
         read(unit_nml,nml=random_number_controls,iostat=ios)
         close(unit_nml)
         if(ios .ne. 0) normal_generator='random_d'
         write(emp_o,*) 'random_number_controls: normal_generator read as: &
              &',normal_generator
      else
         normal_generator = 'random_d'
      end if
    end subroutine set_random_number_controls
end module random_number_controls

!> generate one dimension of uniform random numbers 
Subroutine UniformRandomNumbers1D(minv, maxv, n,phi)
!use random
implicit none
integer, parameter :: rk = kind(1.0D0)
integer, intent(in) :: n !< @param[in] n size of output vector
real(kind=rk), intent(in) :: minv!< @param[in] minv minimum value of uniform distribution
real(kind=rk), intent(in) :: maxv !< @param[in] maxv maximum value of uniform distribution
real(kind=rk), dimension(n), intent(out) :: phi !<@param[out] phi n dimensional uniform random numbers

call random_number(phi)

phi = minv + (maxv-minv)*phi
end Subroutine UniformRandomNumbers1D

!> generate one dimension of Normal random numbers 
Subroutine NormalRandomNumbers1D(mean,stdev,n,phi)
  use output_empire, only : emp_e
use random_number_controls
use random
use ziggurat
IMPLICIT NONE
integer, parameter :: rk = kind(1.0D0)
integer, intent(in) :: n !< @param[in] n size of output vector
real(kind=rk), INTENT(IN) :: mean  !< @param[in] mean mean of normal distribution
real(kind=rk), INTENT(IN) :: stdev !< @param[in] stdev Standard Deviation of normal distribution
real(kind=rk), dimension(n), INTENT(OUT) :: phi  !<@param[out] phi n dimensional normal random numbers
integer :: i

select case(normal_generator)
case('random_d')
   do i = 1,n
      phi(i) = mean+stdev*random_normal()
   end do
case('ziggurat')
   do i = 1,n
      phi(i) = mean+stdev*rnor()
   end do
case default
   write(emp_e,*) 'EMPIRE ERROR: wrong normal_generator selected in Nor&
        &malRandomNumbers1D'
   write(emp_e,*) 'EMPIRE ERROR: normal_generator = ',normal_generator&
        &,' is unknown'
   stop '-12'
end select



End Subroutine NormalRandomNumbers1D

!> generate two dimensional Normal random numbers
Subroutine NormalRandomNumbers2D(mean,stdev,n,k,phi)
  use output_empire, only : emp_e
use random_number_controls
use ziggurat
use random
IMPLICIT NONE
integer, parameter :: rk = kind(1.0D0)
integer, intent(in) :: n !< @param[in] n first dimension of output vector
integer, intent(in) :: k !< @param[in] k second dimension of output vector
real(kind=rk), INTENT(IN) :: mean  !< @param[in] mean mean of normal distribution
real(kind=rk), INTENT(IN) :: stdev !< @param[in] stdev Standard Deviation of normal distribution
real(kind=rk), dimension(n,k), INTENT(OUT) :: phi !< @param[out] phi n,k dimensional normal random numbers                                      
integer :: i,j

select case(normal_generator)
case('random_d')
   do j = 1,k
      do i = 1,n
         phi(i,j) = mean+stdev*random_normal()
      end do
   end do
case('ziggurat')
   do j = 1,k
      do i = 1,n
         phi(i,j) = mean+stdev*rnor()
      end do
   end do
case default
   write(emp_e,*) 'EMPIRE ERROR: wrong normal_generator selected in Nor&
        &malRandomNumbers2D'
   write(emp_e,*) 'EMPIRE ERROR: normal_generator = ',normal_generator&
        &,' is unknown'
   stop '-12'
end select

End Subroutine NormalRandomNumbers2D


!> generate one dimensional vector drawn from mixture density
!! @param[in] mean Mean of normal distribution
!! @param[in] stdev Standard deviation of normal distribution
!! @param[in] ufac half-width of uniform distribution
!! that is centered on the mean
!! @param[in] epsi Proportion controlling mixture draw.
!! if random_number > epsi then draw from uniform, else normal
!! @param[in] n size of output vector
!! @param[out] phi n dimensional mixture random numbers
!! @param[out] uniform True if mixture drawn from uniform. False if
!! drawn from normal
subroutine MixtureRandomNumbers1D(mean,stdev,ufac,epsi,n,phi,uniform)
use comms
use random
implicit none
include 'mpif.h'
real(kind=kind(1.0D0)), intent(in) :: mean,stdev,ufac,epsi
integer, intent(in) :: n
real(kind=kind(1.0D0)), dimension(n), intent(out) :: phi
logical, intent(out) :: uniform
real(kind=kind(1.0D0)) :: draw
integer :: mpi_err

if(comm_version .eq. 1 .or. comm_version .eq. 2) then
   call random_number(draw)
elseif(comm_version .eq. 3) then
   if(pf_member_rank .eq. 0) then
      call random_number(draw)
   end if
   call mpi_scatter(draw,1,MPI_DOUBLE_PRECISION,draw,1&
        &,MPI_DOUBLE_PRECISION,0,pf_member_rank,mpi_err)
else
    print*,'EMPIRE VERSION ',comm_version,' NOT SUPPORTED IN gen_rand'
    print*,'THIS IS AN ERROR. STOPPING'
    stop '-24'
end if

if(draw .gt. epsi) then
   call UniformRandomNumbers1D(mean-ufac,mean+ufac, n,phi)
   uniform = .true.
else
   call NormalRandomNumbers1D(mean,stdev,n,phi)
   uniform = .false.
end if

end subroutine MixtureRandomNumbers1D

!> generate two dimensional vector, each drawn from mixture density
!! @param[in] mean Mean of normal distribution
!! @param[in] stdev Standard deviation of normal distribution
!! @param[in] ufac half-width of uniform distribution
!! that is centered on the mean
!! @param[in] epsi Proportion controlling mixture draw.
!! if random_number > epsi then draw from uniform, else normal
!! @param[in] n first dimension of output vector
!! @param[in] k second dimension of output vector
!! @param[out] phi n,k dimensional mixture random numbers
!! @param[out] uniform k dimensional logical with uniform(i) True if
!! phi(:,i) drawn from uniform. False if
!! drawn from normal
subroutine MixtureRandomNumbers2D(mean,stdev,ufac,epsi,n,k,phi,uniform)
use random
implicit none
real(kind=kind(1.0D0)), intent(in) :: mean,stdev,ufac,epsi
integer, intent(in) :: n,k
real(kind=kind(1.0D0)), dimension(n,k), intent(out) :: phi
logical, dimension(k), intent(out) :: uniform
real(kind=kind(1.0D0)) :: draw
integer :: i

do  i = 1,k
   call random_number(draw)

   if(draw .gt. epsi) then
      call UniformRandomNumbers1D(mean-ufac,mean+ufac,n,phi(:,i))
      uniform(i) = .true.
   else
      call NormalRandomNumbers1D(mean,stdev,n,phi(:,i))
      uniform(i) = .false.
   end if
end do
end subroutine MixtureRandomNumbers2D


!> Subroutine to set the random seed across MPI threads
!! @param[in] pfid The process identifier of the MPI process
subroutine random_seed_mpi(pfid)
  use random_number_controls
  use pf_control
  use ziggurat
  integer, intent(in) :: pfid

  integer :: n
  integer, allocatable, dimension(:) :: seed

  call set_random_number_controls

  call random_seed(SIZE=n)
  allocate(seed(n))
  call random_seed(GET=seed)
  !add the particle filter id to the seed to make it              
  !independent on each process              
  if(.not. pf%gen_data) then
     seed = seed+pfid
  else
     seed = seed-1
  end if
  call random_seed(PUT=seed)
  deallocate(seed)

  call zigset(pfid+1)

  

end subroutine random_seed_mpi

