!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-01-26 19:32:26 pbrowne>
!!!
!!!    Program to implement 4dEnVar
!!!    Copyright (C) 2015  Philip A. Browne
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

!> the main program to run 4DEnVar
program FourDEnVar
  use comms
  use var_data
  implicit none


  integer :: n !size of optimization state vector

  
  real(kind=kind(1.0d0)), allocatable, dimension(:) :: x0
  

  write(6,'(A)') 'PF: Starting PF code'
  call flush(6)
  !> set up EMPIRE coupling
  call initialise_mpi
  print*,'PF: setting controls'
  !> read in controlling data
  call set_var_controls
  print*,'PF: configuring model'
  !> call user specific routine for initialisation
  call configure_model
  print*,'allocating pf'
  !> allocate space for the filter
  call allocate_var


  !get the initial ensemble

  !get the initial guess

  !choose method and particulars
  

  select case (vardata%opt_method)
  case('cg')

     
     n = 2 !rosenbrock test function
     

     call subroutine_cg(vardata%cg_method,vardata%n,x0)

  case('lbfgs')
     call lbfgs_sub(vardata%n,x0)

  case('lbfgsb')
     call read_lbfgs_bounds
     call lbfgsb_sub(vardata%n,x0,vardata%nbd,vardata%l,vardata%u)

  case default

  end select


  !finish model by sending state with tag = 0


end program FourDEnVar
