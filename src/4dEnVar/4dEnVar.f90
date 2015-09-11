!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-09-11 15:35:55 pbrowne>
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
  use sizes
  use comms
  use var_data
  use output_empire
  use fourdenvardata
!  use pf_control
  implicit none
  include 'mpif.h'

  integer :: message,mpi_err
  real(kind=kind(1.0d0)), dimension(3) :: xbar
  
  print*,'Starting 4DEnVar'
  call flush(6)
  !> set up EMPIRE coupling
  call initialise_mpi

  !> define output files
  call open_emp_o(pfrank)


  call random_seed_mpi(pfrank)
  !> read in controlling data
  call set_var_controls

  !> call user specific routine for initialisation
  call configure_model

  !> allocate space for the filter


  call allocate_vardata

  call allocate4denvardata

  !get initial states from ensemble members
  call recv_all_models(state_dim,cnt,xt)


  !get the initial ensemble perturbation matrix
  call read_ensemble_perturbation_matrix


  !get the initial guess
  call read_background_term


  vardata%x0 = 0.0d0


  call convert_control_to_state(vardata%n,vardata%x0,state_dim,xbar)
  print*,'init guess = ',xbar

  
  call read_observation_numbers

!  vardata%ny = 0
!  call NormalRandomNumbers1d(0.0d0,1.0d0,vardata%n,vardata%x0)
  !run optimization method

  if(pfrank == 0) then
     print*,'vardata%n = ',vardata%n

     select case (vardata%opt_method)
     case('cg')
        
        call subroutine_cg(vardata%cg_method,vardata%n,vardata%cg_eps,vardata%x0)
        
     case('lbfgs')
        call  lbfgs_sub(vardata%n,vardata%lbfgs_factr,vardata%lbfgs_pgtol,vardata%x0)
        
     case('lbfgsb')
        call read_lbfgsb_bounds
        call lbfgsb_sub(vardata%n,vardata%lbfgs_factr,vardata%lbfgs_pgtol,&
             vardata%x0,vardata%nbd,vardata%l,vardata%u)
        
     case default
        
     end select
     

     !finish other mpi loops with a broadcast
     message = 0
     call mpi_bcast(message,1,mpi_integer,0,pf_mpi_comm,mpi_err)

  else
     call fourdenvar_fcn
  end if

  print*,'optimization solution is:'
  print*,vardata%x0

  print*,'that means optimal state is:'
  call convert_control_to_state(vardata%n,vardata%x0,state_dim,xt(:&
       &,1))
  print*,xt(:,1)
  print*,'background state was:'
  print*,xb

  !finish model by sending state with tag = 3
  call send_all_models(state_dim,cnt,xt,3)

  call mpi_finalize(mpi_err)


end program FourDEnVar
