!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-05-31 10:28:16 pbrowne>
!!!
!!!    Program to implement 3DVar
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

!> @todo make work with empire version 3
subroutine three_d_var(x)
  use threedvar_data
  use sizes
  use var_data
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(state_dim), intent(inout) :: x

  vardata%n = state_dim
  
  if(allocated(xb)) deallocate(xb)
  allocate(xb(state_dim))
  xb = x
  vardata%total_timesteps = 1
  call allocate_vardata
  
  vardata%x0 = x
  
  
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
     print*,'three_d_var ERROR: vardata%opt_method incorrect. Stopping&
          &'
     stop '-78'
  end select
  
  x = vardata%x0
  call deallocate_vardata

end subroutine three_d_var
