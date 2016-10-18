!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:34:20 pbrowne>
!!!
!!!    Subroutine to allocate space for the filtering code
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
!> subroutine to allocate space for the filtering code        
subroutine allocate_pf
  use output_empire, only : emp_e
  use pf_control
  use sizes
  use histogram_data
  integer :: st
  allocate(pf%weight(pf%nens),stat=st)
  if(st .ne. 0) then
     write(emp_e,*) 'Error in allocating pf%weight'
     stop
  end if
  pf%weight = -log(1.0D0/pf%nens)
  allocate(pf%psi(state_dim,pf%count),stat=st)
  if(st .ne. 0) then
     write(emp_e,*) 'Error in allocating pf%psi'
     stop
  end if
  
  if(pf%use_talagrand) then
     allocate(pf%talagrand(rhn_n,pf%nens+1),stat=st)
     if(st .ne. 0) then
        write(emp_e,*) 'Error in allocating pf%talagrand'
        stop
     end if
     pf%talagrand = 0
  end if
  
end subroutine allocate_pf
