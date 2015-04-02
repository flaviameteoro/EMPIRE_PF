!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-04-02 14:19:04 pbrowne>
!!!
!!!    minimal code to setup and test empire comms
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

!> the main program
!!
program minimal_empire_comms
  use comms
  use pf_control
  use sizes
  implicit none
  include 'mpif.h'
  integer :: mpi_err

  print*,'RUNNING MINIMAL_EMPIRE_COMMS'
  print*,'EMPIRE COUPLING VERSION ',empire_version
  !> set up EMPIRE coupling
  call initialise_mpi

  call mpi_finalize(mpi_err)
  print*,'MINIMAL_EMPIRE_COMMS got to the end without breaking!'

end program minimal_empire_comms


