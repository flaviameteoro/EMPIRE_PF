!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-03-25 13:54:08 pbrowne>
!!!
!!!    A module to store data about the observation error covariance
!!!    matrix
!!!    The user should change according to their own needs
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
!> @brief Module to hold user supplied data for \f$R\f$
!> observation error covariance matrix
module Rdata
implicit none
contains
  !> Subroutine to load data for R
  subroutine loadR
  end subroutine loadR

  !> SUbroutine to deallocate R data
  subroutine killR
  end subroutine killR

end module Rdata

module hqht_plus_r

  implicit none

contains

  subroutine load_HQHTR
    call HQHTR_factor
  end subroutine load_HQHTR

  subroutine HQHTR_factor

    
  end subroutine HQHTR_factor
  
  subroutine kill_HQHTR

  end subroutine kill_HQHTR




end module hqht_plus_r
