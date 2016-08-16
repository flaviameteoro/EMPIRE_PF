!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-08-16 17:26:10 pbrowne>
!!!
!!!    Module to store the variable comm_version that the user wishes
!!!    to use
!!!    Copyright (C) 2016  Philip A. Browne
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
!> module to store the parameter \link comms::comm_version
!> comm_version \endlink
!> to control the communication pattern that
!> empire will use.
!>
!> this should be set by the user before compilation
!> so that the correct communicator version is used.
!> see @ref comms_comm_version for an up-to-date description of the
!> options for this.
!>
!> This file is not tracked by git, so any changes that the user
!> makes here will not be updated by a *git pull* command
module communicator_version
  integer, parameter :: comm_version=1
  !< \memberof comms
  !< The style of communication
  !! between the model and empire.
  !! See @ref comms_comm_version for an up-to-date description of the
  !! options implemented
end module communicator_version
  
