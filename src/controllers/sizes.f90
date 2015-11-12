!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-10-19 10:12:14 pbrowne>
!!!
!!!    Module that stores the dimension of observation and state spaces
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

!> @brief Module that stores the dimension of observation and state spaces
Module Sizes
  
  integer :: obs_dim     !< size of the observations held on this
                         !! process. For empire versions 1 and 2,
                         !! this is the total number of observations
  
  integer :: state_dim   !< size of the state held on this process. 
                         !! For empire versions 1 and 2,
                         !! this is the total size of the state vector
  
  integer :: obs_dim_g   !< global size of  obs  dim over all processes
  integer :: state_dim_g !< global size of state dim over all processes

End Module Sizes



