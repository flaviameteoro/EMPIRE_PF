!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-11-23 11:33:01 pbrowne>
!!!
!!!    Subroutine to compute the relaxation strength
!!!    Copyright (C) 2016 Philip A. Browne
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
!> subroutine to compute the relaxation strength
!>
!>
!> there are currrently 4 different relaxation profiles implemented:
!>  - 'zero_linear'
!>  - 'power_law'
!>  - 'exponential'
!>  - 'constant'
!>
!> plus
!>
!>  - 'user'
!>
!> Details:
!>
!>   - 'zero_linear'
!>          relaxation strength is zero until relaxation_freetime.
!>          relaxation strength grows linearly from 0 at
!>          tau=relaxation_freetime until nudgefac at tau=1
!>
!>   - 'power_law'
!>          relaxation strength is zero until relaxation_freetime.
!>          relaxation strength = nudgefac *
!>          ((tau-relaxation_freetime)/(1-relaxation_freetime))^p
!>
!>   - 'exponential'
!>          relaxation strength is zero until relaxation_freetime.
!>          relaxation strength = nudgefac *
!>          (exp((tau-relaxation_freetime)/(1-relaxation_freetime))-1)/exp(1)
!>
!>   - 'constant'
!>          relaxation strength is zero until relaxation_freetime.
!>          relaxation strength = nudgefac from
!>          tau=relaxation_freetime until nudgefac at tau=1
!>
!>   - 'user'
!>          Calls a user defined function to calculate the relaxation strength
!>
!> @todo produce a plot of implemented relaxation profiles
subroutine relaxation_profile(tau,p,zero)
  use output_empire, only : emp_e
  use pf_control
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), intent(in) :: tau !< the pseudotime between
  !! observations.  \$f 0 \le \tau \le 1 \$f
  real(kind=rk), intent(out) :: p !< the relaxation strength at time tau
  logical, intent(out) :: zero !< switch if false is no relaxation

  select case (pf%relaxation_type)
  case ('zero_linear')
     if ( tau < pf%relaxation_freetime ) then
        p = 0.0d0
        zero = .true.
     else
        p = pf%nudgefac*(tau-pf%relaxation_freetime)/(1.0d0-pf&
             &%relaxation_freetime)
        zero = .false.
     end if
  case ('power_law')
     if ( tau < pf%relaxation_freetime ) then
        p = 0.0d0
        zero = .true.
     else
        p = pf%nudgefac * ((tau-pf%relaxation_freetime)/(1.0d0-pf&
             &%relaxation_freetime))**pf%power_law_p
     end if
  case ('exponential')
     if ( tau < pf%relaxation_freetime ) then
        p = 0.0d0
        zero = .true.
     else
        p = pf%nudgefac*(exp( (tau-pf%relaxation_freetime)/(1.0d0-pf&
             &%relaxation_freetime) ) - 1 )/exp(1.0d0)
        zero = .false.
     end if
  case ('constant')
     if ( tau < pf%relaxation_freetime ) then
        p = 0.0d0
        zero = .true.
     else
        p = pf%nudgefac
        zero = .false.
     end if
  case ('user')
     call user_relaxation_profile(tau,p,zero)
  case default
     write(emp_e,*) 'Wrong relaxation profile'
  end select
end subroutine relaxation_profile
