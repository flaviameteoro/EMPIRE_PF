!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-10-13 13:51:31 pbrowne>
!!!
!!!    module to store data for 4DEnVar
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


!> module holding data specific for 4denvar, not var itself.
!> this is necessary because of the difference in x in optimization
!> and in the model state.
module fourdenvardata
  implicit none
  integer :: m !< the number of perturbations, or nens-1
  real(kind=kind(1.0d0)), allocatable, dimension(:) :: xb !< the background
  !<guess
  real(kind=kind(1.0d0)), allocatable, dimension(:,:) :: x0 !< the
  !<initial ensemble perturbation matrix
  !!
  !! THIS IS ***NOT*** SCALED!!
  real(kind=kind(1.0d0)), allocatable, dimension(:,:) :: xt !< the
  !<current ensemble
contains
  subroutine allocate4denvardata
    use sizes, only : state_dim
    use comms, only : cnt,nens,pfrank
    m = nens-1
    allocate(xb(state_dim))
    if(pfrank .ne. 0) then
       allocate(x0(state_dim,cnt))
       allocate(xt(state_dim,cnt)) 
    else
       allocate(x0(state_dim,cnt-1)) !no perturbation associated
       !with optimization solution
       allocate(xt(state_dim,cnt))!here we include the optimization
       !current solution as well as the perturbations
    end if

  end subroutine allocate4denvardata

  !> subroutine to read xb from file
  subroutine read_background_term()
    real(kind=kind(1.0d0)), dimension(3) :: tempxb
    !read xb from file
    print*,'READ xb from file not implemented'
    !stop 10      
    !type in a normally distributed random vector
    call NormalRandomNumbers1d(0.0d0,1.0d0,3,tempxb)
    !      tempxb = (/1.8368174447696750d-1,-1.3102740999288161d-2,  &
    !           &-5.9318375014562030d-1/)
    call Bhalf(1,tempxb,xb)
    print*,'background guess pert = '
    print*,xb

    xb = xb + (/-3.12346395, -3.12529803, 20.69823159/)
    print*,'background guess xb = '
    print*,xb
  end subroutine read_background_term

  subroutine deallocate4denvardata
    deallocate(xb)
  end subroutine deallocate4denvardata

  !> subroutine to read in the ensemble perturbation matrix
  !!
  !! we need to fill in the entries of x0 here
  !!
  !! 
  subroutine read_ensemble_perturbation_matrix
    use comms, only : cnt,pfrank
    use sizes
    real(kind=kind(1.0d0)),allocatable, dimension(:,:) :: tempx0

    if(pfrank .ne. 0) then
       ! READ x0(particle) from file
       print*,'READ FROM FILE NOT IMPLEMENTED'
       !stop 7
       allocate(tempx0(state_dim,cnt))
       call NormalRandomNumbers2D(0.0d0,1.0d0,state_dim,cnt,tempx0)
       !         print*,'sup tempx0 = ',tempx0
       call Bhalf(cnt,tempx0,x0)
       !         print*,'arse x0 = ',x0
       deallocate(tempx0)
    else!no perturbation                  
       !associated with optimization solution

       ! READ x0(particle-1) from file
       print*,'READ FROM FILE NOT IMPLEMENTED'
       !stop 8
       allocate(tempx0(state_dim,cnt-1))
       call NormalRandomNumbers2D(0.0d0,1.0d0,state_dim,cnt-1,tempx0)
       !         print*,'sup tempx0 = ',tempx0
       call Bhalf(cnt-1,tempx0,x0)
       !         print*,'arse x0 = ',x0
       deallocate(tempx0)
    end if
    x0 = x0!/(real(m-1,kind(1.0d0))**0.5d0)

    print*,'read_ensemble_perturbation_matrix'
    print*,x0
    print*,'read_ensemble_perturbation_matrix'
    print*,'x0(1,3) = ',x0(1,3)
    print*,'x0(3,2) = ',x0(3,2)
  end subroutine read_ensemble_perturbation_matrix



end module fourdenvardata
