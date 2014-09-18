!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-18 10:09:54 pbrowne>
!!!
!!!    {one line to give the program's name and a brief idea of what it does.}
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
!!! Time-stamp: <2014-07-25 11:12:12 pbrowne>

module pf_control
  implicit none
  type, public :: pf_control_type
     integer :: nens !the number of ensemble members
     real(kind=kind(1.0D0)), allocatable, dimension(:) :: weight !stores the weights of the particles
     integer :: time_obs !the number of observations we will assimilate
     integer :: time_bwn_obs !the number of model timesteps between observations
     real(kind=kind(1.0D0)) :: nudgefac !the nudging factor
     logical :: gen_data,gen_Q,human_readable
     integer :: timestep=0
     real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: psi
     real(kind=kind(1.0D0)), allocatable, dimension(:) :: mean
     real(kind=kind(1.0D0)) :: nfac                !standard deviation of normal distribution in mixture density
     real(kind=kind(1.0D0)) :: ufac                !half width of the uniform distribution in mixture density
     real(kind=kind(1.0D0)) :: efac
     real(kind=kind(1.0D0)) :: keep,time
     real(kind=kind(1.0D0)) :: Qscale
     integer :: couple_root
     logical :: use_talagrand,use_weak,use_mean,use_var,use_traj,use_rmse
     integer, dimension(:,:), allocatable :: talagrand
     integer :: count
     integer,allocatable, dimension(:) :: particles
     character(2) :: type
     character(1) :: init
  end type pf_control_type
  type(pf_control_type) :: pf
  contains
    subroutine set_pf_controls
      integer :: ios
      write(6,'(A)') 'Opening pf_parameters.dat'
      open(32,file='pf_parameters.dat',iostat=ios,action='read',status='old')
      if(ios .ne. 0) stop 'Cannot open pf_parameters.dat'

      read(32,*) pf%time_obs
      read(32,*) pf%time_bwn_obs
      read(32,*) pf%nudgefac
      read(32,*) pf%gen_data
      read(32,*) pf%nfac
      read(32,*) pf%ufac
      read(32,*) pf%keep
      read(32,*) pf%Qscale
      read(32,*) pf%human_readable
      read(32,*) pf%use_talagrand
      read(32,*) pf%use_weak
      read(32,*) pf%use_mean
      read(32,*) pf%use_var
      read(32,*) pf%use_rmse
      read(32,*) pf%gen_Q
      read(32,*) pf%use_traj
      read(32,*) pf%type
      read(32,*) pf%init
      close(32)
      pf%efac = 0.001/pf%nens
      write(6,'(A)') 'pf_parameters.dat successfully read to control pf code.'
      call flush(6)
      if(pf%human_readable .and. pf%gen_data) then
         open(64,file='pf_data',iostat=ios,action='read',status='replace')
         if(ios .ne. 0) stop 'Error checking pf_data'
         close(64)
      end if

      !let us verify pf%type
      if(    pf%type .eq. 'EW') then
         print*,'Running the equivalent weights particle filter'
      elseif(pf%type .eq. 'SE') then
         print*,'Running a stochastic ensemble'
      elseif(pf%type .eq. 'SI') then
         print*,'Running the SIR particle filter'
      elseif(pf%type .eq. 'ET') then
         print*,'Running the Ensemble Transform Kalman Filter'
         print*,'Error: The ETKF is not implemented here'
         stop
      elseif(pf%type .eq. 'EA') then
         print*,'Running the Ensemble Adjustment Kalman Filter'
         print*,'Error: The EAKF is not implemented here yet'
         stop
      else
         print*,'Error: Incorrect filter type selected'
         print*,'Please ensure that pf%type in pf_parameters.dat is either:'
         print*,'EW                  the equivalent weights particle filter'
         print*,'SE                  a stochastic ensemble'
         print*,'SI                  the SIR particle filter'
         print*,'ET                  the Ensemble Transform Kalman Filter'
         print*,'EA                  the Ensemble Adjustment Kalman Filter'
         stop
      end if
         

    end subroutine set_pf_controls

    subroutine allocate_pf
      use sizes
      use histogram_data
      integer :: st
      allocate(pf%weight(pf%nens),stat=st)
      if(st .ne. 0) stop 'Error in allocating pf%weight'
      pf%weight = -log(1.0D0/pf%nens)
      allocate(pf%psi(state_dim,pf%count),stat=st)
      if(st .ne. 0) stop 'Error in allocating pf%psi'

      if(pf%use_talagrand) then
         allocate(pf%talagrand(rhn_n,pf%nens+1),stat=st)
         if(st .ne. 0) stop 'Error in allocating pf%talagrand'
         pf%talagrand = 0
      end if
      
!      allocate(pf%particles(pf%count),stat=st)
!      if(st .ne. 0) stop 'Error in allocating pf%particles'

    end subroutine allocate_pf

    subroutine deallocate_pf
      deallocate(pf%weight)
      deallocate(pf%psi)
      if(allocated(pf%talagrand)) deallocate(pf%talagrand)
      deallocate(pf%particles)
    end subroutine deallocate_pf


end module pf_control
