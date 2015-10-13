!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-10-13 15:07:23 pbrowne>
!!!
!!!    module to store data for variational methods
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

!> module holding data for variational problems
module var_data

  implicit none
  type, public :: var_control_type
     character(6) :: opt_method  !< which optimization method to use
     !< currently this has a number of
     !< options: \n
     !< - 'cg' \n
     !< - 'lbfgs' \n
     !< - 'lbfgsb' \n
     integer :: cg_method !< which type of nonlinear CG method to use
     !< options are: \n
     !< 1  FLETCHER-REEVES \n
     !< 2  POLAK-RIBIERE (DEFAULT) \n
     !< 3  POSITIVE POLAK-RIBIERE (BETA=MAX{BETA,0} )
     real(kind=kind(1.0d0)) :: cg_eps !< convergence tolerance for CG
     !< method \n
     !< DEFAULT = 1.0d-5

     real(kind=kind(1.0d0)) :: lbfgs_factr

     !<    factr is a DOUBLE PRECISION variable that must be set by the user.\n
     !!       It is a tolerance in the termination test for the algorithm.
     !!       The iteration will stop when
     !!
     !!        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
     !!
     !!       where epsmch is the machine precision which is automatically
     !!       generated by the code. Typical values for factr on a computer
     !!       with 15 digits of accuracy in double precision are:\n
     !!       factr=1.d+12 for low accuracy;\n
     !!             1.d+7  for moderate accuracy;\n 
     !!             1.d+1  for extremely high accuracy.\n
     !!       The user can suppress this termination test by setting factr=0.
     !!
     !!   DEFAULT = 1.0d7

     real(kind=kind(1.0d0)) :: lbfgs_pgtol
     !<     pgtol is a double precision variable.\n
     !!       On entry pgtol >= 0 is specified by the user.  The iteration
     !!         will stop when
     !!
     !!                 max{|proj g_i | i = 1, ..., n} <= pgtol
     !!
     !!         where pg_i is the ith component of the projected gradient.\n
     !!       The user can suppress this termination test by setting pgtol=0.
     !!
     !!     DEFAULT = 1.0d-5

     real(kind=kind(1.0d0)), dimension(:), allocatable :: l,u,x0
     integer, dimension(:), allocatable :: nbd

     integer :: n !< the size of the state vector
     integer :: total_timesteps !< the total number of timesteps in
     !! the assimilation window

     integer, allocatable, dimension(:) :: ny !< array containing the number of
     !! observations. 
     !!
     !! ny(t) contains the number of observations at time t
     !!
     !! if no observations at time t then ny(t) = 0

  end type var_control_type
  type(var_control_type), save :: vardata !< the derived data type holding all controlling data


contains
  !> subroutine to ensure vardata is ok
  subroutine set_var_controls
    !    integer :: ios
    write(6,'(A)') 'Opening vardata.nml'

    call parse_vardata

    write(6,'(A)') 'vardata.nml successfully read to control pf code.'

  end subroutine set_var_controls


  !>subroutine to read the namelist file and save it to vardata datatype
  !!Here we read vardata.nml
  !!
  !! vardata.nml is a fortran namelist file. As such, within
  !! it there must be a line beginning
  !!
  !! &var_params
  !!
  !! To make it (probably) work, ensure there is a forward slash on
  !! the penultimate line and a blank line to end the file
  !!
  !! This is just the fortran standard for namelists though.
  !!
  !!
  !! On to the content...in any order, the vardata.nml may
  !! contain the following things:
  !! 
  !! 
  !!
  !! Integers:
  !! - \link var_data::var_control_type::cg_method cg_method\endlink
  !!
  !! Reals, double precision:
  !! - \link var_data::var_control_type::lbfgs_factr lbfgs_factr\endlink
  !! - \link var_data::var_control_type::lbfgs_pgtol lbfgs_pgtol\endlink
  !!
  !! 6 Characters:
  !! - \link var_data::var_control_type::opt_method opt_method\endlink
  !!
  subroutine parse_vardata
    implicit none
    character(*), parameter :: filename='vardata.nml'
    integer :: ios

    character(6) :: opt_method='CG'
    integer :: cg_method=2
    real(kind=kind(1.0d0)) :: cg_eps=1.0d-5
    real(kind=kind(1.0d0)) :: lbfgs_factr=1.0d+7
    real(kind=kind(1.0d0)) :: lbfgs_pgtol=1.0d-5 
    integer :: n=-1
    integer :: total_timesteps = -1

    namelist/var_params/opt_method,&
         &cg_method,&
         &cg_eps,&
         &lbfgs_factr,&
         &lbfgs_pgtol,&
         &total_timesteps,&
         &n

    !set defaults:
    vardata%opt_method = opt_method
    vardata%cg_method = cg_method
    vardata%cg_eps = cg_eps
    vardata%lbfgs_factr = lbfgs_factr
    vardata%lbfgs_pgtol = lbfgs_pgtol
    vardata%n = n
    vardata%total_timesteps = total_timesteps



    open(32,file=filename,iostat=ios,action='read'&
         &,status='old')
    if(ios .ne. 0) then
       write(*,*) 'Cannot open ',filename
       stop 2
    end if
    read(32,nml=var_params) 
    !      print*,time_obs,time_bwn_obs,nudgefac,gen_data
    close(32)

    select case (opt_method)
    case('cg')
       write(*,*) 'VAR_DATA: Nonlinear Conjugate Gradient Method selected'
       vardata%opt_method = opt_method
    case('lbfgs')
       write(*,*) 'VAR_DATA: Unconstrained L-BFGS Method selected'
       vardata%opt_method = opt_method
    case('lbfgsb')
       write(*,*) 'VAR_DATA: Bound constrained L-BFGS Method selected'
       vardata%opt_method = opt_method
    case default
       write(*,*) 'VAR_DATA ERROR: opt_method in ',filename,' incorrect&
            &ly given as ',opt_method
       write(*,*) 'VAR_DATA ERROR: Correct inputs are:'
       write(*,*) 'VAR_DATA ERROR: cg'
       write(*,*) 'VAR_DATA ERROR: lbfgs'
       write(*,*) 'VAR_DATA ERROR: lbfgsb'
       stop 3
    end select

    if(cg_method .lt. 1 .or. cg_method .gt. 3) then
       write(*,*) 'VAR_DATA ERROR: cg_method in ',filename, ' incorrect&
            &ly given as ',cg_method
       write(*,*) 'VAR_DATA ERROR: Correct inputs are:'
       write(*,*) 'VAR_DATA ERROR: 1  FLETCHER-REEVES '
       write(*,*) 'VAR_DATA ERROR: 2  POLAK-RIBIERE (DEFAULT)'
       write(*,*) 'VAR_DATA ERROR: 3  POSITIVE POLAK-RIBIERE ( BETA=MAX{BETA,0} )'
       stop 4
    else
       vardata%cg_method = cg_method
       if(vardata%opt_method .eq. 'cg') then
          select case (vardata%cg_method)
          case(1) 
             write(*,*) 'CG METHOD 1:  FLETCHER-REEVES '
          case(2)
             write(*,*) 'CG METHOD 2:  POLAK-RIBIERE (DEFAULT)'
          case(3)
             write(*,*) 'CG METHOD 3:  POSITIVE POLAK-RIBIERE'
          case default
          end select
       end if
    end if


    if(cg_eps .ne. 1.0d-5) then
       if(cg_eps .lt. 0.0d0) then
          write(*,*) 'VAR_DATA ERROR: cg_eps read as negative: ',cg_eps
          write(*,*) 'VAR_DATA ERROR: Please make cg_eps positive (small).'
          stop 5
       elseif(cg_eps .ge. 0.5d0) then
          write(*,*) 'VAR_DATA WARNING: cg_eps read as "large": '&
               &,cg_eps
          write(*,*) 'VAR_DATA WARNING: cg_eps default is 1.0d-5'
       else
          write(*,*) 'VAR_DATA: CG tolerance read in as ',cg_eps
          vardata%cg_eps = cg_eps
       end if
    else
       vardata%cg_eps = 1.0d-5
    end if

    if(lbfgs_factr .ne. 1.0d+7) then
       write(*,*) 'VAR_DATA: LBFGS factr read in as ',lbfgs_factr
       vardata%lbfgs_factr = lbfgs_factr    
    end if


    if(lbfgs_pgtol .ne. 1.0d-5) then
       write(*,*) 'VAR_DATA: LBFGS pgtol read in as ',lbfgs_pgtol
       vardata%lbfgs_pgtol = lbfgs_pgtol
    end if


    !    if(n .lt. 1) then
    !       write(*,*) 'VAR_DATA ERROR: n < 1'
    !       stop 2
    !    else
    !       vardata%n = n
    !    end if
    !print*,'nens = ',nens
    !vardata%n = nens-1
    
    vardata%n = -2
    

    if(total_timesteps .lt. 1) then
       write(*,*) 'VAR_DATA ERROR: total_timesteps < 1'
       stop 6
    else
       vardata%total_timesteps = total_timesteps
    end if



  end subroutine parse_vardata


  !> subroutine to allocate space for 4denvar
  subroutine allocate_vardata
    allocate(vardata%x0(vardata%n))

    allocate(vardata%ny(vardata%total_timesteps))


    if(vardata%opt_method == 'lbfgsb') then
       allocate(vardata%nbd(vardata%n))
       allocate(vardata%l(vardata%n))
       allocate(vardata%u(vardata%n))
    end if

  end subroutine allocate_vardata

  !> subroutine to deallocate space for 4denvar
  subroutine deallocate_vardata
  end subroutine deallocate_vardata


  !> subroutine to somehow read in bounds data
  subroutine read_lbfgsb_bounds
  end subroutine read_lbfgsb_bounds

  !> subroutine to somehow read in observation numbers
  subroutine read_observation_numbers
    vardata%ny = 0
    vardata%ny(8) = 3
    vardata%ny(16) = 3
    vardata%ny(24) = 3
    print*,'vardata%ny = ',vardata%ny
  end subroutine read_observation_numbers


end module var_data



