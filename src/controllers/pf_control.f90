!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-09-30 12:24:29 pbrowne>
!!!
!!!    module to hold all the information to control the the main program
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

!> module pf_control holds all the information to control the the main program
module pf_control
  implicit none
  type, public :: pf_control_type
     integer :: nens !<the total number of ensemble members
     real(kind=kind(1.0D0)), allocatable, dimension(:) :: weight !< the negative log of the weights of the particles
     integer :: time_obs !< the number of observations we will assimilate
     integer :: time_bwn_obs !< the number of model timesteps between observations
     real(kind=kind(1.0D0)) :: nudgefac !< the nudging factor
     logical :: gen_data !< true generates synthetic obs for a twin experiment
     logical :: gen_Q    !< true attempts to build up \f$Q\f$ from
     !<long model run. UNUSED. DOES NOTHING!
     integer :: timestep=0     !< the current timestep as the model progresses
     real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: psi !< state vector of ensemble members on this mpi process
     real(kind=kind(1.0D0)), allocatable, dimension(:) :: mean !< mean state vector
     real(kind=kind(1.0D0)) :: nfac                !< standard deviation of normal distribution in mixture density
     real(kind=kind(1.0D0)) :: ufac                !< half width of the uniform distribution in mixture density
     real(kind=kind(1.0D0)) :: efac
     real(kind=kind(1.0D0)) :: keep                !< proportion of particles to keep in EWPF EW step
     real(kind=kind(1.0D0)) :: time                !< dunno
     real(kind=kind(1.0D0)) :: Qscale              !< scalar to
                                                   !!multiply Q by
     real(kind=kind(1.0d0)) :: rho                 !< enkf inflation factor
                                                   !! so that \f$P_f =
                                                   !! (1+\rho)^2P_f\f$
     real(kind=kind(1.0d0)) :: len                 !< R localisation
                                                   !! length scale
     !! The entries in the observation error covariance matrix
     !! \f$R\f$ 
     !! are multiplied by the function
     !! \f$ \exp\left(\frac{\mathrm{dist}^2}{2 \mathrm{len}^2}\right) \f$
     integer :: couple_root                        !< empire master processor
     logical :: use_talagrand !< switch if true outputs rank
     !!histograms. See \link histogram_data::load_histogram_data
     !!load_histogram_data \endlink for details.
     logical :: output_weights      !< switch if true outputs
     !!ensemble weights
     logical :: output_forecast      !< switch if true outputs
     !!the entire forecast ensemble
     logical :: use_mean      !< switch if true outputs ensemble mean
     logical :: use_variance  !< switch if true outputs ensemble variance
     logical :: use_traj      !< switch if true outputs trajectories
     logical :: use_spatial_rmse      !< switch if true outputs Root Mean
     !!Square Errors \n See @ref rmse for more information
     logical :: use_ens_rmse  !< switich if true outputs the field of
     !!root mean squared errors where \f$ rmse(j) = \sqrt{
     !!\frac{1}{N_e}\sum_{i=1}^{N_e} (x_i(j)-x^t(j))^2}\f$
     character(250) :: rmse_filename !< string to hold the name of
     !< the file to output rmse to
     character(250) :: forecast_path !< string to hold the path to
     !!output the forecast ensemble
     
     integer, dimension(:,:), allocatable :: talagrand !< storage for rank histograms
     integer :: count         !< number of ensemble members associated with this MPI process
     integer,allocatable, dimension(:) :: particles !< particles associates with this MPI process
     character(2) :: filter   !< which filter to use
                              !< currently this has a number of
                              !<options:
                              !< - SE -- a stochastic ensemble
                              !< - DE -- a deterministic ensemble
                              !< - SI -- the SIR filter
                              !< - LE -- the L-ETKF with noise
                              !< - LD -- the L-ETKF without noise
                              !< - EW -- the Equivalent Weights filter
                              !< - EZ -- the Zhu equal weights filter
                              !< particle filter
                              !< - LS -- the L-ETKS with noise
                              !< - 3D -- 3DVar
                              !< @todo change these to a longer string
     character(1) :: init     !< which method to initialise ensemble
                              !< currently this has a number of
                              !< options:
                              !< - N -- perturb around the model
                              !< initial conditions with random noise
                              !< distributed \f$\mathcal{N}(0,I)\f$
                              !< - P -- perturb around the model
                              !< initial conditions with random noise
                              !< distributed \f$\mathcal{N}(0,Q)\f$
                              !< - B -- perturb around the model
                              !< initial conditions with random noise
                              !< distributed \f$\mathcal{N}(0,B)\f$
                              !< - R -- read model states from
                              !< rstrt folder where each ensemble member
                              !< is stored in the file rstrt/##.state
                              !< - S -- read model states from
                              !< start folder where each ensemble member
                              !< is stored in the file start/##.state
                              !< - U -- call user defined
                              !< perturbation routine. This assumes
                              !< the user has implemented their own
                              !< perturbation in \link
                              !< user_perturb_particle \endlink
                              !< - Z -- do not perturb particles.
                              !< This will assume each model is 
                              !< received with initial spread
     
  end type pf_control_type
  type(pf_control_type), save :: pf !< the derived data type holding all controlling data

  
contains
  !> subroutine to ensure pf_control data is ok
  subroutine set_pf_controls
    use output_empire, only : emp_o
    
    write(emp_o,'(A)') 'Opening empire namelist file:'
    
    call parse_pf_parameters
    
    pf%efac = 0.001/pf%nens
    write(emp_o,'(A)') 'empire namelist successfully read to control pf code.'
    call flush(emp_o)
    
    
  end subroutine set_pf_controls


    !>subroutine to read the namelist file and save it to pf datatype
    !!Here we read pf_parameters.dat or empire.nml
    !!
    !! pf_parameters.dat or empire.nml is a fortran namelist file. As such, within
    !! it there must be a line beginning
    !!
    !! &pf_params
    !!
    !! To make it (probably) work, ensure there is a forward slash on
    !! the penultimate line and a blank line to end the file
    !!
    !! This is just the fortran standard for namelists though.
    !!
    !!
    !! On to the content...in any order, the pf_parameters.dat (or
    !! empire.nml) file may
    !! contain the following things:
    !! 
    !! Integers:
    !! - \link pf_control::pf_control_type::time_obs time_obs \endlink
    !! - \link pf_control::pf_control_type::time_bwn_obs
    !! time_bwn_obs\endlink
    !!
    !! Reals, double precision:
    !! - \link pf_control::pf_control_type::nudgefac nudgefac\endlink
    !! - \link pf_control::pf_control_type::nfac nfac\endlink
    !! - \link pf_control::pf_control_type::ufac ufac\endlink
    !! - \link pf_control::pf_control_type::qscale Qscale \endlink
    !! - \link pf_control::pf_control_type::keep keep  \endlink
    !! - \link pf_control::pf_control_type::rho rho  \endlink
    !! - \link pf_control::pf_control_type::len len  \endlink
    !!
    !! 2 Characters:
    !! - \link pf_control::pf_control_type::filter filter\endlink
    !!
    !! 1 Character:
    !! - \link pf_control::pf_control_type::init init\endlink
    !!
    !! Logicals:
    !! - \link pf_control::pf_control_type::gen_q gen_Q\endlink
    !! - \link pf_control::pf_control_type::gen_data gen_data\endlink
    !! - \link pf_control::pf_control_type::use_talagrand use_talagrand\endlink
    !! - \link pf_control::pf_control_type::use_mean use_mean\endlink
    !! - \link pf_control::pf_control_type::use_variance use_variance\endlink
    !! - \link pf_control::pf_control_type::use_traj use_traj\endlink
    !! - \link pf_control::pf_control_type::use_spatial_rmse use_spatial_rmse\endlink
    !! - \link pf_control::pf_control_type::output_forecast output_forecast\endlink
    !! - \link pf_control::pf_control_type::use_ens_rmse use_ens_rmse\endlink
    !! - \link pf_control::pf_control_type::output_weights output_weights\endlink
    !!
    !! 250 Character string:
  !! - \link pf_control::pf_control_type::rmse_filename rmse_filename\endlink
      !! - \link pf_control::pf_control_type::forecast_path forecast_path\endlink
    !!

    subroutine parse_pf_parameters
      use var_data
      use output_empire, only : unit_nml
      implicit none
      integer :: ios

      integer :: time_obs=-1
      integer :: time_bwn_obs=-1 
      real(kind=kind(1.0D0)) :: nudgefac=-1.0d0
      logical :: gen_data,gen_Q
      real(kind=kind(1.0D0)) :: nfac=-1.0d0                
      real(kind=kind(1.0D0)) :: ufac=-1.0d0
      real(kind=kind(1.0D0)) :: Qscale=-1.0d0
      real(kind=kind(1.0D0)) :: rho=0.0d0
      real(kind=kind(1.0d0)) :: len=-1.0d0
      real(kind=kind(1.0D0)) :: keep
      logical :: use_talagrand,use_mean,use_variance,use_traj&
           &,use_spatial_rmse,use_ens_rmse,output_weights,output_forecast
      character(2) :: filter='++'
      character(1) :: init='+'
      character(250) :: rmse_filename='rmse'
      character(250) :: forecast_path='forecast/'
      
      logical :: file_exists

      namelist/pf_params/time_obs,time_bwn_obs,&
      &nudgefac,& 
      &gen_data,gen_Q,&
      &nfac,&
      &keep,&
      &ufac,&                
      &Qscale,&
      &rho,&
      &len,&
      &use_talagrand,use_mean,use_variance,use_traj,use_spatial_rmse,use_ens_rmse,&
      &output_weights,&
      &output_forecast,&
      &filter,&
      &init,&
      &rmse_filename,&
      &forecast_path


      gen_data = .false.
      gen_Q = .false.
      use_talagrand = .false.
      use_mean = .false.
      use_variance = .false.
      use_traj = .false.
      use_spatial_rmse = .false.
      use_ens_rmse = .false.
      output_weights=.false.
      output_forecast=.false.


      inquire(file='pf_parameters.dat',exist=file_exists)
      if(file_exists) then
         open(unit_nml,file='pf_parameters.dat',iostat=ios,action='read'&
              &,status='old')
         if(ios .ne. 0) stop 'Cannot open pf_parameters.dat'
      else
         inquire(file='empire.nml',exist=file_exists)
         if(file_exists) then
            open(unit_nml,file='empire.nml',iostat=ios,action='read'&
                 &,status='old')
            if(ios .ne. 0) stop 'Cannot open empire.nml'
         else
            print*,'ERROR: cannot find pf_parameters.dat or empire.nml'
            stop '-1'
         end if
      end if
         
      read(unit_nml,nml=pf_params) 
      close(unit_nml)

      if(time_obs .gt. -1) then
         print*,'read time_obs = ',time_obs
         pf%time_obs = time_obs
      end if
      if(time_bwn_obs .gt. -1) then
         print*,'read time_bwn_obs = ',time_bwn_obs
         pf%time_bwn_obs = time_bwn_obs
      end if
      if(nudgefac .gt. -1.0d0) then
         print*,'read nudgefac = ',nudgefac
         pf%nudgefac = nudgefac
      end if

      if(nfac .gt. -1.0d0) then
         print*,'read nfac = ',nfac
         pf%nfac = nfac
      end if

      if(keep .gt. -1.0d0) then
         print*,'read keep = ',keep
         pf%keep = keep
      end if

      if(ufac .gt. -1.0d0) then
         print*,'read ufac = ',ufac
         pf%ufac = ufac
      end if

      !real(kind=kind(1.0D0)) :: ufac=-1.0d0
      if(Qscale .gt. -1.0d0) then
         print*,'read Qscale = ',Qscale
         pf%Qscale = Qscale
      end if
      

      if(rho .gt. 0.0d0) then
         print*,'read rho = ',rho
         pf%rho = rho
      elseif(rho .lt. 0.0d0) then
         print*,'read rho = ',rho,' WARNING ABOUT THAT ONE! rho is nor&
              &mally positive'
         pf%rho = rho
      else
         pf%rho = rho
      end if


      if(len .ge. 0.0d0) then
         print*,'read len = ',len
         pf%len = len
       else
         pf%len = len
      end if


      !logical ::
      !use_talagrand,use_mean,use_variance,use_traj,use_spatial_rmse
      
      if(use_spatial_rmse) then
         print*,'going to output spatial Root mean squared errors'
      end if

      if(use_ens_rmse) then
         print*,'going to output field of Root mean squared errors'
      end if
      
      if(use_mean) then
         print*,'going to output ensemble mean'
      end if

      if(use_variance) then
         print*,'going to output ensemble variance'
      end if
      
      if(output_weights) then
         print*,'going to output ensemble weights'
      end if

      if(output_forecast) then
         print*,'going to output ensemble forecasts'
      end if

      
      if(use_traj) then
         print*,'going to output trajectories'
      end if


      if(rmse_filename .ne. 'rmse') then
         print*,'read rmse_filename = ',rmse_filename
      end if
      pf%rmse_filename = rmse_filename

      if(forecast_path .ne. 'forecast/') then
         print*,'read forecast_path = ',forecast_path
      end if
      pf%forecast_path = forecast_path
      
      
      if(filter .ne. '++') then
         print*,'read filter = ',filter
         pf%filter = filter
      end if


      !ensure that if we are generating the data then SE is selected
      if(gen_data) then
         pf%filter = 'SE'         
      end if



      !let us verify pf%filter
      select case(pf%filter)
      case('EW')
         print*,'Running the equivalent weights particle filter'
      case('EZ')
         print*,'Running the Zhu equal weights particle filter'
      case('SE')
         print*,'Running a stochastic ensemble'
      case('DE')
         print*,'Running a deterministic ensemble'
      case('SI')
         print*,'Running the SIR particle filter'
      case('ET')
         print*,'filter read as ET. This is depreciated. Changing to L&
              &E'
         pf%filter = 'LE'
         print*,'Running the Local Ensemble Transform Kalman Filter'
         print*,'With random noise'
         print*,'For LETKF without random noise set filter="LD"'
      case('EA')
         print*,'Running the Ensemble Adjustment Kalman Filter'
         print*,'Error: The EAKF is not implemented here yet'
         stop '-557'
      case('LE')
         print*,'Running the Local Ensemble Transform Kalman Filter'
         print*,'With random noise'
      case('LS')
         print*,'Running the Local Ensemble Transform Kalman Smoother'
         print*,'With random noise'
      case('LD')
         print*,'Running the Local Ensemble Transform Kalman Filter'
         print*,'With NO random noise'
      case('3D')
         print*,'Running a stochastic ensemble and 3DVar at observatio&
              &n times'
         call set_var_controls
      case default
         print*,'Error: Incorrect filter type selected:', pf%filter
         print*,'Please ensure that pf%filter in empire.nml is either:'
         print*,'EW        the equivalent weights particle f&
              &ilter'
         print*,'EZ        the Zhu equal weights particle f&
              &ilter'
         print*,'SE        a stochastic ensemble'
         print*,'DE        a deterministic ensemble'
         print*,'SI        the SIR particle filter'
         print*,'LE        the Local Ensemble Transform Kalman Filter'
         print*,'          with random noise'
         print*,'LD        the Local Ensemble Transform Kalman Filter'
         print*,'          without random noise'
         print*,'LS        the Local Ensemble Transform Kalman Smoother'
         print*,'          with random noise'
         print*,'3D        3DVar'
         stop
      end select

      
      if(init .ne. '+') then
         print*,'read init = ',init
         pf%init = init
      end if

      pf%gen_data = gen_data
      pf%gen_Q = gen_Q
      pf%use_talagrand = use_talagrand
      pf%use_mean = use_mean
      pf%use_variance = use_variance
      pf%use_traj = use_traj
      pf%use_spatial_rmse = use_spatial_rmse
      pf%use_ens_rmse = use_ens_rmse
      pf%output_weights = output_weights
      pf%output_forecast = output_forecast
      
      !check for specific things we can't do if running the truth
      if(pf%gen_data) then
         pf%gen_Q = .false.
         pf%use_talagrand=.false.
         pf%use_spatial_rmse = .false.
         pf%use_ens_rmse = .false.
         pf%output_forecast = .false.
      end if

      
    end subroutine parse_pf_parameters


    
    !> subroutine to deallocate space for the filtering code
    subroutine deallocate_pf
      deallocate(pf%weight)
      deallocate(pf%psi)
      if(allocated(pf%talagrand)) deallocate(pf%talagrand)
      deallocate(pf%particles)
    end subroutine deallocate_pf


end module pf_control
