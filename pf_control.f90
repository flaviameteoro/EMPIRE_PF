module pf_control
  implicit none
  type, public :: pf_control_type
     integer :: ngrand !the number of ensemble members
     real(kind=kind(1.0D0)), allocatable, dimension(:) :: weight !stores the weights of the particles
     integer :: time_obs !the number of observations we will assimilate
     integer :: time_bwn_obs !the number of model timesteps between observations
     real(kind=kind(1.0D0)) :: nudgefac !the nudging factor
     logical :: gen_data,human_readable
     integer :: timestep=0
     real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: psi
     real(kind=kind(1.0D0)) :: nfac                !standard deviation of normal distribution in mixture density
     real(kind=kind(1.0D0)) :: ufac                !half width of the uniform distribution in mixture density
     real(kind=kind(1.0D0)) :: efac
     real(kind=kind(1.0D0)) :: keep
     integer :: couple_root
     logical :: use_talagrand,use_weak,use_mean
     integer, dimension(:,:), allocatable :: talagrand
     integer :: tala_obs_num
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
      read(32,*) pf%human_readable
      read(32,*) pf%use_talagrand
      read(32,*) pf%use_weak
      read(32,*) pf%tala_obs_num
      read(32,*) pf%use_mean
      close(32)
      pf%efac = 0.001/pf%ngrand
      write(6,'(A)') 'pf_parameters.dat successfully read to control pf code.'
      call flush(6)
      if(pf%human_readable .and. pf%gen_data) then
         open(64,file='pf_data',iostat=ios,action='read',status='replace')
         if(ios .ne. 0) stop 'Error checking pf_data'
         close(64)
      end if

    end subroutine set_pf_controls

    subroutine allocate_pf
      use sizes
      integer :: st
      allocate(pf%weight(pf%ngrand),stat=st)
      if(st .ne. 0) stop 'Error in allocating pf%weight'
      pf%weight = -log(1.0D0/pf%ngrand)
      allocate(pf%psi(state_dim,pf%ngrand),stat=st)
      if(st .ne. 0) stop 'Error in allocating pf%psi'

      if(pf%use_talagrand) then
         allocate(pf%talagrand(pf%ngrand+1,1),stat=st)
         if(st .ne. 0) stop 'Error in allocating pf%talagrand'
         do st = 1,pf%ngrand+1
            pf%talagrand(st,:) = 0
         end do
      end if

    end subroutine allocate_pf

    subroutine deallocate_pf
      deallocate(pf%weight)
      deallocate(pf%psi)
      if(allocated(pf%talagrand)) deallocate(pf%talagrand)
    end subroutine deallocate_pf


end module pf_control
