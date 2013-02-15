module pf_control
  implicit none
  type, public :: pf_control_type
     integer :: ngrand !the number of ensemble members
     real(kind=kind(1.0D0)), allocatable, dimension(:) :: weight !stores the weights of the particles
     integer :: time_obs !the number of observations we will assimilate
     integer :: time_bwn_obs !the number of model timesteps between observations
     real(kind=kind(1.0D0)) :: nudgefac !the nudging factor
     logical :: gen_data
     integer :: timestep=0
     real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: psi
     real(kind=kind(1.0D0)) :: nfac                !standard deviation of normal distribution in mixture density
     real(kind=kind(1.0D0)) :: ufac                !half width of the uniform distribution in mixture density
     real(kind=kind(1.0D0)) :: efac
  end type pf_control_type
  type(pf_control_type) :: pf
  contains
    subroutine set_pf_controls
      integer, parameter :: rk=kind(1.0D0)
      integer :: ios
      write(6,'(A)') 'Opening pf_parameters.dat'
      open(32,file='pf_parameters.dat',iostat=ios,action='read',status='old')
      if(ios .ne. 0) stop 'Cannot open pf_parameters.dat'

      read(32,*) pf%ngrand
      read(32,*) pf%time_obs
      read(32,*) pf%time_bwn_obs
      read(32,*) pf%nudgefac
      read(32,*) pf%gen_data
!      read(32,*) pf%nfac
!      read(32,*) pf%ufac
      close(32)
      pf%efac = 0.001/pf%ngrand
      write(6,'(A)') 'pf_parameters.dat successfully read to control pf code.'
      call flush(6)
    end subroutine set_pf_controls

    subroutine allocate_pf
      use sizes
      integer :: st
      allocate(pf%weight(pf%ngrand),stat=st)
      if(st .ne. 0) stop 'Error in allocating pf%weight'
      allocate(pf%psi(state_dim,pf%ngrand),stat=st)
      if(st .ne. 0) stop 'Error in allocating pf%weight'
    end subroutine allocate_pf

    subroutine deallocate_pf
      deallocate(pf%weight)
    end subroutine deallocate_pf


end module pf_control
