!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-08-16 14:57:26 pbrowne>
!!!
!!!    Module and subroutine to intitalise EMPIRE coupling to models
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


!> Module containing EMPIRE coupling data
!! @todo Need to see what happens if some process has no observations
!! in comms_v3
module comms
  integer :: CPL_MPI_COMM !< the communicator between the empire
  !< codes and the model master nodes
  integer :: world_rank !< the rank of this process on MPI_COMM_WORLD
  integer :: cpl_rank !< the rank of this process on CPL_MPI_COMM
  integer :: nProc !< the total number of processes
  integer :: pf_mpi_comm !< the communicator between DA processes
  integer :: pfrank      !< the rank of this process on PF_MPI_COMM
  integer :: npfs        !< the total number of DA processes
  integer, allocatable, dimension(:) :: gblcount !< the number of
  !< ensemble members associated with each DA process
  integer, allocatable, dimension(:) :: gbldisp !< the displacements
  !< of each each ensemble member relative to pfrank=0. VERY useful
  !< for mpi_gatherv and mpi_scatterv on pf_mpi_comm
  integer :: nens !< the total number of ensemble members
  integer :: cnt !< the number of ensemble members associated with
  !<this process
  integer, allocatable, dimension(:) :: particles !< the ensemble
  !!members associated with this process
  integer, allocatable, dimension(:) :: cpl_mpi_comms !<communicators
  !! for if we are using empire v2 or v3
  integer, allocatable, dimension(:) :: state_dims !<state dimensions
  !!on each model process for empire v2
  integer, allocatable, dimension(:) :: state_displacements
                                !< displacements of the various parts 
                                !! of the state vector for empire v2
  integer, allocatable, dimension(:) :: obs_dims !<obs dimensions
                                !! on each model process for empire v3
  integer, allocatable, dimension(:) :: obs_displacements
                                !< displacements of the various parts 
                                !! of the obs vector for empire v3
  integer :: mdl_num_proc !< number of processes of each ensemble
                          !! member
  integer :: pf_member_comm !< communicator for empire v3 which
                            !! contains all processes of individual
                            !! ensemble members
  integer :: pf_ens_comm    !< communicator for empire v3 which
                            !! contains all ensemble members for that
                            !! specific part of the state vector
  integer :: pf_ens_rank    !< rank of the process on pf_ens_comm
  integer :: pf_ens_size    !< size of pf_ens_comm for comms v3
  integer :: pf_member_rank !< rank of the process on pf_member_comm
                            !! for empire v3
  integer :: pf_member_size !< size of pf_member_comm
                            !! for empire v3
  integer, parameter :: comm_version=1 !< The style of communication
  !! between the model and empire.
  !! 
  !! - 1 = MPI SEND/RECV pairs between a single model process (single
  !! EMPIRE process per ensemble member)
  !!
  !! - 2 = MPI GATHERV/SCATTERV between (possibly) multiple model processes
  !! (single EMPIRE process per ensemble member)
  !!
  !! - 3 = MPI SEND/RECV pairs between multiple model processes and
  !! the same parallel process disribution in EMPIRE
  !! 
  !! - 4 = MODEL AS A SUBROUTINE OF EMPIRE @todo Fully document how
  !! to specify the model_as_subroutine calls in src/user/model 
  !!
  !! - 5 = Similar to 2, but with multiple ensemble members for each
  !! model process (TOMCAT CASE) 
contains

  subroutine allocate_data

    implicit none

  end subroutine allocate_data

  subroutine deallocate_data
    implicit none

  end subroutine deallocate_data

  !> subroutine to select which mpi comms to use
  subroutine initialise_mpi
    implicit none

    select case(comm_version)
    case(0)
       call user_initialise_mpi
    case(1)
       call initialise_mpi_v1
    case(2)
       call initialise_mpi_v2
    case(3)
       call initialise_mpi_v3
    case(4)
       call initialise_mpi_v4
    case(5)
       call initialise_mpi_v5
    case default
       print*,'ERROR: comm_version ',comm_version,' not implemente&
            &d.'
       print*,'STOPPING.'
       stop '-6'
    end select

  end subroutine initialise_mpi

  !> subroutine to make EMPIRE connections and saves details into
  !! pf_control module
  subroutine initialise_mpi_v1

    use pf_control
    implicit none
    include 'mpif.h'

    integer :: mpi_err
    integer :: couple_colour
    integer :: myrank
    integer :: i
    integer :: da,count
    integer :: pf_colour
    integer :: world_size

    if(comm_version .eq. 2) then
       call initialise_mpi_v2
       return
    end if

    pf_colour = 10000
    couple_colour=9999
    call MPI_INIT (mpi_err)

    da = 1
    call MPI_COMM_RANK (MPI_COMM_WORLD,world_rank,   mpi_err)
    call mpi_comm_size (mpi_comm_world,world_size,   mpi_err)
    call mpi_comm_split(mpi_comm_world,da,           world_rank,  pf_mpi_comm, mpi_err)
    call mpi_comm_rank (pf_mpi_comm,   pfrank,       mpi_err)
    call mpi_comm_size (pf_mpi_comm,   npfs,         mpi_err)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,couple_colour,world_size,CPL_MPI_COMM,mpi_err)
    call MPI_COMM_RANK (CPL_MPI_COMM,  myRank,       mpi_err)
    call MPI_COMM_SIZE (CPL_MPI_COMM,  nens,         mpi_err)

    nens = nens-npfs
    print*,'DA'
    print*,'nens = ',nens
    print*,'npfs = ',npfs


    !lets find the particles:

    count = ceiling(real((myrank-nens+1)*nens)/real(npfs)) -&
         & ceiling(real((myrank-nens)*nens)/real(npfs))

    allocate(pf%particles(count))
    allocate(   particles(count))

    pf%particles = (/ (i, i = ceiling(real((myrank-nens)*nens)&
         &/real(npfs)),&
         ceiling(real((myrank-nens+1)*nens)/real(npfs))-1) /)+1
    particles = pf%particles

    allocate(gblcount(npfs))
    allocate(gbldisp(npfs))


    call mpi_allgather(count,1,MPI_INTEGER,gblcount,1,MPI_INTEGER&
         &,pf_mpi_comm,mpi_err)
    

    gbldisp = 0
    if(npfs .gt. 1) then
       do i = 2,npfs
          gbldisp(i) = gbldisp(i-1) + gblcount(i-1)
       end do
    end if
    pf%count = count
    cnt = count

    pf%nens = nens
    PRINT*,'PF_rank = ',pfrank,' and I own particles ',pf%particles

    pf_ens_comm=pf_mpi_comm
    pf_ens_size=npfs
    pf_ens_rank=pfrank
    pf_member_rank=0
    pf_member_size=1
  end subroutine initialise_mpi_v1


  !> subroutine to initialise new version of empire
  subroutine initialise_mpi_v2
    use pf_control
    use sizes
    implicit none
    include 'mpif.h'

    integer :: mpi_err
    integer :: world_size
    integer,parameter :: da=1
    integer, parameter :: rk = kind(1.0d0)
    integer :: i
    integer :: mdl_procs
    integer :: first_ptcl
    integer :: final_ptcl
    integer :: tmp_cpl_comm


    call mpi_init(mpi_err)
    call mpi_comm_rank(mpi_comm_world,world_rank,mpi_err)
    call mpi_comm_size(mpi_comm_world,world_size,mpi_err)
    print*,'EMPIRE:  rank = ',world_rank,' on mpi_comm_world which has size ',world_size


    !get the number of processes per model
    call mpi_allreduce(0,mdl_num_proc,1,MPI_INTEGER,MPI_MAX&
         &,MPI_COMM_WORLD,mpi_err)

    if(mdl_num_proc .lt. 1) then
       print*,'EMPIRE COMMS v2 ERROR: mdl_num_proc < 1'
       print*,'mdl_num_proc = ',mdl_num_proc
       print*,'THIS SUGGESTS YOU HAVE NOT LINKED TO A MODEL. STOP.'
       stop
    else
       print*,'mdl_num_proc = ',mdl_num_proc
    end if

    !split into models and da processes. create pf_mpi_comm
    call mpi_comm_split(MPI_COMM_WORLD,da,world_rank,pf_mpi_comm,mpi_err)
    call mpi_comm_size(pf_mpi_comm,npfs,mpi_err)
    call mpi_comm_rank(pf_mpi_comm,pfrank,mpi_err)
    
    !compute number of model processes
    mdl_procs = world_size-npfs
    print*,'npfs = ',npfs

    print*,'mdl_procs = ',mdl_procs


    !compute number of ensemble members
    nens = mdl_procs/mdl_num_proc
    print*,'nens = ',nens


    !compute range of particles that this mpi process communicates with
    first_ptcl = ceiling(real(pfrank)*real(nens)/real(npfs))
    final_ptcl = ceiling(real(pfrank+1)*real(nens)/real(npfs))-1
    particles = (/ (i, i = first_ptcl,final_ptcl) /)
    print*,'range of particles = ',first_ptcl,final_ptcl


    !create a temporary communicator with all associated model processes
    call mpi_comm_split(mpi_comm_world,pfrank,world_size+pfrank&
         &,tmp_cpl_comm,mpi_err)
    print*,'split and created tmp_cpl_comm'


    ! count the number of particles associated with this process
    cnt = final_ptcl-first_ptcl+1
    if(cnt .lt. 1) then
       print*,'EMPIRE ERROR: YOU HAVE LAUNCHED MORE EMPIRE DA PROCESSES'
       print*,'EMPIRE ERROR: THAN MODELS. I AM REDUDANT AND STOPPING.'
       print*,'EMPIRE ERROR: RECONSIDER HOW YOU EXECUTE NEXT TIME. xx'
       stop
    end if


    !allocate a communicator for each ensemble member
    allocate(cpl_mpi_comms(cnt))


    !split the temporary communicator into individual ones for each
    !ensemble member
    do i = 1,cnt
       call mpi_comm_split(tmp_cpl_comm,1,world_size,cpl_mpi_comms(i)&
            &,mpi_err)
       write(*,'(A,i3.3,A)') 'created cpl_mpi_comms(',i,')'
    end do


    !the rank of this mpi process each of cpl_mpi_comms(:) is the highest
    cpl_rank = mdl_num_proc

    
    !free up the temporary communicator
    call mpi_comm_free(tmp_cpl_comm,mpi_err)


    !allocate space to get the sizes from each model process
    allocate(state_dims(mdl_num_proc+1))
    allocate(state_displacements(mdl_num_proc+1))

    
    !set the sizes to zero, specifically the final entry in the arrays
    state_dims = 0
    state_dim = 0

    
    print*,'doing a gather on cpl_mpi_comm'
    !for each ensemble member, gather the number of variables stored
    !on each process
    do i = 1,cnt
       call mpi_gather(state_dim,1,MPI_INTEGER,state_dims&
            &,1,MPI_INTEGER,cpl_rank,cpl_mpi_comms(i),mpi_err)
    end do
    print*,'state_dims = ',state_dims

    
    !compute the relevant displacements for the gathers
    state_displacements = 0
    do i = 2,mdl_num_proc
       state_displacements(i:) = state_displacements(i:) + state_dims(i-1)
    end do
    print*,'state_displacements = ',state_displacements

    
    !compute the total size of the state
    state_dim = sum(state_dims)
    print*,'total state_dim = ',state_dim


    !compute counts and displacements of particles associated with da
    !processes
    allocate(gblcount(npfs))
    allocate(gbldisp(npfs))


    call mpi_allgather(cnt,1,MPI_INTEGER,gblcount,1,MPI_INTEGER&
         &,pf_mpi_comm,mpi_err)
    if(mpi_err .eq. 0) then
       print*,'mpi_allgather successful: gblcount known on all da proc&
            &esses'
       print*,'gblcount = ',gblcount
    else
       print*,'mpi_allgather unsucessful'
    end if
    

    gbldisp = 0
    if(npfs .gt. 1) then
       do i = 2,npfs
          gbldisp(i) = gbldisp(i-1) + gblcount(i-1)
       end do
    end if    





    pf_ens_comm=pf_mpi_comm
    pf_ens_size=npfs
    pf_ens_rank=pfrank
    pf_member_rank=0
    pf_member_size=1

    
    pf%particles = particles+1
    pf%count = cnt
    pf%nens = nens
  end subroutine initialise_mpi_v2

  !> subroutine to initialise even newer version of empire
  subroutine initialise_mpi_v3
    use sizes
    use pf_control
    implicit none
    include 'mpif.h'

    integer :: mpi_err
    integer :: world_size
    integer,parameter :: da=1
    integer, parameter :: rk = kind(1.0d0)
    integer :: i
    integer :: mdl_procs
    integer :: first_ptcl
    integer :: final_ptcl
    integer :: tmp_cpl_comm
    integer :: tmp_cpl_comm2
    integer :: tmp_cpl_rank
    integer :: tmp_cpl_colour2
    integer :: pf_member_colour
    integer :: pf_ens_colour
    integer :: status(MPI_STATUS_SIZE)


    call mpi_init(mpi_err)
    call mpi_comm_rank(mpi_comm_world,world_rank,mpi_err)
    call mpi_comm_size(mpi_comm_world,world_size,mpi_err)
    print*,'EMPIRE:  rank = ',world_rank,' on mpi_comm_world which has size ',world_size


    !get the number of processes per model
    call mpi_allreduce(0,mdl_num_proc,1,MPI_INTEGER,MPI_MAX&
         &,MPI_COMM_WORLD,mpi_err)

    if(mdl_num_proc .lt. 1) then
       print*,'EMPIRE COMMS v3 ERROR: mdl_num_proc < 1'
       print*,'mdl_num_proc = ',mdl_num_proc
       print*,'THIS SUGGESTS YOU HAVE NOT LINKED TO A MODEL. STOP.'
       stop
    else
       print*,'mdl_num_proc = ',mdl_num_proc
    end if
!!!==================================================!!!
!!!================ BRUCE LEE =======================!!!
!!!==================================================!!!
    

    !split into models and da processes. create pf_mpi_comm
    call mpi_comm_split(MPI_COMM_WORLD,da,world_rank,pf_mpi_comm,mpi_err)
    call mpi_comm_size(pf_mpi_comm,npfs,mpi_err)
    call mpi_comm_rank(pf_mpi_comm,pfrank,mpi_err)
    
!!!==================================================!!!
!!!================ JEAN CLAUDE VAN DAMME ===========!!!
!!!==================================================!!!
    

    !split pf_mpi_comm into communicator for ensemble member
    pf_member_colour= pfrank/mdl_num_proc
    call mpi_comm_split(pf_mpi_comm,pf_member_colour,pfrank&
         &,pf_member_comm,mpi_err)
    call mpi_comm_rank(pf_member_comm,pf_member_rank,mpi_err)
    
    !split pf_mpi_comm into communicator for ensemble (local to this
    !part of the state vector)
    pf_ens_colour = mod(pfrank,mdl_num_proc)
    call mpi_comm_split(pf_mpi_comm,pf_ens_colour,pfrank,pf_ens_comm&
         &,mpi_err)
    call mpi_comm_size(pf_ens_comm,pf_ens_size,mpi_err)
    call mpi_comm_rank(pf_ens_comm,pf_ens_rank,mpi_err)

    
    

    !compute number of model processes
    mdl_procs = world_size-npfs
    npfs = pf_ens_size

    print*,'npfs = ',npfs
    print*,'mdl_procs = ',mdl_procs


    !compute number of ensemble members
    nens = mdl_procs/mdl_num_proc
    print*,'nens = ',nens


    !compute range of particles that this mpi process communicates with
    first_ptcl = ceiling(real(pf_member_colour)*real(nens)/real(pf_ens_size))
    final_ptcl = ceiling(real(pf_member_colour+1)*real(nens)/real(pf_ens_size))-1
    particles = (/ (i, i = first_ptcl,final_ptcl) /)
    print*,'range of particles = ',first_ptcl,final_ptcl
   

    ! count the number of particles associated with this process
    cnt = final_ptcl-first_ptcl+1
    if(cnt .lt. 1) then
       print*,'EMPIRE ERROR: YOU HAVE LAUNCHED MORE EMPIRE DA PROCESSES'
       print*,'EMPIRE ERROR: THAN MODELS. I AM REDUDANT AND STOPPING.'
       print*,'EMPIRE ERROR: RECONSIDER HOW YOU EXECUTE NEXT TIME. xx'
       stop
    end if

    !!!==================================================!!!
    !!!============== STEVEN SEGAL ======================!!!
    !!!==================================================!!!

    !create a temporary communicator with all associated model processes
    call mpi_comm_split(mpi_comm_world,pf_member_colour,world_size+pf_member_colour&
         &,tmp_cpl_comm,mpi_err)
    print*,'split and created tmp_cpl_comm'
    
!!!==================================================!!!
!!!============== CHUCK NORRIS ======================!!!
!!!==================================================!!! 

    !get the rank and set the colour across the model mpi processes
    call mpi_comm_rank(tmp_cpl_comm,tmp_cpl_rank,mpi_err)
    tmp_cpl_colour2 = mod(tmp_cpl_rank,mdl_num_proc)
    
    !split this temp communicator into a new temporary one
    call mpi_comm_split(tmp_cpl_comm,tmp_cpl_colour2,tmp_cpl_rank&
         &,tmp_cpl_comm2,mpi_err)
  
!!!==================================================!!!
!!!=============== JACKIE CHAN ======================!!!
!!!==================================================!!!

    !allocate a communicator for each ensemble member
    allocate(cpl_mpi_comms(cnt))


    !split the second temporary communicator into individual ones for each
    !ensemble member
    do i = 1,cnt
       call mpi_comm_split(tmp_cpl_comm2,1,world_size,cpl_mpi_comms(i)&
            &,mpi_err)
       write(*,'(A,i3.3,A,i0)') 'created cpl_mpi_comms(',i,') ',cpl_mpi_comms(i)
    end do

    !the rank of this mpi process each of cpl_mpi_comms(:) is the highest
    cpl_rank = 1

    
!!!==================================================!!!
!!!=============== MICHELLE YEOH ====================!!!
!!!==================================================!!! 
    
    !free up the temporary communicators
    call mpi_comm_free(tmp_cpl_comm ,mpi_err)
    call mpi_comm_free(tmp_cpl_comm2,mpi_err)
    
    state_dim = 0

    
    !for each ensemble member, gather the number of variables stored
    !on each process
    do i = 1,cnt
       call mpi_recv(state_dim,1,MPI_INTEGER,0,MPI_ANY_TAG&
            &,cpl_mpi_comms(i),status,mpi_err)
    end do
    print*,'state_dim = ',state_dim
    
!!!==================================================!!!
!!!=============== SCARLETT JOHANSSON ===============!!!
!!!==================================================!!! 

    !compute the total state dimension:
    call mpi_allreduce(state_dim,state_dim_g,1,MPI_INTEGER,MPI_SUM&
         &,pf_member_comm,mpi_err)
    print*,'total state vector size = ',state_dim_g

   
    !compute counts and displacements of particles associated with da
    !processes
    allocate(gblcount(pf_ens_size))
    allocate(gbldisp(pf_ens_size))


    call mpi_allgather(cnt,1,MPI_INTEGER,gblcount,1,MPI_INTEGER&
         &,pf_ens_comm,mpi_err)
    if(mpi_err .eq. 0) then
       print*,'mpi_allgather successful: gblcount known on all da proc&
            &esses'
       print*,'gblcount = ',gblcount
    else
       print*,'mpi_allgather unsucessful'
    end if
    

    gbldisp = 0
    if(npfs .gt. 1) then
       do i = 2,pf_ens_size
          gbldisp(i) = gbldisp(i-1) + gblcount(i-1)
       end do
    end if    






    pf%particles = particles+1
    pf%count = cnt
    pf%nens = nens
  end subroutine initialise_mpi_v3


  !> subroutine to initialise empire communicators when the model is
  !! to be a subroutine itself
  subroutine initialise_mpi_v4
    use pf_control
    use sizes
    use output_empire, only : unit_nml
    implicit none
    include 'mpif.h'

    integer :: mpi_err
    integer :: world_size
    integer, parameter :: rk = kind(1.0d0)
    integer :: i
    integer :: first_ptcl
    integer :: final_ptcl
    integer :: ios
    logical :: file_exists
    namelist/comms_v4/nens

    
    call mpi_init(mpi_err)
    call mpi_comm_rank(mpi_comm_world,world_rank,mpi_err)
    call mpi_comm_size(mpi_comm_world,world_size,mpi_err)
    print*,'EMPIRE:  rank = ',world_rank,' on mpi_comm_world which has size ',world_size


    npfs = world_size
    pfrank = world_rank
    
    !now need to get the total number of ensemble members that will
    !be used. let us read it in from empire.nml!
    
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
    ! set nens to be negative as a default value
    nens = -1
    !now read it in
    read(unit_nml,nml=comms_v4) 
    close(unit_nml)

    if( nens .lt. 1 ) then
       print*,'EMPIRE ERROR: __________initialise_mpi_v4_____________'
       print*,'EMPIRE ERROR: nens is less than 1... nens = ',nens
       print*,'EMPIRE ERROR: please correctly specify this in empire.n&
            &ml'
       stop '-1'
    end if

    if (npfs .gt. nens) then
       print*,'EMPIRE ERROR: __________initialise_mpi_v4_____________'
       print*,'EMPIRE ERROR: npfs is great than nens...'
       print*,'EMPIRE ERROR: npfs = ',npfs,' nens = ',nens
       stop '-1'
    end if

    !compute range of particles that this mpi process communicates with
    first_ptcl = ceiling(real(pfrank)*real(nens)/real(npfs))
    final_ptcl = ceiling(real(pfrank+1)*real(nens)/real(npfs))-1
    particles = (/ (i, i = first_ptcl,final_ptcl) /)
    print*,'range of particles = ',first_ptcl,final_ptcl

    !set the da communicator:
    pf_mpi_comm = MPI_COMM_WORLD
    pf_ens_comm = MPI_COMM_WORLD
    pf_ens_size=npfs
    pf_ens_rank=pfrank

    pf_member_rank=0
    pf_member_size=1
    
    ! count the number of particles associated with this process
    cnt = final_ptcl-first_ptcl+1

    !compute counts and displacements of particles associated with da
    !processes
    allocate(gblcount(npfs))
    allocate(gbldisp(npfs))

    do i = 1,npfs
       gblcount(i) = ceiling(real(i)*real(nens)/real(npfs)) &
            &- ceiling(real(i-1)*real(nens)/real(npfs))
    end do
    
    gbldisp = 0
    if(npfs .gt. 1) then
       do i = 2,npfs
          gbldisp(i) = gbldisp(i-1) + gblcount(i-1)
       end do
    end if    
    
    pf%particles = particles+1
    pf%count = cnt
    pf%nens = nens
    
  end subroutine initialise_mpi_v4


  !> subroutine to initialise empire communication pattern similarly
  !! to v2 but with multiple ensemble members per model process
  subroutine initialise_mpi_v5
    use pf_control
    use sizes
    implicit none
    include 'mpif.h'

    integer :: mpi_err
    integer :: world_size
    integer, parameter :: da=1
    integer, parameter :: rk = kind(1.0d0)
    integer :: i,j
    integer :: mdl_procs
    integer :: first_ptcl
    integer :: final_ptcl
    integer :: tmp_cpl_comm
    integer :: n_mdl_instances
    integer :: nens_per_instance

    call mpi_init(mpi_err)
    call mpi_comm_rank(mpi_comm_world,world_rank,mpi_err)
    call mpi_comm_size(mpi_comm_world,world_size,mpi_err)
    print*,'EMPIRE:  rank = ',world_rank,' on mpi_comm_world which has size ',world_size


    !get the number of processes per model
    call mpi_allreduce(0,mdl_num_proc,1,MPI_INTEGER,MPI_MAX&
         &,MPI_COMM_WORLD,mpi_err)

    if(mdl_num_proc .lt. 1) then
       print*,'EMPIRE COMMS v5 ERROR: mdl_num_proc < 1'
       print*,'mdl_num_proc = ',mdl_num_proc
       print*,'THIS SUGGESTS YOU HAVE NOT LINKED TO A MODEL. STOP.'
       stop
    else
       print*,'mdl_num_proc = ',mdl_num_proc
    end if

    !get the number of ensemble members per model process
    call mpi_allreduce(0,nens_per_instance,1,MPI_INTEGER,MPI_MAX&
         &,MPI_COMM_WORLD,mpi_err)
    if(nens_per_instance .lt. 1) then
       print*,'EMPIRE COMMS v5 ERROR: nens_per_instance < 1'
       print*,'nens_per_instance = ',nens_per_instance
       print*,'THIS SUGGESTS YOU HAVE NOT LINKED TO A MODEL. STOP.'
       stop
    else
       print*,'nens_per_instance = ',nens_per_instance
    end if
    call  flush(6)


    !split into models and da processes. create pf_mpi_comm
    call mpi_comm_split(MPI_COMM_WORLD,da,world_rank,pf_mpi_comm,mpi_err)
    call mpi_comm_size(pf_mpi_comm,npfs,mpi_err)
    call mpi_comm_rank(pf_mpi_comm,pfrank,mpi_err)
    
    !compute number of model processes
    mdl_procs = world_size-npfs
    print*,'npfs = ',npfs

    print*,'mdl_procs = ',mdl_procs


    !compute number of ensemble members
    n_mdl_instances = mdl_procs/mdl_num_proc
    print*,'n_mdl_instances = ',n_mdl_instances
    
    nens = n_mdl_instances*nens_per_instance
    print*,'nens = ',nens


    !compute range of particles that this mpi process communicates with
    first_ptcl = ceiling(real(pfrank)*real(nens)/real(npfs))
    final_ptcl = ceiling(real(pfrank+1)*real(nens)/real(npfs))-1
    particles = (/ (i, i = first_ptcl,final_ptcl) /)
    print*,'range of particles = ',first_ptcl,final_ptcl


    ! count the number of particles associated with this process
    cnt = final_ptcl-first_ptcl+1
    if(cnt .lt. 1) then
       print*,'EMPIRE ERROR: YOU HAVE LAUNCHED MORE EMPIRE DA PROCESSES'
       print*,'EMPIRE ERROR: THAN MODELS. I AM REDUDANT AND STOPPING.'
       print*,'EMPIRE ERROR: RECONSIDER HOW YOU EXECUTE NEXT TIME. xx'
       stop
    end if


    !allocate a communicator for each ensemble member
    allocate(cpl_mpi_comms(cnt))

    print*,'you have got to this point:',mod(n_mdl_instances,npfs)

    if(mod(n_mdl_instances,npfs) .ne. 0 ) then !number of instances
       !of the model executable doesnt divide the number of empire
       !processes launched. this is difficult to get an efficient
       !splitting so let's just do this sequentially for each
       !ensemble member...
       j = 0
       print*,'EMPIRE V5: SEQUENTIAL SPLITTING ON MPI_COMM_WORLD...'
       print*,'nens = ',nens
       print*,'first_ptcl = ',first_ptcl
       print*,'final_ptcl = ',final_ptcl
       do i = 0,nens-1
          if( i .ge. first_ptcl .and. i .le. final_ptcl) then
             !we should be in this communicator
             j = j + 1
             call mpi_comm_split(MPI_COMM_WORLD,1,world_size,cpl_mpi_comms(j)&
                  &,mpi_err)
             write(*,'(A,i3.3,A)') 'created cpl_mpi_comms(',j,')'
          else
             ! we should not be part of this communicator
             call mpi_comm_split(MPI_COMM_WORLD,0,world_size,tmp_cpl_comm,mpi_err)
             call mpi_comm_free(tmp_cpl_comm,mpi_err)
          end if
       end do
       
       
    else !number of model executables does divide the number of
       !empire processes: that way we can split easily...
       
       !first split based on pfrank
       print*,'doing first split based on pfrank',pfrank
       call mpi_comm_split(MPI_COMM_WORLD,pfrank,world_size,tmp_cpl_comm,mpi_err)
       print*,'finished first split mpi_err = ',mpi_err
       call flush(6)
       do i = 1,cnt
          print*,'i = ',i
          call mpi_comm_split(tmp_cpl_comm,1,world_size,cpl_mpi_comms(i)&
               &,mpi_err)
          write(*,'(A,i3.3,A)') 'created cpl_mpi_comms(',i,')'
       end do
       call mpi_comm_free(tmp_cpl_comm,mpi_err)
    end if
       
    print*,'EMPIRE: all commiunicators made'

    !the rank of this mpi process each of cpl_mpi_comms(:) is the highest
    cpl_rank = mdl_num_proc

    !allocate space to get the sizes from each model process
    allocate(state_dims(mdl_num_proc+1))
    allocate(state_displacements(mdl_num_proc+1))

   
    !set the sizes to zero, specifically the final entry in the arrays
    state_dims = 0
    state_dim = 0

    
    print*,'doing a gather on cpl_mpi_comm'
    print*,'cnt = ',cnt
    !for each ensemble member, gather the number of variables stored
    !on each process
    do i = 1,cnt
       print*,'cpl_mpi_comms(i) = ',cpl_mpi_comms(i),cpl_rank
       call mpi_gather(state_dim,1,MPI_INTEGER,state_dims&
            &,1,MPI_INTEGER,cpl_rank,cpl_mpi_comms(i),mpi_err)
    end do
    print*,'state_dims = ',state_dims

    
    !compute the relevant displacements for the gathers
    state_displacements = 0
    do i = 2,mdl_num_proc
       state_displacements(i:) = state_displacements(i:) + state_dims(i-1)
    end do
    print*,'state_displacements = ',state_displacements

    
    !compute the total size of the state
    state_dim = sum(state_dims)
    print*,'total state_dim = ',state_dim


    !compute counts and displacements of particles associated with da
    !processes
    allocate(gblcount(npfs))
    allocate(gbldisp(npfs))


    call mpi_allgather(cnt,1,MPI_INTEGER,gblcount,1,MPI_INTEGER&
         &,pf_mpi_comm,mpi_err)
    if(mpi_err .eq. 0) then
       print*,'mpi_allgather successful: gblcount known on all da proc&
            &esses'
       print*,'gblcount = ',gblcount
    else
       print*,'mpi_allgather unsucessful'
    end if
    

    gbldisp = 0
    if(npfs .gt. 1) then
       do i = 2,npfs
          gbldisp(i) = gbldisp(i-1) + gblcount(i-1)
       end do
    end if    

    pf_ens_comm = pf_mpi_comm
    pf_ens_size=npfs
    pf_ens_rank=pfrank

    pf_member_rank=0
    pf_member_size=1
    
    pf%particles = particles+1
    pf%count = cnt
    pf%nens = nens
  end subroutine initialise_mpi_v5










  !> subroutine to send all the model states to the models
  subroutine send_all_models(stateDim,nrhs,x,tag)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: stateDim
    integer, intent(in) :: nrhs
    real(kind=kind(1.0d0)), intent(in), dimension(stateDim,nrhs) :: x
    integer, intent(in) :: tag
    integer :: k
    integer :: mpi_err
    integer :: particle
    real(kind=kind(1.0d0)), dimension(0) :: send_null
    select case(comm_version)
       case(0)
          call user_mpi_send(statedim,nrhs,x,tag)
       case(1)
          do k =1,cnt
             particle = particles(k)
             call mpi_send(x(:,k),stateDim,MPI_DOUBLE_PRECISION&
                  &,particle-1,tag,CPL_MPI_COMM,mpi_err)
          end do
       case(2,5)
          do k = 1,cnt
             call mpi_scatterv(x(:,k),state_dims,state_displacements&
                  &,MPI_DOUBLE_PRECISION,send_null,0,MPI_DOUBLE_PRECISION&
               &,cpl_rank,cpl_mpi_comms(k),mpi_err)
             call mpi_bcast(tag,1,MPI_INTEGER,cpl_rank,cpl_mpi_comms(k)&
                  &,mpi_err)
          end do
       case(3)
          do k = 1,cnt
             call mpi_send(x(:,k),stateDim,MPI_DOUBLE_PRECISION&
                  &,0,tag,cpl_mpi_comms(k),mpi_err)
          end do
       case(4)
          do k = 1,cnt
             particle = particles(k)
             call model_as_subroutine_start(x(:,k),particle,tag)
          end do
       case default
          print*,'EMPIRE ERROR: THIS ISNT BACK TO THE FUTURE. empire_v&
               &ersion not yet implemented'
          stop
       end select
  end subroutine send_all_models


  !> subroutine to receive all the model states from the models after
  !it has updated them one timestep
  subroutine recv_all_models(stateDim,nrhs,x)
    use timestep_data
    implicit none
    include 'mpif.h'
    integer, intent(in) :: stateDim
    integer, intent(in) :: nrhs
    real(kind=kind(1.0d0)), intent(out), dimension(stateDim,nrhs) :: x
    integer :: k
    integer, dimension(MPI_STATUS_SIZE) :: mpi_status
    integer :: mpi_err
    integer :: particle
    real(kind=kind(1.0d0)), dimension(0) :: send_null
    select case(comm_version)
    case(0)
       call user_mpi_recv(statedim,nrhs,x)
    case(1)
       DO k = 1,cnt
          particle = particles(k)
          CALL MPI_RECV(x(:,k), stateDim, MPI_DOUBLE_PRECISION, &
               particle-1, MPI_ANY_TAG, CPL_MPI_COMM,mpi_status, mpi_err)
       END DO
    case(2,5)
       do k = 1,cnt
          call mpi_gatherv(send_null,0,MPI_DOUBLE_PRECISION,x(:,k)&
               &,state_dims,state_displacements,MPI_DOUBLE_PRECISION,cpl_rank&
               &,cpl_mpi_comms(k),mpi_err)
       end do
    case(3)
       DO k = 1,cnt
          CALL MPI_RECV(x(:,k), stateDim, MPI_DOUBLE_PRECISION, &
               0, MPI_ANY_TAG, cpl_mpi_comms(k),mpi_status, mpi_err)
       END DO
    case(4)
       do k = 1,cnt
          particle = particles(k)
          call model_as_subroutine_return(x(:,k),particle)
       end do
    case default
       print*,'EMPIRE ERROR: THIS ISNT BACK TO THE FUTURE PART 2. empire_v&
            &ersion not yet implemented'
       stop
    end select
    
    
    !at this point the ensemble has been updated by the model so is
    !no longer an analysis
    call timestep_data_set_no_analysis
  end subroutine recv_all_models


  !> subroutine to receive all the model states from the models after
  !it has updated them one timestep in a non-blocking manner
  subroutine irecv_all_models(stateDim,nrhs,x,requests)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: stateDim
    integer, intent(in) :: nrhs
    real(kind=kind(1.0d0)), intent(out), dimension(stateDim,nrhs) :: x
    integer, dimension(nrhs), intent(inout) :: requests
    integer :: k
    integer :: mpi_err
    integer :: particle
    real(kind=kind(1.0d0)), dimension(0) :: send_null
    integer, dimension(MPI_STATUS_SIZE) :: mpi_status
    select case(comm_version)
       case(0)
          call user_mpi_irecv(statedim,nrhs,x,requests)
       case(1)
          DO k = 1,cnt
             particle = particles(k)
             CALL MPI_IRECV(x(:,k), stateDim, MPI_DOUBLE_PRECISION, &
                  particle-1, MPI_ANY_TAG, CPL_MPI_COMM,requests(k), mpi_err)
          end DO
       case(2,5)
          do k = 1,cnt
             !I DONT THINK MY VERSION OF MPI HAS THE MPI_IGATHERV
             !call mpi_igatherv(send_null,0,MPI_DOUBLE_PRECISION,x(:,k)&
             !     &,state_dims,state_displacements,MPI_DOUBLE_PRECISION,cpl_rank&
             !     &,cpl_mpi_comms(k),requests(k),mpi_err)

             !replace with standard blocking gather
             call mpi_gatherv(send_null,0,MPI_DOUBLE_PRECISION,x(:,k)&
                  &,state_dims,state_displacements,MPI_DOUBLE_PRECISION,cpl_rank&
                  &,cpl_mpi_comms(k),mpi_err)

             !get the requests to be null
             call mpi_isend(send_null,0,MPI_DOUBLE_PRECISION,cpl_rank&
                  &,1,cpl_mpi_comms(k),requests(k),mpi_err)
             CALL MPI_RECV (send_null,0,MPI_DOUBLE_PRECISION,cpl_rank&
                  &,1,cpl_mpi_comms(k),mpi_status, mpi_err)

          end do
       case(3)
          DO k = 1,cnt
             CALL MPI_IRECV(x(:,k), stateDim, MPI_DOUBLE_PRECISION, &
                  0, MPI_ANY_TAG, cpl_mpi_comms(k),requests(k), mpi_err)
          end DO
       case(4)
          do k = 1,cnt
             particle = particles(k)
             call model_as_subroutine_return(x(:,k),particle)
             ! the following may be totally wrong and need to change
             ! to something like the v2 case! see if it breaks if
             ! this is ever called.
             requests(k) = MPI_REQUEST_NULL
          end do
       case default
          print*,'EMPIRE ERROR: THIS ISNT BACK TO THE FUTURE PART 3. empire_v&
               &ersion not yet implemented'
          stop
       end select

  end subroutine irecv_all_models


  subroutine verify_sizes
    use sizes
    implicit none
    include 'mpif.h'
    integer :: mpi_err,i

    select case(comm_version)
    case(1,2,4,5)
       
       state_dim_g = state_dim
       obs_dim_g = obs_dim
       pf_ens_rank = pfrank

    case(3)
       if(allocated(obs_dims)) deallocate(obs_dims)
       allocate(obs_dims(pf_member_size))
       if(allocated(obs_displacements)) deallocate(obs_displacements)
       allocate(obs_displacements(pf_member_size))
       
       call mpi_allgather(obs_dim,1,MPI_INTEGER,obs_dims,1&
            &,MPI_INTEGER,pf_member_comm,mpi_err)
       
       obs_displacements(1) = 0
       if(pf_member_size .gt. 1) then
          do i = 2,pf_member_size
             obs_displacements(i) = obs_displacements(i-1) + obs_dims(i&
                  &-1)
          end do
       end if


       if(.not. allocated(state_dims)) then
          allocate(state_dims(pf_member_size))
          allocate(state_displacements(pf_member_size))

          call mpi_allgather(state_dim,1,MPI_INTEGER,state_dims,1&
               &,MPI_INTEGER,pf_member_comm,mpi_err)
          
          state_displacements(1) = 0
          if(pf_member_size .gt. 1) then
             do i = 2,pf_member_size
                state_displacements(i) = state_displacements(i-1) +&
                     & state_dims(i-1)
             end do
          end if
       end if
    case default
       print*,'EMPIRE ERROR: COMM VERSION IN VERIFY SIZES NOT IMPLEMENTED'
    end select

  end subroutine verify_sizes

end module comms
