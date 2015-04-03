!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-04-03 15:04:44 pbrowne>
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
  !! for if we are using empire v2
  integer, allocatable, dimension(:) :: state_dims !<state dimensions
  !!on each model process
  integer, allocatable, dimension(:) :: state_displacements&
       & !<displacements of the various parts of the state vector
  integer :: mdl_num_proc !< number of processes of each ensemble
  !!member
  integer, parameter :: empire_version=1

contains

  subroutine allocate_data

    implicit none

  end subroutine allocate_data

  subroutine deallocate_data
    implicit none

  end subroutine deallocate_data

  !> subroutine to make EMPIRE connections and saves details into
  !! pf_control module
  subroutine initialise_mpi

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

    if(empire_version .eq. 2) then
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


    call mpi_allgather(count,1,mpi_integer,gblcount,1,mpi_integer&
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


  end subroutine initialise_mpi


  !> subroutine to initialise new version of empire
  subroutine initialise_mpi_v2
    use pf_control
    implicit none
    include 'mpif.h'

    integer :: mpi_err
    integer :: world_size
    integer,parameter :: da=1
    integer :: state_dim
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


    call mpi_allgather(cnt,1,mpi_integer,gblcount,1,mpi_integer&
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






    pf%particles = particles
    pf%count = cnt
    pf%nens = nens
  end subroutine initialise_mpi_v2


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
    select case(empire_version)
       case(1)
          do k =1,cnt
             particle = particles(k)
             call mpi_send(x(:,k),stateDim,MPI_DOUBLE_PRECISION&
                  &,particle-1,tag,CPL_MPI_COMM,mpi_err)
          end do
       case(2)
          do k = 1,cnt
             call mpi_scatterv(x(:,k),state_dims,state_displacements&
                  &,MPI_DOUBLE_PRECISION,send_null,0,MPI_DOUBLE_PRECISION&
               &,cpl_rank,cpl_mpi_comms(k),mpi_err)
             call mpi_bcast(tag,1,mpi_integer,cpl_rank,cpl_mpi_comms(k)&
                  &,mpi_err)
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
    select case(empire_version)
       case(1)
          DO k = 1,cnt
             particle = particles(k)
             CALL MPI_RECV(x(:,k), stateDim, MPI_DOUBLE_PRECISION, &
                  particle-1, MPI_ANY_TAG, CPL_MPI_COMM,mpi_status, mpi_err)
          END DO
       case(2)
          do k = 1,cnt
             call mpi_gatherv(send_null,0,MPI_DOUBLE_PRECISION,x(:,k)&
                  &,state_dims,state_displacements,MPI_DOUBLE_PRECISION,cpl_rank&
                  &,cpl_mpi_comms(k),mpi_err)
          end do
       case default
          print*,'EMPIRE ERROR: THIS ISNT BACK TO THE FUTURE PART 2. empire_v&
               &ersion not yet implemented'
          stop
       end select
       
  end subroutine recv_all_models


  !> subroutine to receive all the model states from the models after
  !it has updated them one timestep in a non-blocking manner
  subroutine irecv_all_models(stateDim,nrhs,x,requests)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: stateDim
    integer, intent(in) :: nrhs
    real(kind=kind(1.0d0)), intent(out), dimension(stateDim,nrhs) :: x
    integer, dimension(nrhs) :: requests
    integer :: k
    integer :: mpi_err
    integer :: particle
    real(kind=kind(1.0d0)), dimension(0) :: send_null
    integer, dimension(MPI_STATUS_SIZE) :: mpi_status
    select case(empire_version)
       case(1)
          DO k = 1,cnt
             particle = particles(k)
             CALL MPI_IRECV(x(:,k), stateDim, MPI_DOUBLE_PRECISION, &
                  particle-1, MPI_ANY_TAG, CPL_MPI_COMM,requests(k), mpi_err)
          end DO
       case(2)
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
       case default
          print*,'EMPIRE ERROR: THIS ISNT BACK TO THE FUTURE PART 3. empire_v&
               &ersion not yet implemented'
          stop
       end select

  end subroutine irecv_all_models




end module comms
