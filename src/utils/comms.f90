!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-04-01 22:46:41 pbrowne>
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
  integer :: mype_id !< the rank of this process on MPI_COMM_WORLD
  integer :: myRank !< the rank of this process on CPL_MPI_COMM
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
  integer, allocatable, dimension(:) :: particles !< the ensemble members associated with this process
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
    integer :: world_id
    integer :: myrank
    integer :: i
    integer :: da,count
    integer :: pf_colour
    integer :: world_size

    pf_colour = 10000
    couple_colour=9999
    call MPI_INIT (mpi_err)

    da = 1
    call MPI_COMM_RANK (MPI_COMM_WORLD,world_id,     mpi_err)
    call mpi_comm_size (mpi_comm_world,world_size,   mpi_err)
    call mpi_comm_split(mpi_comm_world,da,           world_id,  pf_mpi_comm, mpi_err)
    call mpi_comm_rank (pf_mpi_comm,   pfrank,       mpi_err)
    call mpi_comm_size (pf_mpi_comm,   npfs,          mpi_err)
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

    gblcount=count

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
    do k =1,cnt
       particle = particles(k)
       call mpi_send(x(:,k),stateDim,MPI_DOUBLE_PRECISION&
            &,particle-1,tag,CPL_MPI_COMM,mpi_err)
    end do    
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
    DO k = 1,cnt
       particle = particles(k)
       CALL MPI_RECV(x(:,k), stateDim, MPI_DOUBLE_PRECISION, &
            particle-1, MPI_ANY_TAG, CPL_MPI_COMM,mpi_status, mpi_err)
    END DO
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
    DO k = 1,cnt
       particle = particles(k)
       CALL MPI_IRECV(x(:,k), stateDim, MPI_DOUBLE_PRECISION, &
            particle-1, MPI_ANY_TAG, CPL_MPI_COMM,requests(k), mpi_err)
    end DO
  end subroutine irecv_all_models


  
  
end module comms
