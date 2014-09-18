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
module comms
!  use hadcm3_config
!  use hadcm3_data

  integer :: CPL_MPI_COMM,mype_id,myRank,nProc
  integer :: pf_mpi_comm,pfrank
  integer*8 :: npfs
  integer, allocatable, dimension(:) :: gblcount,gbldisp
contains
  
  subroutine allocate_data

    implicit none
    
  end subroutine allocate_data
  
  subroutine deallocate_data
    implicit none

  end subroutine deallocate_data
  
  subroutine initialise_mpi
    use pf_control
    implicit none
    include 'mpif.h'
    
    integer :: mpi_err!dummy_colour,mpi_err
    integer :: couple_colour !DUMMY_MPI_COMMUNICATOR,couple_colour
    !integer :: couple_mype_id,couple_root
    integer :: rtmp!,ctmp

 !   integer :: tag!,state_dim!,iter
 !   integer :: num_iters
    integer :: particle,mype_id
    integer :: myrank !nproc,myrank
!    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: nens,i
    integer*8 :: da
    integer :: count,pf_colour,pf_id!,pf_mpi_comm

    
    pf_colour = 10000
    couple_colour=9999
    CALL MPI_INIT (mpi_err)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,mype_id,mpi_err)
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,pf_colour,pf_id,pf_mpi_comm&
         &,mpi_err)
    
    CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,couple_colour,mype_id&
         &,CPL_MPI_COMM,mpi_err)
    CALL MPI_COMM_RANK (CPL_MPI_COMM, myRank, mpi_err)
    CALL MPI_COMM_SIZE (CPL_MPI_COMM, nens, mpi_err)
    
    da = 1
    CALL MPI_ALLREDUCE(da,npfs,1,mpi_integer8,mpi_sum,cpl_mpi_comm&
         &,mpi_err)
    nens = nens-npfs
    
    pfrank = myrank-nens
    
    !lets find the particles:
    count = 0
    do particle = 1,nens
       if( real(particle-1) .ge. real(nens*(pfrank))/real(npfs) .and.&
            & real(particle-1) .lt. real(nens*(pfrank+1))/real(npfs)) then
          count = count + 1
       end if
    end do
    
    allocate(pf%particles(count))
    rtmp = 0
    do particle = 1,nens
       if(real(particle-1) .ge. real(nens*(pfrank))/real(npfs) .and.&
            & real(particle-1) .lt. real(nens*(pfrank+1))/real(npfs))&
            & then
          rtmp = rtmp + 1
          pf%particles(rtmp) = particle
       end if
    end do
    
    allocate(gblcount(npfs))
    allocate(gbldisp(npfs))
!    print*,'woohoo allgather'
!    print*,count
!    print*,gblcount
!    print*,pf_mpi_comm
    call mpi_allgather(count,1,mpi_integer,gblcount,1,mpi_integer&
         &,pf_mpi_comm,mpi_err)
!    print*,'allgather did not break'
    gbldisp = 0
    if(npfs .gt. 1) then
       do i = 2,npfs
          gbldisp(i) = gbldisp(i-1) + gblcount(i-1)
       end do
    end if
    pf%count = count

    pf%nens = nens
    PRINT*,'PF_rank = ',pfrank,' and I own particles ',pf%particles

    
  end subroutine initialise_mpi

end module comms