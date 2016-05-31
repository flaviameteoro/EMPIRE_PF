!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    linear_empire_vader.f90 Implements a linear model with EMPIRE
!    coupling extended to include reverse communication (VADER)
!
!
!The MIT License (MIT)
!
!Copyright (c) 2015 Philip A. Browne
!
!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:
!
!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.
!
!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.
!
!Email: p.browne@reading.ac.uk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> program to implement a simple linear model of no use to anyone
!! but for testing and debugging purposes :)
!!
!! NOTE: THIS PROGRAM ***MUST*** RECIEVE A COUPLET OF INTEGERS
!!       FROM THE DATA ASSIMILATION CODE CONTAINING 
!!
!!       FIRST : THE SIZE OF THE DIMENSION OF THE MODEL
!!
!!       SECOND: THE NUMBER OF TIMESTEPS THE MODEL SHOULD DO
!!
!!       THIS IS A BIT WEIRD, AS NORMALLY THE MODEL DICTATES
!!       SUCH THINGS. BUT THIS IS A USELESS TOY MODEL.
!!       SO WE MIGHT AS WELL MAKE IT EASY TO USE TO TEST DA.
program linear
  implicit none
  include 'mpif.h'
  real(kind=kind(1.0D0)), dimension(:), allocatable :: x
  integer :: i,n,maxt
  integer, dimension(2) :: data
  integer :: mpi_err,mdl_id,cpl_root,cpl_mpi_comm
  integer,  dimension(MPI_STATUS_SIZE) :: mpi_status
  
  ! SET UP EMPIRE COMMUNICATORS
  call initialise_mpi(mdl_id,cpl_root,cpl_mpi_comm)

  
  ! GET THE COUPLET FROM THE DATA ASSIMILATION CODE
  call mpi_recv(data,2,MPI_INTEGER,cpl_root,MPI_ANY_TAG,cpl_mpi_comm&
       &,mpi_status,mpi_err)

  ! PUT THE DATA INTO THE RIGHT VARIABLES AND ALLOCATE
  n = data(1)
  maxt = data(2)
  allocate(x(n))
  x = 1.0d0

  
  ! DO THE INITIAL SEND AND RECIEVES FROM EMPIRE
  call mpi_send(x,n,MPI_DOUBLE_PRECISION,cpl_root&
       &,1,cpl_mpi_comm,mpi_err)
  call mpi_recv(x,n,MPI_DOUBLE_PRECISION,cpl_root&
       &,MPI_ANY_TAG,cpl_mpi_comm,mpi_status,mpi_err)
2 continue

  write(6,*) 0,x
  ! START THE TIMESTEP LOOP
  do i = 1,maxt

     ! CALL THE SIMPLE LINEAR MODEL AND UPDATE MODEL STATE
     x = f(n,x)

     ! DO THE TIMESTEP LOOP SEND AND RECIEVES WITH EMPIRE
     call mpi_send(x,n,MPI_DOUBLE_PRECISION,cpl_root&
          &,1,cpl_mpi_comm,mpi_err)
     call mpi_recv(x,n,MPI_DOUBLE_PRECISION,cpl_root&
          &,MPI_ANY_TAG,cpl_mpi_comm,mpi_status,mpi_err)

     write(6,*) i,x
     if(mpi_status(MPI_TAG) .eq. 1) then
        go to 1 !simply continue code
     elseif(mpi_status(MPI_TAG) .eq. 2) then
        go to 2 !restart code
     elseif(mpi_status(MPI_TAG) .eq. 3) then
        go to 3 !end code
     else
        print*,'Linear model error: unknown MPI_TAG: ',mpi_status(MPI_TAG)
        stop '-1'
     end if
1    continue


     !END TIMESTEP LOOP
  end do
3 continue
  
  ! END CODE
  deallocate(x)
  call mpi_finalize(mpi_err)
contains

  ! DEFINE A SIMPLE LINEAR MODEL
  function f (n,x)
    implicit none
    integer, intent(in) :: n
    real(kind=kind(1.0D0)),intent(in),dimension (n) :: x
    real(kind=kind(1.0D0)),dimension(n) :: f
    f = x
  end function f

  ! STANDARD EMPIRE INITIALISATION SUBROUTINE
  subroutine initialise_mpi(mdl_id,cpl_root,cpl_mpi_comm)
    implicit none
    include 'mpif.h'
    integer, intent(out) :: mdl_id,cpl_root,cpl_mpi_comm
    integer :: mdl_num_proc=1
    integer :: mpi_err,world_size,world_id
    integer :: cpl_colour
    integer :: particle_id,nens, da, nda
    integer :: mdl_mpi_comm,mdlcolour
    integer :: tmp_mdls_comm,models_id,models_size
    call mpi_init(mpi_err)
    da = 0
    call mpi_comm_rank (MPI_COMM_WORLD,world_id,mpi_err)
    call mpi_comm_size (MPI_COMM_WORLD,world_size,mpi_err)
    call mpi_comm_split(MPI_COMM_WORLD,da,world_id,tmp_mdls_comm,mpi_err)
    call mpi_comm_size (tmp_mdls_comm,models_size,mpi_err)
    call mpi_comm_rank (tmp_mdls_comm,models_id,  mpi_err)
    mdlcolour = models_id/mdl_num_proc
    call mpi_comm_split(tmp_mdls_comm,mdlcolour,models_id,mdl_mpi_comm,mpi_err)
    call mpi_comm_rank (mdl_mpi_comm,mdl_id,mpi_err)
    if(mdl_id .eq. 0) then
       cpl_colour = 9999
    else
       cpl_colour = MPI_UNDEFINED
    end if
    call mpi_comm_split(MPI_COMM_WORLD,cpl_colour,mdlcolour,cpl_mpi_comm,mpi_err)
    if(mdl_id .eq. 0) then
       call mpi_comm_size(cpl_mpi_comm,nens,mpi_err)
       call mpi_comm_rank(cpl_mpi_comm,particle_id,mpi_err)
       nda = world_size-models_size;nens = nens - nda
       cpl_root = ((nda*particle_id)/nens)+nens
       if(nda ==0) cpl_root = particle_id
    else
       cpl_root = -1
    end if
  end subroutine initialise_mpi

end program linear
