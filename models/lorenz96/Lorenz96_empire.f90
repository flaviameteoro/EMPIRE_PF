!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Lorenz96_empire.f90 Implements Lorenz 1996 model with EMPIRE coupling
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

program lorenz96
  implicit none
  include 'mpif.h'
  real(kind=kind(1.0D0)) :: dt=1.0d-2
  real(kind=kind(1.0D0)), allocatable, dimension(:) :: x,k1,k2,k3,k4
  real(kind=kind(1.0d0)) :: F=8.0d0
  integer :: N=40
  integer :: total_timesteps=100
  integer :: t
  integer :: mpi_err,mdl_id,cpl_root,cpl_mpi_comm
  integer, dimension(MPI_STATUS_SIZE) :: mpi_status
  integer :: ios
  logical :: l96_exists
  
  namelist/l96/N,&
       &total_timesteps,&
       &F,&
       &dt

  inquire(file='l96.nml', exist=l96_exists)
  if(l96_exists) then
     open(32,file='l96.nml',iostat=ios,action='read'&
          &,status='old')
     if(ios .ne. 0) stop 'Cannot open l96.nml'
     read(32,nml=l96) 
     close(32)
  end if








  call initialise_mpi(mdl_id,cpl_root,cpl_mpi_comm)


  allocate(x(N),k1(N),k2(N),k3(N),k4(N))
  x = F
  x(N/2) = F+0.05


  if(mdl_id .eq. 0) then
     call mpi_send(x,N,MPI_DOUBLE_PRECISION,cpl_root&
          &,1,cpl_mpi_comm,mpi_err)
     call mpi_recv(x,N,MPI_DOUBLE_PRECISION,cpl_root&
          &,MPI_ANY_TAG,cpl_mpi_comm,mpi_status,mpi_err)
  end if
2 continue
  print*,x


  do t = 1,total_timesteps
     k1 = g (x                  , N , F )
     k2 = g (x +0.5D0 * dt * k1 , N , F )
     k3 = g (x +0.5D0 * dt * k2 , N , F )
     k4 = g (x +        dt * k3 , N , F )
     x = x + dt *( k1 + 2.0D0 *( k2 + k3 ) + k4 )/6.0D0

     if(mdl_id .eq. 0) then
        call mpi_send(x,N,MPI_DOUBLE_PRECISION,cpl_root&
             &,1,cpl_mpi_comm,mpi_err)
        call mpi_recv(x,N,MPI_DOUBLE_PRECISION,cpl_root&
             &,MPI_ANY_TAG,cpl_mpi_comm,mpi_status,mpi_err)
     end if
     print*,x
     
     select case(mpi_status(MPI_TAG))
     case(2)
        go to 2
     case(3)
        go to 3
     case default
     end select

  end do
3 continue
  call mpi_finalize(mpi_err)

contains
  function g (x , N, F )
    implicit none
    real(kind=kind(1.0D0)),intent(in),dimension(0:N-1) :: x
    integer, intent(in) :: N
    real(kind=kind(1.0D0)),dimension(0:N-1) :: g
    real(kind=kind(1.0D0)),intent(in) :: F
    integer :: j
    do j = 0,N-1
       g(j) = ( x(modulo(j+1,N)) -x( modulo(j-2,N)) )*&
            &x(modulo(j-1,N)) - x(j) + F
    end do
  end function g
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
    else
       cpl_root = -1
    end if
  end subroutine initialise_mpi
end program lorenz96
