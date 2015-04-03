!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Lorenz63_empire.f90 Implements Lorenz 1963 model with EMPIRE coupling
!
!The MIT License (MIT)
!
!Copyright (c) 2014 Philip A. Browne
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

program lorenz63
  implicit none
  include 'mpif.h'
  real(kind=kind(1.0D0)) :: t,sigma,rho,beta,dt,tstart,tstop
  real(kind=kind(1.0D0)), dimension(3) :: x,k1,k2,k3,k4
  integer :: mpi_err,mdl_id,cpl_root,cpl_mpi_comm
  integer, dimension(MPI_STATUS_SIZE) :: mpi_status
  call initialise_mpi(mdl_id,cpl_root,cpl_mpi_comm)
  tstart =0.0D0 ; dt = 0.01D0 ; tstop = real(40*100)*dt
  sigma = 10.0D0 ; rho = 28.0D0 ; beta = 8.0D0 /3.0D0
  x = (/ 1.508870D0, -1.531271D0 , 25.46091D0 /)
  call mpi_send(x,3,MPI_DOUBLE_PRECISION,cpl_root&
       &,1,cpl_mpi_comm,mpi_err)
  call mpi_recv(x,3,MPI_DOUBLE_PRECISION,cpl_root&
       &,MPI_ANY_TAG,cpl_mpi_comm,mpi_status,mpi_err)
2 continue
  t = tstart
  do; if ( t .ge. tstop -1.0D-10) exit
     k1 = f (x                  , sigma , rho , beta )
     k2 = f (x +0.5D0 * dt * k1 , sigma , rho , beta )
     k3 = f (x +0.5D0 * dt * k2 , sigma , rho , beta )
     k4 = f (x +        dt * k3 , sigma , rho , beta )
     x = x + dt *( k1 + 2.0D0 *( k2 + k3 ) + k4 )/6.0D0
     call mpi_send(x,3,MPI_DOUBLE_PRECISION,cpl_root&
          &,1,cpl_mpi_comm,mpi_err)
     call mpi_recv(x,3,MPI_DOUBLE_PRECISION,cpl_root&
          &,MPI_ANY_TAG,cpl_mpi_comm,mpi_status,mpi_err)
     t = t + dt
     print*,x(1),x(2),x(3)
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
  function f (x , sigma , rho , beta )
    implicit none
    real(kind=kind(1.0D0)),intent(in),dimension (3) :: x
    real(kind=kind(1.0D0)),dimension(3) :: f
    real(kind=kind(1.0D0)),intent(in) :: sigma , rho , beta
    f = (/sigma *(x(2)-x(1)),x(1)*(rho-x(3)) -x(2),x(1)*x(2)-beta*x(3)/)
  end function f
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
end program lorenz63
