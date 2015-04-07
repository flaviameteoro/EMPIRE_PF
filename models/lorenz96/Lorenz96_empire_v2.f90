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

program lorenz96_v2
  implicit none
  include 'mpif.h'
  real(kind=kind(1.0D0)) :: dt=1.0d-2
  real(kind=kind(1.0D0)), allocatable, dimension(:) :: x,k1,k2,k3,k4
  real(kind=kind(1.0d0)) :: F=8.0d0
  integer :: N=40
  integer :: total_timesteps=100
  integer :: t
  integer :: mpi_err,mdl_rank,cpl_root,cpl_mpi_comm
  real(kind=kind(1.0d0)), dimension(0) :: send_null
  integer :: tag
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


  call initialise_mpi_v2(mdl_rank,cpl_root,cpl_mpi_comm)

  call empire_process_dimensions(N,cpl_root,cpl_mpi_comm)
  allocate(x(N),k1(N),k2(N),k3(N),k4(N))
  x = F
  x(N/2) = F+0.05
  print*,x


  !send the state to da code with mpi_gatherv
  call mpi_gatherv(x,N,MPI_DOUBLE_PRECISION,x&
       &,N,N,MPI_DOUBLE_PRECISION,cpl_root&
       &,cpl_mpi_comm,mpi_err)


  !get the state back from da code with mpi_gatherv
  call mpi_scatterv(send_null,0,0,MPI_DOUBLE_PRECISION,x&
       &,N,MPI_DOUBLE_PRECISION,cpl_root,cpl_mpi_comm,mpi_err)
  !get the tag from the da code
  call mpi_bcast(tag,1,mpi_integer,cpl_root,cpl_mpi_comm,mpi_err)
2 continue


  do t = 1,total_timesteps
     k1 = g (x                  , N , F )
     k2 = g (x +0.5D0 * dt * k1 , N , F )
     k3 = g (x +0.5D0 * dt * k2 , N , F )
     k4 = g (x +        dt * k3 , N , F )
     x = x + dt *( k1 + 2.0D0 *( k2 + k3 ) + k4 )/6.0D0


     !send the state to da code with mpi_gatherv
     call mpi_gatherv(x,N,MPI_DOUBLE_PRECISION,x&
          &,N,N,MPI_DOUBLE_PRECISION,cpl_root&
          &,cpl_mpi_comm,mpi_err)

     !get the state back from da code with mpi_gatherv
     call mpi_scatterv(send_null,0,0,MPI_DOUBLE_PRECISION,x&
          &,N,MPI_DOUBLE_PRECISION,cpl_root,cpl_mpi_comm,mpi_err)
     !get the tag from the da code
     call mpi_bcast(tag,1,mpi_integer,cpl_root,cpl_mpi_comm,mpi_err)

     print*,x

     select case(tag)
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
  subroutine initialise_mpi_v2(mdl_rank,cpl_root,cpl_mpi_comm)


    implicit none

    include 'mpif.h'
    integer, intent(out) :: mdl_rank
    integer, intent(out) :: cpl_root    
    integer, intent(out) :: cpl_mpi_comm

    integer, parameter :: mdl_num_proc=1
    integer :: mdl_mpi_comm

    integer :: mpi_err
    integer :: world_rank
    integer :: world_size

    integer :: temp_mdls_size
    integer :: temp_cpl_comm
    integer :: temp_mdls_comm
    integer :: temp_mdls_rank
    integer :: da
    integer :: i
    integer :: particle_rank
    integer :: nda
    integer :: nens
    integer :: first_ptcl
    integer :: final_ptcl
    integer :: null_mpi_comm

    logical :: msg=.true.


    call mpi_init(mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_init successful'
    else
       print*,'mpi_init unsuccessful'
    end if


    call mpi_comm_rank(mpi_comm_world,world_rank,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_comm_rank successful'
       if(msg) print*,'world_rank = ',world_rank
    else
       print*,'mpi_comm_rank unsuccessful'
    end if


    call mpi_comm_size(mpi_comm_world,world_size,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_comm_size successful'
       if(msg) print*,'world_size = ',world_size
    else
       print*,'mpi_comm_size unsuccessful'
    end if


    cpl_root = world_size-1
    if(msg) then
       print*,'rank = ',world_rank,' on mpi_comm_world which has &
            &size ',world_size
    end if

    da = 0


    call mpi_allreduce(mdl_num_proc,i,1,MPI_INTEGER,MPI_MAX&
         &,MPI_COMM_WORLD,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_allreduce successful'
       if(msg) print*,'i = ',i
    else
       print*,'mpi_allreduce unsuccessful'
    end if



    call mpi_comm_split(MPI_COMM_WORLD,da,world_rank,temp_mdls_comm&
         &,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_comm_split successful: temp_mdls_comm created'
    else
       print*,'mpi_comm_split unsuccessful: temp_mdls_comm not created'
    end if

    call mpi_comm_size(temp_mdls_comm,temp_mdls_size,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_comm_size successful'
       if(msg) print*,'temp_mdls_size = ',temp_mdls_size
    else
       print*,'mpi_comm_size unsuccessful'
    end if


    if(mod(temp_mdls_size,mdl_num_proc) .ne. 0) then
       print*,'MINIMAL MODEL LAUNCH ERROR.'
       print*,'MUST LAUNCH A MULTIPLE OF ',mdl_num_proc,' copies of the &
            &model'
       stop
    end if


    nda = world_size-temp_mdls_size
    if(nda .lt. 1) then
       print*,'MINIMAL MODEL COMMS v2 ERROR: nda is less than 1.'
       print*,'Make sure you launch with a DA CODE'
       stop
    end if



    nens = temp_mdls_size/mdl_num_proc
    call mpi_comm_rank(temp_mdls_comm,temp_mdls_rank,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_comm_rank successful'
       if(msg) print*,'temp_mdls_rank = ',temp_mdls_rank
    else
       print*,'mpi_comm_rank unsuccessful'
    end if


    particle_rank = temp_mdls_rank/mdl_num_proc

    call mpi_comm_split(temp_mdls_comm,particle_rank,temp_mdls_rank&
         &,mdl_mpi_comm,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_comm_split successful: mdl_mpi_comm created'
    else
       print*,'mpi_comm_split unsuccessful: mdl_mpi_comm not created'
    end if



    call mpi_comm_rank(mdl_mpi_comm,mdl_rank,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_comm_rank successful'
       if(msg) print*,'mdl_rank = ',mdl_rank
    else
       print*,'mpi_comm_rank unsuccessful'
    end if


    cpl_root = nda*particle_rank/nens
    if(msg) print*,'cpl_root = ',cpl_root

    if(cpl_root .lt. 0) then
       print*,'MINIMAL MODEL COMMS v2 ERROR: cpl_root is less than 0.'
       print*,'Make sure you launch with a DA CODE'
       stop
    end if

    call mpi_comm_split(mpi_comm_world,cpl_root,temp_mdls_rank,temp_cpl_comm,mpi_err)
    if(mpi_err .eq. 0) then    
       if(msg) print*,'mpi_comm_split successful: temp_cpl_comm created'
    else
       print*,'mpi_comm_split unsuccessful: temp_cpl_comm not created'
    end if






    first_ptcl = ceiling(real(cpl_root)*real(nens)/real(nda))
    final_ptcl = ceiling(real(cpl_root+1)*real(nens)/real(nda))-1


    if(msg) print*,'range of particles = ',first_ptcl,final_ptcl



    do i = first_ptcl,final_ptcl
       if(msg) print*,'i = ',i,' particle_rank = ',particle_rank
       if(i .eq. particle_rank) then
          call mpi_comm_split(temp_cpl_comm,1,temp_mdls_rank&
               &,cpl_mpi_comm,mpi_err)
          if(msg) print*,'created cpl_mpi_comm'
       else
          if(msg) print*,'doing null splitting'
          call mpi_comm_split(temp_cpl_comm,0,temp_mdls_rank&
               &,null_mpi_comm,mpi_err)
          if(msg) print*,'created mpi_comm_null'
          call mpi_comm_free(null_mpi_comm,mpi_err)
          if(msg) print*,'freed up null_mpi_comm'
       end if


    end do

    cpl_root = mdl_num_proc




  end subroutine initialise_mpi_v2

  subroutine empire_process_dimensions(N,cpl_root,cpl_mpi_comm)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: N
    integer, intent(in) :: cpl_root
    integer, intent(in) :: cpl_mpi_comm
    integer :: mpi_err
    logical :: msg = .true.
    
    if(msg) print*,'called empire_process_dimensions'
    call mpi_gather(N,1,MPI_INTEGER,N&
         &,1,MPI_INTEGER,cpl_root,cpl_mpi_comm,mpi_err)
    if(msg) print*,'finished the gather on cpl_mpi_comm for empire_process_dimensions'  
  end subroutine empire_process_dimensions



end program lorenz96_v2
