!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    minimal_model_comms just sets up communicators for empire for
!    testing purposes
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

program minimal_model_comms_v2
  implicit none

  include 'mpif.h'
  integer :: mpi_err
  integer :: world_rank
  integer :: world_size
  integer, parameter :: mdl_num_proc=5
  integer :: mdl_mpi_comm
  integer :: temp_mdls_comm
  integer :: temp_mdls_rank
  integer :: mdl_rank
  integer :: da
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), allocatable, dimension(:) :: state_vector
  integer :: state_dim
  integer :: i
  integer :: cpl_root
  integer :: particle_rank
  integer :: nda
  integer :: nens
  integer :: temp_mdls_size
  integer :: temp_cpl_comm
  integer :: first_ptcl
  integer :: final_ptcl
  integer :: cpl_mpi_comm
  integer :: null_mpi_comm
  real(kind=rk), dimension(0) :: send_null
  integer :: total_timesteps
  integer :: tag

  call mpi_init(mpi_err)
      if(mpi_err .eq. 0) then    
       print*,'mpi_init successful'
    else
       print*,'mpi_init unsuccessful'
    end if


  call mpi_comm_rank(mpi_comm_world,world_rank,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_rank successful'
       print*,'world_rank = ',world_rank
    else
       print*,'mpi_comm_rank unsuccessful'
    end if


  call mpi_comm_size(mpi_comm_world,world_size,mpi_err)
      if(mpi_err .eq. 0) then    
       print*,'mpi_comm_size successful'
       print*,'world_size = ',world_size
    else
       print*,'mpi_comm_size unsuccessful'
    end if


  cpl_root = world_size-1
  print*,'rank = ',world_rank,' on mpi_comm_world which has size ',world_size
  da = 0


  call mpi_allreduce(mdl_num_proc,i,1,MPI_INTEGER,MPI_MAX&
       &,MPI_COMM_WORLD,mpi_err)
  if(mpi_err .eq. 0) then    
       print*,'mpi_allreduce successful'
       print*,'i = ',i
    else
       print*,'mpi_allreduce unsuccessful'
    end if



  call mpi_comm_split(MPI_COMM_WORLD,da,world_rank,temp_mdls_comm&
       &,mpi_err)
      if(mpi_err .eq. 0) then    
       print*,'mpi_comm_split successful: temp_mdls_comm created'
    else
       print*,'mpi_comm_split unsuccessful: temp_mdls_comm not created'
    end if
  
    call mpi_comm_size(temp_mdls_comm,temp_mdls_size,mpi_err)
    if(mpi_err .eq. 0) then    
     print*,'mpi_comm_size successful'
     print*,'temp_mdls_size = ',temp_mdls_size
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
       print*,'mpi_comm_rank successful'
       print*,'temp_mdls_rank = ',temp_mdls_rank
    else
       print*,'mpi_comm_rank unsuccessful'
    end if


  particle_rank = temp_mdls_rank/mdl_num_proc

  call mpi_comm_split(temp_mdls_comm,particle_rank,temp_mdls_rank&
       &,mdl_mpi_comm,mpi_err)
  if(mpi_err .eq. 0) then    
       print*,'mpi_comm_split successful: mdl_mpi_comm created'
    else
       print*,'mpi_comm_split unsuccessful: mdl_mpi_comm not created'
    end if
    


  call mpi_comm_rank(mdl_mpi_comm,mdl_rank,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_rank successful'
       print*,'mdl_rank = ',mdl_rank
    else
       print*,'mpi_comm_rank unsuccessful'
    end if


  cpl_root = nda*particle_rank/nens
  print*,'cpl_root = ',cpl_root

  if(cpl_root .lt. 0) then
     print*,'MINIMAL MODEL COMMS v2 ERROR: cpl_root is less than 0.'
     print*,'Make sure you launch with a DA CODE'
     stop
  end if

  call mpi_comm_split(mpi_comm_world,cpl_root,temp_mdls_rank,temp_cpl_comm,mpi_err)
  if(mpi_err .eq. 0) then    
     print*,'mpi_comm_split successful: temp_cpl_comm created'
  else
     print*,'mpi_comm_split unsuccessful: temp_cpl_comm not created'
  end if






  first_ptcl = -((-cpl_root)*nens/nda)
  final_ptcl = -((-cpl_root-1)*nens/nda)-1

  first_ptcl = ceiling(real(cpl_root)*real(nens)/real(nda))
  final_ptcl = ceiling(real(cpl_root+1)*real(nens)/real(nda))-1


  print*,'range of particles = ',first_ptcl,final_ptcl

 

  do i = first_ptcl,final_ptcl
     print*,'i = ',i,' particle_rank = ',particle_rank
     if(i .eq. particle_rank) then
        call mpi_comm_split(temp_cpl_comm,1,temp_mdls_rank&
             &,cpl_mpi_comm,mpi_err)
        print*,'created cpl_mpi_comm'
     else
        print*,'doing null splitting'
        call mpi_comm_split(temp_cpl_comm,0,temp_mdls_rank&
             &,null_mpi_comm,mpi_err)
        print*,'created mpi_comm_null'
        call mpi_comm_free(null_mpi_comm,mpi_err)
        print*,'freed up null_mpi_comm'
     end if


  end do

  cpl_root = mdl_num_proc








  select case(mdl_rank)
  case(0)
     state_dim = 1
  case(1)
     state_dim = 3
  case(2)
     state_dim = 2
  case(3)
     state_dim = 5
  case(4)
     state_dim = 1
  case default
     print*,'it was at this point, model realised, he fucked up'
     stop
  end select
  allocate(state_vector(state_dim))




  state_vector = 10*mdl_rank + (/ (real(i,rk), i = 1,state_dim) /)

  print*,'state_vector = '
  print*,state_vector

  print*,'doing a gather on cpl_mpi_comm'
  call mpi_gather(state_dim,1,MPI_INTEGER,state_dim&
       &,1,MPI_INTEGER,cpl_root,cpl_mpi_comm,mpi_err)
  print*,'finished the gather on cpl_mpi_comm'





















  print*,'Reading total_timesteps from file timesteps: '
  open(11,file='timesteps',action='read',status='old')
  read(11,*) total_timesteps
  close(11)

  !send the state to da code with mpi_gatherv
  call mpi_gatherv(state_vector,state_dim,MPI_DOUBLE_PRECISION,state_vector&
       &,state_dim,state_dim,MPI_DOUBLE_PRECISION,cpl_root&
       &,cpl_mpi_comm,mpi_err)


  !get the state back from da code with mpi_gatherv
  call mpi_scatterv(send_null,0,0,MPI_DOUBLE_PRECISION,state_vector&
       &,state_dim,MPI_DOUBLE_PRECISION,cpl_root,cpl_mpi_comm,mpi_err)
  !get the tag from the da code
  call mpi_bcast(tag,1,mpi_integer,cpl_root,cpl_mpi_comm,mpi_err)
  print*,'Received tag = ',tag


  do i = 1,total_timesteps
     print*,'Timestep = ',i

     !send the state to da code with mpi_gatherv
     call mpi_gatherv(state_vector,state_dim,MPI_DOUBLE_PRECISION,state_vector&
          &,state_dim,state_dim,MPI_DOUBLE_PRECISION,cpl_root&
          &,cpl_mpi_comm,mpi_err)
     
     !get the state back from da code with mpi_gatherv
     call mpi_scatterv(send_null,0,0,MPI_DOUBLE_PRECISION,state_vector&
          &,state_dim,MPI_DOUBLE_PRECISION,cpl_root,cpl_mpi_comm,mpi_err)
     !get the tag from the da code
     call mpi_bcast(tag,1,mpi_integer,cpl_root,cpl_mpi_comm,mpi_err)
     print*,'Received tag = ',tag
  end do
  


















  call mpi_finalize(mpi_err)
  print*,'MINIMAL_MODEL_COMMS_v2 finished nicely'
end program minimal_model_comms_v2
