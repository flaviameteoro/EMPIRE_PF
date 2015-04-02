!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    minimal_model_comms sets up communicators for empire for
!    testing purposes and sends the state back and forward
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

program minimal_model_comms
implicit none
include 'mpif.h'
integer :: mpi_err,mdl_id,cpl_root,cpl_mpi_comm
integer :: total_timesteps,i
call initialise_mpi(mdl_id,cpl_root,cpl_mpi_comm)

print*,'Reading total_timesteps from standard input: '
read*,total_timesteps


do i = 1,total_timesteps

end do



call mpi_finalize(mpi_err)
contains
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
    if(mpi_err .eq. 0) then    
       print*,'mpi_init successful'
    else
       print*,'mpi_init unsuccessful'
    end if

    da = 0
    call mpi_comm_rank (MPI_COMM_WORLD,world_id,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_rank successful'
       print*,'world_id = ',world_id
    else
       print*,'mpi_comm_rank unsuccessful'
    end if


    call mpi_comm_size (MPI_COMM_WORLD,world_size,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_size successful'
       print*,'world_size = ',world_size
    else
       print*,'mpi_comm_size unsuccessful'
    end if


    call mpi_comm_split(MPI_COMM_WORLD,da,world_id,tmp_mdls_comm,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_split successful'
    else
       print*,'mpi_comm_split unsuccessful'
    end if


    call mpi_comm_size (tmp_mdls_comm,models_size,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_size successful'
       print*,'models_size = ',models_size
    else
       print*,'mpi_comm_size unsuccessful'
    end if


    call mpi_comm_rank (tmp_mdls_comm,models_id,  mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_rank successful'
       print*,'models_id = ',models_id
    else
       print*,'mpi_comm_rank unsuccessful'
    end if

    mdlcolour = models_id/mdl_num_proc
    call mpi_comm_split(tmp_mdls_comm,mdlcolour,models_id,mdl_mpi_comm,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_split successful'
    else
       print*,'mpi_comm_split unsuccessful'
    end if


    call mpi_comm_rank (mdl_mpi_comm,mdl_id,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_rank successful'
       print*,'mdl_id = ',mdl_id
    else
       print*,'mpi_comm_rank unsuccessful'
    end if

    if(mdl_id .eq. 0) then
       cpl_colour = 9999
    else
       cpl_colour = MPI_UNDEFINED
    end if
    call mpi_comm_split(MPI_COMM_WORLD,cpl_colour,mdlcolour,cpl_mpi_comm,mpi_err)
    if(mpi_err .eq. 0) then    
       print*,'mpi_comm_split successful'
    else
       print*,'mpi_comm_split unsuccessful'
    end if

    if(mdl_id .eq. 0) then
       call mpi_comm_size(cpl_mpi_comm,nens,mpi_err)
       if(mpi_err .eq. 0) then    
          print*,'mpi_comm_size successful'
          print*,'nens = ',nens
       else
          print*,'mpi_comm_size unsuccessful'
       end if
    
       call mpi_comm_rank(cpl_mpi_comm,particle_id,mpi_err)
       if(mpi_err .eq. 0) then    
          print*,'mpi_comm_rank successful'
          print*,'particle_id = ',particle_id
       else
          print*,'mpi_comm_rank unsuccessful'
       end if
       
       nda = world_size-models_size;nens = nens - nda
       cpl_root = ((nda*particle_id)/nens)+nens
    else
       cpl_root = -1
    end if
    print*,'cpl_root = ',cpl_root
  end subroutine initialise_mpi
end program minimal_model_comms
