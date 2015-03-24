!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-03-24 11:16:35 pbrowne>
!!!
!!!    subroutine to provide objective function and gradient for var
!!!    Copyright (C) 2015  Philip A. Browne
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


!> This is the subroutine which the optimization routines call
!! to get the objective function value and its gradient
subroutine fcn( n, x, f, g )
  implicit none
  integer,intent(in) :: n !< the dimension of the optimzation problem
  real(kind=kind(1.0d0)), dimension(n), intent(in) :: x !< the
  !!current optimization state
  real(kind=kind(1.0d0)), intent(out) :: f !< the objective function value
  real(kind=kind(1.0d0)), dimension(n), intent(out) :: g !< the
  !!gradient of the objective function
  
  call fourdenvar_fcn(n,x,f,g)
  print*,'function = ',f
  print*,'gradient = ',g
end subroutine fcn


!> subroutine to provide the objective function and gradient for
!! 4dEnVar.
!!
!!
!! Let \f$x\f$ be the state we wish to find using Var.
!!
!! Then we define \f$X_k := \{ x_1(k)-x(k);\ldots;x_m(k)-x(k)\}\f$
!!
!! to be the ensemble perturbation matrix, 
!!
!! where \f$x_j(k)\f$ is the jth
!! ensemble member at time \f$k\f$
!!
!! and \f$x(k)\f$ is the 
!! optimization solution integrated forward in time to timestep
!! \f$k\f$.
!!
!! The objective function considered is
!!
!! \f$J(x) = \frac{1}{2}(x-x_b)^TB^{-1}(x-x_b) +
!! \frac{1}{2}\sum_i(y_i-H_i(M_i(x)))^T R_i^{-1} (y_i - H_i(M_i(x)) )\f$
!!
!! where \f$x_b\f$ is a background guess, \f$B\f$ the background
!! error covariance matrix,
!!
!! \f$y_i\f$ observations at a timestep \f$i\f$, \f$H_i\f$ the
!! corresponding observation operator with associated observation
!! error covariance matrix \f$R_i\f$ and \f$M_i\f$ the model which
!! propogates a state from time \f$0\f$ to the observation timestep
!! \f$i\f$.
!!
!! In this code, \f$B:= \frac{1}{m-1}X_0X_0^T\f$.
!!
!! We make the following control variable transform:
!!
!! \f$x = X_0 v + x_b\f$ where \f$v \in \mathbb{R}^m\f$.
!!
!! Then the objective function can be re-written as a function of v,
!!
!! \f$f=J(x)=J(v) = \frac{1}{2}(m-1)v^Tv +
!! \frac{1}{2}\sum_i(y_i-H_i(M_i(X_0v+x_b)))^T R_i^{-1} (y_i -
!! H_i(M_i(X_0v+x_b)) )\f$.
!!
!! The gradient of the objective function can then be written
!!
!! \f$g = \nabla_vJ(v) \approx (m-1)v -
!! \sum_i (H_i(X_i))^T R_i^{-1} (y_i -    
!! H_i(M_i(X_0v+x_b))
!! \f$
!!
!! which is exact if \f$H_i\f$ and \f$M_i\f$ are linear and \f$X_0\f$
!! is invertible (at least I assume there has to be this condition on
!! \f$X_0\f$...).
!!
!! Note this is not exactly the 4dEnVar algorithm as given by Liu,
!! Xian and Wang (2008) as the ensemble perturbation matrix here are
!! perturbations around the current var solution, not the
!! perturbations around \f$M_i(x_b)\f$. 
!! i.e. \f$X_i = X_i(x)=X_i(v)\f$.
subroutine fourdenvar_fcn(n, v, f, g )
  use comms
  use var_data
  use fourdenvardata
  use sizes, only : state_dim
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: n !< this is the dimension of the
  !!optimization control variable, so the number of ensemble members-1
  real(kind=rk), dimension(n), intent(in) :: v!< this is the
  !!optimization control variable
  real(kind=rk), intent(out) :: f !< the 4dvar objective function
  real(kind=kind(1.0d0)), dimension(n), intent(out) :: g !< the
  !<gradient of the 4dvar objective function
  real(kind=kind(1.0d0)), dimension(n) :: tempg 
  integer :: t
  integer :: tag
  real(kind=rk), dimension(:), allocatable :: y
  real(kind=rk), dimension(:,:), allocatable :: hxt
  real(kind=rk), dimension(state_dim) :: xdet
  real(kind=rk), dimension(:), allocatable :: RY_HMX
  real(kind=rk) :: ft
  integer :: message,mpi_err,particle,k
  integer, dimension(mpi_status_size) :: mpi_status


  if(pfrank == 0) then
     !tell the other mpi processes to continue in this function...
     message = -1
     call mpi_bcast(message,1,mpi_integer,0,pf_mpi_comm,mpi_err)     
     
     !initialise f as background term
     f = 0.5d0*sum(v*v)*real(m-1,rk)


     !initialise tempg, the gradient, to zero
     tempg = 0.0d0


     !convert control variable v into model state:
     call convert_control_to_state(n,v,state_dim,xdet)


     !now make xt centered on xdet, the optimization state
     xt(:,1) = xdet(:)
     if(cnt .gt. 1) then
        do particle = 2,cnt
           xt(:,particle) = x0(:,particle-1) + xdet(:)
        end do
     end if


     !enter the model timestep loop
     do t = 1,vardata%total_timesteps

        !advance model to timestep t
        if(t == 1) then
           tag = 2
        else
           tag = 1
        end if


        do k =1,cnt
           particle = particles(k)
           call mpi_send(xt(:,k),state_dim,MPI_DOUBLE_PRECISION&
                &,particle-1,tag,CPL_MPI_COMM,mpi_err)
        end do
        
        
        DO k = 1,cnt
           particle = particles(k)
           CALL MPI_RECV(xt(:,k),state_dim, MPI_DOUBLE_PRECISION, &
                particle-1, MPI_ANY_TAG, CPL_MPI_COMM,mpi_status, mpi_err)
        END DO
        

        !check if we have observations this timestep
        if(vardata%ny(t) .gt. 0) then

           
           !if observations exist then allocate data
           allocate(     y(vardata%ny(t)))
           allocate(   hxt(vardata%ny(t),cnt))
           allocate(RY_HMX(vardata%ny(t)))

           !now let us compute the ensemble perturbation matrix

           !first we get the optimization soln to all other processes
           xdet = xt(:,1)
           call mpi_bcast(xdet,state_dim,mpi_double_precision,&
                0,pf_mpi_comm,mpi_err)              


           !subtract the optimization solution to form the perturbation matrix
           if(cnt .gt. 1) then
              do particle = 2,cnt
                 xt(:,particle) = (xt(:,particle) - xdet(:))
              end do
           end if
                      
           
           !get model equivalent of observations, store
           !in variable hxt
           call H(vardata%ny(t),cnt,xt,hxt,t)


           !add xdet back to the perturbation matrix to regain full ensemble
           if(cnt .gt. 1) then
              do particle = 2,cnt
                 xt(:,particle) = xt(:,particle) + xdet(:)
              end do
           end if
           

           !read in the data
           call get_observation_data(y,t)

           
           !compute difference in obs and model obs
           y = y - hxt(:,1)

           
           !compute R_i^{-1}( y_i - H_i M_i X)
           call solve_r(vardata%ny(t),1,y,RY_HMX,t)

           
           !pass this data to all other processes
           call mpi_bcast(RY_HMX,vardata%ny(t),mpi_double_precision,&
                0,pf_mpi_comm,mpi_err)

           
           !compute part of objective function corresponding
           !to the observations
           ft = sum(y*RY_HMX)


           !now compute the gradient part
           !-[H_i(X_i)]^TR_i^{-1}[Y_HMX]
           if(cnt .gt. 1) then
              do particle = 2,cnt
                 tempg(particles(particle)-1) =&
                      & tempg(particles(particle)-1) &
                      &- sum(hxt(:,particle)*RY_HMX)
              end do
           end if


           !deallocate data
           deallocate(y)
           deallocate(hxt)
           deallocate(RY_HMX)
           

           !update objective function with observation component
           !at this timestep
           f = f + 0.5d0*ft
        end if
        
     end do !end of the model timestep loop


     !reduce tempg to the master processor
     call mpi_reduce(tempg,g,n,MPI_DOUBLE_PRECISION,MPI_SUM&
          &,0,pf_mpi_comm,mpi_err)


     !add gradient of background term to g
     g = g + v*real(m-1,rk)
     
  else !(pfrank != 0)

     do
        call mpi_bcast(message,1,mpi_integer,0,pf_mpi_comm,mpi_err)
        if(message .eq. -1) then !compute the gradient etc


           !convert control variable to into model state
           !necessary here as some info only stored on this
           !process
           call convert_control_to_state(n,v,state_dim,xdet)

           !now make xt centered on xdet, the optimization state
           do particle = 1,cnt
              xt(:,particle) = x0(:,particle) + xdet(:)
           end do
           

           !enter the model timestep loop
           do t = 1,vardata%total_timesteps

              !advance model to timestep t
              if(t == 1) then
                 tag = 2
              else
                 tag = 1
              end if


              do k =1,cnt
                 particle = particles(k)
                 call mpi_send(xt(:,k),state_dim,MPI_DOUBLE_PRECISION&
                      &,particle-1,tag,CPL_MPI_COMM,mpi_err)
              end do

              
              DO k = 1,cnt
                 particle = particles(k)
                 CALL MPI_RECV(xt(:,k),state_dim, MPI_DOUBLE_PRECISION, &
                      particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
              END DO


              !check if we have observations this timestep
              if(vardata%ny(t) .gt. 0) then


                 !if observations exist then allocate data
                 allocate(   hxt(vardata%ny(t),cnt))
                 allocate(RY_HMX(vardata%ny(t)))
                 
                 !now let us compute the ensemble perturbation matrix

                 !first we get the optimization soln from process 0
                 call mpi_bcast(xdet,state_dim,mpi_double_precision,&
                      0,pf_mpi_comm,mpi_err)

                 !subtract the optimization solution to form the
                 !perturbation matrix
                 do particle = 1,cnt
                    xt(:,particle) = xt(:,particle) - xdet(:)
                 end do
                 
                 
                 !get model equivalent of observations, store
                 !in variable hxt
                 call H(vardata%ny(t),cnt,xt,hxt,t)


                 !add xdet back to the perturbation matrix to regain
                 !full ensemble         
                 do particle = 1,cnt
                    xt(:,particle) = xt(:,particle) + xdet(:)
                 end do
                                  

                 !get RY_HMX from master process
                 call mpi_bcast(RY_HMX,vardata%ny(t),mpi_double_precision,&
                      0,pf_mpi_comm,mpi_err)              
                 
                 
                 !now compute the gradient part
                 ![H_i(X_i)]^TR_i^{-1}[RY_HMX]
                 do particle = 2,cnt
                    tempg(particles(particle)) = tempg(particles(particle)) &
                         &- sum(hxt(:,particle)*RY_HMX)
                 end do
                 
                 
                 !deallocate data
                 deallocate(hxt)
                 deallocate(RY_HMX)
                 
              end if
           end do
           !reduce tempg to the master processor
           call mpi_reduce(tempg,g,n,MPI_DOUBLE_PRECISION,MPI_SUM&
                &,0,pf_mpi_comm,mpi_err)

        elseif(message .eq. 0) then
           exit !leave this function
        else !this shouldnt happen so is an error
           print*,'4DEnVar ERROR: mpi_bcast had unknown integer&
                & value',message
           stop 9
        end if

     end do
  end if
end subroutine fourdenvar_fcn


!> a subroutine to convert the optimization control variable to a model
!! state vector
!! 
!! this must be called by all processes on pf_mpi_comm
!! 
!! and the result x is known to all processes
subroutine convert_control_to_state(n,v,stateDim,x)
  use comms
  use fourdenvardata
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: n !< the dimension of the control variable
  real(kind=rk), dimension(n), intent(in) :: v !< the optimization control
  !<variable
  integer, intent(in) :: stateDim !< the dimension of the model state 
  real(kind=rk), dimension(stateDim), intent(out) :: x !<the
  !<resulting model state
  real(kind=rk), dimension(stateDim) :: tempx
  integer :: particle,mpi_err
  

  !communicate v to all the mpi processes
  call mpi_bcast(v,n,mpi_double_precision,&
       0,pf_mpi_comm,mpi_err)              
  

  !now do x = X0*v
  tempx = 0.0d0
  do particle = 1,cnt
     if(pfrank .ne. 0) then
        tempx = tempx + x0(:,particle  )*v(particles(particle  ))
     elseif(particle .gt. 1) then
        tempx = tempx + x0(:,particle-1)*v(particles(particle-1))
     end if
  end do


  !reduce this to all processors
  call mpi_allreduce(tempx,x,stateDim,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &pf_mpi_comm,mpi_err)


  ! now add the background term
  x = x+xb

end subroutine convert_control_to_state

