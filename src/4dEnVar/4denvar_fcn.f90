subroutine fcn( n, x, f, g )
implicit none
integer :: n
real(kind=kind(1.0d0)), dimension(n), intent(in) :: x
real(kind=kind(1.0d0)), intent(out) :: f
real(kind=kind(1.0d0)), dimension(n), intent(out) :: g


call fourdenvar_fcn(n,x,f,g)

end subroutine fcn



subroutine fourdenvar_fcn(n, v, f, g )
  use comms
  use var_data
  use fourdenvardata
  use sizes, only : state_dim
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: n !< this is the dimension of the
  !optimization control variable, so the number of ensemble members-1
  real(kind=rk), dimension(n), intent(in) :: v!< this is the
  !optimization control variable
!  real(kind=rk), dimension(state_dim) :: xt
  real(kind=rk), intent(out) :: f !< the 4dvar objective function
  real(kind=kind(1.0d0)), dimension(n), intent(out) :: g !< the
  !<gradient of the 4dvar objective function
  real(kind=kind(1.0d0)), dimension(n) :: tempg 
  integer :: t
  integer :: tag
  real(kind=rk), dimension(:), allocatable :: y
  real(kind=rk), dimension(:,:), allocatable :: hxt
  real(kind=rk), dimension(state_dim) :: xbar
  real(kind=rk), dimension(state_dim) :: tempx
  real(kind=rk), dimension(state_dim) :: xdet
  real(kind=rk), dimension(:), allocatable :: RY_HMX
  real(kind=rk) :: ft
  integer :: message,mpi_err,particle,k
  integer, dimension(mpi_status_size) :: mpi_status
!  print*,'   v = ',v
  if(pfrank == 0) then
     !tell the other mpi processes to continue in this function...
     message = -1
     call mpi_bcast(message,1,mpi_integer,0,pf_mpi_comm,mpi_err)     
     
     !initialise f as background term
     f = 0.5d0*sum(v*v)

 !    print*,0.5d0*sum(v*v)
     !initialise tempg, the gradient, to zero
     tempg = 0.0d0


     !  call inner_b_1(x-vardata%x0,f)
     
     !convert control variable v into model state:
     call convert_control_to_state(n,v,state_dim,xdet)
!     print*,'xdet = ',xdet

     !now make xt centered on xdet, the optimization state
     xt(:,1) = xdet(:)
     if(cnt .gt. 1) then
        do particle = 2,cnt
           xt(:,particle) = x0(:,particle-1) + xdet(:)
        end do
     end if

     
     do t = 1,vardata%total_timesteps
!        print*,'timestep = ',t
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
           !
           !first we get the mean, stored in xbar
           if(pfrank .ne. 0) then
              tempx = sum(xt,dim=2)
           else
              tempx = sum(xt(:,2:cnt),dim=2)
           end if
           
           call mpi_allreduce(tempx,xbar,state_dim,MPI_DOUBLE_PRECISION,MPI_SUM,&
                &pf_mpi_comm,mpi_err)

           
           xbar = xbar/real(m,rk)
!           print*,'xbar=',xbar

           !subtract the mean to form the perturbation matrix (not
           ! from the optimization solution
!           xt(:,1) = xdet(:)
           if(cnt .gt. 1) then
              do particle = 2,cnt
                 xt(:,particle) = (xt(:,particle) - xbar(:)) /(real(m-1,rk)**0.5d0)
              end do
           end if
           
           
           
           !get model equivalent of observations, store
           !in variable hxt
!           print*,'x=',xt(:,1)
           call H(vardata%ny(t),cnt,xt,hxt,t)
!           print*,'hx',hxt(:,1)

           !add the mean back to the perturbation matrix (not
           ! from the optimization solution to regain full ensemble
           !           xt(:,1) = xdet(:)
           if(cnt .gt. 1) then
              do particle = 2,cnt
                 xt(:,particle) = (real(m-1,rk)**0.5d0)*xt(:,particle) + xbar(:)
              end do
           end if
           

           !read in the data
           call get_observation_data(y,t)
!           print*,'y=',y
           
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
!           print*,ft

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
        
        
        
     end do
     !reduce tempg to the master processor
     call mpi_reduce(tempg,g,n,MPI_DOUBLE_PRECISION,MPI_SUM&
          &,0,pf_mpi_comm,mpi_err)

     !add gradient of background term to tempg
!     print*,'v = ',v
!     print*,'g = ',g
     g = g + v
!     print*,'n = ',g
     
     
  else !(pfrank != 0)

     do
        call mpi_bcast(message,1,mpi_integer,0,pf_mpi_comm,mpi_err)
        if(message .eq. -1) then !compute the gradient etc

           call convert_control_to_state(n,v,state_dim,xdet)


           !now make xt centered on xdet, the optimization state
           do particle = 1,cnt
              xt(:,particle) = x0(:,particle) + xdet(:)
           end do
           


           do t = 1,vardata%total_timesteps
              !advance model to timestep t
              if(t == 1) then
                 tag = -1
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
                 !
                 !first we get the mean, stored in xbar
                 if(pfrank .ne. 0) then
                    tempx = sum(xt,dim=2)
                 else
                    tempx = sum(xt(:,2:cnt),dim=2)
                 end if
                 
                 call mpi_allreduce(tempx,xbar,state_dim&
                      &,MPI_DOUBLE_PRECISION,MPI_SUM,&
                      &pf_mpi_comm,mpi_err)

                 xbar = xbar/real(m,rk)
                 
                 !now subtract xbar from xt...
                 do particle = 1,cnt
                    xt(:,particle) = xt(:,particle) - xbar(:)
                 end do
                 
                 
                 !get model equivalent of observations, store
                 !in variable hxt
                 call H(vardata%ny(t),cnt,xt,hxt,t)
                 
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
!print*,'f=',f
!call convert_control_to_state(n,g,state_Dim,xbar)
!print*,xbar
end subroutine fourdenvar_fcn

subroutine convert_control_to_state(n,v,stateDim,x)
  use comms
  use fourdenvardata
  implicit none
  include 'mpif.h'
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: n !< the dimension of the control variable
  real(kind=rk), dimension(n) :: v !< the optimization control
  !<variable
  integer, intent(in) :: stateDim !< the dimension of the model state 
  real(kind=rk), dimension(stateDim), intent(out) :: x !<the
  !<resulting model state
  real(kind=rk), dimension(stateDim) :: tempx
  integer :: particle,mpi_err
  
!  print*,'  v1 = ',v

  !communicate v to all the mpi processes
  call mpi_bcast(v,n,mpi_double_precision,&
       0,pf_mpi_comm,mpi_err)              
  
!  print*,'  v2 = ',v
!  print*,'  x0 = ',x0
  !now do x = X0*v
  tempx = 0.0d0
  do particle = 1,cnt
     if(pfrank .ne. 0) then
        tempx = tempx + x0(:,particle  )*v(particles(particle  ))
     elseif(particle .gt. 1) then
        tempx = tempx + x0(:,particle-1)*v(particles(particle-1))
     end if
  end do

!  print*,'temp = ',tempx

  call mpi_allreduce(tempx,x,stateDim,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &pf_mpi_comm,mpi_err)

!  print*,'xxxx = ',x

  ! now add the background term
  x = x+xb
!  print*,'xnew = ',x
end subroutine convert_control_to_state

