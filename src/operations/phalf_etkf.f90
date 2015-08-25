!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-08-25 10:26:50 pbrowne>
!!!
!!!    Routine to change an ensemble N(0,I) to N(0,P)
!!!    Copyright (C) 2015 Philip A. Browne
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


!> Subroutine to go from N(0,Q) to N(0,P)
subroutine phalf_etkf(nrhs,x,px)
  use pf_control
  use sizes
  use comms
  implicit none
  include 'mpif.h'

  integer, parameter :: rk=kind(1.0D0)

  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x !< the
  !!input vector, assumed to be columns N(0,I)
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: px !< the
  !!resulting vector where Px \f$= P^{\frac{1}{2}}x\f$

  real(kind=rk), dimension(state_dim,pf%count) :: Xp_loc
  real(kind=rk), dimension(obs_dim,pf%count) :: yf_loc
  real(kind=rk), dimension(obs_dim,pf%nens) :: yf,Ysf
  real(kind=rk), dimension(obs_dim) :: mean_yf
  integer, dimension(npfs) :: start_var,stop_var
  !real(kind=rk), dimension(stateDimension,N) :: Xp
  real(kind=rk), dimension(:,:), allocatable :: Xp
  real(kind=rk), dimension(:,:), allocatable :: Xa
  integer :: mpi_err
  integer :: i,j,number_gridpoints

  integer :: red_obsdim
  real(kind=rk), allocatable, dimension(:,:) :: Ysf_red

  ! Local variables for the SVD
  integer :: r
  real(kind=rk), dimension(:,:), allocatable :: V
  real(kind=rk), dimension(:), allocatable :: S
  real(kind=rk), dimension(pf%nens,pf%nens) :: UT
  integer :: LWORK,INFO
  real(kind=rk), dimension(:), allocatable :: WORK
  integer :: stateDim
  !!> Forecast ensemble on entry, analysis ensemble on exit
  !real(kind=rk), dimension(stateDimension,N), intent(inout) :: x
  real(kind=rk), dimension(state_dim,pf%count) :: x_loc
  
  !first lets make something N(0,Q):
  call Qhalf(nrhs,x,Xp_loc)


  !now go into the LETKF code:

  x_loc = Xp_loc


  ! Calculate forecast observations, split into mean and ensemble
  ! perturbation matrix, scale perturbations by inverse square root of
  ! observation covariance

  ! first apply observation operator only to local state vectors
  ! on this mpi process
  call H(obs_dim,pf%count,x_loc,yf_loc,pf%timestep)


  ! as yf_loc should be much smaller than x_loc, send this to mpi processes
  ! need to send round all yf_loc and store in yf on all processes
  call mpi_allgatherv(yf_loc,pf%count*obs_dim,MPI_DOUBLE_PRECISION,yf&
       &,gblcount*obs_dim,gbldisp*obs_dim,MPI_DOUBLE_PRECISION&
       &,pf_mpi_comm,mpi_err)

  ! compute the mean of yf
  mean_yf = sum(yf,dim=2)/real(pf%nens,rk)
  do i = 1,pf%nens
     yf(:,i) = yf(:,i) - mean_yf
  end do

  ! now make yf the forecast perturbation matrix
  yf = yf/sqrt(real(pf%nens-1,rk))
  ! scale yf to become ysf
  call solve_rhalf(obs_dim,pf%nens,yf,Ysf,pf%timestep)


  ! now let us compute which state variables will be analysed on each
  ! MPI process:
  do i = 1,npfs
     start_var(i) = (i-1)*ceiling( real(state_dim,rk)/real(npfs,rk) )+1
     stop_var(i) = min( i*ceiling(real(state_dim,rk)/real(npfs,rk)) ,state_dim)
  end do

  !allocate space for Xp and Xa now that we know how many grid points we consider
  allocate(Xp(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
  allocate(Xa(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))

  !now we have to get Xp filled with all the values from each process
  do i = 1,npfs
     number_gridpoints = stop_var(i)-start_var(i)+1

     call mpi_gatherv(Xp_loc(start_var(i):stop_var(i),:),& !sendbuf
          pf%count*number_gridpoints,&                     !sendcount
          MPI_DOUBLE_PRECISION,&                           !sendtype
          Xp,&                                             !recvbuf
          gblcount*number_gridpoints,&                     !recvcounts
          gbldisp*number_gridpoints,&                      !displs
          MPI_DOUBLE_PRECISION,&                           !recvtype
          i-1,&                                            !root
          pf_mpi_comm,&                                    !comm
          mpi_err)                                         !ierror
  end do


  !PAB  d = y - mean_yf
  !PAB  call solve_rhalf(obs_dim,1,d,dd,pf%timestep)



  !this is for serial processing
  number_gridpoints = stop_var(pfrank+1)-start_var(pfrank+1)+1
  !$OMP PARALLEL DO &
  !$OMP& PRIVATE(stateDim,i), &
  !$OMP& PRIVATE(Ysf_red,red_obsDim,r,V,S,LWORK,WORK,INFO), &
  !$OMP& PRIVATE(UT)
  do j = 1,number_gridpoints
     stateDim = 1

     ! count the total number of observations we shall consider for this
     ! state variable
     red_obsdim = obs_dim


     allocate(Ysf_red(red_obsdim,pf%nens))
     ! reduce the forecast ensemble to only the observations in range
     ! and scale by the distance function
     do i = 1,pf%nens
        Ysf_red(:,i) = Ysf(:,i)
     end do

     ! Compute the SVD
     r = min(red_obsDim,pf%nens)
     allocate(V(red_obsDim,r))
     allocate(S(r))
     LWORK = 2*max( 3*r+max(red_obsDim,pf%nens), 5*r )

     allocate(WORK(LWORK))

     call dgesvd('S','A',red_obsDim,pf%nens,Ysf_red,red_obsDim,S,V,red_obsDim,UT,pf%nens,WORK,LWORK,INFO)
     if(INFO .ne. 0) then
        print*,'SVD failed with INFO = ',INFO
        print*,'FYI WORK(1) = ',WORK(1)
        stop
     end if
     deallocate(WORK)

     ! Compute product of forecast ensemble perturbation matrix and U.  We
     ! store the result in xa, which we here call Xa to distinguish this
     ! secondary use of the variable.
     ! pick out the correct row of Xa and Xp here:
     call dgemm('N','T',stateDim,pf%nens,pf%nens,1.0d0,Xp(j,:),stateDim,UT,pf%nens,0.0d0,Xa(j,:),stateDim)


     ! Build up analysis ensemble perturbation matrix
     do i = 1,r
        Xa(j,i) = Xa(j,i) / sqrt(1.0_rk + S(i)**2)
     end do
     call dgemm('N','N',stateDim,pf%nens,pf%nens,1.0d0,Xa(j,:),stateDim,UT,pf%nens,0.0d0,Xp(j,:),stateDim)


     !put this back into the full state vector
     !wait it is already
     deallocate(Ysf_red)
     deallocate(V)
     deallocate(S)


  end do
  !$OMP END PARALLEL DO



  ! now we must scatter xa back into px
  do i = 1,npfs

     number_gridpoints = stop_var(i)-start_var(i)+1

     call mpi_scatterv(Xa,&                                !sendbuf  
          gblcount*number_gridpoints,&                     !sendcounts  
          gbldisp*number_gridpoints,&                      !displs
          MPI_DOUBLE_PRECISION,&                           !sendtype    
          px(start_var(i):stop_var(i),:),&             !recvbuf  
          pf%count*number_gridpoints,&                     !recvcount  
          MPI_DOUBLE_PRECISION,&                           !recvtype    
          i-1,&                                            !root     
          pf_mpi_comm,&                                    !comm     
          mpi_err)                                         !ierror 
  end do



  deallocate(Xp)
  deallocate(Xa)



end subroutine phalf_etkf
