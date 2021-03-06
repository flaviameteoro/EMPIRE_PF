!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-12-08 11:28:35 pbrowne>
!!!
!!!    Ensemble transform Kalman filter
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

!> subroutine to perform the ensemble transform Kalman filter as part
!! of L-ETKF
!>
!> @todo update to allow for non-diagonal R matrices to be used. 
subroutine letkf_analysis
  use output_empire, only : emp_e
  use timestep_data
  use comms
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)

  !!> Forecast ensemble on entry, analysis ensemble on exit
  !real(kind=rk), dimension(stateDimension,N), intent(inout) :: x
  real(kind=rk), dimension(state_dim,pf%count) :: x_p

  !> The observation
  real(kind=rk), dimension(obs_dim) :: y

  ! Local variables for the SVD
  integer :: r
  real(kind=rk), dimension(:,:), allocatable :: V
  real(kind=rk), dimension(:), allocatable :: S
  real(kind=rk), dimension(pf%nens,pf%nens) :: UT
  integer :: LWORK,INFO
  real(kind=rk), dimension(:), allocatable :: WORK

  ! Miscellaneous local variables
  real(kind=rk), dimension(state_dim) :: mean_x
  real(kind=rk), dimension(:), allocatable :: mean_xa
  real(kind=rk), dimension(state_dim,pf%count) :: Xp_p
  !real(kind=rk), dimension(stateDimension,N) :: Xp
  real(kind=rk), dimension(:,:), allocatable :: Xp
  real(kind=rk), dimension(:,:), allocatable :: Xa
  real(kind=rk), dimension(obs_dim) :: mean_yf,d,dd
  real(kind=rk), dimension(obs_dim_g) :: dd_g
  real(kind=rk), dimension(obs_dim,pf%nens) :: yf,Ysf
  real(kind=rk), dimension(obs_dim_g,pf%nens) :: Ysf_g
  real(kind=rk), dimension(obs_dim,pf%count) :: yf_p
  integer :: i,j,number_gridpoints

  !variables for localisation
  real(kind=rk) :: dist
  real(kind=rk), parameter :: maxscal=exp(8.0d0)
  real(kind=rk),dimension(obs_dim_g) :: scal
  logical, dimension(obs_dim_g) :: yes
  integer :: red_obsdim
  real(kind=rk), allocatable, dimension(:,:) :: Ysf_red
  real(kind=rk), allocatable, dimension(:) :: dd_red
  integer :: stateDim !generally 1 for the letkf

  !variables for mpi
  integer :: mpi_err
  integer, dimension(npfs) :: start_var,stop_var
  integer :: ensemble_comm
  include 'mpif.h'

  if(comm_version .eq. 1 .or. comm_version .eq. 2) then
     ensemble_comm = pf_mpi_comm
  elseif(comm_version .eq. 3) then
     ensemble_comm = pf_ens_comm
  else
     write(emp_e,*) 'EMPIRE VERSION ',comm_version,' NOT SUPPORTED IN letkf_analysis'
     write(emp_e,*) 'THIS IS AN ERROR. STOPPING'
     stop '-24'
  end if


  call get_observation_data(y,pf%timestep)

  ! Split forecast ensemble into mean and perturbation matrix, inflating
  ! if necessary
  ! mean_x will only be the sum of state vectors on this mpi process
  mean_x = sum(pf%psi,dim=2)

  !send mean_x to all processes and add up to get global sum
  call mpi_allreduce(MPI_IN_PLACE,mean_x,state_dim,MPI_DOUBLE_PRECISION&
       &,MPI_SUM,ensemble_comm,mpi_err)

  !now divide by the total number of ensemble members to make it a true
  !mean
  mean_x = mean_x/real(pf%nens,rk)

  ! compute the ensemble perturbation matrix for those ensemble members
  ! stored on this local mpi process
  do i = 1,pf%count
     Xp_p(:,i) = pf%psi(:,i) - mean_x
  end do

  ! inflate the ensemble perturbation matrix
  Xp_p = (1.0_rk + pf%rho) * Xp_p
  ! store the local state vectors back in x_p
  do i = 1,pf%count
     x_p(:,i) = mean_x + Xp_p(:,i)
  end do

  ! make the local ensemble perturbation matrix the correct scale
  Xp_p = Xp_p/sqrt(real(pf%nens-1,rk))

  ! Calculate forecast observations, split into mean and ensemble
  ! perturbation matrix, scale perturbations by inverse square root of
  ! observation covariance

  ! first apply observation operator only to local state vectors
  ! on this mpi process
  call H(obs_dim,pf%count,x_p,yf_p,pf%timestep)


  ! as yf_p should be much smaller than x_p, send this to mpi processes
  ! need to send round all yf_p and store in yf on all processes
  call mpi_allgatherv(yf_p,pf%count*obs_dim,MPI_DOUBLE_PRECISION,yf&
       &,gblcount*obs_dim,gbldisp*obs_dim,MPI_DOUBLE_PRECISION&
       &,ensemble_comm,mpi_err)

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
     start_var(i) = ((i-1)*state_dim)/npfs + 1 
     stop_var(i) = (i*state_dim)/npfs
  end do

  !allocate space for Xp and Xa now that we know how many grid points we consider
  allocate(Xp(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
  allocate(Xa(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
  allocate(mean_xa(stop_var(pfrank+1)-start_var(pfrank+1)+1))
  mean_xa = mean_x(start_var(pfrank+1):stop_var(pfrank+1))

  !now we have to get Xp filled with all the values from each process
  do i = 1,npfs
     number_gridpoints = stop_var(i)-start_var(i)+1

     call mpi_gatherv(Xp_p(start_var(i):stop_var(i),:),& !sendbuf
          pf%count*number_gridpoints,&                     !sendcount
          MPI_DOUBLE_PRECISION,&                           !sendtype
          Xp,&                                             !recvbuf
          gblcount*number_gridpoints,&                     !recvcounts
          gbldisp*number_gridpoints,&                      !displs
          MPI_DOUBLE_PRECISION,&                           !recvtype
          i-1,&                                            !root
          ensemble_comm,&                                    !comm
          mpi_err)                                         !ierror
  end do


  d = y - mean_yf
  call solve_rhalf(obs_dim,1,d,dd,pf%timestep)

  if(comm_version .eq. 1 .or. comm_version .eq. 2) then
     dd_g = dd
     Ysf_g = Ysf
  elseif(comm_version .eq. 3) then
     call mpi_allgatherv(dd,obs_dim,MPI_DOUBLE_PRECISION,dd_g,obs_dims&
          &,obs_displacements,MPI_DOUBLE_PRECISION,pf_member_comm&
          &,mpi_err)
     call mpi_allgatherv(Ysf,obs_dim*pf%nens,MPI_DOUBLE_PRECISION&
          &,Ysf_g,obs_dims*pf%nens,obs_displacements*pf%nens&
          &,MPI_DOUBLE_PRECISION,pf_member_comm,mpi_err)
  else
     print*,'should not enter here. error. stopping :('
  end if

  !this is for serial processing
  number_gridpoints = stop_var(pfrank+1)-start_var(pfrank+1)+1
  !$OMP PARALLEL DO &
  !$OMP& PRIVATE(stateDim,i,dist,scal,yes), &
  !$OMP& PRIVATE(Ysf_red,red_obsDim,r,V,S,LWORK,WORK,INFO), &
  !$OMP& PRIVATE(dd_red,UT,d)
  do j = 1,number_gridpoints
     stateDim = 1
     yes = .false.
     !let us process the observations here:
     do i = 1,obs_dim_g
        ! get the distance between the current state variable
        ! j+start_var(pfrank+1)-1
        ! and the observation i
        ! store it as distance
        call dist_st_ob(j+start_var(pfrank+1)-1,i,dist,pf%timestep)
        ! compute the scaling factor based in this distance and the
        ! length scale
        call loc_function(1,dist,scal(i),yes(i))
     end do

     ! count the total number of observations we shall consider for this
     ! state variable
     red_obsdim = count(yes)


     ! if there are no observations in range, treat this as a special case
     if(red_obsdim .gt. 0) then

        allocate(Ysf_red(red_obsdim,pf%nens))
        !multiply by the distance matrix
        !this line only works for diagonal R matrix...
        ! reduce the forecast ensemble to only the observations in range
        ! and scale by the distance function
        do i = 1,pf%nens
           Ysf_red(:,i) = pack(Ysf_g(:,i),yes)*sqrt(pack(scal,yes))
        end do

        ! Compute the SVD
        r = min(red_obsDim,pf%nens)
        allocate(V(red_obsDim,r))
        allocate(S(r))
        LWORK = 2*max( 3*r+max(red_obsDim,pf%nens), 5*r )

        allocate(WORK(LWORK))

        call dgesvd('S','A',red_obsDim,pf%nens,Ysf_red,red_obsDim,S,V,red_obsDim,UT,pf%nens,WORK,LWORK,INFO)
        if(INFO .ne. 0) then
           write(emp_e,*) 'EMPIRE ERROR IN LETKF WITH THE SVD'
           write(emp_e,*) 'SVD failed with INFO = ',INFO
           write(emp_e,*) 'FYI WORK(1) = ',WORK(1)
           write(emp_e,*) 'STOPPING'
           stop
        end if
        deallocate(WORK)

        ! Compute product of forecast ensemble perturbation matrix and U.  We
        ! store the result in xa, which we here call Xa to distinguish this
        ! secondary use of the variable.
        ! pick out the correct row of Xa and Xp here:
        call dgemm('N','T',stateDim,pf%nens,pf%nens,1.0d0,Xp(j,:),stateDim,UT,pf%nens,0.0d0,Xa(j,:),stateDim)


        ! Build up analysis ensemble mean (to be stored in mean_x); done now
        ! because we shall be reusing X shortly

        allocate(dd_red(red_obsDim))
        dd_red = pack(dd_g,yes)*sqrt(pack(scal,yes))

        ! Only the first r elements of d are in use from here on
        call dgemv('T',red_obsDim,r,1.0d0,V,red_obsDim,dd_red,1,0.0d0,d,1)

        d(1:r) = S * d(1:r) / (1.0_rk + S**2)

        ! Only the first r columns of Xa are used in the following
        call dgemv('N',stateDim,r,1.0d0,Xa(j,:),stateDim,d,1,1.0d0,mean_xa(j),1)

        ! Build up analysis ensemble perturbation matrix
        do i = 1,r
           Xa(j,i) = Xa(j,i) / sqrt(1.0_rk + S(i)**2)
        end do
        call dgemm('N','N',stateDim,pf%nens,pf%nens,1.0d0,Xa(j,:),stateDim,UT,pf%nens,0.0d0,Xp(j,:),stateDim)

        ! Put ensemble mean and perturbation matrix back together
        xa(j,:) = sqrt(real(pf%nens-1,rk))*Xp(j,:)
        do i = 1,pf%nens
           xa(j,i) = mean_xa(j) + xa(j,i)
        end do

        !put this back into the full state vector
        !wait it is already
        deallocate(Ysf_red)
        deallocate(V)
        deallocate(S)
        deallocate(dd_red)

     else !if there are no observations near, the analysis is
        ! just the forecast
        do i = 1,pf%nens
           xa(j,i) =  mean_xa(j) + xp(j,i)*sqrt(real(pf%nens-1,rk))
        end do

     end if
  end do
  !$OMP END PARALLEL DO



  ! now we must scatter xa back into pf%psi
  do i = 1,npfs

     number_gridpoints = stop_var(i)-start_var(i)+1

     call mpi_scatterv(Xa,&                                !sendbuf  
          gblcount*number_gridpoints,&                     !sendcounts  
          gbldisp*number_gridpoints,&                      !displs
          MPI_DOUBLE_PRECISION,&                           !sendtype    
          pf%psi(start_var(i):stop_var(i),:),&             !recvbuf  
          pf%count*number_gridpoints,&                     !recvcount  
          MPI_DOUBLE_PRECISION,&                           !recvtype    
          i-1,&                                            !root     
          ensemble_comm,&                                    !comm     
          mpi_err)                                         !ierror 
  end do


  deallocate(mean_xa)
  deallocate(Xp)
  deallocate(Xa)

  call timestep_data_set_is_analysis
end subroutine letkf_analysis
