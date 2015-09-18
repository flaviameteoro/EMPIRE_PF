!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2015-09-18 15:40:32 pbrowne>
!!!
!!!    Module for doing things related to the LETKS: Local Ensemble
!!!    Kalman Transform Smoother
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

!> module for doing things related to the LETKS:
!>
module LETKS_data
  implicit none
  type LETKS_local
     integer :: red_obsdim
     real(kind=kind(1.0d0)), allocatable, dimension(:,:) :: USIUT
     real(kind=kind(1.0d0)), allocatable, dimension(:) :: Ud
  end type LETKS_local

  type (LETKS_local), allocatable, dimension(:) :: LSD

contains
  subroutine allocate_LETKS(N)
    integer, intent(in) :: N
    allocate(LSD(N))
  end subroutine allocate_LETKS

  subroutine deallocate_LETKS()
    deallocate(LSD)
  end subroutine deallocate_LETKS

  !> subroutine to compute the data for the LETKS, so that the
  !! increments can subsquently be computed
  !> @todo update to allow for non-diagonal R matrices to be used. 
  subroutine LETKS_filter_stage
    use comms
    use pf_control
    use sizes
    implicit none
    integer, parameter :: rk = kind(1.0d0)

    !!> Forecast ensemble on entry, analysis ensemble on exit
    !real(kind=rk), dimension(stateDimension,N), intent(inout) :: x
    real(kind=rk), dimension(state_dim,pf%count) :: x_loc

    !> The observation
    real(kind=rk), dimension(obs_dim) :: y

    ! Local variables for the SVD
    integer :: r
    real(kind=rk), dimension(:,:), allocatable :: V
    real(kind=rk), dimension(:), allocatable :: S
    real(kind=rk), dimension(pf%nens,pf%nens) :: UT
    !SMOOTHER ONLY:
    real(kind=rk), dimension(pf%nens,pf%nens) :: U
    integer :: LWORK,INFO
    real(kind=rk), dimension(:), allocatable :: WORK

    ! Miscellaneous local variables
    real(kind=rk), dimension(state_dim) :: mean_x
    real(kind=rk), dimension(:), allocatable :: mean_xa
    real(kind=rk), dimension(state_dim,pf%count) :: Xp_loc
    !real(kind=rk), dimension(stateDimension,N) :: Xp
    real(kind=rk), dimension(:,:), allocatable :: Xp
    real(kind=rk), dimension(:,:), allocatable :: Xa
    real(kind=rk), dimension(obs_dim) :: mean_yf,d,dd
    real(kind=rk), dimension(obs_dim_g) :: dd_g
    real(kind=rk), dimension(obs_dim,pf%nens) :: yf,Ysf
    real(kind=rk), dimension(obs_dim_g,pf%nens) :: Ysf_g
    real(kind=rk), dimension(obs_dim,pf%count) :: yf_loc
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

    if(empire_version .eq. 1 .or. empire_version .eq. 2) then
       ensemble_comm = pf_mpi_comm
    elseif(empire_version .eq. 3) then
       ensemble_comm = pf_ens_comm
    else
       print*,'EMPIRE VERSION ',empire_version,' NOT SUPPORTED IN letkf_analysis'
       print*,'THIS IS AN ERROR. STOPPING'
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
       Xp_loc(:,i) = pf%psi(:,i) - mean_x
    end do

    ! inflate the ensemble perturbation matrix
    Xp_loc = (1.0_rk + pf%rho) * Xp_loc
    ! store the local state vectors back in x_loc
    do i = 1,pf%count
       x_loc(:,i) = mean_x + Xp_loc(:,i)
    end do

    ! make the local ensemble perturbation matrix the correct scale
    Xp_loc = Xp_loc/sqrt(real(pf%nens-1,rk))

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
       start_var(i) = (i-1)*ceiling( real(state_dim,rk)/real(npfs,rk) )+1
       stop_var(i) = min( i*ceiling(real(state_dim,rk)/real(npfs,rk)) ,state_dim)
    end do

    !allocate space for Xp and Xa now that we know how many grid points we consider
    allocate(Xp(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
    allocate(Xa(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
    allocate(mean_xa(stop_var(pfrank+1)-start_var(pfrank+1)+1))
    mean_xa = mean_x(start_var(pfrank+1):stop_var(pfrank+1))

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
            ensemble_comm,&                                    !comm
            mpi_err)                                         !ierror
    end do


    d = y - mean_yf
    call solve_rhalf(obs_dim,1,d,dd,pf%timestep)

    if(empire_version .eq. 1 .or. empire_version .eq. 2) then
       dd_g = dd
       Ysf_g = Ysf
    elseif(empire_version .eq. 3) then
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


    !SMOOTHER ONLY:
    call allocate_LETKS(number_gridpoints)

    !$OMP PARALLEL DO &
    !$OMP& PRIVATE(stateDim,i,dist,scal,yes), &
    !$OMP& PRIVATE(Ysf_red,red_obsDim,r,V,S,LWORK,WORK,INFO), &
    !$OMP& PRIVATE(dd_red,UT,d)
    do j = 1,number_gridpoints
!              print*,'j = ',j,' Xp(j,1) = ',Xp(j,1)
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


       !SMOOTHER ONLY
       LSD(j)%red_obsdim = red_obsdim

       ! if there are no observations in range, treat this as a special case
       if(red_obsdim .gt. 0) then

          allocate(Ysf_red(red_obsdim,pf%nens))
          !multiply by the distance matrix
          !this line only works for diagonal R matrix...
          ! reduce the forecast ensemble to only the observations in range
          ! and scale by the distance function
          do i = 1,pf%nens
             Ysf_red(:,i) = pack(Ysf_g(:,i),yes)*pack(scal,yes)
          end do

          ! Compute the SVD
          r = min(red_obsDim,pf%nens)

          !SMOOTHER ONLY:
          allocate(LSD(j)%Ud(pf%nens))
          allocate(LSD(j)%USIUT(pf%nens,pf%nens))


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


          ! Build up analysis ensemble mean (to be stored in mean_x); done now
          ! because we shall be reusing X shortly

          allocate(dd_red(red_obsDim))
          dd_red = pack(dd_g,yes)*pack(scal,yes)


          ! Only the first r elements of d are in use from here on
          call dgemv('T',red_obsDim,r,1.0d0,V,red_obsDim,dd_red,1,0.0d0,d,1)

          d(1:r) = S * d(1:r) / (1.0_rk + S**2)

          !SMOOTHER ONLY: calculate U*d
          call dgemv('T',r,pf%nens,1.0d0,UT,r,d(1:r),1,0.0d0&
               &,LSD(j)%Ud,1)

!          print*,LSD(j)%Ud,matmul(transpose(UT),d(1:r))-LSD(j)%Ud
!          print*,dot_product(Xa(j,:),d(1:r)) + mean_xa(j)
!          print*,dot_product(Xp(j,:),LSD(j)%Ud)


          ! Only the first r columns of Xa are used in the following
          call dgemv('N',stateDim,r,1.0d0,Xa(j,:),stateDim,d,1,1.0d0,mean_xa(j),1)
!          print*,'j = ',j,' mean_xa(j) = ',mean_xa(j)

          ! Build up analysis ensemble perturbation matrix
          do i = 1,r
             Xa(j,i) = Xa(j,i) / sqrt(1.0_rk + S(i)**2)
          end do


          !SMOOTHER ONLY: calculate U*(S^2+I)^{-1/2}U^T
          U = transpose(UT)
          do i = 1,r
             U(:,i) = U(:,i) / sqrt(1.0_rk + S(i)**2)
          end do
          call dgemm('N','N',pf%nens,pf%nens,pf%nens,1.0d0,U,pf%nens&
               &,UT,pf%nens,0.0d0,LSD(j)%USIUT,pf%nens)
!          print*,'j = ',j,' XpUSLUT = ',matmul(Xp(j,:),LSD(j)%USIUT)

          call dgemm('N','N',stateDim,pf%nens,pf%nens,1.0d0,Xa(j,:),stateDim,UT,pf%nens,0.0d0,Xp(j,:),stateDim)

!          print*,'j = ',j,' Xp(j,:) = ',Xp(j,:)

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
             xa(j,i) =  mean_xa(j) + xp(j,i)
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
  end subroutine LETKS_filter_stage


  !> subroutine to compute the LETKS increments
  subroutine LETKS_increment(psi,inc)
    use comms
    use pf_control
    use sizes
    implicit none
    integer, parameter :: rk = kind(1.0d0)


    real(kind=rk), intent(in), dimension(state_dim,pf%count) :: psi
    !< input ensemble
    real(kind=rk), intent(out), dimension(state_dim,pf%count) :: inc
    !< LETKS increment


    !!> Forecast ensemble on entry, analysis ensemble on exit
    !real(kind=rk), dimension(stateDimension,N), intent(inout) :: x
    real(kind=rk), dimension(state_dim,pf%count) :: x_loc


    ! Local variables for the SVD
    integer :: r

    ! Miscellaneous local variables
    real(kind=rk), dimension(state_dim) :: mean_x
    real(kind=rk), dimension(:), allocatable :: mean_xa
    real(kind=rk), dimension(state_dim,pf%count) :: Xp_loc
    !real(kind=rk), dimension(stateDimension,N) :: Xp
    real(kind=rk), dimension(:,:), allocatable :: Xp
    real(kind=rk), dimension(:,:), allocatable :: Xa
    integer :: i,j,number_gridpoints

    !variables for localisation
    integer :: red_obsdim
    integer :: stateDim !generally 1 for the letkf

    !variables for mpi
    integer :: mpi_err
    integer, dimension(npfs) :: start_var,stop_var
    integer :: ensemble_comm
    include 'mpif.h'

    if(empire_version .eq. 1 .or. empire_version .eq. 2) then
       ensemble_comm = pf_mpi_comm
    elseif(empire_version .eq. 3) then
       ensemble_comm = pf_ens_comm
    else
       print*,'EMPIRE VERSION ',empire_version,' NOT SUPPORTED IN letks_increment'
       print*,'THIS IS AN ERROR. STOPPING'
       stop '-24'
    end if


    ! Split forecast ensemble into mean and perturbation matrix, inflating
    ! if necessary
    ! mean_x will only be the sum of state vectors on this mpi process
    mean_x = sum(psi,dim=2)

    !send mean_x to all processes and add up to get global sum
    call mpi_allreduce(MPI_IN_PLACE,mean_x,state_dim,MPI_DOUBLE_PRECISION&
         &,MPI_SUM,ensemble_comm,mpi_err)

    !now divide by the total number of ensemble members to make it a true
    !mean
    mean_x = mean_x/real(pf%nens,rk)

    ! compute the ensemble perturbation matrix for those ensemble members
    ! stored on this local mpi process
    do i = 1,pf%count
       Xp_loc(:,i) = psi(:,i) - mean_x
    end do

    ! inflate the ensemble perturbation matrix
    Xp_loc = (1.0_rk + pf%rho) * Xp_loc
    ! store the local state vectors back in x_loc
    do i = 1,pf%count
       x_loc(:,i) = mean_x + Xp_loc(:,i)
    end do

    ! make the local ensemble perturbation matrix the correct scale
    Xp_loc = Xp_loc/sqrt(real(pf%nens-1,rk))


    ! now let us compute which state variables will be analysed on each
    ! MPI process:
    do i = 1,npfs
       start_var(i) = (i-1)*ceiling( real(state_dim,rk)/real(npfs,rk) )+1
       stop_var(i) = min( i*ceiling(real(state_dim,rk)/real(npfs,rk)) ,state_dim)
    end do

    !allocate space for Xp and Xa now that we know how many grid points we consider
    allocate(Xp(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
    allocate(Xa(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
    allocate(mean_xa(stop_var(pfrank+1)-start_var(pfrank+1)+1))


    mean_xa = mean_x(start_var(pfrank+1):stop_var(pfrank+1))

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
            ensemble_comm,&                                    !comm
            mpi_err)                                         !ierror
    end do




    number_gridpoints = stop_var(pfrank+1)-start_var(pfrank+1)+1
    !$OMP PARALLEL DO &
    !$OMP& PRIVATE(stateDim,i), &
    !$OMP& PRIVATE(red_obsDim,r)
    do j = 1,number_gridpoints
!       print*,'j = ',j,' Xp(j,1) = ',Xp(j,1)
       stateDim = 1

       ! count the total number of observations we shall consider for this
       ! state variable
       red_obsdim = LSD(j)%red_obsdim


       ! if there are no observations in range, treat this as a special case
       if(red_obsdim .gt. 0) then


          ! compute
          r = min(red_obsDim,pf%nens)


          !LOOK UP U*d(1:r), call it Ud. It has dimension pf%nens         
          ! Only the first r columns of Xp are used in the following
          call dgemv('N',stateDim,pf%nens,1.0d0,Xp(j,:),stateDim,LSD(j)%Ud,1,1.0d0,mean_xa(j),1)

!          print*,'j = ',j,' mean_xa(j) = ',mean_xa(j)

          !LOOK UP U*(S^2+I)^{-1/2}U^T, call it USIUT.
          !IT has dimension (nens,nens)
          call dgemm('N','N',stateDim,pf%nens,pf%nens,1.0d0,Xp(j,:),stateDim,LSD(j)%USIUT,pf%nens,0.0d0,Xa(j,:),stateDim)

          ! Put ensemble mean and perturbation matrix back together
          xp(j,:) = sqrt(real(pf%nens-1,rk))*Xa(j,:)
          do i = 1,pf%nens
             xa(j,i) = mean_xa(j) + xp(j,i)
          end do


       else !if there are no observations near, the analysis is
          ! just the forecast
          do i = 1,pf%nens
             xa(j,i) =  mean_xa(j) + xp(j,i)
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
            inc(start_var(i):stop_var(i),:),&             !recvbuf  
            pf%count*number_gridpoints,&                     !recvcount  
            MPI_DOUBLE_PRECISION,&                           !recvtype    
            i-1,&                                            !root     
            ensemble_comm,&                                    !comm     
            mpi_err)                                         !ierror 
    end do

    inc = inc-psi

    deallocate(mean_xa)
    deallocate(Xp)
    deallocate(Xa)



  end subroutine LETKS_increment
end module LETKS_DATA

