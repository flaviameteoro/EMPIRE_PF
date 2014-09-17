subroutine configure_model
  use pf_control
  use sizes
  use Qdata
  use Rdata
!  use hqht_plus_r
  implicit none
  include 'mpif.h'
  
  real(kind=kind(1.0d0)) :: t1
!this is for Lorenz 63
!  state_dim=3
!  obs_dim = 1

  !this is for hadcm3
  state_dim = 2314430
  obs_dim = 27370

  print*,'#################################'
  print*,'######### SANITY CHECK ##########'
  print*,'#################################'
  print*,'## STATE DIMENSION = ',state_dim
  print*,'##  OBS  DIMENSION = ',obs_dim
  print*,'#################################'
  if(.not. pf%gen_Q) then
     t1 = mpi_wtime()
     call loadQ
     print*,'load Q     took ',mpi_wtime()-t1,' seconds'
     t1 = mpi_wtime()
     call loadR
     print*,'load R     took ',mpi_wtime()-t1,' seconds'
     t1 = mpi_wtime()
!     call load_HQHTR
     print*,'load HQHTR took ',mpi_wtime()-t1,' seconds'
     
  end if
end subroutine configure_model



subroutine solve_r(y,v)
  !subroutine to take an observation vector y and return v
  !in observation space.
  use sizes
  use pf_control
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim,pf%count), intent(in) :: y
  real(kind=rk), dimension(obs_dim,pf%count), intent(out) :: v

  !R = diag(0.1)
  !v = y/0.3_rk
  v = 0.0_rk
  call daxpy(obs_dim*pf%count,1.0D0/0.3_rk,y,1,v,1)

  
end subroutine solve_r

subroutine solve_rhalf(y,v)
  !subroutine to take an observation vector y and return v
  !in observation space.
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim,pf%count), intent(in) :: y
  real(kind=rk), dimension(obs_dim,pf%count), intent(out) :: v

  !R^1/2 is diag(sqrt(0.1))
!  v = y/sqrt(0.3_rk)
  v = 0.0_rk
  call daxpy(obs_dim*pf%count,1.0D0/sqrt(0.3_rk),y,1,v,1)



end subroutine solve_rhalf


subroutine solve_hqht_plus_r(y,v)
  !subroutine to take an observation vector y and return v
  !in observation space.
!  use hsl_ma87_double
!  use hqht_plus_r
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim), intent(out) :: v

  !HQHt = (sqrt(2)*0.01)**2 = 0.0002
  !R = 2
  !thus (HQHt + R) = 2.0002

!  v = y/2.0002_rk
  v = y
!  call ma87_solve(v,HQHTR_order, HQHTR_keep, HQHTR_control, HQHTR_info)

!!$use sizes
!!$use hqht_plus_r
!!$implicit none
!!$integer, parameter :: rk=kind(1.0D+0)
!!$real(kind=rk), dimension(obs_dim), intent(in) :: y
!!$real(kind=rk), dimension(obs_dim), intent(out) :: v
!!$include 'mpif.h'
!!$real(kind=kind(1.0d0)), dimension(HQHTRnev) :: s
!!$real(kind=kind(1.d00)) :: start
!!$start = mpi_wtime()
!!$call dgemv('T',obs_dim,HQHTRnev,1.0D0,HQHTRU,obs_dim,y,1,0.0D0,s,1)
!!$s = s / HQHTRD
!!$call dgemv('N',obs_dim,HQHTRnev,1.0D0,HQHTRU,obs_dim,s,1,0.0D0,v,1)
!!$print*,'time spent in solve_hqhtr_plus_r = ',mpi_wtime()-start



end subroutine solve_hqht_plus_r

subroutine Q(nrhs,x,Qx)
  !subroutine to take a full state vector x and return Qx
  !in state space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: Qx
  real(kind=rk), dimension(state_dim,nrhs) :: temp

  call Qhalf(nrhs,x,temp)

  call Qhalf(nrhs,temp,Qx)
  

end subroutine Q
!!$
!!$subroutine Q(x,Qx)
!!$  use sizes
!!$  use Qevd
!!$  implicit none
!!$  include 'mpif.h'
!!$  real(kind=kind(1.0d0)), dimension(state_dim), intent(in) :: x
!!$  real(kind=kind(1.0d0)), dimension(state_dim), intent(out) :: qx
!!$  real(kind=kind(1.0d0)), dimension(Qnev) :: y
!!$  real(kind=kind(1.d00)) :: start
!!$  start = mpi_wtime()
!!$  call dgemv('T',state_dim,Qnev,1.0D0,QU,state_dim,x,1,0.0D0,y,1)
!!$  y = y * QD
!!$  call dgemv('N',state_dim,Qnev,1.0D0,QU,state_dim,y,1,0.0D0,Qx,1)
!!$  print*,'time spent in Q = ',mpi_wtime()-start
!!$end subroutine Q

subroutine Qhalf(nrhs,x,Qx)
  !subroutine to take a full state vector x and return Q^{1/2}x
  !in state space.
  use sizes
  use Qdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: qx

!  call coord_matmul_t(state_dim,Qne,Qrow,Qcol,Qval,x,qx)
  call openmp_matmul_t(state_dim,Qne,Qrow,Qcol,Qval,Qdiag,nrhs,x,qx)

  
end subroutine Qhalf


!!$subroutine Qhalf(x,Qx)
!!$  use sizes
!!$  use Qevd
!!$  implicit none
!!$  include 'mpif.h'
!!$  real(kind=kind(1.0d0)), dimension(Qnev) :: y
!!$  real(kind=kind(1.0d0)), dimension(state_dim), intent(in) :: x
!!$  real(kind=kind(1.0d0)), dimension(state_dim), intent(out) :: qx
!!$  real(kind=kind(1.d00)) :: start
!!$  start = mpi_wtime()
!!$  call dgemv('T',state_dim,Qnev,1.0D0,QU,state_dim,x,1,0.0D0,y,1)
!!$  y = y * sqrt(QD)
!!$  call dgemv('N',state_dim,Qnev,1.0D0,QU,state_dim,y,1,0.0D0,Qx,1)
!!$  print*,'time spent in Qhalf = ',mpi_wtime()-start
!!$end subroutine QHALF


subroutine R(nrhs,y,Ry)
  !subroutine to take an observation vector x and return Rx
  !in observation space.
  use sizes
  use Rdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs
  real(kind=rk), dimension(obs_dim,nrhs), intent(in) :: y
  real(kind=rk), dimension(obs_dim,nrhs), intent(out) :: Ry

!  call openmp_matmul_t(obs_dim,Rne,Rrow,Rcol,Rval,Rdiag,y,Ry)
  call openmp_matmul_t(obs_dim,Rne,Rrow,Rcol,Rval,Rdiag,nrhs,y,Ry)


end subroutine R



subroutine H(x,hx)
  !subroutine to take a full state vector x and return H(x)
  !in observation space.
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim,pf%count), intent(in) :: x
  real(kind=rk), dimension(obs_dim,pf%count), intent(out) :: hx

!  include 'H_output.f90'
  hx(:,:) = x(539617:566986,:)
  !now convert from celcius to Kelvin
!  hx = hx + 273.15_rk
end subroutine H

subroutine HT(y,x)
  !subroutine to take an observation vector y and return x = H^T(y)
  !in full state space.
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim,pf%count), intent(out) :: x
  real(kind=rk), dimension(obs_dim,pf%count), intent(in) :: y

!  x=273.15_rk
  x = 0.0_rk
  x(539617:566986,:) = y(:,:)
!  include 'HT_output.f90'
  !convert from Kelvin to celcius
!  x = x - 273.15_rk
end subroutine HT


!!$subroutine coord_matmul_t(n,ne,row,col,val,x,b)
!!$  integer, parameter :: rk = kind(1.0D+0)
!!$  integer, intent(in) :: n,ne
!!$  integer, intent(in), dimension(ne) :: row,col
!!$  real(kind = rk), dimension(ne), intent(in) :: val
!!$  real(kind = rk), dimension(n), intent(in) :: x
!!$  real(kind = rk), dimension(n), intent(out) :: b
!!$  integer :: k
!!$  b = 0.0_rk
!!$  !$omp parallel do
!!$  do k = 1,ne
!!$     if(row(k) .eq. col(k)) then
!!$!        !$OMP ATOMIC
!!$        b(row(k)) = b(row(k)) + val(k)*x(col(k))
!!$     else
!!$!        !$OMP ATOMIC
!!$        b(row(k)) = b(row(k)) + val(k)*x(col(k))
!!$!        !$OMP ATOMIC
!!$        b(col(k)) = b(col(k)) + val(k)*x(row(k))
!!$     end if
!!$  end do
!!$  !$omp end parallel do
!!$end subroutine coord_matmul_t

subroutine openmp_matmul_t(n,ne,row,col,val,diag,nrhs,x,b)
use omp_lib
integer, parameter :: rk = kind(1.0D+0)
integer, intent(in) :: n,ne,nrhs
integer, intent(in), dimension(ne) :: row,col
real(kind = rk), dimension(ne), intent(in) :: val
real(kind = rk), dimension(n), intent(in) :: diag
real(kind = rk), dimension(n,nrhs), intent(in) :: x
real(kind = rk), dimension(n,nrhs), intent(out) :: b
real(kind = rk), allocatable, dimension(:,:,:) :: temp
integer :: k,i,num_thr,thr,iam
include 'mpif.h'
real(kind=rk) :: t
real(kind=rk), dimension(5) :: ti
logical, parameter :: time = .false.

!print*,'in openmp'
!print*,'n = ',n
!print*,'ne = ',ne
!print*,'nrhs = ',nrhs
if(time) t = mpi_wtime()
b = 0.0_rk
if(time) ti(1) = mpi_wtime()
!first do the diagonal part
!$omp parallel do
do k = 1,n
   b(k,:) = b(k,:) + diag(k)*x(k,:)
end do
!$omp end parallel do
if(time) ti(2) = mpi_wtime()
!print*,'finished computing b'


!now do the off diagonal terms:
num_thr = omp_get_max_threads()
!print*,'num_thr = ',num_thr
allocate(temp(n,nrhs,num_thr))
temp = 0.0D0
if(time) ti(3) = mpi_wtime()
!print*,'before paralelel'
!$omp parallel private(iam)
iam = omp_get_thread_num()+1
!$omp do
do i = 1, ne
!   temp(row(i),:,omp_get_thread_num()+1) = temp(row(i),:&
!        &,omp_get_thread_num()+1) + val(i)*x(col(i),:)
!   temp(col(i),:,omp_get_thread_num()+1) = temp(col(i),:,omp_get_thread_num()&
!        &+1) + val(i)*x(row(i),:)

   temp(row(i),:,iam) = temp(row(i),:,iam) + val(i)*x(col(i),:)
   temp(col(i),:,iam) = temp(col(i),:,iam) + val(i)*x(row(i),:)

end do
!$omp end do
!$omp end parallel
if(time) ti(4) = mpi_wtime()

!print*,'before reduction'

!$omp parallel do private(thr)
do i = 1, n
   do thr = 1, num_thr
      b(i,:) = b(i,:) + temp(i,:,thr)
   end do
end do
!$omp end parallel do 
deallocate(temp)
if(time) then
   ti(5) = mpi_wtime()
   ti(5) = ti(5)-ti(4)
   ti(4) = ti(4)-ti(3)
   ti(3) = ti(3)-ti(2)
   ti(2) = ti(2)-ti(1)
   ti(1) = ti(1)-t
   write(6,'(4(es7.1,x),es7.1)') ti(1),ti(2),ti(3),ti(4),ti(5)
   print*,'###########'
end if
!print*,'finished ompenmp_matmul_y'
end subroutine openmp_matmul_t
