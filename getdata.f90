!program to run the particle filter on the model HadCM3.
!this shall hopefully have minimal changes specific to the model.
!Simon Wilson and Philip Browne 2013
!----------------------------------------------------------------

program couple_pf
  use comms
  use pf_control
  use sizes
  use hadcm3_config
  implicit none
  include 'mpif.h'
  integer :: i,j,k,l,count,line
  integer :: mpi_err,particle,tag
  INTEGER, DIMENSION(:), ALLOCATABLE  :: requests
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: mpi_statuses
  logical :: mpi_flag
  logical, dimension(:), ALLOCATABLE :: received
!  real(kind=kind(1.0D0)), dimension(o_nxn*o_nyn*o_levels) :: ocean_temp
  real(kind=kind(1.0D0)), dimension(o_nxn,o_nyn,o_levels) :: o_u,o_v,o_theta&
       &,o_sal
  real(kind=kind(1.0D0)), dimension(a_nxn,a_nyn,a_levels) :: a_u,a_v,a_theta,a_q
  real(kind=kind(1.0D0)), dimension(a_nxn,a_nyn) :: a_pstar
  integer, dimension(o_nxn,o_nyn) :: bath,q_depth,p_depth
  integer, dimension(o_nyn*o_nxn) :: bog
  real(kind=kind(1.0)) :: aa1,bb1,cc1,dd1
  character(7) :: a1,b1,c1,d1
  

  write(6,'(A)') 'PF: Starting PF code'
  call flush(6)
  call initialise_mpi
  print*,'PF: configuring model'
  call configure_model
  print*,'PF: seting controls'
  call set_pf_controls
  print*,'allocating pf'
  call allocate_pf

  call random_seed_mpi(pfrank)
  write(6,*) 'PF: starting to receive from model'
! 1st call to model to get psi
!  do k =1,pf%count
!     call receive_from_model(pf%psi(:,k),pf%particles(k))

     !lets add some random noise to the initial conditions
!     call perturb_particle(pf%psi(:,k))

!  enddo
  print*,'launching the receives'
 ! print*,pf%count
 ! print*,pf%particles
 ! print*,allocated(requests)
  tag = 1        
  allocate(requests(pf%count))
  allocate(mpi_statuses(mpi_status_size,pf%count))
  allocate(received(pf%count))


do i = 1,50
  DO k = 1,pf%count
     particle = pf%particles(k)
     CALL MPI_IRECV(pf%psi(:,k), state_dim, MPI_DOUBLE_PRECISION, &
          particle-1, tag, CPL_MPI_COMM,requests(k), mpi_err)
  end DO
  print*,'receives launched'
  k = 0
  received = .false.
  do
     k = mod(k,pf%count)+1
!     print*,k
     if(.not. received(k)) then
        particle = pf%particles(k)
!        print*,particle ,'not received so testing it'
        call MPI_TEST(requests(k), mpi_flag, mpi_statuses(:,k), mpi_err)
        
        if(mpi_flag) then
!           PRINT*,'Particle filter ',pfrank,'has received initial state_v&
!                &ector over mpi from ensemble member ',particle
           received(k) = .true.
           call perturb_particle(pf%psi(:,k))
        end if
     end if
     if(all(received)) exit
  end do
  write(6,*) 'PF: All models received in pf couple' 
  call flush(6)

  do k=1,pf%count
     particle = pf%particles(k)
     call mpi_send(pf%psi(:,k), state_dim, MPI_DOUBLE_PRECISION,particle-1, 1, CPL_MPI_COMM, mpi_err)
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! lets work in here now.


open(4,file='bath.txt',action='read',status='old')
i = 0 
do line = 1,o_nyn*o_nxn/4
read(4,'(3(A7,f7.0,4x),A7,f7.0)') a1,aa1,b1,bb1,c1,cc1,d1,dd1
!print*,a,aa1,b,bb1,c,cc1,d,dd1
bog(i+1) = nint(aa1)
bog(i+2) = nint(bb1)
bog(i+3) = nint(cc1)
bog(i+4) = nint(dd1)
i = i + 4
end do
close(4)
line = 0
do j = 1,o_nyn
   do i = 1,o_nxn
      line = line + 1
      bath(i,j) = bog(line) 
   end do
end do

DO j = 1,o_nyn
   DO i = 1,o_nxn
      q_depth(I,J)=0
      p_depth(I,J)=bath(i,j)
   ENDDO
ENDDO

!  2ND, COMPUTE NUMBER OF VERTICAL LEVELS AT EACH U,V POINT
DO j=1,o_nyn-1
   DO i=1,o_nxn-1
      q_depth(i,j)=MIN(p_depth(i,j),&
     &          p_depth(i+1,j),p_depth(i,j+1),&
     &          p_depth(i+1,j+1))
  ENDDO
ENDDO


!first we count how many blank entries there are:
tag = 0
do i = 1,state_dim
   if(pf%psi(i,1) .lt. -2.0**29) tag = tag+1
end do
print*,'total tag count which are not active = ',tag
print*,'compare this with the state dimension ',state_dim


!let us put everything in its place:
count = 0
!first is PSTAR (2d)
do j = 1,a_nyn
   do i = 1,a_nxn
      count = count + 1
      a_pstar(i,j) = pf%psi(count,1)
   end do
end do
print*,'a_pstar finishes on the ',count,' element'
!second is  U
do k = 1,a_levels
   do j = 1,a_nyn
      do i = 1,a_nxn
         count = count + 1
         a_u(i,j,k) = pf%psi(count,1)
      end do
   end do
end do
print*,'a_u finishes on the ',count,' element'
!third is V
do k = 1,a_levels
   do j = 1,a_nyn
      do i = 1,a_nxn
         count = count + 1
         a_v(i,j,k) = pf%psi(count,1)
      end do
   end do
end do
print*,'a_v finishes on the ',count,' element'
!fouth is THETA
do k = 1,a_levels
   do j = 1,a_nyn
      do i = 1,a_nxn
         count = count + 1
         a_theta(i,j,k) = pf%psi(count,1)
      end do
   end do
end do
print*,'a_theta finishes on the ',count,' element'
!fifth is Q
do k = 1,a_levels
   do j = 1,a_nyn
      do i = 1,a_nxn
         count = count + 1
         a_q(i,j,k) = pf%psi(count,1)
      end do
   end do
end do
print*,'a_q finishes on the ',count,' element'
!now we are onto the ocean:
!sixth is THETA
print*,'first entry of SST is then:', count+1
do k = 1,o_levels
   do j = 1,o_nyn
      do i = 1,o_nxn
         if(k .lt. p_depth(i,j)) then
            count = count + 1
            o_theta(i,j,k) = pf%psi(count,1)
         end if
      end do
   end do
   if(k .eq. 1) print*,'final entry of SST is then:',count
end do
print*,'o_theta finished on the ',count,' element'
!seventh is SALINITY
do k = 1,o_levels
   do j = 1,o_nyn
      do i = 1,o_nxn
         if(k .lt. p_depth(i,j)) then
            count = count + 1
            o_sal(i,j,k) = pf%psi(count,1)
         end if
      end do
   end do
end do
print*,'o_sal finished on the ',count,' element'
!eighth is U
do k = 1,o_levels
   do j = 1,o_nyn
      do i = 1,o_nxn
         if(k .lt. q_depth(i,j)) then
            count = count + 1
            o_u(i,j,k) = pf%psi(count,1)
         end if
      end do
   end do
end do
print*,'o_u finished on the ',count,' element'
!ninth is V
do k = 1,o_levels
   do j = 1,o_nyn
      do i = 1,o_nxn
         if(k .lt. q_depth(i,j)) then
            count = count + 1
            o_v(i,j,k) = pf%psi(count,1)
         end if
      end do
   end do
end do
print*,'o_v finished on the ',count,' element'
























!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call deallocate_data

  write(6,*) 'PF: finished deallocate_data - off to mpi_finalize'
  call flush(6)

  STOP
  call MPI_Finalize(mpi_err)
  deallocate(requests)
  deallocate(mpi_statuses)
  deallocate(received)
  write(*,*) 'Program couple_pf terminated successfully.'









end program couple_pf


