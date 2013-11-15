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
  integer :: i,j,k,l
  integer :: mpi_err,particle,tag
  INTEGER, DIMENSION(:), ALLOCATABLE  :: requests
  INTEGER, DIMENSION(:,:), ALLOCATABLE  :: mpi_statuses
  logical :: mpi_flag
  logical, dimension(:), ALLOCATABLE :: received
  real(kind=kind(1.0D0)), dimension(o_nxn*o_nyn*o_levels) :: ocean_temp,sali,uu,vv

  logical, dimension(o_nxn,o_nyn,o_levels) :: tmask,smask,umask,vmask
  logical, dimension(o_nxn*o_nyn*o_levels) :: t1,s1,u1,v1

!  real(kind=kind(1.0D0)), dimension(o_nxn,o_nyn,o_levels) :: ocean_levels




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

!first we count how many blank entries there are:
tag = 0
do i = 1,state_dim
   if(pf%psi(i,1) .lt. -2.0**29) tag = tag+1
end do
print*,'total tag count which are not active = ',tag
print*,'compare this with the state dimension ',state_dim



open(7,file='tmask',status='replace')
open(8,file='smask',status='replace')
open(9,file='umask',status='replace')
open(10,file='vmask',status='replace')


ocean_temp = pf%psi(a_nxn*a_nyn*(1+4*a_levels)+1:a_nxn*a_nyn*(1+4&
     &*a_levels)+o_nxn*o_nyn*o_levels,1)
sali = pf%psi(a_nxn*a_nyn*(1+4*a_levels)+o_nxn*o_nyn*o_levels + 1:a_nxn*a_nyn&
     &*(1+4*a_levels)+2*o_nxn*o_nyn*o_levels,1)
uu = pf%psi(a_nxn*a_nyn*(1+4*a_levels)+2*o_nxn*o_nyn*o_levels + 1:a_nxn*a_nyn&
     &*(1+4*a_levels)+3*o_nxn*o_nyn*o_levels,1)
vv = pf%psi(a_nxn*a_nyn*(1+4*a_levels)+3*o_nxn*o_nyn*o_levels + 1:a_nxn*a_nyn&
     &*(1+4*a_levels)+4*o_nxn*o_nyn*o_levels,1)


do j = 1,o_nxn*o_nyn*o_levels
   if(ocean_temp(j) .lt. -2.0**29) then
      t1(j) = .false.
   else
      t1(j) = .true.
   end if

    if(sali(j) .lt. -2.0**29) then
      s1(j) = .false.
   else
      s1(j) = .true.
   end if
   if(uu(j) .lt. -2.0**29) then
      u1(j) = .false.
   else
      u1(j) = .true.
   end if
   if(vv(j) .lt. -2.0**29) then
      v1(j) = .false.
   else
      v1(j) = .true.
   end if
end do

tmask = reshape(t1,(/o_nxn,o_nyn,o_levels/) )
smask = reshape(s1,(/o_nxn,o_nyn,o_levels/) )
umask = reshape(u1,(/o_nxn,o_nyn,o_levels/) )
vmask = reshape(v1,(/o_nxn,o_nyn,o_levels/) )


write(7,*) tmask
write(8,*) smask
write(9,*) umask
write(10,*) vmask





close(7)
close(8)
close(9)
close(10)

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


