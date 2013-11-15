subroutine genQ
  use sizes
  use pf_control
  use hadcm3_config
  use comms
  use Qdata
  use hqht_plus_r
  implicit none
  include 'mpif.h'
  integer, parameter :: rk=(kind(1.0d0))
  !real(kind=rk), dimension(a_nxn,a_nyn,a_levels) :: a_u,a_v,a_theta,a_q
  !real(kind=rk), dimension(a_nxn,a_nyn) :: a_pstar
  !real(kind=rk), dimension(o_nxn,o_nyn,o_levels) :: o_u,o_v,o_theta,o_sal
  integer, dimension(a_nxn,a_nyn,a_levels) :: a_u_vec,a_v_vec,a_theta_vec,a_q_vec
  integer, dimension(a_nxn,a_nyn) :: a_pstar_vec
  integer, dimension(o_nxn,o_nyn,o_levels) :: o_u_vec,o_v_vec,o_theta_vec,o_sal_vec
  integer :: i,j,k,count,radius,nnz
  integer, parameter :: n = 426655238
  integer, allocatable, dimension(:) :: row,col
  real(kind=rk), allocatable, dimension(:) :: val
  integer :: ne,iter,day,tag,mpi_err,days
  integer :: mpi_status(MPI_STATUS_SIZE)
  real(kind=rk), dimension(state_dim) :: x 
  real(kind=rk) :: start_t,end_t

  allocate(pf%mean(state_dim),row(n),col(n),val(n))
  print*,'Going to generate Q ~ the model error covariance matrix'
  start_t = mpi_wtime()
  call load_bathimitry
  end_t = mpi_wtime()
  print*,'Load bathimitry took ',end_t - start_t,' seconds'

  start_t = end_t
  !let us calculate the position vectors: put everything in its place:
  count = 0
  !first is PSTAR (2d)
  do j = 1,a_nyn
     do i = 1,a_nxn
        count = count + 1
        a_pstar_vec(i,j) = count
     end do
  end do
  print*,'a_pstar finishes on the ',count,' element'
  !second is  U
  do k = 1,a_levels
     do j = 1,a_nyn
        do i = 1,a_nxn
           count = count + 1
           a_u_vec(i,j,k) = count
        end do
     end do
  end do
  print*,'a_u finishes on the ',count,' element'
  !third is V
  do k = 1,a_levels
     do j = 1,a_nyn
        do i = 1,a_nxn
           count = count + 1
           a_v_vec(i,j,k) = count
        end do
     end do
  end do
  print*,'a_v finishes on the ',count,' element'
  !fouth is THETA
  do k = 1,a_levels
     do j = 1,a_nyn
        do i = 1,a_nxn
           count = count + 1
           a_theta_vec(i,j,k) = count
        end do
     end do
  end do
  print*,'a_theta finishes on the ',count,' element'
  !fifth is Q
  do k = 1,a_levels
     do j = 1,a_nyn
        do i = 1,a_nxn
           count = count + 1
           a_q_vec(i,j,k) = count
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
           if(k .le. p_depth(i,j)) then
              count = count + 1
              o_theta_vec(i,j,k) = count
           else
              o_theta_vec(i,j,k) = 0
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
           if(k .le. p_depth(i,j)) then
              count = count + 1
              o_sal_vec(i,j,k) = count
           else
              o_sal_vec(i,j,k) = 0
           end if
        end do
     end do
  end do
  print*,'o_sal finished on the ',count,' element'
  !eighth is U
  do k = 1,o_levels
     do j = 1,o_nyn
        do i = 1,o_nxn
           if(k .le. q_depth(i,j)) then
              count = count + 1
              o_u_vec(i,j,k) = count
           else
              o_u_vec(i,j,k) = 0
           end if
        end do
     end do
  end do
  print*,'o_u finished on the ',count,' element'
  !ninth is V
  do k = 1,o_levels
     do j = 1,o_nyn
        do i = 1,o_nxn
           if(k .le. q_depth(i,j)) then
              count = count + 1
              o_v_vec(i,j,k) = count
           else
              o_v_vec(i,j,k) = 0
           end if
        end do
     end do
  end do
  print*,'o_v finished on the ',count,' element'
  !FINISHED COMPUTING THE VECTOR POSITIONS OF EACH COMPONENT
  print*,'FINISHED COMPUTING THE VECTOR POSITIONS OF EACH COMPONENT'
  end_t = mpi_wtime()
  print*,'and it took ',end_t-start_t,' seconds'
!  stop
  !initialise the value thingy
  val = 0.0_rk

  !initialise the mean
  pf%mean = 0.0_rk
  tag = 1
  print*,'first recv...'
  call mpi_recv(pf%psi(:,1), state_dim, MPI_DOUBLE_PRECISION, &
       0, tag, CPL_MPI_COMM,mpi_status, mpi_err)
  print*,'.............recvd'
  
  !loop for 5 years:
  days = 360*5
  do day = 1,days
     start_t = mpi_wtime()
     print*,'day = ',day
     !update by a day, or 72 iterations:
     do iter = 1,72
        call mpi_send(pf%psi(:,1), state_dim , MPI_DOUBLE_PRECISION, &
             0, tag, CPL_MPI_COMM, mpi_err)
        print*,'iter = ',iter
        call mpi_recv(pf%psi(:,1), state_dim, MPI_DOUBLE_PRECISION, &
             0, tag, CPL_MPI_COMM,mpi_status, mpi_err)
     end do
     end_t = mpi_wtime()
     print*,'72 model runs on ',day,' took ',end_t-start_t,' seconds'
     !update the mean:
     pf%mean = pf%mean + pf%psi(:,1)

     x = pf%psi(:,1)

     radius = 2
     start_t = mpi_wtime()
     print*,'going into genQ_order'
     include 'genQ_order.f90'
     print*,'finished genQ_order at the end of the day with ne = ',ne
     end_t = mpi_wtime()
     print*,'generating Q on ',day,' took ',end_t-start_t,' seconds.'

  end do !end of the daily loop



  !now we should divide the value by N-1 = 360*5-1...
  val = val/real(days-1,rk)


  !now subtract off the mean components:
  pf%mean = pf%mean/real(days,rk)


  val = -1.0_rk*val
  x =  pf%mean
  include 'genQ_order.f90'
  val = -1.0_rk*val


  print*,'MAXIMUM COVARIANCE VALUE = ',maxval(val)
  print*,'MIMINUM COVARIANCE VALUE = ',minval(val)
  print*,'MINIMUM ABSOLUTE C VALUE = ',minval(abs(val))







  !do the final mpi send to end model cleanly
  call mpi_send(pf%psi(:,1), state_dim , MPI_DOUBLE_PRECISION, &
       0, tag, CPL_MPI_COMM, mpi_err)


  !let us check the diagonal terms:
!  print*,'checking diagonal entries'
!  do i = 1,ne
!     if(col(i) .eq. row(i)) then
!        if(abs(val(i)) .lt. 1.0D-13) then
!           print*,'diagonal term at i = ',i,' in row ',row(i),' is ',val(i)
!        end if
!     end if
!  end do
!  print*,'diagonal entries checked'

  

  !now let us scale Q
  val = val/1.0D4







  nnz = 0
  do i = 1,ne
     if(row(i) .eq. col(i) .or. abs(val(i)) .gt. 1.0D-13) nnz = nnz + 1
  end do

  open(2,file='Qdata.dat',action='write',status='replace',form='unformatted')
  write(2) state_dim
  write(2) nnz

  do i = 1,ne
     if(row(i) .eq. col(i) .or. abs(val(i)) .gt. 1.0D-13) write(2) val(i)
  end do
  do i = 1,ne
     if(row(i) .eq. col(i) .or. abs(val(i)) .gt. 1.0D-13) write(2) row(i)
  end do
  do i = 1,ne
     if(row(i) .eq. col(i) .or. abs(val(i)) .gt. 1.0D-13) write(2) col(i)
  end do
  close(2)
  print*,'finished generating Q'
  !put it into the data structure

!!$
!!$  call loadQ
!!$  call inflate
!!$
!!$  open(2,file='QdataNEW.dat',action='write',status='replace',form&
!!$         &='unformatted')
!!$  write(2) state_dim
!!$  write(2) Qne
!!$
!!$  do i = 1,Qne
!!$     write(2) Qval(i)
!!$  end do
!!$  do i = 1,Qne
!!$     write(2) Qrow(i)
!!$  end do
!!$  do i = 1,Qne
!!$     write(2) Qcol(i)
!!$  end do
!!$  close(2)
!!$  print*,'finished inflating Q'
!!$
!!$  stop
!!$
!!$  call calc_HQHTR
!!$  
!!$  print*,'factored (HQHt+R) and wrote HQHTRdata.dat'

  deallocate(pf%mean,row,col,val)



  


end subroutine genQ


subroutine genQ_at2d_at2d(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn), intent(in) :: field1,field2
  integer :: i,j,ii,jj

  do j = 1,a_nyn
     do i = 1,a_nxn
        do jj = max(1,j-radius),min(a_nyn,j+radius)
           do ii = max(1,i-radius),min(a_nxn,i+radius)
              if(field1(i,j) .le. field2(ii,jj)) then
                 ne = ne + 1
                 row(ne) = field1(i,j)
                 col(ne) = field2(ii,jj)
                 val(ne) = val(ne)+x(row(ne))*x(col(ne))
              end if
           end do
           if(i .le. radius) then
              do ii = a_nxn-radius+i,a_nxn
                 if(field1(i,j) .le. field2(ii,jj)) then
                    ne = ne + 1
                    row(ne) = field1(i,j)
                    col(ne) = field2(ii,jj)
                    val(ne) = val(ne)+x(row(ne))*x(col(ne))
                 end if
              end do
           end if
           if(i+radius .gt. a_nxn) then
              do ii = 1,i+radius-a_nxn
                 if(field1(i,j) .le. field2(ii,jj)) then
                    ne = ne + 1
                    row(ne) = field1(i,j)
                    col(ne) = field2(ii,jj)
                    val(ne) = val(ne)+x(row(ne))*x(col(ne))
                 end if
              end do
           end if
        end do
     end do
  end do
end subroutine genQ_at2d_at2d


subroutine genQ_at2d_atq(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn), intent(in) :: field1
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field2
  integer :: i,j,ii,jj,kk

  do j = 1,a_nyn
     do i = 1,a_nxn
        do kk = 1,1+radius
           do jj = max(1,j-radius),min(a_nyn,j+radius)
              do ii = max(1,i-radius),min(a_nxn,i+radius)
                 if(field1(i,j) .le. field2(ii,jj,kk)) then
                    ne = ne + 1
                    row(ne) = field1(i,j)
                    col(ne) = field2(ii,jj,kk)
                    val(ne) = val(ne)+x(row(ne))*x(col(ne))
                 end if
              end do

              if(i .le. radius) then
                 do ii = a_nxn-radius+i,a_nxn
                    if(field1(i,j) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end do
              end if
              if(i+radius .gt. a_nxn) then
                 do ii = 1,i+radius-a_nxn
                    if(field1(i,j) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end do
              end if

           end do
        end do
     end do
  end do
end subroutine genQ_at2d_atq

subroutine genQ_at2d_atu(ne,row,col,val,n,radius,field1,field2,x)
  !field1 is offset of field2

  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn), intent(in) :: field1
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field2
  integer :: i,j,ii,jj,kk

  do j = 1,a_nyn
     do i = 1,a_nxn
        do kk = 1,1+radius
           !jj is now offset: so go from j-radius
           !up to j+radius-1
           !stop at a_nyn-1 as there are only that many in field2
           do jj = max(1,j-radius),min(a_nyn-1,j+radius-1)
              !ii is offset so go from i-radius to i+radius-1
              !stop at a_nxn as there are enough in that direction
              do ii = max(1,i-radius),min(a_nxn,i+radius-1)
                 if(field1(i,j) .le. field2(ii,jj,kk)) then
                    ne = ne + 1
                    row(ne) = field1(i,j)
                    col(ne) = field2(ii,jj,kk)
                    val(ne) = val(ne)+x(row(ne))*x(col(ne))
                 end if
              end do

              if(i .le. radius) then
                 do ii = a_nxn-radius+i,a_nxn
                    if(field1(i,j) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end do
              end if
              if(i+radius-1 .gt. a_nxn) then
                 do ii = 1,i+radius-a_nxn-1
                    if(field1(i,j) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end do
              end if

           end do
        end do
     end do
  end do
end subroutine genQ_at2d_atu

subroutine genQ_atq_atq(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field1
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,a_levels
     do j = 1,a_nyn
        do i = 1,a_nxn
           do kk = max(1,k-radius),min(a_levels,k+radius)
              do jj = max(1,j-radius),min(a_nyn,j+radius)
                 do ii = max(1,i-radius),min(a_nxn,i+radius)
                    if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j,k)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end do

                 if(i .lt. radius) then
                    do ii = a_nxn+1-radius+i,a_nxn
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end do
                 end if
                 if(i+radius .gt. a_nxn) then
                    do ii = 1,i+radius-a_nxn-1
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end do
                 end if

              end do
           end do
        end do
     end do
  end do
end subroutine genQ_atq_atq



subroutine genQ_atu_atu(ne,row,col,val,n,radius,field1,field2,x)
!both fields on u grid
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field1
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,a_levels
     do j = 1,a_nyn-1
        do i = 1,a_nxn
           do kk = max(1,k-radius),min(a_levels,k+radius)
              do jj = max(1,j-radius),min(a_nyn-1,j+radius)
                 do ii = max(1,i-radius),min(a_nxn,i+radius)
                    if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j,k)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end do

                 if(i .le. radius) then
                    do ii = a_nxn-radius+i,a_nxn
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end do
                 end if
                 if(i+radius .gt. a_nxn) then
                    do ii = 1,i+radius-a_nxn
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end do
                 end if

              end do
           end do
        end do
     end do
     !now do the top polar values
     !these are blank, so lets set them to be 1 each time we pass...
     j = a_nyn
     do i = 1,a_nxn
        if(field1(i,j,k) .eq. field2(i,j,k)) then
           ne = ne + 1
           row(ne) = field1(i,j,k)
           col(ne) = field2(i,j,k)
           val(ne) = 1.0D0
        end if
     end do
  end do
end subroutine genQ_atu_atu

subroutine genQ_atu_atq(ne,row,col,val,n,radius,field1,field2,x)
!both fields on u grid
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field1
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,a_levels
     !first entry is on u grid so only a_nyn-1 points in that direction
     do j = 1,a_nyn-1
        do i = 1,a_nxn
           do kk = max(1,k-radius),min(a_levels,k+radius)
              !jj are offset so plus 1 on the lower bound
              do jj = max(1,j-radius+1),min(a_nyn,j+radius)
                 !ii are offset so plus 1 on the lower bound
                 do ii = max(1,i-radius+1),min(a_nxn,i+radius)
                    if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j,k)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end do

                 !use lt here because of offset
                 if(i .lt. radius) then
                    do ii = a_nxn+1-radius+i,a_nxn
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end do
                 end if
                 if(i+radius .gt. a_nxn) then
                    do ii = 1,i+radius-a_nxn
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end do
                 end if

              end do
           end do
        end do
     end do
     !as we are using different fields, dont worry about the
     !diagonal term being included here
  end do
end subroutine genQ_atu_atq

subroutine genQ_ocq_ocq(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,o_levels
     do j = 1,o_nyn
        do i = 1,o_nxn-2
           do kk = max(1,k-radius),min(o_levels,k+radius)
              do jj = max(1,j-radius),min(o_nyn,j+radius)
                 do ii = max(1,i-radius),min(o_nxn-2,i+radius)
                    if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                         & .gt. 0) then
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if
                 end do

                 if(i .lt. radius) then
                    do ii = o_nxn-2+1-radius+i,o_nxn-2
                       if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                            & .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if
                 if(i+radius .gt. o_nxn-2) then
                    do ii = 1,i+radius-(o_nxn-2)-1
                       if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                            & .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if

              end do
           end do
        end do
        do i = o_nxn-1,o_nxn
           if(field1(i,j,k) .gt. 0) then
              if(field1(i,j,k) .eq. field2(i,j,k)) then
                 ne = ne + 1
                 row(ne) = field1(i,j,k)
                 col(ne) = field2(i,j,k)
                 val(ne) = val(ne)+x(row(ne))*x(col(ne))
              end if
           end if
        end do
     end do
  end do
end subroutine genQ_ocq_ocq

subroutine genQ_ocu_ocu(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,o_levels
     do j = 1,o_nyn-1
        do i = 1,o_nxn-2
           do kk = max(1,k-radius),min(o_levels,k+radius)
              do jj = max(1,j-radius),min(o_nyn-1,j+radius)
                 do ii = max(1,i-radius),min(o_nxn-2,i+radius)
                    if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                         & .gt. 0) then
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if
                 end do

                 if(i .lt. radius) then
                    do ii = o_nxn-2+1-radius+i,o_nxn-2
                       if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                            & .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if
                 if(i+radius .gt. o_nxn-2) then
                    do ii = 1,i+radius-(o_nxn-2)-1
                       if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                            & .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if

              end do
           end do
        end do
        do i = o_nxn-1,o_nxn
           if(field1(i,j,k) .gt. 0) then
              if(field1(i,j,k) .eq. field2(i,j,k)) then
                 ne = ne + 1
                 row(ne) = field1(i,j,k)
                 col(ne) = field2(i,j,k)
                 val(ne) = 1.0d0
              end if
           end if
        end do
     end do

     j = o_nyn
     do i = 1,o_nxn
        if(field1(i,j,k) .gt. 0) then
           if(field1(i,j,k) .eq. field2(i,j,k)) then
              ne = ne + 1
              row(ne) = field1(i,j,k)
              col(ne) = field2(i,j,k)
              val(ne) = val(ne)+x(row(ne))*x(col(ne))
           end if
        end if
     end do
  end do

end subroutine genQ_ocu_ocu


subroutine genQ_ocq_ocu(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,o_levels
     do j = 1,o_nyn
        do i = 1,o_nxn-2
           do kk = max(1,k-radius),min(o_levels,k+radius)
              do jj = max(1,j-radius),min(o_nyn,j+radius-1)
                 do ii = max(1,i-radius),min(o_nxn-2,i+radius-1)
                    if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                         & .gt. 0) then
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if
                 end do

                 if(i .lt. radius) then
                    do ii = o_nxn-2-radius+i,o_nxn-2
                       if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                            & .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if
                 if(i+radius .gt. o_nxn-2) then
                    do ii = 1,i+radius-(o_nxn-2)-1
                       if(field1(i,j,k) .gt. 0 .and. field2(ii,jj,kk)&
                            & .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if

              end do
           end do
        end do
     end do
  end do
end subroutine genQ_ocq_ocu


subroutine genQ_atq_ocq(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,radius
     do j = 1,a_nyn
        do i = 1,a_nxn
           do kk = 1,(radius+1-k)
              do jj = max(2*j-1-radius,1),min(2*j-1+radius-1,o_nyn)
                 
                 do ii = max(3*i-2-radius,1),min(3*i-2+radius,o_nxn)
                    if(field2(ii,jj,kk) .gt. 0) then
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if                       
                 end do
                 
                 if(i .eq. 1) then
                    do ii = o_nxn-2+1-radius,o_nxn-2-1
                       if(field2(ii,jj,kk) .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if

                 if(i .eq. a_nxn) then
                    do ii = 1,radius
                       if(field2(ii,jj,kk) .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if
                 
              end do
           end do
        end do
     end do
  end do
end subroutine genQ_atq_ocq


subroutine genQ_atq_ocu(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,radius
     do j = 1,a_nyn
        do i = 1,a_nxn
           do kk = 1,(radius+1-k)
              do jj = max(2*j-2-radius,1),min(2*j-2+radius,o_nyn)
                 
                 do ii = max(3*i-3-radius+1,1),min(3*i-2+radius-1,o_nxn)
                    if(field2(ii,jj,kk) .gt. 0) then
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if                       
                 end do
                 
                 if(i .eq. 1) then
                    do ii = o_nxn-2+1-radius,o_nxn-2-1
                       if(field2(ii,jj,kk) .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if

                 if(i .eq. a_nxn) then
                    do ii = 1,radius
                       if(field2(ii,jj,kk) .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if
                 
              end do
           end do
        end do
     end do
  end do
end subroutine genQ_atq_ocu

subroutine genQ_atu_ocq(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,radius
     do j = 1,a_nyn-1
        do i = 1,a_nxn
           do kk = 1,(radius+1-k)
              do jj = max(2*j-1-radius+1,1),min(2*j+radius-1,o_nyn)
                 
                 do ii = max(3*i-1-radius+1,1),min(3*i+radius-1,o_nxn)
                    if(field2(ii,jj,kk) .gt. 0) then
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if                       
                 end do
                 
                 if(i .eq. 1) then
                    do ii = o_nxn-2+1-radius,o_nxn-2-1
                       if(field2(ii,jj,kk) .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if

                 if(i .eq. a_nxn) then
                    do ii = 1,radius
                       if(field2(ii,jj,kk) .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if
                 
              end do
           end do
        end do
     end do
  end do
end subroutine genQ_atu_ocq


subroutine genQ_atu_ocu(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn,a_levels), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,k,ii,jj,kk

  do k = 1,radius
     do j = 1,a_nyn-1
        do i = 1,a_nxn
           do kk = 1,(radius+1-k)
              do jj = max(2*j-1-radius,1),min(2*j-1+radius,o_nyn)
                 
                 do ii = max(3*i-1-radius,1),min(3*i-1+radius,o_nxn)
                    if(field2(ii,jj,kk) .gt. 0) then
                       if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j,k)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if                       
                 end do
                 
                 if(i .eq. 1) then
                    do ii = o_nxn-2+1-radius,o_nxn-2-1
                       if(field2(ii,jj,kk) .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if

                 if(i .eq. a_nxn) then
                    do ii = 1,radius
                       if(field2(ii,jj,kk) .gt. 0) then
                          if(field1(i,j,k) .le. field2(ii,jj,kk)) then
                             ne = ne + 1
                             row(ne) = field1(i,j,k)
                             col(ne) = field2(ii,jj,kk)
                             val(ne) = val(ne)+x(row(ne))*x(col(ne))
                          end if
                       end if
                    end do
                 end if
                 
              end do
           end do
        end do
     end do
  end do
end subroutine genQ_atu_ocu


subroutine genQ_at2d_ocq(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,ii,jj,kk

  do j = 1,a_nyn
     do i = 1,a_nxn
        do kk = 1,radius
           do jj = max(2*j-1-radius,1),min(2*j-1+radius-1,o_nyn)
              
              do ii = max(3*i-2-radius,1),min(3*i-2+radius,o_nxn)
                 if(field2(ii,jj,kk) .gt. 0) then
                    if(field1(i,j) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end if
              end do
              
              if(i .eq. 1) then
                 do ii = o_nxn-2+1-radius,o_nxn-2-1
                    if(field2(ii,jj,kk) .gt. 0) then
                       if(field1(i,j) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if
                 end do
              end if
              
              if(i .eq. a_nxn) then
                 do ii = 1,radius
                    if(field2(ii,jj,kk) .gt. 0) then
                       if(field1(i,j) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if
                 end do
              end if
              
           end do
        end do
     end do
  end do

end subroutine genQ_at2d_ocq


subroutine genQ_at2d_ocu(ne,row,col,val,n,radius,field1,field2,x)
  use hadcm3_config
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(inout) :: ne
  integer, intent(in) :: n,radius
  integer, dimension(n), intent(inout) :: row,col
  real(kind=rk), dimension(n), intent(inout) :: val
  real(kind=rk), dimension(state_dim), intent(in) :: x
  integer, dimension(a_nxn,a_nyn), intent(in) :: field1
  integer, dimension(o_nxn,o_nyn,o_levels), intent(in) :: field2
  integer :: i,j,ii,jj,kk
  
  do j = 1,a_nyn
     do i = 1,a_nxn
        do kk = 1,radius
           do jj = max(2*j-2-radius,1),min(2*j-2+radius,o_nyn)
              
              do ii = max(3*i-3-radius+1,1),min(3*i-2+radius-1,o_nxn)
                 if(field2(ii,jj,kk) .gt. 0) then
                    if(field1(i,j) .le. field2(ii,jj,kk)) then
                       ne = ne + 1
                       row(ne) = field1(i,j)
                       col(ne) = field2(ii,jj,kk)
                       val(ne) = val(ne)+x(row(ne))*x(col(ne))
                    end if
                 end if
              end do
              
              if(i .eq. 1) then
                 do ii = o_nxn-2+1-radius,o_nxn-2-1
                    if(field2(ii,jj,kk) .gt. 0) then
                       if(field1(i,j) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if
                 end do
              end if
              
              if(i .eq. a_nxn) then
                 do ii = 1,radius
                    if(field2(ii,jj,kk) .gt. 0) then
                       if(field1(i,j) .le. field2(ii,jj,kk)) then
                          ne = ne + 1
                          row(ne) = field1(i,j)
                          col(ne) = field2(ii,jj,kk)
                          val(ne) = val(ne)+x(row(ne))*x(col(ne))
                       end if
                    end if
                 end do
              end if
              
           end do
        end do
     end do
  end do

end subroutine genQ_at2d_ocu
