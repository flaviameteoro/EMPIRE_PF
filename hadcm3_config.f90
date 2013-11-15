module hadcm3_config
  Integer, parameter :: a_nxn=96
  integer, parameter :: a_nyn=73
  integer, parameter :: a_levels=19
  integer, parameter :: a_ntimesteps=48

  Integer, parameter :: o_nxn=290
  integer, parameter :: o_nyn=144
  integer, parameter :: o_levels=20
  integer, parameter :: o_ntimesteps=24

  integer, dimension(o_nxn,o_nyn) :: q_depth,p_depth
  integer, dimension(:), allocatable :: ocn_grp

  contains
    subroutine load_bathimitry
      integer, dimension(o_nyn*o_nxn) :: bog
      integer, dimension(o_nxn,o_nyn) :: bath
      real(kind=kind(1.0)) :: aa1,bb1,cc1,dd1
      character(7) :: a1,b1,c1,d1
      integer :: i,j,line
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

      q_depth(o_nxn,:) = q_depth(2,:)
      
    end subroutine load_bathimitry

    subroutine split_oceans
      use sizes
      integer :: count,i,j
      call load_bathimitry
      allocate(ocn_grp(obs_dim))
      count = 0
      ocn_grp = 0
      do j = 1,o_nyn
         do i = 1,o_nxn-1
            if(p_depth(i,j) .ne. 0) then
               count = count + 1
               if(mod(j,5) .eq. 4 .and. mod(i,5) .eq. 1) then
                  if(j .le. 24) then
                     ocn_grp(count) = 1
                  elseif(j .ge. 125) then
                     ocn_grp(count) = 2
                  elseif(i .gt. 30 .and. i .le. 110 .and. &
                       &real(j) .gt. 51.5-0.09*real(i) .and. &
                       &real(j) .le. 170.0-1.1*real(i) .and. &
                       &real(j) .le. 40.0 + real(i) ) then
                     !                  write(*,'(A)',advance='no') '3'
                     ocn_grp(count) = 3
                  elseif(j .ge. 64 .and. j .le. 80 .and. i .le. 15) then
                     ocn_grp(count) = 4
                  elseif(j .ge. 64 .and. j .le. 80 .and. i .gt. 240) then
                     ocn_grp(count) = 4
                  elseif(j .ge. 64 .and. j .le. 80) then
                     ocn_grp(count) = 5
                  elseif(j .ge. 81 .and. j .le. 124 .and. i .ge. 215 ) then
                     ocn_grp(count) = 6
                  elseif(j .ge. 81 .and. j .le. 124 .and. i .ge. 80 .and. i&
                       & .le. 214) then
                     ocn_grp(count) = 7
                  elseif(j .le. 80 .and. real(j) .ge.  51.5&
                       &-0.09*real(i) .and. i .ge. 120 .and. i .le. 240) then
                     ocn_grp(count) = 8
                  elseif(i .gt. 30 .and. i .lt. 240 .and. real(j) .le. 51.5&
                       &-0.09*real(i)) then
                     ocn_grp(count) = 1
                  elseif(i.gt. 240 .and. real(j) .lt. 0.27*real(i) - 38.0) then
                     ocn_grp(count) = 1
                  elseif(i .lt. 30 .and. real(j) .lt. 0.33*real(i) + 40.0) then
                     ocn_grp(count) = 1
                  elseif(i .gt. 240 .and. real(j) .ge. 0.27*real(i) - 38.0)&
                       & then
                     ocn_grp(count) = 9
                  elseif(i .le. 15 .and. j .le. 63 .and. j .ge. 25) then
                     ocn_grp(count) = 9
                  else
                     ocn_grp(count) = 0
                  end if
               else
                  ocn_grp(count) = 0
               end if
            end if
         end do
         i = o_nxn
         if(p_depth(i,j) .ne. 0) then
            count = count + 1
            if(mod(j,5) .eq. 4) then
               if(j .le. 24) then
                  ocn_grp(count) = 1
               elseif(j .ge. 125) then
                  ocn_grp(count) = 2
               elseif(j .ge. 64 .and. j .le. 80 .and. i .gt. 240) then
                  ocn_grp(count) = 4
               elseif(j .ge. 64 .and. j .le. 80) then
                  ocn_grp(count) = 5
               elseif(j .ge. 81 .and. j .le. 124 .and. i .ge. 215 ) then
                  ocn_grp(count) = 6
               elseif(i.gt. 240 .and. real(j) .lt. 0.27*real(i) - 38.0) then
                  ocn_grp(count) = 1
               elseif(i .gt. 240 .and. real(j) .ge. 0.27*real(i) - 38.0)&
                    & then
                  ocn_grp(count) = 9
               else
                  ocn_grp(count) = 0
               end if
            else
               ocn_grp(count) = 0
            end if
         end if
      end do
      print*,'at end of splitting oceans, maxval = ',maxval(ocn_grp)
    end subroutine split_oceans

end module hadcm3_config

module hadcm3_data
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:,:) :: u,v,qt,thetal
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:) :: pstar
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:,:) :: t_o,sal,b_u,b_v
  real(kind=kind(1.0D0)), allocatable, dimension(:,:,:) :: mld
end module hadcm3_data

  
