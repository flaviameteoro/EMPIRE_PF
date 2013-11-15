subroutine perturb_particle(x)
use sizes
integer, parameter :: rk=kind(1.0D0)
real(kind=rk), dimension(state_dim), intent(inout) :: x
real(kind=rk), dimension(state_dim) :: rdom,y

call NormalRandomNumbers1D(0.0D0,1.0D0,state_dim,rdom)
call Qhalf(1,rdom,y)
x = x + y

end subroutine perturb_particle


subroutine update_state(state,fpsi,kgain,betan)
use sizes
integer, parameter :: rk=kind(1.0D0)
real(kind=rk), dimension(state_dim), intent(out) :: state
real(kind=rk), dimension(state_dim), intent(in) :: fpsi,kgain
real(kind=rk), dimension(state_dim), intent(inout) :: betan
real(kind=rk), dimension(state_dim) :: diff
integer :: k,start
real(kind=rk) :: dnrm2
!logical, parameter :: debug=.false.
logical, parameter :: norms=.false.


!!atmosphere humidity = state_vector(406465:539616)
!do k = 406465,539616
!   ! fpsi(i)+kgain(i)+betan(i) has to lie between 0 and 1   
!   betan(k) = max(-fpsi(k)-kgain(k),min(betan(k),1.0D0-fpsi(k)-kgain(k)))
!end do


!!ocean salinity = state_vector(997128:1454638)
!do k = 997128,1454638
!   ! fpsi(i)+kgain(i)+betan(i) has to lie between 0 and 1
!   betan(k) = max(-fpsi(k)-kgain(k),min(betan(k),1.0D0-fpsi(k)-kgain(k)))
!end do


!do the addition
state = fpsi+kgain+betan
!if(debug) print*,'|fpsi|=',dnrm2(state_dim,fpsi,1),' |kgain|= ',dnrm2(state_dim,kgain&
!     &,1),' |betan| = ',dnrm2(state_dim,betan,1)
if(norms) print*,' |kgain|= ',dnrm2(state_dim,kgain,1),' |betan| = ',dnrm2(state_dim,betan,1)


!now fix the atmospheric polar values:
! PSTAR     U     V     THETA    Q

!index = 0
!do k=1,levels
!  do j=1,columns
!    do i=1,rows
!       index = index + 1
!       data_out(index)=data(i,j,k) == data(longitude,latitude,level)
!    enddo
!  enddo
!enddo

!at the poles, all the longitude values are equal
!first lets do Pstar
!state(1:a_nxn) = state(1)
!state((a_nyn-1)*a_nxn+1:a_nyn*a_nxn) = state((a_nyn-1)*a_nxn+1)

!start = a_nxn*a_nyn   !for U
!do k = 1,a_levels
   !j = 1
!   state( start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn )  =&
!        & state( start+(k-1)*a_nyn*a_nxn + 1 )
!if(debug)   print*,'fpsi test U k:',k,&
!        maxval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : &
!        start+(k-1)*a_nyn*a_nxn + a_nxn))-&
!        minval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : &
!        start+(k-1)*a_nyn*a_nxn + a_nxn))
   !j = a_nyn
!   state( start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn) = state(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)&
!        &*a_nxn + 1)
!if(debug)   print*,'fpsi test U k:',k,maxval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))-minval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))
!end do
!start = a_nxn*a_nyn + 1*a_levels*a_nyn*a_nxn !for V
!do k = 1,a_levels
   !j = 1
   !        
!   state( start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn )  =&
!        & state( start+(k-1)*a_nyn*a_nxn + 1 )
!if(debug)   print*,'fpsi test V k:',k,&
!        maxval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & a_nxn))&
!        -minval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn))
   !j = a_nyn
   !        
!   state( start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn) = state(start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)&
!        &*a_nxn + 1)
!if(debug)   print*,'fpsi test V k:',k,maxval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))-minval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))
!end do
!start = a_nxn*a_nyn + 2*a_levels*a_nyn*a_nxn !for Theta
!do k = 1,a_levels
   !j = 1
   !        
!   state( start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn )  =&
!        & state( start+(k-1)*a_nyn*a_nxn + 1 )
!if(debug)   print*,'fpsi test Theta k:',k,&
!        maxval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & a_nxn))&
!        -minval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn))
   !j = a_nyn
   !        
!   state( start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn) = state(start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)&
!        &*a_nxn + 1)
!if(debug)   print*,'fpsi test Theta k:',k,maxval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))-minval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))
!end do
!start = a_nxn*a_nyn + 3*a_levels*a_nyn*a_nxn !for Q
!do k = 1,a_levels
   !j = 1
   !        
!   state( start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn )  =&
!        & state( start+(k-1)*a_nyn*a_nxn + 1 )
!if(debug)   print*,'fpsi test Q:',k,':',maxval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)&
!        &*a_nyn*a_nxn + a_nxn))-minval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k&
!        &-1)*a_nyn*a_nxn + a_nxn ))
   !j = a_nyn
   !        
!   state( start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn) = state(start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)&
!        &*a_nxn + 1)
!if(debug)   print*,'fpsi test Q:',k,':',&
!        maxval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)&
!        &*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + a_nxn))&
!        -minval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)&
!        &*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + a_nxn))
!end do

!if(debug)print*,'PSTAR: min fpsi= ',minval(fpsi(1:7008))
!if(debug)print*,'PSTAR: max fpsi= ',maxval(fpsi(1:7008))
!if(debug)print*,'PSTAR: min state= ',minval(state(1:7008))
!if(debug)print*,'PSTAR: max state= ',maxval(state(1:7008))




!if(debug)diff = abs(state-fpsi)
!if(debug)print*,'max perturbation is ',maxval(diff),' at ',maxloc(diff),' of '&
!     &,state(maxloc(diff))
!if(debug)diff = diff/(max(abs(fpsi),1.0d-8))
!if(debug)print*,'max prop perturb is ',maxval(diff),' at ',maxloc(diff),' of '&
!     &,fpsi(maxloc(diff))




end subroutine update_state
