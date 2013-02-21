Module FullQ

Contains
!!$   subroutine solve_r(y,v)
!!$    !subroutine to take an observation vector y and return v
!!$    !in observation space.
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(obs_dim), intent(in) :: y
!!$    real(kind=rk), dimension(obs_dim), intent(out) :: v
!!$    
!!$    v = y
!!$  end subroutine solve_r
!!$
!!$   subroutine solve_q(x,v)
!!$    !subroutine to take a full state vector x and return v
!!$    !in state space.
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(state_dim), intent(in) :: x
!!$    real(kind=rk), dimension(state_dim), intent(out) :: v
!!$    
!!$    v = x
!!$  end subroutine solve_q
!!$
!!$  subroutine solve_hqht_plus_r(y,v)
!!$    !subroutine to take an observation vector y and return v
!!$    !in observation space.
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(obs_dim), intent(in) :: y
!!$    real(kind=rk), dimension(obs_dim), intent(out) :: v
!!$    
!!$    v = y
!!$  end subroutine solve_hqht_plus_r
!!$  
!!$  subroutine Qhalf(x,Qx)
!!$    !subroutine to take a full state vector x and return Q^(1/2)x
!!$    !in state space.
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(state_dim), intent(in) :: x
!!$    real(kind=rk), dimension(state_dim), intent(out) :: qx
!!$    
!!$    qx = sqrt(x)
!!$  end subroutine QHALF
!!$
!!$  subroutine Q(x,Qx)
!!$    !subroutine to take a full state vector x and return Qx
!!$    !in state space.
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(state_dim), intent(in) :: x
!!$    real(kind=rk), dimension(state_dim), intent(out) :: qx
!!$    
!!$    qx = x
!!$  end subroutine Q

!!$  subroutine innerHQHt_plus_R_1(y,w)
!!$    !subroutine to take an observation vector y and return w = y^T (HQH^T+R)^(-1) y
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(obs_dim), intent(in) :: y
!!$    real(kind=rk), dimension(obs_dim) :: v
!!$    real(kind=rk), intent(out) :: w
!!$
!!$    call solve_hqht_plus_r(y,v)
!!$    
!!$    !this can defo be done better using BLAS PAB...
!!$    w = sum(y*v)
!!$
!!$  end subroutine innerHQHt_plus_R_1
!!$
!!$  subroutine innerR_1(y,w)
!!$    !subroutine to take an observation vector y and return w = y^T R^(-1) y
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(obs_dim), intent(in) :: y
!!$    real(kind=rk), dimension(obs_dim) :: v
!!$    real(kind=rk), intent(out) :: w
!!$
!!$    call solve_r(y,v)
!!$    
!!$    !this can defo be done better using BLAS PAB...
!!$    w = sum(y*v)
!!$
!!$  end subroutine innerR_1
!!$
!!$  subroutine innerQ_1(x,w)
!!$    !subroutine to take a full state vector x and return w = x^T Q^(-1) x
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(state_dim), intent(in) :: x
!!$    real(kind=rk), dimension(state_dim) :: v
!!$    real(kind=rk), intent(out) :: w
!!$
!!$    call solve_q(x,v)
!!$    
!!$    !this can defo be done better using BLAS PAB...
!!$    w = sum(x*v)
!!$
!!$  end subroutine innerQ_1



!!$
!!$  subroutine K(y,x)
!!$    !subroutine to apply the operator K to a vector y in obs space and return
!!$    !the vector x in full state space.
!!$    use pf_control
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(state_dim), intent(out) :: x
!!$    real(kind=rk), dimension(obs_dim), intent(in) :: y
!!$    
!!$    real(kind=rk), dimension(obs_dim) :: v
!!$    real(kind=rk), dimension(state_dim) :: vv
!!$    real(kind=rk) :: temp
!!$    
!!$    call solve_hqht_plus_r(y,v)
!!$
!!$    call HT(v,vv)
!!$
!!$    call Q(vv,x)
!!$
!!$    !now apply the scaling of the K operator
!!$
!!$    temp = pf%nudgeFac*real(modulo(pf%timestep,pf%time_bwn_obs),rk)/real(pf%time_bwn_obs,rk)
!!$
!!$    x = temp*x
!!$
!!$  end subroutine K


!!$  subroutine H(x,hx)
!!$    !subroutine to take a full state vector x and return H(x)
!!$    !in observation space.
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(state_dim), intent(in) :: x
!!$    real(kind=rk), dimension(obs_dim), intent(out) :: hx
!!$    
!!$    hx = x(1:obs_dim)
!!$  end subroutine H
!!$
!!$  subroutine HT(y,x)
!!$    !subroutine to take an observation vector y and return x = H^T(y)
!!$    !in full state space.
!!$    use sizes
!!$    implicit none
!!$    integer, parameter :: rk=kind(1.0D+0)
!!$    real(kind=rk), dimension(state_dim), intent(out) :: x
!!$    real(kind=rk), dimension(obs_dim), intent(in) :: y
!!$    
!!$    x=0.0_rk
!!$    x(1:obs_dim) = y
!!$  end subroutine HT


!!$  Subroutine SolveQ (inVector,resultQ)
!!$
!!$    use Sizes
!!$    use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT) :: resultQ
!!$
!!$    integer :: status, i
!!$    real(kind=kind(1.0D0)), dimension(:), allocatable :: xVector, psi, psiSin
!!$
!!$    write(6,*) 'Subroutine SolveQ not yet implemented - returning 0 anyway'
!!$    call flush(6)
!!$    resultQ = 0.0

!!$    allocate(xVector(3*nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of storage failed"
!!$    allocate(psi(nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of storage2 failed"
!!$    allocate(psiSin(nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of help failed"         
!!$
!!$    xVector=inVector
!!$
!!$    call uInverse(xVector,psi)
!!$    call uHatInverse(psi,psiSin)
!!$
!!$    resultQ=0
!!$    do i=1,nxn*nyn
!!$       if (qv(i) .gt. 0) then
!!$          resultQ=resultQ+psiSin(i)*psiSin(i)*1./qv(i)
!!$       endif
!!$    enddo
!!$
!!$    resultQ=1./(qo*qo*qrel*qrel)*resultQ

!!$  End Subroutine SolveQ

!!$  Subroutine SolveQR (inVector,resVector)
!!$
!!$    use Sizes
!!$    use ErrorsAndVariances
!!$    use ObservationField
!!$    use minresModule
!!$    use minresDataModule
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    integer :: n, itnlim, nout, istop, itn, status, i
!!$    logical :: checkA, precon       
!!$    real(dp), allocatable, dimension(:) :: b, x
!!$    real(dp) :: shift, rtol, Anorm, Acond, rnorm, Arnorm, ynorm
!!$
!!$    write(6,*) 'Subroutine multiplySqrtQ not yet implemented - returning 0 anyway'
!!$    call flush(6)
!!$    resVector = 0.0

!!$    n=nObs
!!$    allocate(b(n),STAT=status)
!!$    if (status /= 0) stop " Allocation of b failed"
!!$    allocate(x(n),STAT=status)
!!$    if (status /= 0) stop " Allocation of x failed"
!!$
!!$    checkA=.true.
!!$    precon=.false.
!!$    shift=-1.
!!$    nout=0
!!$    itnlim=500 !Mel-change this?????
!!$    rtol=1.e-16 !Mel-change this?????
!!$
!!$    do i=1, nObs
!!$       b(i) = 1./qdData(i)*inVector(i)
!!$    enddo
!!$
!!$    call MINRES(n,Aprod,Aprod,b,shift,checkA,precon,x,itnlim,nout,rtol,istop,itn,Anorm,Acond,rnorm,Arnorm,ynorm)
!!$
!!$    do i=1,nObs
!!$       resVector(i) = 1./qdData(i)*x(i)
!!$    enddo
!!$
!!$    write(*,*) 'istop: ', istop
!!$    write(*,*) 'itn: ', itn

!!$  End Subroutine SolveQR


!!$  Subroutine Aprod(n,x,y)
!!$
!!$    use minresDataModule
!!$    use ErrorsAndVariances
!!$    use ObservationField
!!$    use Sizes
!!$
!!$    IMPLICIT NONE
!!$
!!$    integer, intent(in) :: n
!!$    real(dp), intent(in) :: x(n)
!!$    real(dp), intent(out) :: y(n)
!!$
!!$    real(kind=kind(1.0D0)), dimension(n) :: inVector, resVector
!!$    real(kind=kind(1.0D0)), dimension(3*nxn*nyn) :: xVector
!!$    real(kind=kind(1.0D0)), dimension(nxn*nyn) :: psi, psiSin
!!$    integer :: i
!!$
!!$    inVector=x
!!$
!!$    do i=1,nObs
!!$       inVector(i)=1./qdData(i)*inVector(i)
!!$    enddo
!!$
!!$    inVector=(qo*qo*qrel*qrel)*inVector
!!$
!!$    call expand(inVector,xVector)
!!$    call uAdjoint(xVector,psi)
!!$    call uHatAdjoint(psi,psiSin)
!!$
!!$    do i=1,nxn*nyn
!!$       psiSin(i)=psiSin(i)*qv(i)
!!$    enddo
!!$
!!$    call uHatTransform(psiSin,psi)
!!$    call uTransform(psi,xVector)
!!$    call measureVec(xVector,resVector)
!!$
!!$    do i=1,nObs
!!$       resVector(i)=1./qdData(i)*resVector(i)
!!$    enddo
!!$
!!$    y=resVector
!!$
!!$  End Subroutine Aprod

!!$  Subroutine CheckQRinverse (inVector, resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$    Use ObservationField
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    real(kind=kind(1.0D0)), dimension(nObs) :: obsVector
!!$    real(kind=kind(1.0D0)), dimension(3*nxn*nyn) :: xVector, QxVector
!!$    integer :: i
!!$
!!$    obsVector=inVector
!!$    call expand(obsVector,xVector)
!!$    call MultiplyQ(xVector,QxVector)
!!$    call measureVec(QxVector,obsVector)
!!$
!!$    do i=1,nObs
!!$       resVector(i)=obsVector(i)+qdData(i)*qdData(i)*inVector(i)
!!$    enddo
!!$
!!$  End Subroutine CheckQRInverse


!!$  Subroutine MultiplyQ (inVector, resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$
!!$    real(kind=kind(1.0D0)), dimension(:), allocatable :: xVector, psi, psiSin
!!$    integer :: status, i
!!$
!!$    write(6,*) 'Subroutine multiplyQ not yet implemented - returning 0 anyway'
!!$    call flush(6)
!!$    resVector = 0.0
!!$    allocate(xVector(3*nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of storage failed"
!!$    allocate(psi(nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of storage2 failed"
!!$    allocate(psiSin(nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of help failed"  
!!$
!!$    xVector=inVector
!!$
!!$    call uAdjoint(xVector,psi)
!!$    call uHatAdjoint(psi,psiSin)
!!$
!!$    do i=1,nxn*nyn
!!$       psiSin(i)=qv(i)*psiSin(i)
!!$    enddo
!!$
!!$    call uHatTransform(psiSin,psi)
!!$    call uTransform(psi,xVector)
!!$
!!$    resVector=qo*qo*qrel*qrel*xVector

!!$  End Subroutine MultiplyQ


!!$  Subroutine MultiplySqrtQ (inVector, resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$
!!$    real(kind=kind(1.0D0)), dimension(:), allocatable :: psiSin, psi, xVector
!!$
!!$    integer :: status, i
!!$    write(6,*) 'Subroutine multiplySqrtQ not yet implemented - returning 0 anyway'
!!$    call flush(6)
!!$    resVector = 0.0
!!$    allocate(xVector(3*nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of storage failed"
!!$    allocate(psi(nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of storage2 failed"
!!$    allocate(psiSin(nxn*nyn), STAT=status)
!!$    if (status /= 0) stop " Allocation of help failed"           
!!$
!!$    do i=1,nxn*nyn
!!$       psiSin(i)=sqrtQv(i)*inVector(i)
!!$    enddo
!!$
!!$    call uHatTransform(psiSin,psi)
!!$    call uTransform(psi,xVector)
!!$
!!$    resVector=qrel*qo*xVector

!!$  End Subroutine MultiplySqrtQ


!!$  Subroutine uTransform(inVector,resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    real(kind=kind(1.0D0)), dimension(nxn,nyn)::psi,u,v,e
!!$    real(kind=kind(1.0D0)) :: helpdy, helpdx, help
!!$    integer :: i,j
!!$
!!$    u(:,:)=0.
!!$    v(:,:)=0.
!!$    e(:,:)=0.
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          psi(i,j)=inVector((i-1)*nyn+j)
!!$       enddo
!!$    enddo
!!$
!!$    helpdy=1./dy
!!$    do i=2,nxn
!!$       do j=2,nyn-1
!!$          !u(i,j)=-(psi(i,j+1)-psi(i,j))*1./dy
!!$          u(i,j)=u(i,j)-helpdy*psi(i,j+1)
!!$          u(i,j)=u(i,j)+helpdy*psi(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    helpdx=1./dx
!!$    do i=2,nxn-1
!!$       do j=2,nyn
!!$          !v(i,j)=(psi(i+1,j)-psi(i,j))*1./dx
!!$          v(i,j)=v(i,j)+helpdx*psi(i+1,j)
!!$          v(i,j)=v(i,j)-helpdx*psi(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    do i=2,nxn-1
!!$       do j=2,nyn-1
!!$          help=f(j)*1./(4*gac)
!!$          !e(i,j)=1./(4*gac)*f(j)*(psi(i,j+1)+psi(i+1,j+1)+psi(i+1,j)+psi(i,j))
!!$          e(i,j)=e(i,j)+help*psi(i,j+1)
!!$          e(i,j)=e(i,j)+help*psi(i+1,j+1)
!!$          e(i,j)=e(i,j)+help*psi(i+1,j)
!!$          e(i,j)=e(i,j)+help*psi(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          resVector((i-1)*nyn+j)=u(i,j)
!!$          resVector(nxn*nyn+(i-1)*nyn+j)=v(i,j)
!!$          resVector(2*nxn*nyn+(i-1)*nyn+j)=e(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    call uBoundary(resVector,resVector)
!!$
!!$  End Subroutine uTransform


!!$  Subroutine uBoundary(inVector,resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    real(kind=kind(1.0D0)), dimension(nxn,nyn)::u,v,e,uF,vF,eF
!!$    integer :: i,j
!!$
!!$    u(:,:)=0.
!!$    v(:,:)=0.
!!$    e(:,:)=0.
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          uF(i,j)=inVector((i-1)*nyn+j)
!!$          vF(i,j)=inVector(nxn*nyn+(i-1)*nyn+j)
!!$          eF(i,j)=inVector(2*nxn*nyn+(i-1)*nyn+j)
!!$       enddo
!!$    enddo
!!$
!!$    do i=2,nxn-1
!!$       do j=2,nyn-1
!!$          u(i,j)=u(i,j)+uF(i,j)
!!$          v(i,j)=v(i,j)+vF(i,j)
!!$          e(i,j)=e(i,j)+eF(i,j)
!!$       enddo
!!$    enddo
!!$    do i=2,nxn
!!$       v(i,1)=v(i,1)+vF(i,1)
!!$    enddo
!!$    do j=2,nyn
!!$       u(1,j)=u(1,j)+uF(1,j)
!!$    enddo
!!$
!!$    !free slip boundaries at north and south sides
!!$    do i=1,nxn
!!$       u(i,1)=u(i,1)-uF(i,2)
!!$       u(i,nyn)=u(i,nyn)-uF(i,nyn-1)
!!$       v(i,2)=0
!!$       v(i,nyn)=0
!!$       !e(i,1)=e(i,1)+eF(i,2)
!!$       e(i,1)=e(i,1)-eF(i,2)
!!$       !e(i,nyn)=e(i,nyn)+eF(i,nyn-1)
!!$       e(i,nyn)=e(i,nyn)-eF(i,nyn-1)
!!$    enddo
!!$
!!$    ! free slip boundaries at west and east sides (u=v_y=0)
!!$    do j=1,nyn
!!$       u(2,j)=0.
!!$       u(nxn,j)=0.
!!$       v(1,j)=v(1,j)-vF(2,j)
!!$       v(nxn,j)=v(nxn,j)-vF(nxn-1,j)
!!$       !e(1,j)=e(1,j)+eF(2,j)
!!$       e(1,j)=e(1,j)-eF(2,j)
!!$       !e(nxn,j)=e(nxn,j)+eF(nxn-1,j)
!!$       e(nxn,j)=e(nxn,j)-eF(nxn-1,j)
!!$    enddo
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          resVector((i-1)*nyn+j)=u(i,j)
!!$          resVector(nxn*nyn+(i-1)*nyn+j)=v(i,j)
!!$          resVector(2*nxn*nyn+(i-1)*nyn+j)=e(i,j)
!!$       enddo
!!$    enddo
!!$
!!$  End Subroutine uBoundary

!!$  Subroutine uTBoundary(inVector,resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    real(kind=kind(1.0D0)), dimension(nxn,nyn)::u,v,e,uhat,vhat,ehat
!!$    integer :: i,j
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          u(i,j)=inVector((i-1)*nyn+j)
!!$          v(i,j)=inVector(nxn*nyn+(i-1)*nyn+j)
!!$          e(i,j)=inVector(2*nxn*nyn+(i-1)*nyn+j)
!!$       enddo
!!$    enddo
!!$
!!$    uhat(:,:)=0.
!!$    vhat(:,:)=0.
!!$    ehat(:,:)=0.
!!$
!!$    ! free slip boundaries at west and east sides (u=v_y=0)
!!$    do j=1,nyn
!!$       vhat(2,j)=vhat(2,j)-v(1,j)
!!$       vhat(nxn-1,j)=vhat(nxn-1,j)-v(nxn,j)
!!$       !ehat(2,j)=ehat(2,j)+e(1,j)
!!$       ehat(2,j)=ehat(2,j)-e(1,j)
!!$       !ehat(nxn-1,j)=ehat(nxn-1,j)+e(nxn,j)
!!$       ehat(nxn-1,j)=ehat(nxn-1,j)-e(nxn,j)
!!$    enddo
!!$
!!$    !free slip boundaries at north and south sides
!!$    do i=1,nxn
!!$       uhat(i,2)=uhat(i,2)-u(i,1)
!!$       uhat(i,nyn-1)=uhat(i,nyn-1)-u(i,nyn)
!!$       !ehat(i,2)=ehat(i,2)+e(i,1)
!!$       ehat(i,2)=ehat(i,2)-e(i,1)
!!$       !ehat(i,nyn-1)=ehat(i,nyn-1)+e(i,nyn)
!!$       ehat(i,nyn-1)=ehat(i,nyn-1)-e(i,nyn)
!!$    enddo
!!$
!!$    do j=2,nyn
!!$       uhat(1,j)=uhat(1,j)+u(1,j)
!!$    enddo
!!$    do i=2,nxn
!!$       vhat(i,1)=vhat(i,1)+v(i,1)
!!$    enddo
!!$    do i=2,nxn-1
!!$       do j=2,nyn-1
!!$          uhat(i,j)=uhat(i,j)+u(i,j)
!!$          vhat(i,j)=vhat(i,j)+v(i,j)
!!$          ehat(i,j)=ehat(i,j)+e(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          resVector((i-1)*nyn+j)=uhat(i,j)
!!$          resVector(nxn*nyn+(i-1)*nyn+j)=vhat(i,j)
!!$          resVector(2*nxn*nyn+(i-1)*nyn+j)=ehat(i,j)
!!$       enddo
!!$    enddo
!!$
!!$  End Subroutine uTBoundary


!!$  Subroutine uInverse(inVector, resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    real(kind=kind(1.0D0)), dimension(nxn,nyn):: psi,u,v,e,psiTempR,psiTempI
!!$    real(kind=kind(1.0D0)), dimension(2*nxn) :: TRIGM
!!$    real(kind=kind(1.0D0)), dimension(2*nyn) :: TRIGN
!!$    real(kind=kind(1.0D0)), dimension(2*nxn*nyn) :: WORK
!!$    real(kind=kind(1.0D0)) :: help
!!$    integer :: ifail=0.
!!$    integer :: i,j,n2
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          u(i,j)=inVector((i-1)*nyn+j)
!!$          v(i,j)=inVector(nxn*nyn+(i-1)*nyn+j)
!!$          e(i,j)=inVector(2*nxn*nyn+(i-1)*nyn+j)
!!$       enddo
!!$    enddo
!!$
!!$    do i=2,nxn
!!$       do j=2,nyn
!!$          psi(i,j)=1./(4*f(j))*gac*(e(i-1,j)+e(i,j)+e(i,j-1)+e(i-1,j-1))
!!$       enddo
!!$    enddo
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          resVector((i-1)*nyn+j)=psi(i,j)
!!$       enddo
!!$    enddo
!!$
!!$  End Subroutine uInverse


!!$  Subroutine uAdjoint(inVector, resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    real(kind=kind(1.0D0)), dimension(nxn,nyn):: psiHat,u,v,e
!!$    real(kind=kind(1.0D0)), dimension(3*nxn*nyn) :: tempVector
!!$    real :: dyHelp, dxHelp, help
!!$    integer :: i,j
!!$! PAB : this subroutine needs replacing!
!!$resVector=0.0
!!$    tempVector=inVector
!!$
!!$    call uTBoundary(tempVector,tempVector)
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          u(i,j)=tempVector((i-1)*nyn+j)
!!$          v(i,j)=tempVector(nxn*nyn+(i-1)*nyn+j)
!!$          e(i,j)=tempVector(2*nxn*nyn+(i-1)*nyn+j)
!!$       enddo
!!$    enddo
!!$
!!$    psiHat(:,:)=0.
!!$
!!$    dyHelp=1./dy
!!$    do i=2,nxn
!!$       do j=2,nyn-1
!!$          psiHat(i,j+1)=psiHat(i,j+1)-dyHelp*u(i,j)
!!$          psiHat(i,j)=psiHat(i,j)+dyHelp*u(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    dxHelp=1./dx
!!$    do i=2,nxn-1
!!$       do j=2,nyn
!!$          psiHat(i+1,j)=psiHat(i+1,j)+dxHelp*v(i,j)
!!$          psiHat(i,j)=psiHat(i,j)-dxHelp*v(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    do i=2,nxn-1
!!$       do j=2,nyn-1
!!$          help=f(j)*1./(4*gac)
!!$          psiHat(i,j+1)=psiHat(i,j+1)+help*e(i,j)
!!$          psiHat(i+1,j+1)=psiHat(i+1,j+1)+help*e(i,j)
!!$          psiHat(i+1,j)=psiHat(i+1,j)+help*e(i,j)
!!$          psiHat(i,j)=psiHat(i,j)+help*e(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          resVector((i-1)*nyn+j)=psiHat(i,j)
!!$       enddo
!!$    enddo
!!$
!!$  End Subroutine uAdjoint
!!$
!!$
!!$  Subroutine uHatAdjoint(inVector,resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    integer :: ifail, nxf, nyf, i, j, k
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:) :: psiC, psiR, work
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: psiSwap, psiFinal
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:) :: trigC, trigR
!!$    character*1 :: init
!!$
!!$! PAB : this subroutine needs replacing!

!!$    nxf=nxn-2
!!$    nyf=nyn-2
!!$
!!$    allocate(psiC(nxf*(nyf-1)), psiR((nxf-1)*nyf), psiSwap(nxf-1,nyf-1), psiFinal(nxn,nyn), work(nxf*nyf))
!!$    allocate(trigC(2*nxf),trigR(2*nyf))
!!$
!!$    !Take in input minus first and second row
!!$    !(first not needed, 2nd zero from boundary conditions and not required by nag routine)
!!$    k=1
!!$    do i=3,nxn
!!$       do j=3,nyn-1
!!$          psiC(k)=inVector((i-1)*nyn+j)
!!$          k=k+1
!!$       enddo
!!$    enddo
!!$
!!$    !sin transform over columns
!!$    init='I'
!!$    ifail=0
!!$    call C06HAF(nyf-1,nxf,psiC,init,trigC,work,ifail)
!!$
!!$    psiC=psiC*1./sqrt(2./nyf)
!!$
!!$    !Put into matrix to swap rows and columns
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiSwap(i,j)=psiC((i-1)*(nyf-1)+j)
!!$       enddo
!!$    enddo
!!$
!!$    psiR(:)=0.
!!$
!!$    !psiR now contains correcting ordering to sin transform over rows
!!$    k=1
!!$    do j=1,nyf-1
!!$       do i=1,nxf-1
!!$          psiR(k)=psiSwap(i,j)
!!$          k=k+1
!!$       enddo
!!$    enddo
!!$
!!$    !sin transform over rows
!!$    init='I'
!!$    ifail=0
!!$    call C06HAF(nxf-1,nyf,psiR,init,trigR,work,ifail)
!!$
!!$    psiR=psiR*1./sqrt(2./nxf)
!!$
!!$    !swap back rows and columns to keep ordering consistent
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiSwap(i,j)=psiR((j-1)*(nxf-1)+i)
!!$       enddo
!!$    enddo
!!$
!!$    psiFinal(:,:)=0.
!!$
!!$    !add back in zero rows and columns
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiFinal(i+2,j+2)=psiSwap(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    !put into vector form to return
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          resVector((i-1)*nyn+j)=psiFinal(i,j)
!!$       enddo
!!$    enddo
!!$
!!$  End Subroutine uHatAdjoint
!!$
!!$
!!$  Subroutine uHatTransform(inVector,resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    integer :: ifail, nxf, nyf, i, j, k
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:) :: psiC, psiR, work
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: psiSwap, psiFinal
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:) :: trigC, trigR
!!$    character*1 :: init
!!$
!!$! PAB: this subroutine needs replacing
!!$
!!$    nxf=nxn-2
!!$    nyf=nyn-2
!!$
!!$    allocate(psiC(nxf*(nyf-1)), psiR((nxf-1)*nyf), psiSwap(nxf-1,nyf-1), psiFinal(nxn,nyn), work(nxf*nyf))
!!$    allocate(trigC(2*nxf),trigR(2*nyf))
!!$
!!$    !Take in input minus first and second row
!!$    !(first not needed, 2nd zero from boundary conditions and not required by nag routine)
!!$    k=1
!!$    do i=3,nxn
!!$       do j=3,nyn-1
!!$          psiC(k)=inVector((i-1)*nyn+j)
!!$          k=k+1
!!$       enddo
!!$    enddo
!!$
!!$    psiC=psiC*1./sqrt(2./nyf)
!!$
!!$    !sin transform over columns
!!$    init='I'
!!$    ifail=0
!!$    call C06HAF(nyf-1,nxf,psiC,init,trigC,work,ifail)
!!$
!!$    !Put into matrix to swap rows and columns
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiSwap(i,j)=psiC((i-1)*(nyf-1)+j)
!!$       enddo
!!$    enddo
!!$
!!$    psiR(:)=0.
!!$
!!$    !psiR now contains correcting ordering to sin transform over rows
!!$    k=1
!!$    do j=1,nyf-1
!!$       do i=1,nxf-1
!!$          psiR(k)=psiSwap(i,j)
!!$          k=k+1
!!$       enddo
!!$    enddo
!!$
!!$    psiR=psiR*1./sqrt(2./nxf)
!!$
!!$    !sin transform over rows
!!$    init='I'
!!$    ifail=0
!!$    call C06HAF(nxf-1,nyf,psiR,init,trigR,work,ifail)
!!$
!!$    !swap back rows and columns to keep ordering consistent
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiSwap(i,j)=psiR((j-1)*(nxf-1)+i)
!!$       enddo
!!$    enddo
!!$
!!$    psiFinal(:,:)=0.
!!$
!!$    !add back in zero rows and columns
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiFinal(i+2,j+2)=psiSwap(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    !put into vector form to return
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          resVector((i-1)*nyn+j)=psiFinal(i,j)
!!$       enddo
!!$    enddo
!!$
!!$  End Subroutine uHatTransform
!!$
!!$
!!$  Subroutine uHatInverse(inVector, resVector)
!!$
!!$    Use Sizes
!!$    Use ErrorsAndVariances
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: inVector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    integer :: ifail, nxf, nyf, i, j, k
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:) :: psiC, psiR, work
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:,:) :: psiSwap, psiFinal
!!$    real(kind=kind(1.0D0)), allocatable, dimension(:) :: trigC, trigR
!!$    character*1 :: init
!!$
!!$
!!$!  PAB: this subroutine needs replacing!
!!$    nxf=nxn-2
!!$    nyf=nyn-2
!!$
!!$    allocate(psiC(nxf*(nyf-1)), psiR((nxf-1)*nyf), psiSwap(nxf-1,nyf-1), psiFinal(nxn,nyn), work(nxf*nyf))
!!$    allocate(trigC(2*nxf),trigR(2*nyf))
!!$
!!$    !Take in input minus first and second row
!!$    !(first not needed, 2nd zero from boundary conditions and not required by nag routine)
!!$    k=1
!!$    do i=3,nxn
!!$       do j=3,nyn-1
!!$          psiC(k)=inVector((i-1)*nyn+j)
!!$          k=k+1
!!$       enddo
!!$    enddo
!!$
!!$    !sin transform over columns
!!$    init='I'
!!$    ifail=0
!!$    call C06HAF(nyf-1,nxf,psiC,init,trigC,work,ifail)
!!$
!!$    psiC=psiC*sqrt(2./nxf)
!!$
!!$    !Put into matrix to swap rows and columns
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiSwap(i,j)=psiC((i-1)*(nyf-1)+j)
!!$       enddo
!!$    enddo
!!$
!!$    psiR(:)=0.
!!$
!!$    !psiR now contains correcting ordering to sin transform over rows
!!$    k=1
!!$    do j=1,nyf-1
!!$       do i=1,nxf-1
!!$          psiR(k)=psiSwap(i,j)
!!$          k=k+1
!!$       enddo
!!$    enddo
!!$
!!$    !sin transform over rows
!!$    init='I'
!!$    ifail=0
!!$    call C06HAF(nxf-1,nyf,psiR,init,trigR,work,ifail)
!!$
!!$    psiR=psiR*sqrt(2./nyf)
!!$
!!$    !swap back rows and columns to keep ordering consistent
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiSwap(i,j)=psiR((j-1)*(nxf-1)+i)
!!$       enddo
!!$    enddo
!!$
!!$    psiFinal(:,:)=0.
!!$
!!$    !add back in zero rows and columns
!!$    do i=1,nxf-1
!!$       do j=1,nyf-1
!!$          psiFinal(i+2,j+2)=psiSwap(i,j)
!!$       enddo
!!$    enddo
!!$
!!$    !put into vector form to return
!!$    do i=1,nxn
!!$       do j=1,nyn
!!$          resVector((i-1)*nyn+j)=psiFinal(i,j)
!!$       enddo
!!$    enddo
!!$
!!$  End Subroutine uHatInverse
!!$
!!$
!!$  Subroutine Measure (psi, lpsi) 
!!$    !I think this computes lpsi=H(psi)
!!$    !PAB
!!$
!!$    use Sizes
!!$    use ObservationField
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN),  dimension(:,:)   :: psi
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:)     :: lpsi
!!$
!!$    integer                             :: i,j,k
!!$
!!$    k=1
!!$    do j=1,nyy
!!$       do i=1,nxx
!!$          if (pseudoH(i+(j-1)*(nxx)) .eq. 1) then
!!$             lpsi(k)=psi(i,j)
!!$             k=k+1
!!$          endif
!!$       enddo
!!$    enddo
!!$
!!$  End Subroutine Measure


!!$  Subroutine MeasureVec (psi, lpsi) 
!!$
!!$    use Sizes
!!$    use ObservationField
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN),  dimension(:)   :: psi
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:)     :: lpsi
!!$
!!$    integer                             :: i,k
!!$
!!$    k=1
!!$    do i=1,3*nxn*nyn
!!$       if (pseudoH(i) .eq. 1) then
!!$          lpsi(k)=psi(i)
!!$          k=k+1
!!$       endif
!!$    enddo
!!$
!!$
!!$  End Subroutine MeasureVec
!!$
!!$
!!$  Subroutine Expand(vector,resVector)
!!$
!!$    Use Sizes
!!$    Use ObservationField
!!$
!!$    IMPLICIT NONE
!!$
!!$    real(kind=kind(1.0D0)), INTENT(IN), dimension(:) :: vector
!!$    real(kind=kind(1.0D0)), INTENT(OUT), dimension(:) :: resVector
!!$
!!$    integer :: i,k
!!$
!!$    resVector(:)=0.0
!!$
!!$    k=1
!!$    do i=1,3*nxn*nyn
!!$       if (pseudoH(i) .eq. 1) then
!!$          resVector(i)=vector(k)
!!$          k=k+1
!!$       endif
!!$    enddo
!!$
!!$  End Subroutine Expand

End Module FullQ
