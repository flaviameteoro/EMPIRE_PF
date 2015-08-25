!> subroutine to take a full state vector x and return \f$P^{1/2}x\f$
!> in state space.
!!
!! Given \f$x\f$ compute \f$P^{\frac{1}{2}}x\f$
!!
!! where \f$P = (Q^{-1} +H^TR^{-1}H)^{-1} = Q^{\frac{1}{2}}(I +
!! Q^{\frac{1}{2}}H^TR^{-1}HQ^{\frac{1}{2}})^{-1}Q^{\frac{1}{2}}\f$ 
!!
!! This is required for the Zhu Equal weights particle filter
!! @ref equivalent_weights_filter_zhu
subroutine Phalf(nrhs,x,Px)
  use sizes
  use Qdata
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  integer, intent(in) :: nrhs !< the number of right hand sides
  real(kind=rk), dimension(state_dim,nrhs), intent(in) :: x !< the
  !!input vector
  real(kind=rk), dimension(state_dim,nrhs), intent(out) :: px !< the
  !!resulting vector where Px \f$= P^{\frac{1}{2}}x\f$


  call phalf_etkf(nrhs,x,Px)
end subroutine Phalf
