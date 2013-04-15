subroutine K(y,x)
  !subroutine to apply the operator K to a vector y in obs space and return
  !the vector x in full state space.
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim), intent(out) :: x
  real(kind=rk), dimension(obs_dim), intent(in) :: y

  real(kind=rk), dimension(obs_dim) :: v
  real(kind=rk), dimension(state_dim) :: vv

  call solve_hqht_plus_r(y,v)

  call HT(v,vv)

  call Q(vv,x)

end subroutine K

subroutine innerR_1(y,w)
  !subroutine to take an observation vector y and return w = y^T R^(-1) y
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim) :: v
  real(kind=rk), intent(out) :: w

  call solve_r(y,v)

  !this can defo be done better using BLAS PAB...
  w = sum(y*v)

end subroutine innerR_1

subroutine innerQ_1(x,w)
  !subroutine to take a full state vector x and return w = x^T Q^(-1) x
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim), intent(in) :: x
  real(kind=rk), dimension(state_dim) :: v
  real(kind=rk), intent(out) :: w

  call solve_q(x,v)

  !this can defo be done better using BLAS PAB...
  w = sum(x*v)

end subroutine innerQ_1

subroutine innerHQHt_plus_R_1(y,w)
  !subroutine to take an observation vector y and return w = y^T (HQH^T+R)^(-1) y
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim) :: v
  real(kind=rk), intent(out) :: w

  call solve_hqht_plus_r(y,v)

  !this can defo be done better using BLAS PAB...
  w = sum(y*v)

end subroutine innerHQHt_plus_R_1


subroutine B(y,x)
use pf_control
use sizes
implicit none
integer, parameter :: rk = kind(1.0D0)
real(kind=rk), dimension(obs_dim), intent(in) :: y
real(kind=rk), dimension(state_dim), intent(out) :: x
real(kind=rk), dimension(obs_dim) :: R_1y
real(kind=rk), dimension(state_dim) :: HtR_1y,QHtR_1y
real(kind=rk) :: freetime,p,tau

freetime = 0.6_rk

tau = real(modulo(pf%timestep,pf%time_bwn_obs),rk)/real(pf&
     &%time_bwn_obs,rk)

if(tau .lt. freetime) then
   x = 0.0_rk
else
   
   call solve_rhalf(y,R_1y)
   
   call HT(R_1y,HtR_1y)
   
   call Qhalf(HtR_1y,QHtR_1y)

   p = pf%nudgefac*(tau-freetime)/(1.0_rk-freetime)

   x = p*QHtR_1y
end if

end subroutine B
