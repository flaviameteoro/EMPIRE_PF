subroutine configure_model
!Set up internal variables for HadCM3
  use hadcm3_config
  use hadcm3_data
  use sizes
  integer, parameter :: a_num_points= &
       4*a_nxn*a_nyn*a_levels+a_nxn*a_nyn
!       2*(a_nxn*(a_nyn-1))*a_levels+2*a_nxn*a_nyn*a_levels+a_nxn*a_nyn
  integer, parameter :: o_num_points= &
       4*o_nxn*o_nyn*o_levels+o_nxn*o_nyn
  integer, parameter :: num_points=a_num_points+o_num_points

  !Atmosphere only for now
  state_dim=a_num_points

end subroutine configure_model



subroutine solve_r(y,v)
  !subroutine to take an observation vector y and return v
  !in observation space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim), intent(out) :: v

  v = y
end subroutine solve_r

subroutine solve_q(x,v)
  !subroutine to take a full state vector x and return v
  !in state space.
  use sizes
 implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim), intent(in) :: x
  real(kind=rk), dimension(state_dim), intent(out) :: v

  v = x
end subroutine solve_q

subroutine solve_hqht_plus_r(y,v)
  !subroutine to take an observation vector y and return v
  !in observation space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  real(kind=rk), dimension(obs_dim), intent(out) :: v

  v = y
end subroutine solve_hqht_plus_r

subroutine Qhalf(x,Qx)
  !subroutine to take a full state vector x and return Q^(1/2)x
  !in state space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim), intent(in) :: x
  real(kind=rk), dimension(state_dim), intent(out) :: qx

  qx = sqrt(abs(x))
end subroutine QHALF

subroutine Q(x,Qx)
  !subroutine to take a full state vector x and return Qx
  !in state space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim), intent(in) :: x
  real(kind=rk), dimension(state_dim), intent(out) :: qx

  qx = x
end subroutine Q

subroutine H(x,hx)
  !subroutine to take a full state vector x and return H(x)
  !in observation space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim), intent(in) :: x
  real(kind=rk), dimension(obs_dim), intent(out) :: hx

  hx = x(1:obs_dim)
end subroutine H

subroutine HT(y,x)
  !subroutine to take an observation vector y and return x = H^T(y)
  !in full state space.
  use sizes
  implicit none
  integer, parameter :: rk=kind(1.0D+0)
  real(kind=rk), dimension(state_dim), intent(out) :: x
  real(kind=rk), dimension(obs_dim), intent(in) :: y

  x=0.0_rk
  x(1:obs_dim) = y
end subroutine HT
