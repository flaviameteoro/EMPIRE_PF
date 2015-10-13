!> subroutine to call 3DVar for each particle
subroutine three_d_var_all_particles
  use comms
  use pf_control
  implicit none
  
  integer :: i
  
  do i = 1,pf%count
     call three_d_var(pf%psi(:,i))
  end do
  
  
end subroutine three_d_var_all_particles
