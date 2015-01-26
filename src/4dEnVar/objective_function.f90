subroutine objective_function(n,x,f)
  use vardata
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: n
  real(kind=rk), dimension(n), intent(in) :: x
  real(kind=rk), intent(out) :: f
  
  !initialise f
  f = 0.0d0

  !add background term
  call inner_b_minus_1(x-vardata%x0,f)
  
  
  do t = 1,timesteps
     !advance model to timestep t
     if(t == 1) then
        tag = -1
     else
        tag = 1
     end if
     call mpi_send
     call mpi_recv
     

     if(observation time) then
        call H(1,xt,hxt,t)
        call get_observation_data(y)
        y_hxt = y - hxt
        call inner_R_1(y_hxt,ft)
        f = f+ft
     end if



  end do







end subroutine objective_function

!> subroutine to compute the inner product of a vector x 
!> in the matrix norm \f$B^{-1}\f$
!> \f$f = x^TB^{-1}x
subroutine inner_b_minus_1(x,f)
  use vardata
  real(kind=kind(1.0d0)), dimension(vardata%n), intent(in) :: x
  real(kind=kind(1.0d0)), intent(out) :: f
  real(kind=kind(1.0d0))`, dimension(vardata%n) :: temp
  real(kind=kind(1.0d0)) :: ddot

  call solve_b(1,x,temp)
  f = ddot(vardata%n,temp,1,x,1)

end subroutine inner_b_minus_1
