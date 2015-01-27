!> subroutine to compute the 4DVar objective function
!!
!! \f$ (x-x_0)^TB^{-1}(x-x_0) + \sum_{i=1}^K [y_i-H_i(x(t_i))]^T R_i
!! [y_i - H_i( x( t_i ) )]  \f$
subroutine objective_function(n,x,f)
  use var_data
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: n
  real(kind=rk), dimension(n), intent(in) :: x
  real(kind=rk), dimension(n) :: xt
  real(kind=rk), intent(out) :: f
  integer :: t
  integer :: tag
  real(kind=rk), dimension(:), allocatable :: y,hxt
  real(kind=rk) :: ft

  !initialise f
  f = 0.0d0

  !add background term
  call inner_b_1(x-vardata%x0,f)
  
  
  do t = 1,vardata%total_timesteps
     !advance model to timestep t
     if(t == 1) then
        tag = -1
     else
        tag = 1
     end if
     call mpi_send
     call mpi_recv
     

     if(vardata%ny(t) .gt. 0) then
        allocate(  y(vardata%ny(t)))
        allocate(hxt(vardata%ny(t)))
        call H(vardata%ny(t),1,xt,hxt,t)
        call get_observation_data(y)
        y = y - hxt
        call innerR_1(vardata%ny(t),1,y,ft,t)
        
        deallocate(y)
        deallocate(hxt)

        f = f+ft
     end if



  end do







end subroutine objective_function

!> subroutine to compute the inner product of a vector x 
!> in the matrix norm \f$B^{-1}\f$
!>
!> \f$f = x^TB^{-1}x\f$
subroutine inner_b_1(x,f)
  use var_data
  !> the vector to make the product with \f$B^{-1}\f$
  real(kind=kind(1.0d0)), dimension(vardata%n), intent(in) :: x
  !> the result \f$f = x^TB^{-1}x\f$ 
  real(kind=kind(1.0d0)), intent(out) :: f
  real(kind=kind(1.0d0)), dimension(vardata%n) :: temp
  real(kind=kind(1.0d0)) :: ddot

  call solve_b(1,x,temp)
  f = ddot(vardata%n,temp,1,x,1)

end subroutine inner_b_1
