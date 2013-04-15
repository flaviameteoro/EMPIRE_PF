program test_model_specific
use sizes
implicit none
integer, parameter :: rk = kind(1.0D0)
integer :: i
real(kind=rk), allocatable, dimension(:) :: x,y,xres,yres,xtemp,ytemp

call configure_model


allocate(x(state_dim),xres(state_dim),xtemp(state_dim))
allocate(y(obs_dim),yres(obs_dim),ytemp(obs_dim))

x = (/ (1.0_rk*i, i=1,state_dim) /)
y = (/ (1.0_rk*i, i=1,obs_dim) /)

print*,'x:'
print*,x
print*,'y:'
print*,y

call Q(x,xres)
print*,'Qx:'
print*,xres

call Qhalf(x,xtemp)
call Qhalf(xtemp,xres)
print*,'QhalfQhalfx:'
print*,xres

call solve_r(y,yres)
print*,'R^-1y:'
print*,yres
call solve_rhalf(y,ytemp)
call solve_rhalf(ytemp,yres)
print*,'R^-1/2R^-1/2y:'
print*,yres


deallocate(x,xres,xtemp,y,yres,ytemp)
end program test_model_specific
