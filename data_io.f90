subroutine get_observation_data(y)
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(obs_dim), intent(out) :: y
  integer :: obs_number,ios
  character(14) :: filename

  obs_number = ((pf%timestep-1)/pf%time_bwn_obs) + 1

  write(filename,'(A,i6.6)') 'obs_num_',obs_number

  open(67,file=filename,iostat=ios,action='read',status='old',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Check it exists. I need a lie down.'
     stop
  end if
  read(67) y
  close(67)
end subroutine get_observation_data

subroutine save_observation_data(y)
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(obs_dim), intent(in) :: y
  integer :: obs_number,ios
  character(14) :: filename

  obs_number = ((pf%timestep-1)/pf%time_bwn_obs) + 1

  write(filename,'(A,i6.6)') 'obs_num_',obs_number

  open(67,file=filename,iostat=ios,action='write',status='replace',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Very strange that I couldn''t open it. I''m going to stop now.'
     stop
  end if
  write(67) y
  close(67)
end subroutine save_observation_data
