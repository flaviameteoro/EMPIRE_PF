subroutine diagnostics
  use pf_control
  use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(pf%ngrand) :: bin_marker
  real(kind=rk), dimension(obs_dim) :: Hpsi,y
  integer :: particle,ios
  logical :: placed


  if(.not. pf%gen_data) then
     if(pf%use_talagrand) then

        do particle = 1,pf%ngrand
           call H(pf%psi(:,particle),Hpsi)
           bin_marker(particle) = Hpsi(pf%tala_obs_num)
        end do

        call quicksort_d(bin_marker,pf%ngrand)

        call get_observation_data(y)

        print*,'observation is:',y(pf%tala_obs_num)
        print*,'bin_marker',bin_marker

        do particle = 1,pf%ngrand
           if(y(pf%tala_obs_num) .lt. bin_marker(particle)) then
              pf%talagrand(particle,1) = pf%talagrand(particle,1) + 1
              placed = .true.
              exit
           end if
        end do

        if(.not. placed) then
           if(y(pf%tala_obs_num) .ge. bin_marker(pf%ngrand)) then
              pf%talagrand(pf%ngrand+1,1) = pf%talagrand(pf%ngrand+1,1) + 1
           else
              stop 'There was an error in the calculation of the placement &
                   &in the rank histogram. Bums.'
           end if
        end if

        !now output the stuff if we are at the last timestep...
        if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
           open(78,file='pf_talagrand',iostat=ios,action='write',status='replace')    
           if(ios .ne. 0)  then
              write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_talagrand'
              write(*,*) 'Very strange that I couldn''t open it. I''m going to stop now.'
              stop
           end if

           do particle = 1,pf%ngrand+1
              write(78,*) pf%talagrand(particle,:)
           end do

           close(78)
           open(78,file='gnuplot_talagrand.cfg',iostat=ios,action='write',status='replace')    
           if(ios .ne. 0)  then
              write(*,*) 'Cannot open file gnuplot_talagrand'
              write(*,*) 'Very strange that I couldn''t open it. I''m going to stop now.'
              stop
           end if
           write(78,'(A)') 'reset'
           !      write(78,'(A)') 'set term png truecolor'
           !      write(78,'(A)') 'set output "pf_talagrand.png"'
           write(78,'(A)') 'set xlabel "Observation rank"'
           write(78,'(A)') 'set ylabel "Frequency"'
           write(78,'(A)') 'set grid'
           write(78,'(A)') 'set boxwidth 0.95 relative'
           write(78,'(A)') 'set style fill transparent solid 0.5 noborder'
           write(78,'(A,i0,A)') 'set xrange [-0.5:',pf%ngrand,'.5]'
           write(78,'(A)') 'plot "pf_talagrand" w boxes notitle'
           close(78)
        end if
     end if
  end if

end subroutine diagnostics
