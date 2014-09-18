!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-18 10:09:54 pbrowne>
!!!
!!!    {one line to give the program's name and a brief idea of what it does.}
!!!    Copyright (C) 2014  Philip A. Browne
!!!
!!!    This program is free software: you can redistribute it and/or modify
!!!    it under the terms of the GNU General Public License as published by
!!!    the Free Software Foundation, either version 3 of the License, or
!!!    (at your option) any later version.
!!!
!!!    This program is distributed in the hope that it will be useful,
!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!!    GNU General Public License for more details.
!!!
!!!    You should have received a copy of the GNU General Public License
!!!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!!
!!!    Email: p.browne @ reading.ac.uk
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagnostics
  use pf_control
  use sizes
  use comms
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(obs_dim,pf%nens) :: HHpsi
  real(kind=rk), dimension(pf%nens) :: bin_marker
  real(kind=rk), dimension(obs_dim) :: y
  real(kind=rk), dimension(obs_dim,pf%count) :: Hpsi
  integer :: particle,time,i,j,mpi_err
  logical :: placed
  character(27) :: filename
  integer, dimension(9,pf%nens+1) :: reduced_talagrand
  include 'mpif.h'
  
  if(.not. pf%gen_data) then
     if(pf%use_talagrand) then
        
        call H(pf%psi,Hpsi)
        
        do particle = 1,pf%count
           write(filename,'(A,i6.6,A,i5.5)') 'timestep',((pf%timestep-1)/pf&
                &%time_bwn_obs) + 1,'particle',pf%particles(particle)
           open(12,file=filename,action='write',status='replace',form='unforma&
                &tted')
!           call H(pf%psi(:,particle),Hpsi)
           write(12) Hpsi(:,particle)
           close(12)
        end do


        if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
           pf%talagrand = 0
           call mpi_barrier(pf_mpi_comm,mpi_err)

           do time = 1,pf%time_obs
!              print*,'time = ',time
              if(mod(time,npfs) .eq. pfrank) then
              pf%timestep = time
              call get_observation_data(y)
!              print*,'got obs data at timestep ',pf%timestep
              do particle = 1,pf%nens
                 write(filename,'(A,i6.6,A,i5.5)') 'timestep',pf%timestep,'particle',particle
                 open(13, file=filename,action='read',status='old',form='unfor&
                      &matted')
                 read(13) HHpsi(:,particle)
                 close(13)
              end do
 !             print*,'read HHpsi'
           
              do i = 1,obs_dim

                    do particle = 1,pf%nens
                       bin_marker(particle) = HHpsi(i,particle)
                    end do
                    
                    call quicksort_d(bin_marker,pf%nens) 
!                    print*,'quicksorted'
                    do particle  = 1,pf%nens
                       if(y(i) .lt. bin_marker(particle)) then
                          pf%talagrand(i,particle) = pf&
                               &%talagrand(i,particle)&
                               & + 1
                          placed = .true.
                          exit
                       end if
                    end do
!                    print*,'did we place?'

                    if(.not. placed) then
                       if(y(i) .ge. bin_marker(pf%nens)) then
                          pf%talagrand(i,pf%nens+1) = pf&
                               &%talagrand(i,pf%nens+1) + 1
                       else
                          stop 'There was an error in the calculation of the placement &
                               &in the rank histogram. Bums.'
                       end if
                    end if
                 
                 end do

           end if !end of mpi splitting by timestep
           end do !end of the timestep
!           print*,'end of all the timesteps'

           !now let us reduce the information to the master processor:
           call mpi_reduce(pf%talagrand,reduced_talagrand,9*(pf%nens+1)&
                &,mpi_integer,mpi_sum,0,pf_mpi_comm,mpi_err)
           
!           print*,'some reduction just happened'

           if(pfrank .eq. 0) then
              do i = 1,1
                 write(filename,'(A,i1.0)') 'histogram_',i
                 open(17,file=filename,action='write',status='replace')
                 do j = 1,pf%nens+1
                    write(17,'(i8.8)') reduced_talagrand(i,j)
                 end do
                 close(17)
              end do
              !now output the image
           end if
        end if !end of if we are the last step in the pf


!        do particle = 1,pf%nens
!           call H(pf%psi(:,particle),Hpsi(:,particle))
!           bin_marker(particle) = Hpsi(i)
!        end do

!        call get_observation_data(y)

!        do i = 1,obs_dim
!           do particle = 1,pf%nens
!              bin_marker(i,particle) = Hpsi(i,particle)
!           end do
           

!           call quicksort_d(bin_marker(i,:),pf%nens)
           
!        call get_observation_data(y)

!        print*,'observation is:',y(pf%tala_obs_num)
!        print*,'bin_marker',bin_marker

!           do particle = 1,pf%nens
!              if(y(i) .lt. bin_marker(i,particle)) then
!                 pf%talagrand(particle,i) = pf%talagrand(particle,i) + 1
!                 placed = .true.
!                 exit
!              end if
!           end do
!           
!           if(.not. placed) then
!              if(y(i) .ge. bin_marker(i,pf%nens)) then
!                 pf%talagrand(pf%nens+1,i) = pf%talagrand(pf%nens+1,i) + 1
!              else
!                 stop 'There was an error in the calculation of the placement &
!                      &in the rank histogram. Bums.'
!              end if
!           end if
!           
!        end do
        !now output the stuff if we are at the last timestep...
!        if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
!           open(78,file='pf_talagrand',iostat=ios,action='write',status='replace')    
!           if(ios .ne. 0)  then
!              write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_talagrand'
!              write(*,*) 'Very strange that I couldn''t open it. I''m going to stop now.'
!              stop
!           end if

!           do particle = 1,pf%nens+1
!              write(78,*) pf%talagrand(particle,:)
!           end do

!           close(78)
!           open(78,file='gnuplot_talagrand.cfg',iostat=ios,action='write',status='replace')    
!           if(ios .ne. 0)  then
!              write(*,*) 'Cannot open file gnuplot_talagrand'
!              write(*,*) 'Very strange that I couldn''t open it. I''m going to stop now.'
!              stop
!           end if
!           write(78,'(A)') 'reset'
!           !      write(78,'(A)') 'set term png truecolor'
!           !      write(78,'(A)') 'set output "pf_talagrand.png"'
!           write(78,'(A)') 'set xlabel "Observation rank"'
!           write(78,'(A)') 'set ylabel "Frequency"'
!           write(78,'(A)') 'set grid'
!           write(78,'(A)') 'set boxwidth 0.95 relative'
!           write(78,'(A)') 'set style fill transparent solid 0.5 noborder'
!           write(78,'(A,i0,A)') 'set xrange [-0.5:',pf%nens,'.5]'
!           write(78,'(A)') 'plot "pf_talagrand" w boxes notitle'
!           close(78)
!        end if
     end if
  end if

end subroutine diagnostics
