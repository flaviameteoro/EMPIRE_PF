!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:36:44 pbrowne>
!!!
!!!    Subroutine to give output diagnositics such as rank histograms
!!!    and trajectories
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

!> Subroutine to give output diagnositics such as rank histograms  
!!    @todo test in anger with empire version 3. will probably segfault
subroutine diagnostics
  use output_empire, only : unit_hist_write,unit_hist_readt&
       &,unit_hist_readp,emp_e
  use timestep_data
  use pf_control
  use sizes
  use comms
  use histogram_data
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(rhl_n,pf%nens) :: HHpsi
  real(kind=rk), dimension(pf%nens) :: bin_marker
  real(kind=rk), dimension(rhl_n) :: y
  integer :: particle,time,i,j,mpi_err,lb,length
  logical :: placed
  character(32) :: filename
  integer, dimension(rhn_n,pf%nens+1) :: reduced_talagrand
  include 'mpif.h'

  if(.not. pf%gen_data) then
     if(pf%use_talagrand) then
        
        inquire(iolength=length) pf%psi(1,1)

        do particle = 1,pf%count
           !write(filename,'(A,i6.6,A,i5.5)') 'hist/timestep',((pf%timestep)/pf&
           !     &%time_bwn_obs) ,'particle',pf%particles(particle)
           if(comm_version .eq. 1 .or. comm_version .eq. 2) then
           write(filename,'(A,i7.7,A,i5.5)') 'hist/timestep',&
                &pf%timestep,'particle',pf%particles(particle)
           else
              write(filename,'(A,i7.7,A,i5.5,A,i0)') 'hist/timestep',&
                &pf%timestep,'particle',pf%particles(particle),'.',pf_member_rank
           end if
           open(unit_hist_write,file=filename,action='write',status='replace',form='unforma&
                &tted',access='direct',recl=length)

           do i = 1,rhl_n
              write(unit_hist_write,rec=i) pf%psi(rank_hist_list(i),particle)
           end do
           close(unit_hist_write)
        end do


!        if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
        if(TSdata%completed_timesteps .eq. TSdata%total_timesteps) then
           pf%talagrand = 0
           call mpi_barrier(pf_mpi_comm,mpi_err)
           
           do time = 1,pf%time_obs
              !              print*,'time = ',time
              if(mod(time,npfs) .eq. pf_ens_rank) then
                 !cock cock cock adjusted the below to make it sensible
                 !pf%timestep = time*pf%time_bwn_obs
                 pf%timestep = TSData%obs_times(time)
!                 print*,'pfrank = ',pfrank,'picking up truth at ',pf%timestep

                 !write(filename,'(A,i6.6,A,i5.5)') 'hist/timestep',((pf%timestep)/pf&
                 !     &%time_bwn_obs) ,'truth'
                 if(comm_version .eq. 1 .or. comm_version .eq. 2) then
                    write(filename,'(A,i7.7,A)') 'hist/timestep',&
                         &pf%timestep,'truth'
                 else
                    write(filename,'(A,i7.7,A,i0)') 'hist/timestep',&
                         &pf%timestep,'truth.',pf_member_rank
                 end if

                 open(unit_hist_readt,file=filename,action='read',status='old',form='unforma&
                      &tted',access='direct',recl=length)

                 do i = 1,rhl_n

                    read(unit_hist_readt,rec=i) y(i)

                 end do
                 close(unit_hist_readt)

                 do particle = 1,pf%nens

                    if(comm_version .eq. 1 .or. comm_version .eq. 2) then
                       write(filename,'(A,i7.7,A,i5.5)') 'hist/timestep',&
                            &pf%timestep,'particle',pf%particles(particle)
                    else
                       write(filename,'(A,i7.7,A,i5.5,A,i0)') 'hist/timestep',&
                            &pf%timestep,'particle',pf&
                            &%particles(particle),'.',pf_member_rank
                    end if
                    
                    open(unit_hist_readp, file=filename,action='read',status='old',&
                    form='unformatted',access='direct',recl=length)
                    do i = 1,rhl_n
                       read(unit_hist_readp,rec=i) HHpsi(i,particle)
                    end do
                    close(unit_hist_readp)
                 end do



                 do j = 1,rhn_n !for each histogram we want to make
                    if( j .eq. 1) then
                       lb = 1
                    else
                       lb = sum(rank_hist_nums(1:j-1)) + 1
                    end if
                 do i = lb,sum(rank_hist_nums(1:j))

                    do particle = 1,pf%nens
                       bin_marker(particle) = HHpsi(i,particle)
                    end do

                    call quicksort_d(bin_marker,pf%nens)
                    
                  
                    !                    print*,'quicksorted'
                    placed = .false.
                    do particle  = 1,pf%nens
                       if(y(i) .lt. bin_marker(particle)) then
                          pf%talagrand(j,particle) = pf&
                               &%talagrand(j,particle)&
                               & + 1
                          placed = .true.
                          exit
                       end if
                    end do


                    if(.not. placed) then
                       if(y(i) .ge. bin_marker(pf%nens)) then
                          pf%talagrand(j,pf%nens+1) = pf&
                               &%talagrand(j,pf%nens+1) + 1
                       else
                          write(emp_e,*) 'There was an error in the calculation of the p&
                               &lacement &
                               &in the rank histogram. Bums.'
                          stop
                       end if
                    end if

                 end do
                 end do

              end if !end of mpi splitting by timestep
           end do !end of the timestep
           !           print*,'end of all the timesteps'

           !now let us reduce the information to the master processor:
           call mpi_reduce(pf%talagrand,reduced_talagrand,rhn_n*(pf%nens+1)&
                &,MPI_INTEGER,MPI_SUM,0,pf_mpi_comm,mpi_err)

           !           print*,'some reduction just happened'

           if(pfrank .eq. 0) then
              do i = 1,rhn_n
                 write(filename,'(A,i0)') 'histogram_',i
                 open(unit_hist_write,file=filename,action='write',status='replace')
                 do j = 1,pf%nens+1
                    write(unit_hist_write,'(i8.8)') reduced_talagrand(i,j)
                 end do
                 close(unit_hist_write)
              end do
              !now output the image
           end if
           !pf%timestep = pf%time_obs*pf%time_bwn_obs
           pf%timestep = TSData%total_timesteps
        end if !end of if we are the last step in the pf

     end if


  else
     if(pf%use_talagrand) then
        !     print*,'in diagnostics at timestep ',pf%timestep
        inquire(iolength=length) pf%psi(1,1)

        do particle = 1,pf%count
           !write(filename,'(A,i6.6,A)') 'hist/timestep',((pf%timestep)/pf&
           !     &%time_bwn_obs) ,'truth'
           if(comm_version .eq. 1 .or. comm_version .eq. 2) then
              write(filename,'(A,i7.7,A)') 'hist/timestep',pf&
                   &%timestep,'truth'
           else
              write(filename,'(A,i7.7,A,i0)') 'hist/timestep',pf&
                   &%timestep,'truth.',pf_member_rank
           end if
           open(unit_hist_write,file=filename,action='write',status='replace',form='unforma&
                &tted',access='direct',recl=length)

           do i = 1,rhl_n
              write(unit_hist_write,rec=i) pf%psi(rank_hist_list(i),particle)
           end do
           close(unit_hist_write)
        end do
     end if
  end if
end subroutine diagnostics

