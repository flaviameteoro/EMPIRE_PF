!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2016-10-18 15:31:47 pbrowne>
!!!
!!!    module to deal with generating and outputting pf matrix
!!!    Copyright (C) 2015 Philip A. Browne
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

!> module to deal with generating and outputting pf matrix
module matrix_pf
  implicit none
  type, public :: matrix_pf_data
     !< Note: this feature is not accessible with empire version 3
     !! communications; the matrices in question are simply too
     !! large to compute the full \f$ P_f \f$ matrix.
     character(30) :: prefix  !< the prefix of the filename to be
     !!                          output
     integer :: k             !< the frequency to output the matrix
     logical :: analysis      !< if true, output at all analysis times
     logical :: frequency     !< if true, output at all timesteps
     !!                          that are 0 mod k
     integer :: output_type   !< output file type. \n
                              !! -  0 - undefined
                              !! -  1 - standard packed format (TP)
                              !! -  2 - rectangular full packed
                              !!                        format (TF)
                              !! \n
                              !! Negative values will be formatted.
                              !! \n
     !! Positive values will be unformatted.
     !! \n
     !! See http://www.netlib.org/lapack/explore-html/db/d37/dtfttp_8f.html for details of the storage stucture formats

  end type matrix_pf_data
  type(matrix_pf_data), save :: matpf !< module holding data for
  !! generating and outputting \f$ P_f  \f$ matrix. \n
  !! Note: this feature is not accessible with empire version 3
  !! communications; the matrices in question are simply too
  !! large to compute the full \f$ P_f \f$ matrix.
contains
  !> subroutine to read namelist to control this output
  subroutine read_matrix_pf_information
    use output_empire, only : unit_nml,emp_e
    implicit none
    character(30) :: prefix
    integer :: k=0
    logical :: analysis=.false.
    logical :: frequency=.false.
    integer :: output_type=0
    logical :: file_exists
    integer :: ios
    namelist/mat_pf/k,analysis,frequency,prefix,output_type

    inquire(file='pf_parameters.dat',exist=file_exists)
    if(file_exists) then
       open(unit_nml,file='pf_parameters.dat',iostat=ios,action='read'&
            &,status='old')
       if(ios .ne. 0) then
          write(emp_e,*) 'Cannot open pf_parameters.dat'
          stop
       end if
    else
       inquire(file='empire.nml',exist=file_exists)
       if(file_exists) then
          open(unit_nml,file='empire.nml',iostat=ios,action='read'&
               &,status='old')
          if(ios .ne. 0) then
             write(emp_e,*) 'Cannot open empire.nml'
             stop
          end if
       else
          write(emp_e,*) 'ERROR: cannot find pf_parameters.dat or empire.nml'
          stop '-1'
       end if
    end if

    read(unit_nml,nml=mat_pf,iostat=ios) 
    close(unit_nml)

    if(ios .ne. 0) then
       print*,'mat_pf not found in namelist file.'
       print*,'no matrix_pf information will be computed'
    else

       matpf%prefix=trim(prefix)
       matpf%analysis = analysis
       matpf%frequency= frequency
       matpf%output_type=output_type

       if(k.gt.0) then
          matpf%k=k
       else
          write(emp_e,*) 'ERROR: k in matrix_pf_data read in less than 1'
          write(emp_e,*) 'ERROR: k should be a natural number'
          write(emp_e,*) 'STOPPING'
          stop '-1'
       end if
    end if

    print*,'matpf=',matpf
    call flush(6)
  end subroutine read_matrix_pf_information
  
  !>subroutine to generate and output matrix Pf
  subroutine matrix_pf_output(root,comm,n,m,x,time,is_analysis)
    implicit none
    include 'mpif.h'
    integer, parameter :: rk = kind(1.0d0)
    integer, intent(in) :: root !< the process to output file
    integer, intent(in) :: comm !< the mpi communicator to build the
    !! matrix on
    integer, intent(in) :: n !< the size of the state vector
    integer, intent(in) :: m !< the number of state vectors on this process
    real(kind=rk), dimension(n,m) :: x !< the local ensemble members
    integer, intent(in) :: time !< the current timestep
    logical, intent(in) :: is_analysis !< true if analysis just performed

    character(40) :: filename
    logical :: actually_output=.false.
    real(kind=rk), dimension(n*(n+1)/2) :: pf
    integer :: rank,mpi_err

    actually_output = .false.

    !check if we should output:
    if(matpf%analysis .and. is_analysis) actually_output = .true.

    if(matpf%frequency .and. mod(time,matpf%k) .eq. 0) actually_output=.true.

    if(matpf%k .gt. 0) print*,'hello: mod(time,matpf%k) = ',mod(time,matpf%k)
    
    if(actually_output) then

       write(filename,'(A,i10.10)') trim(matpf%prefix),time
       
       call generate_pf(n,m,comm,x,pf)
       
       call MPI_COMM_RANK(comm,rank,mpi_err)
       if(rank .eq. root) call output_mat_tri(n,pf,filename,matpf%output_type)

    end if

  end subroutine matrix_pf_output
end module matrix_pf
