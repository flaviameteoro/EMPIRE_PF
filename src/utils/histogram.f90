!!! Time-stamp: <2014-02-26 15:38:58 pbrowne>

module histogram_data
  integer, allocatable, dimension(:) :: rank_hist_list
  integer, allocatable, dimension(:) :: rank_hist_nums
  integer :: rhl_n,rhn_n

contains
  subroutine load_histogram_data
    implicit none
    integer :: i

!    rhn_n = 9
    open(2,file='variables_hist.dat',action='read',status='old')
    read(2,'(i7.7)') rhn_n
    allocate(rank_hist_nums(rhn_n))
    do i = 1,rhn_n
       read(2,'(i7.7)') rank_hist_nums(i)
    end do
    rhl_n = sum(rank_hist_nums)
    allocate(rank_hist_list(rhl_n))
    do i = 1,rhl_n
       read(2,'(i7.7)') rank_hist_list(i)
    end do
    close(2)
  end subroutine load_histogram_data

  subroutine kill_histogram_data
    deallocate(rank_hist_list)
    deallocate(rank_hist_nums)
  end subroutine kill_histogram_data
end module histogram_data
