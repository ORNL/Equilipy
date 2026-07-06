!
!> subroutine to sort using the insertionsort algorithm and return indecies
!! @param[in,out] a, an array of doubles to be sorted
!! @param[in,out] idx_a, an array of integers of sorted indecies
!! @param[in] na, dimension of the array a 
subroutine InsertionSort_d(a,idx_a,na)
!
  implicit none
!
  ! DUMMY ARGUMENTS
  integer,intent(in) :: na
  real(8), dimension(nA), intent(inout) :: a
  integer,dimension(nA), intent(inout) :: idx_a
!
  ! LOCAL VARIABLES
  real(8) :: temp
  integer:: i, j
  integer:: idx_tmp
!
  do i = 2, nA
     j = i - 1
     temp = A(i)
     idx_tmp = idx_a(i)
     do
        if (j == 0) exit
        if (a(j) <= temp) exit
        A(j+1) = A(j)
        idx_a(j+1) = idx_a(j)
        j = j - 1
     end do
     a(j+1) = temp
     idx_a(j+1) = idx_tmp
  end do
!
end subroutine InsertionSort_d
!
!
