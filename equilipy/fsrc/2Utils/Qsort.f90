!
!
recursive subroutine Qsort(a,idx_a,na)
    implicit none 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!    Collection of subroutines to sort and return a one-dimensional array
    !!!    as well as corresponding sorted index of the array a. Original code 
    !!!    (distributed under GNU Free licence 1.2) was taken from 
    !!!    http://rosettacode.org/wiki/Quicksort#Fortran and modified to 
    !!!    also return sorted index of the original array a.
    !!!    Copyright (C) 2015  Sanita Vetra-Carvalho
    !!!
    !!!    This program is distributed under the Lesser General Public License (LGPL) version 3,
    !!!    for more details see <https://www.gnu.org/licenses/lgpl.html>.
    !!!
    !!!    Email: s.vetra-carvalho @ reading.ac.uk
    !!!    Mail:  School of Mathematical and Physical Sciences,
    !!!    	      University of Reading,
    !!!	      Reading, UK
    !!!	      RG6 6BB
    !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> subroutine to sort using the quicksort algorithm
    !! @param[in,out] a, an array of doubles to be sorted
    !! @param[out] idx_a, an array of sorted indecies of original array a
    !! @param[in] na, dimension of the array a
    ! DUMMY ARGUMENTS
    integer, intent(in) :: na ! nr or items to sort
    real(8), dimension(nA), intent(inout) :: a ! vector to be sorted
    integer, dimension(nA), intent(out) :: idx_a ! sorted indecies of a
!
    ! LOCAL VARIABLES
    integer :: left, right, mid, i
    real(8) :: pivot, temp
    integer :: marker, idx_temp
!
!
    idx_a=(/(i, i=1, na, 1)/)
!
    if (nA > 1) then
        ! insertion sort limit of 47 seems best for sorting 10 million
        ! integers on Intel i7-980X CPU.  Derived data types that use
        ! more memory are optimized with smaller values - around 20 for a 16
        ! -byte type.
        if (nA > 47) then
            ! Do quicksort for large groups
            ! Get median of 1st, mid, & last points for pivot (helps reduce
            ! long execution time on some data sets, such as already
            ! sorted data, over simple 1st point pivot)
            mid = (nA+1)/2
            if (a(mid) >= a(1)) then
                if (a(mid) <= a(nA)) then
                    pivot = a(mid)
                else if (a(nA) > a(1)) then
                    pivot = a(nA)
                else
                    pivot = a(1)
                end if
            else if (a(1) <= a(nA)) then
                pivot = a(1)
            else if (a(nA) > a(mid)) then
                pivot = a(nA)
            else
                pivot = a(mid)
            end if
!
            left = 0
            right = nA + 1
!
            do while (left < right)
                right = right - 1
                do while (A(right) > pivot)
                    right = right - 1
                end do
!
                left = left + 1
                do while (A(left) < pivot)
                    left = left + 1
                end do
!
                if (left < right) then
                    temp = A(left)
                    idx_temp = idx_a(left)
                    A(left) = A(right)
                    idx_a(left) = idx_a(right)
                    A(right) = temp
                    idx_a(right) = idx_temp
                end if
            end do
!
            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if
!
            call Qsort(A(:marker-1),idx_A(:marker-1),marker-1)
            call Qsort(A(marker:),idx_A(marker:),nA-marker+1)
!
        else
            call InsertionSort_d(A,idx_a,nA)    ! Insertion sort for small groups is
            !  faster than Quicksort
        end if
    end if
!
end subroutine Qsort
!
!
!
