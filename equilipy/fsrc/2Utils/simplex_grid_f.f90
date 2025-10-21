subroutine simplex_grid_f(L, nVer, nSpacing, res)
    !>
    !> @param[in] L         The number of grid points to generate (rows).
    !> @param[in] nVer      The number of vertices (columns).
    !> @param[in] nSpacing  The integer value to be partitioned.
    !> @param[out] res      The resulting (L, nVer) grid, normalized by nSpacing.
    integer, intent(in) :: L
    integer, intent(in) :: nVer
    integer, intent(in) :: nSpacing
    real, dimension(L, nVer), intent(out) :: res
    
    ! --- Local variables ---
    integer, dimension(nVer) :: vertices  ! Stores the integer partitions
    integer :: h, val, i, j
    real :: nSpacing_real
    
    ! --- Initialization ---
    vertices = 0
    vertices(nVer) = nSpacing
    
    ! --- Copy first row ---
    do j = 1, nVer
        ! Implicit conversion from integer to real
        res(1, j) = vertices(j)
    end do
    
    ! --- Main loop to generate subsequent rows ---
    h = nVer + 1
    
    do i = 2, L
        h = h - 1
        
        val = vertices(h)
        
        vertices(h) = 0
        
        vertices(nVer) = val - 1
        
        vertices(h-1) = vertices(h-1) + 1
        
        ! Copy the new 'vertices' state to the result row
        do j = 1, nVer
            res(i, j) = vertices(j)
        end do
        
        ! (Resets h to one-past-the-end)
        if (val /= 1) then
            h = nVer + 1
        end if
    
    end do
    
    ! --- Final Normalization ---
    
    ! Cast nSpacing to real for floating-point division
    nSpacing_real = real(nSpacing)
    
    ! Fortran can perform element-wise array operations
    res = res / nSpacing_real
    
end subroutine simplex_grid_f