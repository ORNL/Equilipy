!
!
subroutine ShuffleAssemblage(iNewPhase, iPhaseTypeOut)
!
!   PURPOSE: Implement the Minimum Ratio Test for the Simplex method.
!            Reorganizes iAssemblage and dMolesPhase based on the ratio test.
!            Order: smallest positive ratio first, then larger positive ratios,
!            then negative ratios (most negative last).
!
!   The ratio is lambda_i / mu_i where:
!   - lambda = barycentric coordinates of bulk composition w.r.t. current simplex
!   - mu = barycentric coordinates of new phase w.r.t. current simplex
!
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer, intent(in)  :: iNewPhase
    integer, intent(out) :: iPhaseTypeOut
    
    integer :: i, j, INFO
    integer, dimension(nElements) :: IPIV, iSortIndex
    real(8), dimension(nElements) :: dLambda, dMu, dBulkComposition, dNewPhaseComposition
    real(8), dimension(nElements) :: dRatio, dRatioSorted
    real(8), dimension(nElements, 2) :: B
    real(8), dimension(nElements, nElements) :: A, A_copy
    
    ! Temporary storage for reorganization
    integer, dimension(nElements) :: iAssemblageTemp
    real(8), dimension(nElements) :: dMolesPhaseTemp
!
!   Initialize
    iPhaseTypeOut = 1
    iSortIndex = (/(i, i = 1, nElements, 1)/)
!
!   Step 1: Build matrix A from current phase compositions
!           A(j,i) = atom fraction of element j in phase i
    if (allocated(dAtomFractionSpeciesGEM)) then
        do i = 1, nElements
            do j = 1, nElements
                A(j, i) = dAtomFractionSpeciesGEM(i, j)
            end do
        end do
    else
        do i = 1, nElements
            do j = 1, nElements
                A(j, i) = dAtomFractionSpecies(iAssemblage(i), j)
            end do
        end do
    end if
!
!   Step 2: Compute bulk composition (normalized moles)
    dBulkComposition = dMolesElement / sum(dMolesElement)
!
!   Step 3: Compute new phase composition
    dNewPhaseComposition = dAtomFractionSpecies(iNewPhase, :)
!
!   Step 4: Solve A * B = [bulk, new_phase] for barycentric coordinates
!           B(:,1) = lambda (bulk composition)
!           B(:,2) = mu     (new phase composition)
    A_copy = A
    B(:,1) = dBulkComposition
    B(:,2) = dNewPhaseComposition
    INFO = 0
    IPIV = 0
    call DGESV(nElements, 2, A_copy, nElements, IPIV, B, nElements, INFO)
    
    if (INFO /= 0) then
        ! Matrix is singular, don't reorganize
        return
    end if
    dLambda = B(:,1)
    dMu     = B(:,2)
!
!   Step 5: Calculate ratio for each phase
!           For mu_i <= 0, set ratio to large positive value (will be sorted last)
    do i = 1, nElements
        if (dMu(i) > 1D-12) then
            dRatio(i) = dLambda(i) / dMu(i)
        else
            ! mu <= 0: set to large value so it goes to the end
            dRatio(i) = 1D30
        end if
    end do
!
!   Step 6: Sort by ratio (ascending order - smallest positive first)
    dRatioSorted = dRatio
    call Qsort(dRatioSorted, iSortIndex, nElements)
!
!   Step 7: Save current state
    iAssemblageTemp = iAssemblage
    dMolesPhaseTemp = dMolesPhase
!
!   Step 8: Reorganize iAssemblage and dMolesPhase based on sorted order
    do i = 1, nElements
        j = iSortIndex(i)
        iAssemblage(i) = iAssemblageTemp(j)
        dMolesPhase(i) = dMolesPhaseTemp(j)
    end do
!
!   Step 9: Also update iShuffled for use in GetNewAssemblage
    if (allocated(iShuffled)) deallocate(iShuffled)
    allocate(iShuffled(nElements))
    iShuffled = iSortIndex
!
!   The first position now has the phase with smallest positive ratio
!   (the best candidate to be replaced)
    iPhaseTypeOut = 1
!
    return
!
end subroutine ShuffleAssemblage
!
!
