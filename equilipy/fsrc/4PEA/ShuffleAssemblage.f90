!
!
subroutine ShuffleAssemblage(iNewPhase,iPhaseTypeOut)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ShuffleAssemblage.f90
    !> \brief   Shuffle the phase assemblage in the order that is most favorable for phase exchange.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      SortPick.f90
    !> \sa      CheckPureConPhaseAdd.f90
    !> \sa      CheckPureConPhaseRem.f90
    !> \sa      CheckSolnPhaseAdd.f90
    !> \sa      GetNewAssemblage.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   09/22/2011      M.H.A. Piro         Original code
    !   10/21/2011      M.H.A. Piro         Clean up code: modules
    !   01/19/2012      M.H.A. Piro         Added the capability of sorting solution phases.
    !   02/14/2012      M.H.A. Piro         Sort the assemblage in terms of the Euclidean norm.
    !   03/06/2012      M.H.A. Piro         Added the capability to shuffle the assemblage when the new
    !                                       phase and the phase to be removed are both solution phases.
    !   04/26/2012      M.H.A. Piro         Implementing Gibbs energy Minimization algorithm and dOxygen.
    !   09/24/2012      M.H.A. Piro         Sort the entire phase assemblage, not just one type of phase.
    !                                        The advantage of this is that another subroutine can determine
    !                                        whether it is best to swap a pure condensed or solution phase first.
    !   09/29/2012      M.H.A. Piro         If the system contains a pair of miscible phases, then set the
    !                                        Euclidean Norm to an arbitrarily large value to place these
    !                                        phases to the back of the list.
    !   11/4/2012       M.H.A. Piro         Fixed typo in Euclidean norm vector.  Coefficients that do not
    !                                        correspond to a phase should have a value of 1000 instead of 0,
    !                                        otherwise the sorting routine will place it at the top.
    !   04/22/2021      S.Y. Kwon           dEuclideanNorm is modified to dEuclidean vector
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to shuffle the current estimated phase assemblage such that
    !! the order of phases in iAssemblage are in oder of atomic similarity to a new phase that is to be added
    !! to the system.  The Gibbs Phase Rule dictates that the maximum number of phases that can coexist cannot
    !! exceed the number of system components (taken here as chemical elements).  Therefore, if the number of
    !! phases currently expected to be stable is equal to the number of elements, another phase must be
    !! withdrawn from the system to accomodate this new phase.
    !!
    !! The principle of this technique is to quantify the atomic similarity between the new phase and all other
    !! phases in the system.  The phase with the most similar atomic constituency is the most likely best
    !! candidate. The principle of this technique is to compute the Euclidean norm vector in nElements dimensional
    !! space, where one point is the stoichiometry of the new phase to be introduced and each other point
    !! represents the other phases in the current phase assemblage.  The iAssemblage vector is reorganized
    !! in descending order of the Euclidean Norm.
    !!
    !! The principle of this technique are discussed in greater detail in the following literature:
    !! - M.H.A. Piro and S. Simunovic, "Performance Enhancing Algorithms for Computing Thermodynamic
    !!   Equilibria," CALPHAD, 39 (2012) 104-110.
    !!
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   iNewPhase       Index of the new phase to be introduced.  If it is positive, it is a pure
    !!                               condensed phase.  If it is negative, it is a solution phase.
    !> \param[in]   iPhaseTypeOut   An integer idenitifying the type of phase with the smallest Euclidean Norm.
    !!                               This can be used to determine whether a pure condensed phase should be
    !!                               swapped first or a solution phase.
    !!
    !!                               iPhaseTypeOut = 0: pure condensed phase;
    !!                               iPhaseTypeOut = 1: solution phase.
    !
    ! iAssemblage           An integer vector containing the indices of phases in the assemblage
    !                        (1:nConphases represent pure condensed phases and
    !                        (nElements-nSolnPhases:nElements) represent solution phases.
    ! nElements             The number of elements in the system.
    ! nConPhases            The number of pure condensed phases in the assemblage
    ! nSolnPhases           The number of solution phases in the assemblage
    ! nSolnPhasesSys        The number of solution phases in the system
    ! dMolesPhase           The number of moles of a phase.  These are directly mapped to phases in iAssemblage
    ! dEuclideanNorm        A double vector representing the Euclidean norm between the stoichiometry of the new
    !                       phase each existing pure condensed phase in the assemblage.
    ! dAtomFractionSpecies  Atomic fraction of a particular element (column) for a particular species (row).
    ! dEffStoichSolnPhase   The effective stoichiometry of a particular element (column) in a particular solution
    !                       phase (row).
    ! iVec                  An integer vector used for internal operations for the sorting routine.
    ! iTempVec              An integer vector storing the phase assemblage at the beginning of the calculation
    !                        to allow the system to be reverted in case if ShuffleAssemblage fails.
    ! dTempVec              A double real temporary vector storing the number of moles of all phases at the
    !                        beginning of the calculation to allow the system to be reverted in case if
    !                        ShuffleAssemblage fails.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer                      :: i, j, iNewPhase, iPhaseTypeOut
    integer,dimension(nElements) :: iTempVec, iShuffle
    real(8)                      :: dLengthNew, dLengthOld
    real(8),dimension(nElements) :: dTempVec, dAtomFractionBulk
    real(8),dimension(nElements) :: dVecNew,dVecOld, dEuclideanVec
    real(8),dimension(nElements,nElements) :: dAtomFractionSpeciesTemp

    if (allocated(dAtomFractionSpeciesGEM)) then
        dAtomFractionSpeciesTemp=dAtomFractionSpeciesGEM
    else
        do i =1,nElements
            dAtomFractionSpeciesTemp(i,:)=dAtomFractionSpecies(iAssemblage(i),:)
        end do
    end if

    if (allocated(iShuffled))               deallocate(iShuffled)
    allocate(iShuffled(nElements))

    ! Initialize variables:
    iShuffled         = (/(i, i = 1,nElements, 1)/)
    iShuffle          = (/(i, i = 1,nElements, 1)/)
    dEuclideanVec     = -1D0
    dAtomFractionBulk = dMolesElement/sum(dMolesElement) !Bulk(input) composition
    iTempVec          = 0
    dTempVec          = 0D0

    !Conduct calcualtion based on EuclideanVec
    dVecNew        = dAtomFractionSpecies(iNewPhase,:) - dAtomFractionBulk(:)
    dVecNew        = dVecNew/sum(dVecNew**2)
    do i = 1, nElements
        dVecOld=(dAtomFractionSpeciesTemp(i,:) - dAtomFractionBulk(:))
        dVecOld = dVecOld/sum(dVecOld**2)
        dEuclideanVec(i) = dot_product(dVecNew,dVecOld)
    end do
   
!
!
    ! Swap the phase with the lowest Euclidean norm for the first phase in the assemblage:
    IF_Euclid: if (iNewPhase /= 0) then
!
        call Qsort(dEuclideanVec, iShuffle, nElements)
!
        ! Store temporary variables:
        iTempVec      = iAssemblage
        dTempVec      = dMolesPhase
!
        ! Reinitialize variables:
        iAssemblage = 0
        dMolesPhase = 0d0

        ! Shuffle the phase assemblage based on dEuclideanVec
        do i = 1, nElements
            iShuffled(i) = iShuffle(nElements-i+1)

            j = iShuffle(nElements-i+1)
            iAssemblage(i) = iTempVec(j)
            dMolesPhase(i) = dTempVec(j)
        end do
!
        ! Determine what type of phase has the smallest Euclidean norm:
        if (iAssemblage(1) < 0) then
            ! A solution phase has the smallest Euclidean norm.
            iPhaseTypeOut = 1
        else
            ! A pure condensed phase has the smallest Euclidean norm.
            iPhaseTypeOut = 0
        end if
    end if IF_Euclid


!
    return
!
end subroutine ShuffleAssemblage
!
!
