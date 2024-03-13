!
subroutine ShuffleDummySpecies
!
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to shuffle dummy species to the front of
    ! the integer vector iAssemblage.  The motivation for doing this is that
    ! dummy species will be considered first to be exchanged for a new phase.
    ! Note that a unique solution is not guaranteed in Leveling in situations
    ! where zero moles is assigned to one of the dummy phases (degrees of freedom
    ! is non-zero).  Therefore, it is possible for G_sys to be at a minimum
    ! but the phase assemblage contains a dummy phase.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! iPhase        An integer vector dimensioned to the total number of
    !                species in the system.  Each coefficient corresponds to
    !                the phase type corresponding to this species.
    !                iPhase < 0: A dummy species.
    !                iPhase = 0: A pure condensed phase.
    !                iPhase > 0: A solution species, where iPhase corresponds
    !                            to the absolute solution phase index.
    !
    !---------------------------------------------------------------------------
    USE ModuleThermo, ONLY: iAssemblage, nElements, iPhaseLevel, iShuffled
!
    implicit none
!
    integer                     :: i, j, k, nDummy
    integer,dimension(nElements):: iTempVec, iTempVecB
!
!
    ! Initialize variables:
    nDummy   = 0
    iTempVec = 0
    iTempVecB= 0
    ! Count the number of dummy species:
    do i = 1, nElements
        if (iPhaseLevel(iAssemblage(i)) < 0) nDummy = nDummy + 1
    end do
!
    ! Only proceed if a dummy species is currently in the phase assemblage:
    IF_Proceed: if (nDummy > 0) then
!
        ! Shuffle vector:
        j = 0
        k = nDummy
        do i = 1, nElements
            if (iPhaseLevel(iAssemblage(i)) < 0) then
                j = j + 1
                iTempVec(j) = iAssemblage(i)
                iTempVecB(j)= iShuffled(i)
            else
                k = k + 1
                iTempVec(k) = iAssemblage(i)
                iTempVecB(k)= iShuffled(i)
            end if
        end do
!
        ! Update the iAssemblage vector:
        iAssemblage = iTempVec
        iShuffled   = iTempVecB
!
    end if IF_Proceed
!
    return
!
end subroutine ShuffleDummySpecies
!
!
