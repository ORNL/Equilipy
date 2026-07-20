!> \brief Remove the active phase selected by the raw CEF phase-direction audit.
!!
!! \details Drops the phase whose CEF Newton phase-amount direction would make
!! it negative, then compacts the Lagrangian assemblage layout.  Stoichiometric
!! phases remain packed at the front and active solution phases remain packed
!! at the back, which preserves the assumptions used by GEMNewtonCEF.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    RemoveRawNegativeCEFPhase.f90
!> \brief   Remove a raw-negative CEF active phase from the Lagrangian set.
!> \author  S.Y. Kwon
!> \date    Jun. 25, 2026
!> \sa      AuditCEFRawPhaseDirection.f90
!> \sa      RunLagrangianGEM.f90
!> \sa      GEMNewtonCEF.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Added guarded removal of CEF and tiny solution phases that cross the active phase boundary.
    !
    !
! Purpose:
! ========
!
!> \details The purpose of this routine is to execute a clean active-set phase
!! removal when the CEF Newton direction identifies a phase amount boundary.
!! This keeps phase-amount positivity separate from site-fraction positivity:
!! the former changes the active set, while the latter remains the line
!! search/tracer-control problem.
!
!
! Required input variables:
! =========================
!
! iGEMRawNegativePhaseSlot     Lagrangian slot selected by AuditCEFRawPhaseDirection.
! iAssemblage                  Current active assemblage in Lagrangian layout.
! dMolesPhase                  Current active phase amounts.
! dMolFraction                 Current active solution endmember fractions.
!
!
! Output/updated variables:
! =========================
!
!> \param[out]  lRemoved       True when an active phase was removed and the
!!                             assemblage was compacted.
!
! iAssemblage                  Compacted active assemblage.
! dMolesPhase                  Compacted phase amounts.
! dMolesSpecies                Rebuilt species amounts for the reduced active set.
! nConPhases/nSolnPhases       Updated active phase counts.
! lSolnPhases                  Updated active solution-phase flags.
!
!
! Called subroutines/functions:
! =============================
!
! None.
!
!
! Primary callers:
! ================
!
! RunLagrangianGEM            Calls after AuditCEFRawPhaseDirection raises a raw
!                             negative CEF phase event.
!
!
! Numerical assumptions:
! ======================
!
! - The removed phase is already active and has a raw negative full-step phase
!   amount.  This routine does not remove phases merely because they are small.
! - Species amounts for retained solution phases are rebuilt from the retained
!   phase amount and current phase-local mole fractions.
!
!-------------------------------------------------------------------------------------------------------------



subroutine RemoveRawNegativeCEFPhase(lRemoved)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lRemoved

    call RemoveCEFActivePhaseAtSlot(iGEMRawNegativePhaseSlot, lRemoved)

    return

end subroutine RemoveRawNegativeCEFPhase



!> \brief Remove a CEF solution phase pinned at the phase-amount boundary.
!!
!! \details Selects an active CEF solution phase for removal only when both the
!! raw full Newton step and the accepted line-search state put that phase at a
!! numerically zero phase amount.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    RemoveRawNegativeCEFPhase.f90
!> \brief   Remove a CEF active phase pinned at the phase-amount boundary.
!> \author  S.Y. Kwon
!> \date    Jun. 25, 2026
!> \sa      RunLagrangianGEM.f90
!> \sa      GEMLineSearchCEF.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Removed CEF boundary phases only when both raw and accepted phase amounts are non-positive.
    !
    !
! Purpose:
! ========
!
!> \details The purpose of this routine is to avoid wasting Lagrangian
!! iterations on an inactive CEF solution phase whose phase amount is
!! non-positive.  The phase is removed only when line search has
!! failed to make useful progress and the raw Newton step and accepted
!! line-search state agree that the phase is at the boundary.
!
!
! Required input variables:
! =========================
!
! iAssemblage                 Current active assemblage in Lagrangian layout.
! dMolesPhase                 Current active phase amounts.
! dGEMLSRawPhaseMoles         Raw full-step phase amounts from the latest CEF line search.
! dGEMLSFinalPhaseMoles       Accepted phase amounts from the latest CEF line search.
!
!
! Output/updated variables:
! =========================
!
!> \param[out]  lRemoved      True when a boundary CEF phase was removed.
!
! iAssemblage                 Compacted active assemblage.
! dMolesPhase                 Compacted phase amounts.
! dMolesSpecies               Rebuilt species amounts for the reduced active set.
! nConPhases/nSolnPhases      Updated active phase counts.
! lSolnPhases                 Updated active solution-phase flags.
!
!
! Called subroutines/functions:
! =============================
!
! RemoveCEFActivePhaseAtSlot  Compacts the active assemblage after removal.
!
!
! Primary callers:
! ================
!
! RunLagrangianGEM            Calls after GEMLineSearchCEF.
!
!
! Numerical assumptions:
! ======================
!
! - This routine is only for normal minimization; postprocess Cp solves must
!   not change the active assemblage.
! - The active set must contain at least three phases.  This helper cleans an
!   overfull CEF active set; it must not collapse a two-phase assemblage into a
!   single phase.
! - The line search must either fail to find descent or exhaust the
!   backtracking attempts with less than 0.1% relative norm improvement.
! - The phase amount must be non-positive in both raw and accepted line-search
!   phase amounts.  A merely small finite phase is not enough.
!
!-------------------------------------------------------------------------------------------------------------



subroutine RemoveCEFBoundaryPhaseAfterLineSearch(lRemoved)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lRemoved

    integer :: i, iEntry, iFirstSolnSlot, iSelectedSlot
    real(8) :: dAcceptedMoles, dBestAcceptedMoles
    real(8) :: dRawMoles, dNormScale, dRelativeNormImprovement

    lRemoved = .FALSE.
    if (lPostProcess) return
    if (nConPhases + nSolnPhases <= 2) return
    if (.NOT.allocated(dGEMLSRawPhaseMoles)) return
    if (.NOT.allocated(dGEMLSFinalPhaseMoles)) return

    dNormScale = DMAX1(DABS(dGEMLineSearchInitialNorm), 1D0)
    dRelativeNormImprovement = (dGEMLineSearchInitialNorm - dGEMLineSearchBestNorm) / dNormScale
    if (iGEMLineSearchNoDescent <= 0) then
        if (iGEMLineSearchIterationCount < 12) return
        if (dRelativeNormImprovement > 1D-3) return
    end if

    iSelectedSlot = 0
    dBestAcceptedMoles = 0D0
    iFirstSolnSlot = nElements - nSolnPhases + 1
    do i = iFirstSolnSlot, nElements
        if ((i < 1).OR.(i > nElements)) cycle
        iEntry = iAssemblage(i)
        if (iEntry >= 0) cycle

        dRawMoles = dGEMLSRawPhaseMoles(i)
        dAcceptedMoles = dGEMLSFinalPhaseMoles(i)
        if ((dRawMoles <= 0D0).AND.(dAcceptedMoles <= 0D0).AND.&
            ((iSelectedSlot == 0).OR.(dAcceptedMoles < dBestAcceptedMoles))) then
            iSelectedSlot = i
            dBestAcceptedMoles = dAcceptedMoles
        end if
    end do

    if (iSelectedSlot == 0) return

    iGEMBoundaryRemovalSlot = iSelectedSlot
    iGEMBoundaryRemovalPhase = -iAssemblage(iSelectedSlot)

    call RemoveCEFActivePhaseAtSlot(iSelectedSlot, lRemoved)

    return

end subroutine RemoveCEFBoundaryPhaseAfterLineSearch



!> \brief Remove a non-positive active solution phase after Lagrangian stagnation.
!!
!! \details Drops a solution phase whose phase amount is non-positive, then
!! compacts the active assemblage while preserving retained phase compositions.
!! This is intentionally used only after the fixed-assemblage solver has
!! classified stagnation.  Finite-positive phase amounts remain active; active
!! set changes for small but positive phases belong to the phase-assemblage
!! logic, not to a numerical pruning tolerance.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    RemoveRawNegativeCEFPhase.f90
!> \brief   Remove a non-positive active solution phase from a stalled Lagrangian set.
!> \author  S.Y. Kwon
!> \date    Jun. 30, 2026
!> \sa      RunLagrangianGEM.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Restricted tiny-boundary removal to non-positive solution-phase remnants.
    !
!-------------------------------------------------------------------------------------------------------------



subroutine RemoveTinyBoundarySolutionPhase(lRemoved)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lRemoved

    integer :: i, iSlot, iEntry, iSolnPhase, iFirst, iLast, iDest
    integer :: iSelectedSlot, nConPhasesNew, nSolnPhasesNew
    integer, dimension(nElements) :: iAssemblageNew
    real(8) :: dSmallestMoles
    real(8), dimension(nElements) :: dMolesPhaseNew
    real(8), dimension(nSpecies) :: dMolesSpeciesNew

    lRemoved = .FALSE.
    if (lPostProcess) return
    if (nSolnPhases <= 0) return
    if (nConPhases + nSolnPhases <= 1) return

    iSelectedSlot = 0
    dSmallestMoles = 0D0
    do iSlot = nElements - nSolnPhases + 1, nElements
        if ((iSlot < 1).OR.(iSlot > nElements)) cycle
        iEntry = iAssemblage(iSlot)
        if (iEntry >= 0) cycle
        if ((dMolesPhase(iSlot) <= 0D0).AND.&
            ((iSelectedSlot == 0).OR.(dMolesPhase(iSlot) < dSmallestMoles))) then
            iSelectedSlot = iSlot
            dSmallestMoles = dMolesPhase(iSlot)
        end if
    end do

    if (iSelectedSlot == 0) return

    iGEMTinyBoundaryRemovalSlot = iSelectedSlot
    iGEMTinyBoundaryRemovalPhase = -iAssemblage(iSelectedSlot)

    iAssemblageNew = 0
    dMolesPhaseNew = 0D0
    dMolesSpeciesNew = 0D0
    nConPhasesNew = 0
    nSolnPhasesNew = 0

    do iSlot = 1, nElements
        if (iSlot == iSelectedSlot) cycle
        iEntry = iAssemblage(iSlot)
        if (iEntry <= 0) cycle

        nConPhasesNew = nConPhasesNew + 1
        iAssemblageNew(nConPhasesNew) = iEntry
        dMolesPhaseNew(nConPhasesNew) = dMolesPhase(iSlot)
        if ((iEntry > 0).AND.(iEntry <= nSpecies)) then
            dMolesSpeciesNew(iEntry) = dMolesPhaseNew(nConPhasesNew)
        end if
    end do

    do iSlot = nElements, 1, -1
        if (iSlot == iSelectedSlot) cycle
        iEntry = iAssemblage(iSlot)
        if (iEntry >= 0) cycle

        iSolnPhase = -iEntry
        if ((iSolnPhase < 1).OR.(iSolnPhase > nSolnPhasesSys)) cycle

        nSolnPhasesNew = nSolnPhasesNew + 1
        iDest = nElements - nSolnPhasesNew + 1
        iAssemblageNew(iDest) = iEntry
        dMolesPhaseNew(iDest) = dMolesPhase(iSlot)

        iFirst = nSpeciesPhase(iSolnPhase-1) + 1
        iLast = nSpeciesPhase(iSolnPhase)
        dMolesSpeciesNew(iFirst:iLast) = dMolesPhaseNew(iDest) * dMolFraction(iFirst:iLast)
    end do

    if (nConPhasesNew + nSolnPhasesNew <= 0) return

    iAssemblage = iAssemblageNew
    dMolesPhase = dMolesPhaseNew
    dMolesSpecies = dMolesSpeciesNew
    nConPhases = nConPhasesNew
    nSolnPhases = nSolnPhasesNew

    if (allocated(lSolnPhases)) then
        lSolnPhases = .FALSE.
        do iSlot = nElements - nSolnPhases + 1, nElements
            iEntry = iAssemblage(iSlot)
            if (iEntry < 0) lSolnPhases(-iEntry) = .TRUE.
        end do
    end if

    lRemoved = .TRUE.

    return

end subroutine RemoveTinyBoundarySolutionPhase



!> \brief Compact the Lagrangian active set after removing one phase slot.
!!
!! \details Rebuilds the stoichiometric and solution-phase portions of the
!! Lagrangian assemblage after the caller has selected an active phase slot for
!! removal.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    RemoveRawNegativeCEFPhase.f90
!> \brief   Compact the Lagrangian active set after phase removal.
!> \author  S.Y. Kwon
!> \date    Jun. 25, 2026
!> \sa      RemoveRawNegativeCEFPhase.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Compacted the active phase set only when rank and nonnegative mass balance remain valid.
    !
    !
! Purpose:
! ========
!
!> \details The purpose of this routine is to preserve the active-set layout
!! expected by GEMNewton and GEMNewtonCEF after one phase is removed.  Pure
!! condensed phases remain packed from the front of iAssemblage and active
!! solution phases remain packed from the back.
!
!
! Required input variables:
! =========================
!
!> \param[in]   iSlotIn       Active Lagrangian slot to remove.
!
! iAssemblage                 Current active assemblage in Lagrangian layout.
! dMolesPhase                 Current active phase amounts.
! dMolFraction                Current active solution endmember fractions.
!
!
! Output/updated variables:
! =========================
!
!> \param[out]  lRemoved      True when an active phase was removed and the
!!                            assemblage was compacted.
!
! iAssemblage                 Compacted active assemblage.
! dMolesPhase                 Compacted phase amounts.
! dMolesSpecies               Rebuilt species amounts for the reduced active set.
! nConPhases/nSolnPhases      Updated active phase counts.
! lSolnPhases                 Updated active solution-phase flags.
!
!
! Called subroutines/functions:
! =============================
!
! None.
!
!
! Primary callers:
! ================
!
! RemoveRawNegativeCEFPhase   Raw-negative CEF phase-removal path.
! RemoveCEFBoundaryPhaseAfterLineSearch
!                             Phase-boundary CEF phase-removal path.
!
!
! Numerical assumptions:
! ======================
!
! - The caller has already established that the selected phase is a plausible
!   boundary candidate.  This routine only removes it if the remaining active
!   set can still satisfy the bulk mass balance with nonnegative phase amounts.
!
!-------------------------------------------------------------------------------------------------------------



subroutine RemoveCEFActivePhaseAtSlot(iSlotIn, lRemoved)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in) :: iSlotIn
    logical, intent(out) :: lRemoved

    integer :: i, j, iSlot, iEntry, iSolnPhase, iFirst, iLast, INFO
    integer :: nConPhasesNew, nSolnPhasesNew, iDest
    integer, dimension(nElements) :: iAssemblageNew, IPIV
    real(8) :: dActivePhaseScale, dMolesTol
    real(8), dimension(nElements) :: dMolesPhaseNew, dMolesPhaseTrial
    real(8), dimension(nSpecies) :: dMolesSpeciesNew
    real(8), dimension(nElements,nElements) :: A

    lRemoved = .FALSE.
    iSlot = iSlotIn
    if (iSlot < 1 .OR. iSlot > nElements) return
    if (iAssemblage(iSlot) == 0) return

    iAssemblageNew = 0
    dMolesPhaseNew = 0D0
    dMolesSpeciesNew = 0D0
    dMolesPhaseTrial = 0D0
    A = 0D0
    IPIV = 0
    nConPhasesNew = 0
    nSolnPhasesNew = 0

    do i = 1, nElements
        if (i == iSlot) cycle
        iEntry = iAssemblage(i)
        if (iEntry <= 0) cycle

        nConPhasesNew = nConPhasesNew + 1
        iAssemblageNew(nConPhasesNew) = iEntry
        dMolesPhaseNew(nConPhasesNew) = dMolesPhase(i)
    end do

    do i = nElements, 1, -1
        if (i == iSlot) cycle
        iEntry = iAssemblage(i)
        if (iEntry >= 0) cycle

        iSolnPhase = -iEntry
        if (iSolnPhase < 1 .OR. iSolnPhase > nSolnPhasesSys) cycle

        nSolnPhasesNew = nSolnPhasesNew + 1
        iDest = nElements - nSolnPhasesNew + 1
        iAssemblageNew(iDest) = iEntry
        dMolesPhaseNew(iDest) = dMolesPhase(i)
    end do

    if (nConPhasesNew + nSolnPhasesNew < nElements) then
        iGEMBoundaryRankGuardUsed = 1
        iGEMBoundaryRankGuardSlot = iSlot
        iGEMBoundaryRankGuardPhase = 0
        if (iAssemblage(iSlot) < 0) iGEMBoundaryRankGuardPhase = -iAssemblage(iSlot)
        return
    end if

    do i = 1, nElements
        iEntry = iAssemblageNew(i)
        if (iEntry > 0) then
            do j = 1, nElements
                A(j,i) = dStoichSpecies(iEntry,j)
            end do
        else if (iEntry < 0) then
            iSolnPhase = -iEntry
            iFirst = nSpeciesPhase(iSolnPhase-1) + 1
            iLast = nSpeciesPhase(iSolnPhase)
            do j = 1, nElements
                A(j,i) = SUM(dStoichSpecies(iFirst:iLast,j) * dMolFraction(iFirst:iLast))
            end do
        end if
        dMolesPhaseTrial(i) = dMolesElement(i)
    end do

    INFO = 0
    call DGESV(nElements, 1, A, nElements, IPIV, dMolesPhaseTrial, nElements, INFO)
    if (INFO /= 0) then
        iGEMBoundaryRankGuardUsed = 1
        iGEMBoundaryRankGuardSlot = iSlot
        iGEMBoundaryRankGuardPhase = 0
        if (iAssemblage(iSlot) < 0) iGEMBoundaryRankGuardPhase = -iAssemblage(iSlot)
        return
    end if

    dActivePhaseScale = DMAX1(SUM(DABS(dMolesElement)), 1D0)
    dMolesTol = DMAX1(dTolerance(7), 1D-12 * dActivePhaseScale)
    if (MINVAL(dMolesPhaseTrial) < -dMolesTol) then
        iGEMBoundaryRankGuardUsed = 1
        iGEMBoundaryRankGuardSlot = iSlot
        iGEMBoundaryRankGuardPhase = 0
        if (iAssemblage(iSlot) < 0) iGEMBoundaryRankGuardPhase = -iAssemblage(iSlot)
        return
    end if

    do i = 1, nElements
        if (DABS(dMolesPhaseTrial(i)) < dMolesTol) dMolesPhaseTrial(i) = 0D0
    end do

    dMolesPhaseNew = dMolesPhaseTrial
    dMolesSpeciesNew = 0D0

    do i = 1, nElements
        iEntry = iAssemblageNew(i)
        if (iEntry > 0) then
            dMolesSpeciesNew(iEntry) = dMolesPhaseNew(i)
        else if (iEntry < 0) then
            iSolnPhase = -iEntry
            iFirst = nSpeciesPhase(iSolnPhase-1) + 1
            iLast = nSpeciesPhase(iSolnPhase)
            dMolesSpeciesNew(iFirst:iLast) = dMolesPhaseNew(i) * dMolFraction(iFirst:iLast)
        end if
    end do

    iAssemblage = iAssemblageNew
    dMolesPhase = dMolesPhaseNew
    dMolesSpecies = dMolesSpeciesNew
    nConPhases = nConPhasesNew
    nSolnPhases = nSolnPhasesNew

    if (allocated(lSolnPhases)) then
        lSolnPhases = .FALSE.
        do i = nElements - nSolnPhases + 1, nElements
            iEntry = iAssemblage(i)
            if (iEntry < 0) lSolnPhases(-iEntry) = .TRUE.
        end do
    end if

    lRemoved = .TRUE.

    return

end subroutine RemoveCEFActivePhaseAtSlot
