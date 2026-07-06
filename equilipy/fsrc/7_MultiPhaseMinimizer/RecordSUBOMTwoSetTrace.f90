!> \brief Record SUBOM two-set candidate and active-slot constitutions.
!!
!! \details Stores bounded diagnostics for the experimental same-parent SUBOM
!! two-set path.  The trace is passive and is used to compare candidate-pool
!! ordinal-2 rows against the slot-local constitutions evolved by Lagrangian GEM.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    RecordSUBOMTwoSetTrace.f90
!> \brief   Record SUBOM two-set candidate and active-slot constitutions.
!> \author  S.Y. Kwon
!> \date    Jul. 03, 2026
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   07/03/2026      S.Y. Kwon           Original passive trace for same-parent SUBOM candidate debugging.
!
!
! Purpose:
! ========
!
!> \details Record the parent-phase product-fraction vector and corresponding
!! site-fraction projection for important points in the two-set path:
!! registration, refreshed PEA candidate rows, Level2Lagrange handoff, and
!! solver-slot snapshots.  This routine does not affect minimization.
!
!
! Required input variables:
! =========================
!
! iStage      Trace stage code from ModuleGEMSolver.
! iPhase      Parent solution phase id.
! iOrdinal    Composition-set ordinal under the parent phase.
! iSlot       Active Lagrangian slot, or zero for candidate-pool events.
! nFraction  Number of entries in the parent product/endmember fraction vector.
! dFraction  Parent product/endmember fraction vector for the phase.
!
!
! Output/updated variables:
! =========================
!
! nSUBOMTwoSetTrace              Number of stored trace rows.
! iSUBOMTwoSetTrace*             Stage, iteration, parent, ordinal, and slot metadata.
! dSUBOMTwoSetTraceMol           Product-fraction trace rows.
! dSUBOMTwoSetTraceSite          Site-fraction projection trace rows.
! dSUBOMTwoSetTraceAmount        Active phase amount when a valid slot is supplied.
!
!
! Primary callers:
! ================
!
! RegisterSUBOMTwoSetCandidateRows    Candidate-pool registration and refresh.
! Level2Lagrange                      Direct Leveling-to-Lagrangian handoff.
! RecordGEMIterationDiagnostics       Solver-slot snapshots during Lagrangian GEM.
!
!-------------------------------------------------------------------------------------------------------------



subroutine RecordSUBOMTwoSetTrace(iStage, iSolnPhaseIn, iOrdinal, iSlot, nFraction, dFraction)
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in)              :: iStage, iSolnPhaseIn, iOrdinal, iSlot, nFraction
    real(8), dimension(*), intent(in):: dFraction

    integer :: iTrace, iFirst, iLast, nLocal, iLocal, iSpecies
    integer :: iPhaseID, iSub, iCon
    real(8) :: dSumLocal
    real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys) :: dSiteLocal

    if (.NOT.lSUBOMTwoSetCandidateEnabled) return
    if ((iSolnPhaseIn < 1).OR.(iSolnPhaseIn > nSolnPhasesSys)) return
    if (TRIM(cSolnPhaseType(iSolnPhaseIn)) /= 'SUBOM') return
    if (.NOT.allocated(iSUBOMTwoSetTraceStage)) return
    if (.NOT.allocated(dSUBOMTwoSetTraceMol)) return
    if (.NOT.allocated(dSUBOMTwoSetTraceSite)) return
    if (nSUBOMTwoSetTrace >= nSUBOMTwoSetTraceCapacity) return

    iFirst = nSpeciesPhase(iSolnPhaseIn-1) + 1
    iLast  = nSpeciesPhase(iSolnPhaseIn)
    nLocal = iLast - iFirst + 1
    if (nLocal <= 0) return

    iTrace = nSUBOMTwoSetTrace + 1
    nSUBOMTwoSetTrace = iTrace

    iSUBOMTwoSetTraceStage(iTrace) = iStage
    iSUBOMTwoSetTraceIterPEA(iTrace) = iterPEA
    iSUBOMTwoSetTraceIterGlobal(iTrace) = iterGlobal
    iSUBOMTwoSetTracePhase(iTrace) = iSolnPhaseIn
    iSUBOMTwoSetTraceOrdinal(iTrace) = iOrdinal
    iSUBOMTwoSetTraceSlot(iTrace) = iSlot

    dSUBOMTwoSetTraceAmount(iTrace) = 0D0
    if ((iSlot >= 1).AND.(iSlot <= SIZE(dMolesPhase))) then
        dSUBOMTwoSetTraceAmount(iTrace) = dMolesPhase(iSlot)
    end if

    dSUBOMTwoSetTraceMol(iTrace,:) = 0D0
    do iLocal = 1, MIN(nLocal, nFraction)
        iSpecies = iFirst + iLocal - 1
        if ((iSpecies >= 1).AND.(iSpecies <= SIZE(dSUBOMTwoSetTraceMol,2))) then
            dSUBOMTwoSetTraceMol(iTrace,iSpecies) = DMAX1(dFraction(iLocal), 0D0)
        end if
    end do

    dSiteLocal = 0D0
    iPhaseID = iPhaseSublattice(iSolnPhaseIn)
    if (iPhaseID > 0) then
        do iLocal = 1, MIN(nLocal, nFraction)
            do iSub = 1, nSublatticePhase(iPhaseID)
                iCon = iConstituentSublattice(iPhaseID,iSub,iLocal)
                if ((iCon >= 1).AND.(iCon <= nMaxConstituentSys)) then
                    dSiteLocal(iSub,iCon) = dSiteLocal(iSub,iCon) + DMAX1(dFraction(iLocal), 0D0)
                end if
            end do
        end do
        do iSub = 1, nSublatticePhase(iPhaseID)
            dSumLocal = SUM(dSiteLocal(iSub,:))
            if (dSumLocal > 1D-300) dSiteLocal(iSub,:) = dSiteLocal(iSub,:) / dSumLocal
        end do
    end if

    if ((iStage == SUBOM_TRACE_SOLVER_SLOT).OR.(iStage == SUBOM_TRACE_HANDOFF)) then
        if (allocated(dActiveSlotSiteFraction)) then
            if ((iSlot >= 1).AND.(iSlot <= SIZE(dActiveSlotSiteFraction,1))) then
                if (SUM(dActiveSlotSiteFraction(iSlot,:,:)) > 0D0) then
                    dSiteLocal = dActiveSlotSiteFraction(iSlot,:,:)
                end if
            end if
        end if
    end if

    dSUBOMTwoSetTraceSite(iTrace,:,:) = dSiteLocal

    return

end subroutine RecordSUBOMTwoSetTrace
