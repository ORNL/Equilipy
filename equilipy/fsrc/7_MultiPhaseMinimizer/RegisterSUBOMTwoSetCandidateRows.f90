!> \brief Register same-parent SUBOM composition-set candidates.
!!
!! \details Adds optional second PEA pseudo-compound rows for ordered SUBOM
!! parents when the passive ordering-mode diagnostic reports an instability.
!! The rows are Leveling evidence: Leveling may select them directly, and
!! Level2Lagrange can then translate them into independent active slots.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    RegisterSUBOMTwoSetCandidateRows.f90
!> \brief   Register switch-gated same-parent SUBOM composition-set candidates.
!> \author  S.Y. Kwon
!> \date    Jul. 03, 2026
!> \sa      CompInitMinSolnPoint.f90
!> \sa      CompMinSolnPoint.f90
!> \sa      SetLevelingSolutionCandidateRow.f90
!> \sa      ProjectOrderDisorderCompanionFraction.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Registered structurally distinct order/disorder companion constitutions as typed parent composition sets through central grid row layout.
    !
    !
! Purpose:
! ========
!
!> \details This routine creates the second composition-set candidate row for
!! an ordered SUBOM parent from its disordered companion phase.  The row is
!! written as an ordered-parent pseudo-compound with identity ordinal 2, so it
!! can coexist with the ordinary ordered-minimum row during Leveling.
!
!
! Required input variables:
! =========================
!
! lSUBOMTwoSetCandidateEnabled Switch enabling the optional helper-only path.
! lODPartitionUnifiedActive    Per-system activation for parser-declared physical DIS_PART pairs.
! iODCompanionPhase            Active disordered-companion map for ordered SUBOM parents.
! iSubMinCandidateStatusSoln   Subminimization status for each solution candidate row.
! dMolFraction                 Current phase-local fractions after solution candidate refresh.
! dElementPotential            Current PEA elemental potential plane.
!
!
! Output/updated variables:
! =========================
!
! dLevelingCompositionSpecies   Extra row composition for the second parent composition set.
! dStoichSpeciesLevel           Extra row stoichiometry for Leveling mass balance.
! dLevelingChemicalPotential    Extra row Leveling Gibbs potential.
! iPhaseLevel                   Extra row phase id, equal to the ordered SUBOM parent.
! iLevelCandidate*              Candidate-pool identity metadata and stored constitution.
!
!
! Called subroutines/functions:
! =============================
!
! CompOrderingModeSUBOM          Supplies the c1 ordering-instability gate.
! CompExcessGibbsEnergy          Evaluates the ordered parent at the projected disordered-manifold start.
! ProjectOrderDisorderCompanionFraction
!                                Maps a companion constitution into the ordered parent.
! SetLevelingSolutionCandidateRow Writes the actual PEA pseudo-compound row.
!
!
! Primary callers:
! ================
!
! CompInitMinSolnPoint           Registers initial PEA solution minima.
! CompMinSolnPoint               Refreshes PEA solution minima after each PEA iteration.
!
!
! Numerical assumptions:
! ======================
!
! - Physical DIS_PART pairs use the canonical path whenever the active system
!   contains such a pair.  The unrelated optional helper-only path remains
!   controlled by lSUBOMTwoSetCandidateEnabled.
! - The extra row index is nSpecies + nSolnPhasesSys + parent phase id.
! - Canonical partition candidates require typed ordered-parent evidence and an
!   accepted, exactly distinct physical-companion constitution.  A companion
!   ordering-instability fact remains typed; it does not erase the physical
!   composition-set candidate before the coupled active-set solve.
! - The optional helper-only path retains its historical instability and
!   duplicate-distance checks.
! - The first ordinal-2 row is created from the disordered companion.  Once an
!   ordinal-2 active slot exists, refreshes use that slot-local constitution so
!   the candidate pool tracks the accepted PEA-Lagrangian state.
! - A canonical ordinal-2 row is emitted only when no separately declared
!   standalone disordered phase exists.  Otherwise the standalone and the
!   genuinely ordered parent remain independent physical candidates.
! - The ordered parent is evaluated temporarily at the projected fraction and
!   then restored to its original selected constitution.
!
!-------------------------------------------------------------------------------------------------------------



subroutine RegisterSUBOMTwoSetCandidateRows
    USE ModuleThermo
    USE ModuleGEMSolver
    USE GridDiscovery, ONLY: GridTangentRowIndex, GridPartitionRowIndex

    implicit none

    integer :: iOrderedPhase, iHelperPhase, iLevelRow, iFirst, iLast, nLocal
    integer :: iTraceStage, iActiveSlotSource, iProjectionStatus
    integer :: iCand
    logical :: lHasInstability, lHelperValid, lProjected, lDistinct, lCandidateValid
    logical :: lFromActiveSlot, lSourceValid, lCanonicalPair, lStandaloneExists
    logical :: lOrderedCompanionMayEnter
    real(8) :: dDistance
    real(8), dimension(:), allocatable :: dProjectedFraction, dOrderedFractionSave
    logical :: OrderDisorderPhaseIsEligible, IsOrderDisorderHelperAliasPhase

    if ((.NOT.lSUBOMTwoSetCandidateEnabled).AND.(.NOT.lODPartitionUnifiedActive)) return
    if (nSpeciesLevel < nSpecies + 2*nSolnPhasesSys) return
    if (.NOT.allocated(iLevelCandidateFromLevel)) return
    if (.NOT.allocated(iSubMinCandidateStatusSoln)) return

    do iOrderedPhase = 1, nSolnPhasesSys
        if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') cycle
        if (.NOT.OrderDisorderPhaseIsEligible(iOrderedPhase)) cycle
        if (.NOT.lSUBOMTwoSetCandidateEnabled) then
            if (.NOT.allocated(iODTopologyClass)) cycle
            if ((iODTopologyClass(iOrderedPhase) < OD_TOPOLOGY_HELPER_STANDALONE).OR.&
                (iODTopologyClass(iOrderedPhase) > OD_TOPOLOGY_HELPER_ONLY)) cycle
        end if

        iHelperPhase = TwoSetCompanionPhase(iOrderedPhase)
        if (iHelperPhase <= 0) cycle
        if (iHelperPhase == iOrderedPhase) cycle
        lCanonicalPair = lODPartitionUnifiedActive.AND.allocated(iODCompanionPhase).AND.&
            allocated(iODTopologyClass)
        if (lCanonicalPair) then
            lCanonicalPair = iOrderedPhase <= SIZE(iODCompanionPhase)
            if (lCanonicalPair) lCanonicalPair = &
                iODCompanionPhase(iOrderedPhase) == iHelperPhase
            if (lCanonicalPair) lCanonicalPair = &
                (iODTopologyClass(iOrderedPhase) >= OD_TOPOLOGY_HELPER_STANDALONE).AND.&
                (iODTopologyClass(iOrderedPhase) <= OD_TOPOLOGY_HELPER_ONLY)
        end if
        lStandaloneExists = StandaloneDisorderedPhase(iOrderedPhase) > 0

        iLevelRow = GridPartitionRowIndex(iOrderedPhase)
        if ((iLevelRow < 1).OR.(iLevelRow > nSpeciesLevel)) cycle

        iFirst = nSpeciesPhase(iOrderedPhase-1) + 1
        iLast  = nSpeciesPhase(iOrderedPhase)
        nLocal = iLast - iFirst + 1
        if (nLocal <= 0) cycle

        allocate(dProjectedFraction(nLocal), dOrderedFractionSave(nLocal))
        dProjectedFraction = 0D0
        dOrderedFractionSave = dMolFraction(iFirst:iLast)

        lHasInstability = .FALSE.
        iActiveSlotSource = 0
        call GetStoredOrdinalSUBOMFraction(iOrderedPhase, 2, dProjectedFraction, lFromActiveSlot)
        if (.NOT.lFromActiveSlot) then
            call GetActiveOrdinalSUBOMFraction(iOrderedPhase, 2, dProjectedFraction, &
                iActiveSlotSource, lFromActiveSlot)
        end if
        lProjected = lFromActiveSlot
        lHelperValid = CandidatePhaseHasAcceptedSubminResult(iHelperPhase)
        if (.NOT.lFromActiveSlot) then
            call ProjectOrderDisorderCompanionFraction(iHelperPhase, iOrderedPhase, &
                nLocal, dProjectedFraction, lProjected, iProjectionStatus)
        end if

        dDistance = SUM(DABS(dProjectedFraction - dOrderedFractionSave)) / DBLE(MAX(1,nLocal))
        if (lCanonicalPair) then
            ! A parser-declared companion is a second representation of the
            ! same parent model.  Typed branch classification, not a numeric
            ! constitution-distance cutoff, decides whether it supplies a row.
            lDistinct = ANY(dProjectedFraction /= dOrderedFractionSave)
        else
            lDistinct = dDistance > 1D-8
        end if
        lSourceValid = lFromActiveSlot.OR.(lHelperValid.AND.lProjected)
        if ((.NOT.lCanonicalPair).AND.lSourceValid.AND.lProjected.AND.lDistinct) then
            lHasInstability = SUBOMCandidateHasOrderingInstability(iOrderedPhase)
        end if
        if (lCanonicalPair) then
            lOrderedCompanionMayEnter = (.NOT.lStandaloneExists).OR.lGridFrontEndActive
            lCandidateValid = lOrderedCompanionMayEnter.AND.&
                lSourceValid.AND.lProjected.AND.lDistinct.AND.&
                ((iODCandidateClass(iOrderedPhase) == OD_CANDIDATE_ORDERED).OR.&
                (iODCandidateClass(iOrderedPhase) == &
                    OD_CANDIDATE_ORDERED_COMPANION_UNSTABLE))
        else
            lOrderedCompanionMayEnter = .FALSE.
            lCandidateValid = lHasInstability.AND.lSourceValid.AND.lProjected.AND.lDistinct.AND.&
                TwoSetCompanionIsHelperAliasOnly(iOrderedPhase, iHelperPhase)
        end if

        if (lCandidateValid) then
            dMolFraction(iFirst:iLast) = dProjectedFraction
            dGibbsSolnPhase(iOrderedPhase) = 0D0
            call CompExcessGibbsEnergy(iOrderedPhase)
        end if

        call SetLevelingSolutionCandidateRow(iLevelRow, iOrderedPhase, &
            dProjectedFraction, lCandidateValid)
        call MarkTwoSetCandidateIdentity(iLevelRow, iOrderedPhase, iHelperPhase, lCanonicalPair)
        if (lCandidateValid) then
            iTraceStage = SUBOM_TRACE_REFRESH
            if (iterPEA <= 0) iTraceStage = SUBOM_TRACE_REGISTER
            call RecordSUBOMTwoSetTrace(iTraceStage, iOrderedPhase, 2, &
                iActiveSlotSource, nLocal, dProjectedFraction)
        end if
        if (lCandidateValid) then
            dMolFraction(iFirst:iLast) = dOrderedFractionSave
            dGibbsSolnPhase(iOrderedPhase) = 0D0
            call CompExcessGibbsEnergy(iOrderedPhase)
        end if

        deallocate(dProjectedFraction, dOrderedFractionSave)
    end do

    return

contains

    integer function TwoSetCompanionPhase(iOrderedPhaseIn)
        integer, intent(in) :: iOrderedPhaseIn

        TwoSetCompanionPhase = 0
        if ((iOrderedPhaseIn < 1).OR.(iOrderedPhaseIn > nSolnPhasesSys)) return

        if (allocated(iODCompanionPhase)) then
            if (iOrderedPhaseIn <= SIZE(iODCompanionPhase)) then
                TwoSetCompanionPhase = iODCompanionPhase(iOrderedPhaseIn)
            end if
        end if
        if (TwoSetCompanionPhase <= 0) then
            if (allocated(iDisorderedPhase)) then
                if (iOrderedPhaseIn <= SIZE(iDisorderedPhase)) then
                    TwoSetCompanionPhase = iDisorderedPhase(iOrderedPhaseIn)
                end if
            end if
        end if
        if ((TwoSetCompanionPhase < 1).OR.(TwoSetCompanionPhase > nSolnPhasesSys)) then
            TwoSetCompanionPhase = 0
        end if

        return
    end function TwoSetCompanionPhase


    integer function StandaloneDisorderedPhase(iOrderedPhaseIn)
        integer, intent(in) :: iOrderedPhaseIn

        StandaloneDisorderedPhase = 0
        if (.NOT.allocated(iODStandalonePhase)) return
        if ((iOrderedPhaseIn < 1).OR.(iOrderedPhaseIn > SIZE(iODStandalonePhase))) return
        if (iODStandalonePhase(iOrderedPhaseIn) == iOrderedPhaseIn) return
        StandaloneDisorderedPhase = iODStandalonePhase(iOrderedPhaseIn)

        return
    end function StandaloneDisorderedPhase


    logical function CandidatePhaseHasAcceptedSubminResult(iSolnPhaseIn)
        integer, intent(in) :: iSolnPhaseIn

        CandidatePhaseHasAcceptedSubminResult = .FALSE.
        if ((iSolnPhaseIn < 1).OR.(iSolnPhaseIn > SIZE(iSubMinCandidateStatusSoln))) return
        CandidatePhaseHasAcceptedSubminResult = &
            (iSubMinCandidateStatusSoln(iSolnPhaseIn) == SUBMIN_CANDIDATE_CONVERGED).OR.&
            (iSubMinCandidateStatusSoln(iSolnPhaseIn) == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)

        return
    end function CandidatePhaseHasAcceptedSubminResult


    subroutine GetActiveOrdinalSUBOMFraction(iOrderedPhaseIn, iOrdinalIn, dFractionOut, &
        iSlotOut, lFoundOut)
        integer, intent(in) :: iOrderedPhaseIn, iOrdinalIn
        real(8), dimension(:), intent(out) :: dFractionOut
        integer, intent(out) :: iSlotOut
        logical, intent(out) :: lFoundOut

        integer :: iSlotLocal, iFirstLocal, iLastLocal, nLocal
        real(8) :: dNormLocal

        dFractionOut = 0D0
        iSlotOut = 0
        lFoundOut = .FALSE.

        if (.NOT.allocated(iActiveSlotThermoPhase)) return
        if (.NOT.allocated(iActiveSlotIdentityOrdinal)) return
        if (.NOT.allocated(dActiveSlotMolFraction)) return
        if ((iOrderedPhaseIn < 1).OR.(iOrderedPhaseIn > nSolnPhasesSys)) return

        iFirstLocal = nSpeciesPhase(iOrderedPhaseIn-1) + 1
        iLastLocal  = nSpeciesPhase(iOrderedPhaseIn)
        nLocal = iLastLocal - iFirstLocal + 1
        if (nLocal <= 0) return
        if (SIZE(dFractionOut) /= nLocal) return
        if ((iFirstLocal < 1).OR.(iLastLocal > SIZE(dActiveSlotMolFraction,2))) return

        do iSlotLocal = 1, MIN(nElements, SIZE(iActiveSlotThermoPhase))
            if (iActiveSlotThermoPhase(iSlotLocal) /= iOrderedPhaseIn) cycle
            if (iActiveSlotIdentityOrdinal(iSlotLocal) /= iOrdinalIn) cycle
            dFractionOut = DMAX1(dActiveSlotMolFraction(iSlotLocal,iFirstLocal:iLastLocal), 0D0)
            dNormLocal = SUM(dFractionOut)
            if (dNormLocal <= 1D-300) cycle
            dFractionOut = dFractionOut / dNormLocal
            iSlotOut = iSlotLocal
            lFoundOut = .TRUE.
            exit
        end do

        return
    end subroutine GetActiveOrdinalSUBOMFraction


    subroutine GetStoredOrdinalSUBOMFraction(iOrderedPhaseIn, iOrdinalIn, dFractionOut, lFoundOut)
        integer, intent(in) :: iOrderedPhaseIn, iOrdinalIn
        real(8), dimension(:), intent(out) :: dFractionOut
        logical, intent(out) :: lFoundOut

        integer :: iFirstLocal, iLastLocal, nLocal
        real(8) :: dNormLocal

        dFractionOut = 0D0
        lFoundOut = .FALSE.

        if (.NOT.allocated(iSUBOMTwoSetStoredPhase)) return
        if (.NOT.allocated(iSUBOMTwoSetStoredOrdinal)) return
        if (.NOT.allocated(dSUBOMTwoSetStoredMol)) return
        if ((iOrderedPhaseIn < 1).OR.(iOrderedPhaseIn > SIZE(iSUBOMTwoSetStoredPhase))) return
        if (iSUBOMTwoSetStoredPhase(iOrderedPhaseIn) /= iOrderedPhaseIn) return
        if (iSUBOMTwoSetStoredOrdinal(iOrderedPhaseIn) /= iOrdinalIn) return

        iFirstLocal = nSpeciesPhase(iOrderedPhaseIn-1) + 1
        iLastLocal  = nSpeciesPhase(iOrderedPhaseIn)
        nLocal = iLastLocal - iFirstLocal + 1
        if (nLocal <= 0) return
        if (SIZE(dFractionOut) /= nLocal) return
        if ((iFirstLocal < 1).OR.(iLastLocal > SIZE(dSUBOMTwoSetStoredMol,2))) return

        dFractionOut = DMAX1(dSUBOMTwoSetStoredMol(iOrderedPhaseIn,iFirstLocal:iLastLocal), 0D0)
        dNormLocal = SUM(dFractionOut)
        if (dNormLocal <= 1D-300) return
        dFractionOut = dFractionOut / dNormLocal
        lFoundOut = .TRUE.

        return
    end subroutine GetStoredOrdinalSUBOMFraction


    logical function SUBOMCandidateHasOrderingInstability(iOrderedPhaseIn)
        integer, intent(in) :: iOrderedPhaseIn
        integer             :: nModeCapacity, nModeLocal, iInfoLocal
        real(8)             :: dEigenLocal
        real(8), dimension(:), allocatable :: dOrderingGradient
        real(8), dimension(:,:), allocatable :: dOrderingHessian

        SUBOMCandidateHasOrderingInstability = .FALSE.
        if (allocated(iSUBOMOrderingGateIterPEA).AND.&
            allocated(lSUBOMOrderingGateUnstable)) then
            if ((iOrderedPhaseIn >= 1).AND.&
                (iOrderedPhaseIn <= SIZE(iSUBOMOrderingGateIterPEA))) then
                if (iSUBOMOrderingGateIterPEA(iOrderedPhaseIn) == iterPEA) then
                    SUBOMCandidateHasOrderingInstability = &
                        lSUBOMOrderingGateUnstable(iOrderedPhaseIn)
                    return
                end if
            end if
        end if

        nModeCapacity = MAX(1, nMaxSublatticeSys*nMaxConstituentSys)
        allocate(dOrderingGradient(nModeCapacity), dOrderingHessian(nModeCapacity,nModeCapacity))
        call CompOrderingModeSUBOM(iOrderedPhaseIn, nModeCapacity, dOrderingHessian, &
            dOrderingGradient, dEigenLocal, nModeLocal, iInfoLocal)
        nSUBOMOrderingGateEvaluated = nSUBOMOrderingGateEvaluated + 1
        if ((iInfoLocal == 0).AND.(nModeLocal > 0).AND.(dEigenLocal < -1D-8)) then
            SUBOMCandidateHasOrderingInstability = .TRUE.
        end if
        if (allocated(iSUBOMOrderingGateIterPEA).AND.&
            allocated(iSUBOMOrderingGateModeCount).AND.&
            allocated(iSUBOMOrderingGateInfo).AND.&
            allocated(dSUBOMOrderingGateEigenMin).AND.&
            allocated(lSUBOMOrderingGateUnstable)) then
            if ((iOrderedPhaseIn >= 1).AND.&
                (iOrderedPhaseIn <= SIZE(iSUBOMOrderingGateIterPEA))) then
                iSUBOMOrderingGateIterPEA(iOrderedPhaseIn) = iterPEA
                iSUBOMOrderingGateModeCount(iOrderedPhaseIn) = nModeLocal
                iSUBOMOrderingGateInfo(iOrderedPhaseIn) = iInfoLocal
                dSUBOMOrderingGateEigenMin(iOrderedPhaseIn) = dEigenLocal
                lSUBOMOrderingGateUnstable(iOrderedPhaseIn) = &
                    SUBOMCandidateHasOrderingInstability
            end if
        end if
        deallocate(dOrderingGradient, dOrderingHessian)

        return
    end function SUBOMCandidateHasOrderingInstability


    subroutine SuppressHelperOnlyCompanionCandidate(iHelperPhaseIn)
        integer, intent(in) :: iHelperPhaseIn

        integer :: iHelperRow

        if (.NOT.IsOrderDisorderHelperAliasPhase(iHelperPhaseIn)) return
        iHelperRow = GridTangentRowIndex(iHelperPhaseIn)
        if ((iHelperRow < 1).OR.(iHelperRow > nSpeciesLevel)) return

        dChemicalPotential(iHelperRow) = 5D9
        dLevelingChemicalPotential(iHelperRow) = 5D9

        return
    end subroutine SuppressHelperOnlyCompanionCandidate


    logical function TwoSetCompanionIsHelperAliasOnly(iOrderedPhaseIn, iCompanionPhaseIn)
        integer, intent(in) :: iOrderedPhaseIn, iCompanionPhaseIn
        integer             :: iMappedCompanion

        TwoSetCompanionIsHelperAliasOnly = .FALSE.
        if (.NOT.allocated(iODCompanionPhase)) return
        if ((iOrderedPhaseIn < 1).OR.(iOrderedPhaseIn > SIZE(iODCompanionPhase))) return
        if ((iCompanionPhaseIn < 1).OR.(iCompanionPhaseIn > nSolnPhasesSys)) return

        iMappedCompanion = iODCompanionPhase(iOrderedPhaseIn)
        if (iMappedCompanion /= iCompanionPhaseIn) return
        TwoSetCompanionIsHelperAliasOnly = &
            IsOrderDisorderHelperAliasPhase(iCompanionPhaseIn)

        return
    end function TwoSetCompanionIsHelperAliasOnly


    subroutine MarkTwoSetCandidateIdentity(iLevelRowIn, iOrderedPhaseIn, &
        iCompanionPhaseIn, lCanonicalPairIn)
        integer, intent(in) :: iLevelRowIn, iOrderedPhaseIn, iCompanionPhaseIn
        logical, intent(in) :: lCanonicalPairIn

        iCand = 0
        if (allocated(iLevelCandidateFromLevel)) then
            if ((iLevelRowIn >= 1).AND.(iLevelRowIn <= SIZE(iLevelCandidateFromLevel))) then
                iCand = iLevelCandidateFromLevel(iLevelRowIn)
            end if
        end if
        if ((iCand < 1).OR.(iCand > nLevelCandidate)) return

        iLevelCandidateParentPhase(iCand) = iOrderedPhaseIn
        iLevelCandidateDisplayPhase(iCand) = iOrderedPhaseIn
        if (lCanonicalPairIn) iLevelCandidateDisplayPhase(iCand) = iCompanionPhaseIn
        iLevelCandidateIdentityOrdinal(iCand) = 2
        iLevelCandidateSource(iCand) = 2

        return
    end subroutine MarkTwoSetCandidateIdentity


end subroutine RegisterSUBOMTwoSetCandidateRows
