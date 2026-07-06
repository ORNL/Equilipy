!> \brief Register switch-gated same-parent SUBOM composition-set candidates.
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
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   07/03/2026      S.Y. Kwon           Original switch-gated PEA candidate-pool identity persistence helper.
!   07/03/2026      S.Y. Kwon           Suppressed the companion Leveling row when a valid parent ordinal-2
!                                       row is available, so the helper is not an independent competitor.
!   07/03/2026      S.Y. Kwon           Refreshed existing ordinal-2 rows from active slot-local constitutions
!                                       instead of reseeding them from the disordered-manifold start.
!   07/03/2026      S.Y. Kwon           Removed companion suppression so same-parent ordinal rows compete with,
!                                       but never replace, physical disordered phase candidates.
!   07/03/2026      S.Y. Kwon           Cached the ordering-instability gate per PEA iteration and skipped it
!                                       when the candidate source cannot become Leveling evidence.
!   07/03/2026      S.Y. Kwon           Allowed ordinal-2 rows to own known helper-only DIS_PART aliases while
!                                       preserving compete-only behavior for physical disordered phases.
!   07/03/2026      S.Y. Kwon           Kept the companion helper row available when adding ordinal-2 evidence.
!   07/03/2026      S.Y. Kwon           Skipped ordinal-row evidence when the companion helper row is valid.
!   07/03/2026      S.Y. Kwon           Limited ordinal-row evidence to nonphysical helper aliases.
!   07/03/2026      S.Y. Kwon           Disabled ordinal-row evidence when a physical DIS_PART phase is present.
!   07/03/2026      S.Y. Kwon           Replaced name-scanned physical-companion detection with the central
!                                       order/disorder companion map and shared helper-alias classifier.
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
! lSUBOMTwoSetCandidateEnabled Switch enabling the experimental two-set path.
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
! - The second row is created only under lSUBOMTwoSetCandidateEnabled.
! - The extra row index is nSpecies + nSolnPhasesSys + parent phase id.
! - Candidate validity requires an ordering-mode negative eigenvalue, a valid
!   constitution source, and a parent-constitution distance above the duplicate
!   tolerance.
! - The first ordinal-2 row is created from the disordered companion.  Once an
!   ordinal-2 active slot exists, refreshes use that slot-local constitution so
!   the candidate pool tracks the accepted PEA-Lagrangian state.
! - A valid parent ordinal-2 row is helper-only additional Leveling evidence.
!   It must not remove or replace a physical disordered companion row.  When a
!   physical companion is available, ordinal-2 evidence is disabled and the
!   ordinary disordered row owns that branch.
! - The ordered parent is evaluated temporarily at the projected fraction and
!   then restored to its original selected constitution.
!
!-------------------------------------------------------------------------------------------------------------



subroutine RegisterSUBOMTwoSetCandidateRows
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: iOrderedPhase, iHelperPhase, iLevelRow, iFirst, iLast, nLocal
    integer :: iTraceStage, iActiveSlotSource
    integer :: iCand
    logical :: lHasInstability, lHelperValid, lProjected, lDistinct, lCandidateValid
    logical :: lFromActiveSlot, lSourceValid
    real(8) :: dDistance
    real(8), dimension(:), allocatable :: dProjectedFraction, dOrderedFractionSave
    logical :: OrderDisorderPhaseIsEligible, IsOrderDisorderHelperAliasPhase

    if (.NOT.lSUBOMTwoSetCandidateEnabled) return
    if (nSpeciesLevel < nSpecies + 2*nSolnPhasesSys) return
    if (.NOT.allocated(iLevelCandidateFromLevel)) return
    if (.NOT.allocated(iSubMinCandidateStatusSoln)) return

    do iOrderedPhase = 1, nSolnPhasesSys
        if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') cycle
        if (.NOT.OrderDisorderPhaseIsEligible(iOrderedPhase)) cycle

        iHelperPhase = TwoSetCompanionPhase(iOrderedPhase)
        if (iHelperPhase <= 0) cycle
        if (iHelperPhase == iOrderedPhase) cycle

        iLevelRow = nSpecies + nSolnPhasesSys + iOrderedPhase
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
            call ProjectCompanionToOrderedSUBOMFraction(iHelperPhase, iOrderedPhase, &
                dProjectedFraction, lProjected)
        end if

        dDistance = SUM(DABS(dProjectedFraction - dOrderedFractionSave)) / DBLE(MAX(1,nLocal))
        lDistinct = dDistance > 1D-8
        lSourceValid = lFromActiveSlot.OR.(lHelperValid.AND.lProjected)
        if (lSourceValid.AND.lProjected.AND.lDistinct) then
            lHasInstability = SUBOMCandidateHasOrderingInstability(iOrderedPhase)
        end if
        lCandidateValid = lHasInstability.AND.lSourceValid.AND.lProjected.AND.lDistinct.AND.&
            TwoSetCompanionIsHelperAliasOnly(iOrderedPhase, iHelperPhase)

        if (lCandidateValid) then
            dMolFraction(iFirst:iLast) = dProjectedFraction
            dGibbsSolnPhase(iOrderedPhase) = 0D0
            call CompExcessGibbsEnergy(iOrderedPhase)
        end if

        call SetLevelingSolutionCandidateRow(iLevelRow, iOrderedPhase, &
            dProjectedFraction, lCandidateValid)
        call MarkTwoSetCandidateIdentity(iLevelRow, iOrderedPhase)
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
        iHelperRow = nSpecies + iHelperPhaseIn
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


    subroutine ProjectCompanionToOrderedSUBOMFraction(iHelperPhaseIn, iOrderedPhaseIn, &
        dProjectedFractionOut, lProjectedOut)
        integer, intent(in) :: iHelperPhaseIn, iOrderedPhaseIn
        real(8), dimension(:), intent(out) :: dProjectedFractionOut
        logical, intent(out) :: lProjectedOut

        integer :: i, j, s, m, iHelperFirst, iHelperLast, iOrderedFirst, iOrderedLast
        integer :: iOrderedLocal, iSublPhase, iElement
        real(8) :: dHelperTotal, dNorm, dProduct, dX
        real(8), dimension(nElements) :: dHelperElementFraction

        dProjectedFractionOut = 0D0
        lProjectedOut = .FALSE.
        dHelperElementFraction = 0D0
        if ((iHelperPhaseIn < 1).OR.(iHelperPhaseIn > nSolnPhasesSys)) return
        if ((iOrderedPhaseIn < 1).OR.(iOrderedPhaseIn > nSolnPhasesSys)) return

        iHelperFirst = nSpeciesPhase(iHelperPhaseIn-1) + 1
        iHelperLast  = nSpeciesPhase(iHelperPhaseIn)
        do i = iHelperFirst, iHelperLast
            dX = DMAX1(dMolFraction(i), 0D0)
            do j = 1, nElements
                dHelperElementFraction(j) = dHelperElementFraction(j) + dX * DMAX1(dStoichSpecies(i,j), 0D0)
            end do
        end do

        dHelperTotal = SUM(dHelperElementFraction)
        if (dHelperTotal <= 1D-300) return
        dHelperElementFraction = dHelperElementFraction / dHelperTotal

        iSublPhase = iPhaseSublattice(iOrderedPhaseIn)
        iOrderedFirst = nSpeciesPhase(iOrderedPhaseIn-1) + 1
        iOrderedLast  = nSpeciesPhase(iOrderedPhaseIn)
        if (SIZE(dProjectedFractionOut) /= (iOrderedLast - iOrderedFirst + 1)) return

        do i = iOrderedFirst, iOrderedLast
            iOrderedLocal = i - iOrderedFirst + 1
            dProduct = 1D0
            do s = 1, nSublatticePhase(iSublPhase)
                m = iConstituentSublattice(iSublPhase,s,iOrderedLocal)
                if (m <= 0) cycle
                if (IsVacancyName(cConstituentNameSUB(iSublPhase,s,m))) then
                    if (nConstituentSublattice(iSublPhase,s) > 1) dProduct = 0D0
                    cycle
                end if

                iElement = ConstituentElementIndex(cConstituentNameSUB(iSublPhase,s,m))
                if (iElement <= 0) then
                    dProduct = 0D0
                else
                    dProduct = dProduct * DMAX1(dHelperElementFraction(iElement), 0D0)
                end if
            end do
            dProjectedFractionOut(iOrderedLocal) = dProduct
        end do

        dNorm = SUM(dProjectedFractionOut)
        if (dNorm <= 1D-300) return
        dProjectedFractionOut = dProjectedFractionOut / dNorm
        lProjectedOut = .TRUE.

        return
    end subroutine ProjectCompanionToOrderedSUBOMFraction


    subroutine MarkTwoSetCandidateIdentity(iLevelRowIn, iOrderedPhaseIn)
        integer, intent(in) :: iLevelRowIn, iOrderedPhaseIn

        iCand = 0
        if (allocated(iLevelCandidateFromLevel)) then
            if ((iLevelRowIn >= 1).AND.(iLevelRowIn <= SIZE(iLevelCandidateFromLevel))) then
                iCand = iLevelCandidateFromLevel(iLevelRowIn)
            end if
        end if
        if ((iCand < 1).OR.(iCand > nLevelCandidate)) return

        iLevelCandidateParentPhase(iCand) = iOrderedPhaseIn
        iLevelCandidateDisplayPhase(iCand) = iOrderedPhaseIn
        iLevelCandidateIdentityOrdinal(iCand) = 2
        iLevelCandidateSource(iCand) = 2

        return
    end subroutine MarkTwoSetCandidateIdentity


    logical function IsVacancyName(cName)
        character(*), intent(in) :: cName
        character(8)             :: cUpper

        cUpper = UpperName(cName)
        IsVacancyName = (TRIM(cUpper) == 'VA').OR.(TRIM(cUpper) == 'VACANCY')

        return
    end function IsVacancyName


    integer function ConstituentElementIndex(cName)
        character(*), intent(in) :: cName
        integer                  :: iElement
        character(8)             :: cTarget, cElement

        ConstituentElementIndex = 0
        cTarget = UpperName(cName)
        do iElement = 1, nElements
            cElement = UpperName(cElementName(iElement))
            if (TRIM(cTarget) == TRIM(cElement)) then
                ConstituentElementIndex = iElement
                return
            end if
        end do

        return
    end function ConstituentElementIndex


    character(8) function UpperName(cName)
        character(*), intent(in) :: cName
        integer                  :: iChar, iCode

        UpperName = '        '
        do iChar = 1, MIN(LEN(cName), LEN(UpperName))
            iCode = IACHAR(cName(iChar:iChar))
            if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                UpperName(iChar:iChar) = ACHAR(iCode - 32)
            else
                UpperName(iChar:iChar) = cName(iChar:iChar)
            end if
        end do

        return
    end function UpperName

end subroutine RegisterSUBOMTwoSetCandidateRows
