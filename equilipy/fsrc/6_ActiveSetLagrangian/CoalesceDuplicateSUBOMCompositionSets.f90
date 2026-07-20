!> \brief Coalesce duplicate same-parent SUBOM composition-set slots.
!!
!! \details Merges active SUBOM composition-set slots that share one
!! thermodynamic parent and have collapsed to the same slot-local site
!! fractions or the same element composition.  This is representation cleanup
!! only; distinct same-parent composition sets are left untouched.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    CoalesceDuplicateSUBOMCompositionSets.f90
!> \brief   Coalesce duplicate same-parent SUBOM composition-set slots.
!> \author  S.Y. Kwon
!> \date    Jul. 03, 2026
!> \sa      RunLagrangianGEM.f90
!> \sa      Level2Lagrange.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Coalesced same-parent SUBOM composition sets using exact structural and ordering identities while preserving distinct ordered states.
    !
    !
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to remove duplicate active
!! composition-set slots for a SUBOM parent after the Lagrangian solve has
!! collapsed them to the same site-fraction state or to the same element
!! composition with different ordering.  A single global
!! dMolFraction range cannot represent two slots of the same parent, so
!! identical-composition slots must be coalesced before thermodynamic
!! postprocessing.
!
!
! Required input variables:
! =========================
!
! iAssemblage                  Current active Lagrangian assemblage.
! dMolesPhase                  Current active phase amounts.
! dActiveSlotSiteFraction      Slot-local CEF site fractions.
! iActiveSlotThermoPhase       Thermodynamic parent phase id for each slot.
!
!
! Output/updated variables:
! =========================
!
!> \param[out] lChanged        True when one duplicate same-parent pair was merged.
!
! iAssemblage                  Compacted active assemblage.
! dMolesPhase                  Merged phase amounts.
! dMolesSpecies/dMolFraction   Rebuilt from retained slot-local site fractions.
! dActiveSlot*                 Rebuilt to match the compacted active set.
! nSolnPhases/lSolnPhases      Updated active solution phase count and flags.
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
! RunLagrangianGEM             Calls after CEF convergence and before result assembly.
!
!
! Numerical assumptions:
! ======================
!
! - Parser-declared partition companions coalesce only at exact parent-site
!   equality.  The legacy helper-only experimental path retains its historical
!   numerical duplicate rule.
! - Same-composition ordering splits are collapsed to the lower fixed-site
!   parent Gibbs branch; distinct element compositions are not changed.
! - If one duplicate slot is the DIS_PART companion display phase, that
!   physical disordered representation is retained only when the collapsed
!   parent site state lies on the disordered manifold.
! - The retained product/endmember fractions are rebuilt from the authoritative
!   slot-local site fractions, not from the shared global dMolFraction range.
!
!-------------------------------------------------------------------------------------------------------------



subroutine CoalesceDuplicateSUBOMCompositionSets(lChanged)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lChanged

    integer :: iSlotA, iSlotB, iPhaseA, iPhaseB, iKeepSlot, iRemoveSlot, iKeepPhase
    real(8) :: dSiteDifference, dCompositionDifference
    real(8), parameter :: dSameSetTolerance = 1D-6
    logical :: lCanonicalPair, lDisplayChanged

    lChanged = .FALSE.
    if (.NOT.allocated(iActiveSlotThermoPhase)) return
    if (.NOT.allocated(iActiveSlotDisplayPhase)) return
    if (.NOT.allocated(dActiveSlotSiteFraction)) return

    call SyncSUBOMDisplayPhases(lDisplayChanged)
    if (nSolnPhases <= 1) return

    do iSlotA = 1, nElements - 1
        if (iAssemblage(iSlotA) >= 0) cycle
        iPhaseA = iActiveSlotThermoPhase(iSlotA)
        if ((iPhaseA <= 0).OR.(iPhaseA > nSolnPhasesSys)) cycle
        if (TRIM(cSolnPhaseType(iPhaseA)) /= 'SUBOM') cycle

        do iSlotB = iSlotA + 1, nElements
            if (iAssemblage(iSlotB) >= 0) cycle
            iPhaseB = iActiveSlotThermoPhase(iSlotB)
            if (iPhaseB /= iPhaseA) cycle

            dSiteDifference = MAXVAL(DABS(&
                dActiveSlotSiteFraction(iSlotA,:,:) - dActiveSlotSiteFraction(iSlotB,:,:)))
            lCanonicalPair = CanonicalCompanionPhase(iPhaseA) > 0
            if (lCanonicalPair) then
                if (ANY(dActiveSlotSiteFraction(iSlotA,:,:) /= &
                    dActiveSlotSiteFraction(iSlotB,:,:))) cycle
            else
                if (dSiteDifference > dSameSetTolerance) cycle
            end if

            call ChooseRetainedDuplicateSlot(iSlotA, iSlotB, iPhaseA, &
                iKeepSlot, iRemoveSlot, iKeepPhase)
            call MergeDuplicateSlots(iKeepSlot, iRemoveSlot, iPhaseA, iKeepPhase)
            lChanged = .TRUE.
            return
        end do
    end do

    if (lSUBOMTwoSetCandidateEnabled) then
        do iSlotA = 1, nElements - 1
            if (iAssemblage(iSlotA) >= 0) cycle
            iPhaseA = iActiveSlotThermoPhase(iSlotA)
            if ((iPhaseA <= 0).OR.(iPhaseA > nSolnPhasesSys)) cycle
            if (TRIM(cSolnPhaseType(iPhaseA)) /= 'SUBOM') cycle

            do iSlotB = iSlotA + 1, nElements
                if (iAssemblage(iSlotB) >= 0) cycle
                iPhaseB = iActiveSlotThermoPhase(iSlotB)
                if (iPhaseB /= iPhaseA) cycle
                if (CanonicalCompanionPhase(iPhaseA) > 0) cycle

                dCompositionDifference = CompositionDifferenceFromSite(iPhaseA, &
                    dActiveSlotSiteFraction(iSlotA,:,:), dActiveSlotSiteFraction(iSlotB,:,:))
                if (dCompositionDifference > dSameSetTolerance) cycle

                call ChooseLowerGibbsSameCompositionSlot(iSlotA, iSlotB, iPhaseA, &
                    iKeepSlot, iRemoveSlot, iKeepPhase)
                call MergeDuplicateSlots(iKeepSlot, iRemoveSlot, iPhaseA, iKeepPhase)
                lChanged = .TRUE.
                return
            end do
        end do
    end if

    return

contains

    subroutine SyncSUBOMDisplayPhases(lChangedOut)
        logical, intent(out) :: lChangedOut

        integer :: iSlot, iEntry, iDisplayPhase, iParentPhase, iPreferredPhase

        lChangedOut = .FALSE.
        do iSlot = 1, nElements
            iEntry = iAssemblage(iSlot)
            if (iEntry >= 0) cycle

            iDisplayPhase = iActiveSlotDisplayPhase(iSlot)
            if (iDisplayPhase <= 0) iDisplayPhase = -iEntry
            iParentPhase = iActiveSlotThermoPhase(iSlot)
            if ((iParentPhase <= 0).OR.(iParentPhase > nSolnPhasesSys)) cycle
            if (TRIM(cSolnPhaseType(iParentPhase)) /= 'SUBOM') cycle

            call PreferredDisplayForSite(iSlot, iParentPhase, iPreferredPhase)
            if (iPreferredPhase <= 0) cycle
            if (iPreferredPhase == iDisplayPhase) cycle

            iActiveSlotDisplayPhase(iSlot) = iPreferredPhase
            lChangedOut = .TRUE.
        end do

        return
    end subroutine SyncSUBOMDisplayPhases


    subroutine PreferredDisplayForSite(iSlotIn, iParentPhase, iDisplayPhaseOut)
        integer, intent(in)  :: iSlotIn, iParentPhase
        integer, intent(out) :: iDisplayPhaseOut

        integer :: iHelperPhase, iClass, iPriorClass

        iDisplayPhaseOut = iParentPhase
        iHelperPhase = CanonicalCompanionPhase(iParentPhase)
        if (iHelperPhase > 0) then
            iPriorClass = OD_CANDIDATE_NOT_EVALUATED
            if (allocated(iActiveSlotODClass)) iPriorClass = iActiveSlotODClass(iSlotIn)
            if (allocated(iODCandidateClass)) then
                if ((iParentPhase >= 1).AND.&
                    (iParentPhase <= SIZE(iODCandidateClass))) then
                    if (iODCandidateClass(iParentPhase) /= &
                        OD_CANDIDATE_NOT_EVALUATED) then
                        iPriorClass = iODCandidateClass(iParentPhase)
                    end if
                end if
            end if
            call ClassifyOrderDisorderActiveSlot(iParentPhase, &
                dActiveSlotSiteFraction(iSlotIn,:,:), iClass, iHelperPhase)
            ! Structural equality classifies the site geometry; the candidate
            ! sweep's curvature fact decides whether that symmetric point is a
            ! stable disordered chart, an unstable ordering branch, or a node.
            ! Preserve both directions of that typed decision without a
            ! constitution-distance cutoff.
            if (.NOT.RepeatedParentCompositionSetsActive(iParentPhase)) then
                if ((iClass == OD_CANDIDATE_ORDERED).AND.&
                    StableDisorderedCandidateClass(iPriorClass)) then
                    iClass = OD_CANDIDATE_DISORDERED
                end if
                if ((iClass == OD_CANDIDATE_DISORDERED).AND.&
                    (.NOT.StableDisorderedCandidateClass(iPriorClass)).AND.&
                    (iPriorClass /= OD_CANDIDATE_NOT_EVALUATED)) then
                    iClass = iPriorClass
                end if
            end if
            if (allocated(iActiveSlotODClass)) iActiveSlotODClass(iSlotIn) = iClass
            if (iClass == OD_CANDIDATE_DISORDERED) iDisplayPhaseOut = iHelperPhase
            if (iClass == OD_CANDIDATE_ORDERED) iDisplayPhaseOut = iParentPhase
            if ((iClass /= OD_CANDIDATE_DISORDERED).AND.&
                (iClass /= OD_CANDIDATE_ORDERED)) then
                if (iActiveSlotDisplayPhase(iSlotIn) > 0) then
                    iDisplayPhaseOut = iActiveSlotDisplayPhase(iSlotIn)
                end if
            end if
            return
        end if
        if (allocated(iDisorderedPhase)) then
            if ((iParentPhase >= 1).AND.(iParentPhase <= SIZE(iDisorderedPhase))) then
                iHelperPhase = iDisorderedPhase(iParentPhase)
            end if
        end if
        if (IsDisorderedManifold(iParentPhase, dActiveSlotSiteFraction(iSlotIn,:,:))) then
            if (iHelperPhase > 0) iDisplayPhaseOut = iHelperPhase
        end if

        return
    end subroutine PreferredDisplayForSite


    logical function RepeatedParentCompositionSetsActive(iParentPhaseIn)
        integer, intent(in) :: iParentPhaseIn

        integer :: iSlot, nParentSlots

        RepeatedParentCompositionSetsActive = .FALSE.
        if (.NOT.allocated(iActiveSlotThermoPhase)) return

        nParentSlots = 0
        do iSlot = 1, MIN(nElements, SIZE(iActiveSlotThermoPhase))
            if (iAssemblage(iSlot) >= 0) cycle
            if (iActiveSlotThermoPhase(iSlot) /= iParentPhaseIn) cycle
            nParentSlots = nParentSlots + 1
            if (nParentSlots > 1) then
                RepeatedParentCompositionSetsActive = .TRUE.
                return
            end if
        end do

        return
    end function RepeatedParentCompositionSetsActive


    subroutine ChooseRetainedDuplicateSlot(iSlotAIn, iSlotBIn, iParentPhase, &
        iKeepSlotOut, iRemoveSlotOut, iKeepPhaseOut)
        integer, intent(in)  :: iSlotAIn, iSlotBIn, iParentPhase
        integer, intent(out) :: iKeepSlotOut, iRemoveSlotOut, iKeepPhaseOut

        integer :: iDisplayA, iDisplayB, iHelperPhase
        logical :: lDisorderedManifold

        iDisplayA = iActiveSlotDisplayPhase(iSlotAIn)
        iDisplayB = iActiveSlotDisplayPhase(iSlotBIn)
        iHelperPhase = CanonicalCompanionPhase(iParentPhase)
        if ((iHelperPhase == 0).AND.allocated(iDisorderedPhase)) then
            if ((iParentPhase >= 1).AND.(iParentPhase <= SIZE(iDisorderedPhase))) then
                iHelperPhase = iDisorderedPhase(iParentPhase)
            end if
        end if

        iKeepSlotOut = iSlotAIn
        iRemoveSlotOut = iSlotBIn
        iKeepPhaseOut = iParentPhase
        lDisorderedManifold = IsDisorderedManifold(iParentPhase, &
            dActiveSlotSiteFraction(iSlotAIn,:,:))

        if (CanonicalCompanionPhase(iParentPhase) > 0) then
            iKeepPhaseOut = iParentPhase
            return
        end if

        if (lDisorderedManifold.AND.(iHelperPhase > 0).AND.(iDisplayA == iHelperPhase)) then
            iKeepSlotOut = iSlotAIn
            iRemoveSlotOut = iSlotBIn
            iKeepPhaseOut = iHelperPhase
        else if (lDisorderedManifold.AND.(iHelperPhase > 0).AND.(iDisplayB == iHelperPhase)) then
            iKeepSlotOut = iSlotBIn
            iRemoveSlotOut = iSlotAIn
            iKeepPhaseOut = iHelperPhase
        else if (iDisplayA > 0) then
            iKeepPhaseOut = iDisplayA
        end if

        return
    end subroutine ChooseRetainedDuplicateSlot


    subroutine ChooseLowerGibbsSameCompositionSlot(iSlotAIn, iSlotBIn, iParentPhase, &
        iKeepSlotOut, iRemoveSlotOut, iKeepPhaseOut)
        integer, intent(in)  :: iSlotAIn, iSlotBIn, iParentPhase
        integer, intent(out) :: iKeepSlotOut, iRemoveSlotOut, iKeepPhaseOut

        real(8) :: dGibbsA, dGibbsB

        dGibbsA = ScalarGibbsFromSite(iParentPhase, dActiveSlotSiteFraction(iSlotAIn,:,:))
        dGibbsB = ScalarGibbsFromSite(iParentPhase, dActiveSlotSiteFraction(iSlotBIn,:,:))

        iKeepSlotOut = iSlotAIn
        iRemoveSlotOut = iSlotBIn
        if (dGibbsB < dGibbsA) then
            iKeepSlotOut = iSlotBIn
            iRemoveSlotOut = iSlotAIn
        end if

        call DisplayPhaseForSelectedSlot(iKeepSlotOut, iParentPhase, iKeepPhaseOut)

        return
    end subroutine ChooseLowerGibbsSameCompositionSlot


    subroutine DisplayPhaseForSelectedSlot(iSelectedSlot, iParentPhase, iDisplayPhaseOut)
        integer, intent(in)  :: iSelectedSlot, iParentPhase
        integer, intent(out) :: iDisplayPhaseOut

        integer :: iHelperPhase

        iDisplayPhaseOut = iParentPhase
        iHelperPhase = CanonicalCompanionPhase(iParentPhase)
        if (iHelperPhase > 0) then
            call PreferredDisplayForSite(iSelectedSlot, iParentPhase, iDisplayPhaseOut)
            return
        end if
        if ((iHelperPhase == 0).AND.allocated(iDisorderedPhase)) then
            if ((iParentPhase >= 1).AND.(iParentPhase <= SIZE(iDisorderedPhase))) then
                iHelperPhase = iDisorderedPhase(iParentPhase)
            end if
        end if

        if (IsDisorderedManifold(iParentPhase, dActiveSlotSiteFraction(iSelectedSlot,:,:))) then
            if (iHelperPhase > 0) then
                iDisplayPhaseOut = iHelperPhase
            else if (iActiveSlotDisplayPhase(iSelectedSlot) > 0) then
                iDisplayPhaseOut = iActiveSlotDisplayPhase(iSelectedSlot)
            end if
        else if (iActiveSlotDisplayPhase(iSelectedSlot) > 0) then
            if (iActiveSlotDisplayPhase(iSelectedSlot) /= iHelperPhase) then
                iDisplayPhaseOut = iActiveSlotDisplayPhase(iSelectedSlot)
            end if
        end if

        return
    end subroutine DisplayPhaseForSelectedSlot


    logical function StableDisorderedCandidateClass(iClassIn)
        integer, intent(in) :: iClassIn

        StableDisorderedCandidateClass = &
            (iClassIn == OD_CANDIDATE_DISORDERED).OR.&
            (iClassIn == OD_CANDIDATE_DISORDERED_PROJECTED).OR.&
            (iClassIn == OD_CANDIDATE_AMBIGUOUS_ROUNDOFF)

        return

    end function StableDisorderedCandidateClass


    logical function IsDisorderedManifold(iParentPhase, dSiteIn)
        integer, intent(in) :: iParentPhase
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(in) :: dSiteIn

        integer :: iSublPhase, iSubA, iSubB, iCon, nConstituent
        integer :: iClass, iCompanionPhase

        if (CanonicalCompanionPhase(iParentPhase) > 0) then
            call ClassifyOrderDisorderActiveSlot(iParentPhase, dSiteIn, &
                iClass, iCompanionPhase)
            IsDisorderedManifold = iClass == OD_CANDIDATE_DISORDERED
            return
        end if

        IsDisorderedManifold = .TRUE.
        iSublPhase = iPhaseSublattice(iParentPhase)
        if (iSublPhase <= 0) return

        do iSubA = 1, nSublatticePhase(iSublPhase) - 1
            nConstituent = nConstituentSublattice(iSublPhase,iSubA)
            if (nConstituent <= 1) cycle
            do iSubB = iSubA + 1, nSublatticePhase(iSublPhase)
                if (.NOT.EquivalentSublattices(iSublPhase, iSubA, iSubB)) cycle
                do iCon = 1, nConstituent
                    if (DABS(dSiteIn(iSubA,iCon) - dSiteIn(iSubB,iCon)) > dSameSetTolerance) then
                        IsDisorderedManifold = .FALSE.
                        return
                    end if
                end do
            end do
        end do

        return
    end function IsDisorderedManifold


    logical function EquivalentSublattices(iSublPhase, iSubA, iSubB)
        integer, intent(in) :: iSublPhase, iSubA, iSubB

        integer :: iCon, nConstituent

        EquivalentSublattices = .FALSE.
        if (nConstituentSublattice(iSublPhase,iSubA) /= &
            nConstituentSublattice(iSublPhase,iSubB)) return
        if (DABS(dStoichSublattice(iSublPhase,iSubA) - &
            dStoichSublattice(iSublPhase,iSubB)) > 1D-12) return

        nConstituent = nConstituentSublattice(iSublPhase,iSubA)
        do iCon = 1, nConstituent
            if (iConstituentSublattice(iSublPhase,iSubA,iCon) /= &
                iConstituentSublattice(iSublPhase,iSubB,iCon)) return
        end do

        EquivalentSublattices = .TRUE.

        return
    end function EquivalentSublattices


    real(8) function CompositionDifferenceFromSite(iParentPhase, dSiteA, dSiteB)
        integer, intent(in) :: iParentPhase
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(in) :: dSiteA, dSiteB

        real(8), dimension(nElements) :: dCompositionA, dCompositionB

        call CompositionFromSite(iParentPhase, dSiteA, dCompositionA)
        call CompositionFromSite(iParentPhase, dSiteB, dCompositionB)
        CompositionDifferenceFromSite = MAXVAL(DABS(dCompositionA - dCompositionB))

        return
    end function CompositionDifferenceFromSite


    real(8) function ScalarGibbsFromSite(iParentPhase, dSiteIn)
        integer, intent(in) :: iParentPhase
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(in) :: dSiteIn

        integer :: iFirst, iLast, nLocal, nSiteOut, iInfo, nSiteCapacity
        real(8) :: dScalarH, dScalarS, dScalarCp
        real(8), dimension(:), allocatable :: dGradientG, dGradientH, dGradientS, dGradientCp
        real(8), dimension(:), allocatable :: dFractionSave, dProductFraction

        ScalarGibbsFromSite = HUGE(1D0)
        if ((iParentPhase <= 0).OR.(iParentPhase > nSolnPhasesSys)) return
        iFirst = nSpeciesPhase(iParentPhase-1) + 1
        iLast  = nSpeciesPhase(iParentPhase)
        nLocal = iLast - iFirst + 1
        if (nLocal <= 0) return
        nSiteCapacity = nMaxSublatticeSys * nMaxConstituentSys

        allocate(dFractionSave(nLocal), dProductFraction(nLocal))
        allocate(dGradientG(nSiteCapacity), dGradientH(nSiteCapacity), &
            dGradientS(nSiteCapacity), dGradientCp(nSiteCapacity))
        dFractionSave = dMolFraction(iFirst:iLast)
        call ProductFractionFromSite(iParentPhase, dSiteIn, dProductFraction)
        dMolFraction(iFirst:iLast) = dProductFraction

        call CompGradientSUBL(iParentPhase, nSiteCapacity, dGradientG, ScalarGibbsFromSite, &
            dGradientH, dGradientS, dGradientCp, dScalarH, dScalarS, dScalarCp, nSiteOut, iInfo)
        if (iInfo /= 0) ScalarGibbsFromSite = HUGE(1D0)

        dMolFraction(iFirst:iLast) = dFractionSave
        deallocate(dFractionSave, dProductFraction)
        deallocate(dGradientG, dGradientH, dGradientS, dGradientCp)

        return
    end function ScalarGibbsFromSite


    subroutine MergeDuplicateSlots(iKeepSlot, iRemoveSlot, iParentPhase, iKeepPhase)
        integer, intent(in) :: iKeepSlot, iRemoveSlot, iParentPhase, iKeepPhase

        integer :: iFirst, iLast, nLocal
        real(8) :: dMergedAmount
        real(8), dimension(:), allocatable :: dProductFraction

        iFirst = nSpeciesPhase(iKeepPhase-1) + 1
        iLast  = nSpeciesPhase(iKeepPhase)
        nLocal = iLast - iFirst + 1
        if (nLocal <= 0) return

        allocate(dProductFraction(nLocal))
        if (CoalescenceLeavesSinglePhase(iRemoveSlot)) then
            call ProductFractionFromBulkComposition(iKeepPhase, dProductFraction)
        else
            call ProductFractionForRetainedPhase(iParentPhase, iKeepPhase, &
                dActiveSlotSiteFraction(iKeepSlot,:,:), dProductFraction)
        end if

        dMergedAmount = dMolesPhase(iKeepSlot) + dMolesPhase(iRemoveSlot)
        call ClearPhaseSpecies(iParentPhase)
        call ClearPhaseSpecies(iKeepPhase)
        if (iActiveSlotDisplayPhase(iKeepSlot) > 0) call ClearPhaseSpecies(iActiveSlotDisplayPhase(iKeepSlot))
        if (iActiveSlotDisplayPhase(iRemoveSlot) > 0) call ClearPhaseSpecies(iActiveSlotDisplayPhase(iRemoveSlot))

        iAssemblage(iKeepSlot) = -iKeepPhase
        dMolesPhase(iKeepSlot) = dMergedAmount
        dMolFraction(iFirst:iLast) = dProductFraction
        dMolesSpecies(iFirst:iLast) = dMergedAmount * dProductFraction
        dMolesPhase(iRemoveSlot) = 0D0

        call CompactAfterRemovingSlot(iRemoveSlot)
        call RebuildActiveSlotState
        call SyncSUBOMDisplayPhases(lDisplayChanged)

        deallocate(dProductFraction)

        return
    end subroutine MergeDuplicateSlots


    subroutine RebuildGlobalSpeciesFromActiveSlots(lSpeciesChangedOut)
        logical, intent(out) :: lSpeciesChangedOut

        integer :: iSlot, iEntry, iSolnPhase, iParentPhase, iFirst, iLast, nLocal
        integer :: iPrevSlot, iOrdinal
        real(8) :: dSumLocal
        real(8), dimension(nElements) :: dComposition
        real(8), dimension(:), allocatable :: dProductFraction

        lSpeciesChangedOut = .FALSE.
        dMolesSpecies = 0D0
        if (allocated(lSolnPhases)) lSolnPhases = .FALSE.
        if (allocated(iActiveSlotThermoPhase)) iActiveSlotThermoPhase = 0
        if (allocated(iActiveSlotDisplayPhase)) iActiveSlotDisplayPhase = 0
        if (allocated(iActiveSlotIdentityOrdinal)) iActiveSlotIdentityOrdinal = 0
        if (allocated(iActiveSlotODClass)) iActiveSlotODClass = OD_CANDIDATE_NOT_EVALUATED
        if (allocated(dActiveSlotMolFraction)) dActiveSlotMolFraction = 0D0

        do iSlot = 1, nElements
            iEntry = iAssemblage(iSlot)
            if (iEntry == 0) cycle

            if (iEntry > 0) then
                dMolesSpecies(iEntry) = dMolesPhase(iSlot)
                cycle
            end if

            iSolnPhase = -iEntry
            iParentPhase = ParentPhaseForDisplay(iSolnPhase)
            if ((iSolnPhase < 1).OR.(iSolnPhase > nSolnPhasesSys)) cycle
            iFirst = nSpeciesPhase(iSolnPhase-1) + 1
            iLast  = nSpeciesPhase(iSolnPhase)
            nLocal = iLast - iFirst + 1
            if (nLocal <= 0) cycle

            allocate(dProductFraction(nLocal))
            dProductFraction = 0D0
            if ((iParentPhase >= 1).AND.(iParentPhase <= nSolnPhasesSys).AND.&
                (TRIM(cSolnPhaseType(iParentPhase)) == 'SUBOM').AND.&
                (SUM(dActiveSlotSiteFraction(iSlot,:,:)) > 0D0)) then
                if (iSolnPhase == iParentPhase) then
                    call ProductFractionFromSite(iParentPhase, dActiveSlotSiteFraction(iSlot,:,:), &
                        dProductFraction)
                else
                    call CompositionFromSite(iParentPhase, dActiveSlotSiteFraction(iSlot,:,:), dComposition)
                    call ProductFractionFromComposition(iSolnPhase, dComposition, dProductFraction)
                end if
            else if (allocated(dActiveSlotMolFraction)) then
                if (iLast <= SIZE(dActiveSlotMolFraction,2)) then
                    dProductFraction = DMAX1(dActiveSlotMolFraction(iSlot,iFirst:iLast), 0D0)
                    dSumLocal = SUM(dProductFraction)
                    if (dSumLocal > 1D-300) dProductFraction = dProductFraction / dSumLocal
                end if
            end if

            dSumLocal = SUM(dProductFraction)
            if (dSumLocal > 1D-300) then
                dProductFraction = dProductFraction / dSumLocal
                if (MAXVAL(DABS(dMolFraction(iFirst:iLast) - dProductFraction)) > 1D-8) then
                    lSpeciesChangedOut = .TRUE.
                end if
                dMolFraction(iFirst:iLast) = dProductFraction
                dMolesSpecies(iFirst:iLast) = dMolesPhase(iSlot) * dProductFraction
                if (allocated(dActiveSlotMolFraction)) then
                    dActiveSlotMolFraction(iSlot,iFirst:iLast) = dProductFraction
                end if
            end if

            if (allocated(lSolnPhases)) lSolnPhases(iSolnPhase) = .TRUE.
            if (allocated(iActiveSlotThermoPhase)) iActiveSlotThermoPhase(iSlot) = iParentPhase
            if (allocated(iActiveSlotDisplayPhase)) iActiveSlotDisplayPhase(iSlot) = iSolnPhase
            if (allocated(iActiveSlotIdentityOrdinal)) then
                iOrdinal = 1
                do iPrevSlot = 1, iSlot - 1
                    if (iActiveSlotThermoPhase(iPrevSlot) == iParentPhase) iOrdinal = iOrdinal + 1
                end do
                iActiveSlotIdentityOrdinal(iSlot) = iOrdinal
            end if

            deallocate(dProductFraction)
        end do

        return
    end subroutine RebuildGlobalSpeciesFromActiveSlots


    subroutine ProductFractionForRetainedPhase(iParentPhase, iKeepPhase, dParentSiteIn, dProductFraction)
        integer, intent(in) :: iParentPhase, iKeepPhase
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(in) :: dParentSiteIn
        real(8), dimension(:), intent(out) :: dProductFraction

        real(8), dimension(nElements) :: dComposition

        if (iKeepPhase == iParentPhase) then
            call ProductFractionFromSite(iKeepPhase, dParentSiteIn, dProductFraction)
        else
            call CompositionFromSite(iParentPhase, dParentSiteIn, dComposition)
            call ProductFractionFromComposition(iKeepPhase, dComposition, dProductFraction)
        end if

        return
    end subroutine ProductFractionForRetainedPhase


    logical function CoalescenceLeavesSinglePhase(iRemoveSlot)
        integer, intent(in) :: iRemoveSlot

        integer :: iSlot, nRemaining

        nRemaining = 0
        do iSlot = 1, nElements
            if (iSlot == iRemoveSlot) cycle
            if (iAssemblage(iSlot) == 0) cycle
            if (dMolesPhase(iSlot) <= 1D-20) cycle
            nRemaining = nRemaining + 1
        end do

        CoalescenceLeavesSinglePhase = nRemaining == 1

        return
    end function CoalescenceLeavesSinglePhase


    subroutine ProductFractionFromBulkComposition(iKeepPhase, dProductFraction)
        integer, intent(in) :: iKeepPhase
        real(8), dimension(:), intent(out) :: dProductFraction

        real(8) :: dNorm
        real(8), dimension(nElements) :: dComposition

        dComposition = DMAX1(dMolesElement, 0D0)
        dNorm = SUM(dComposition)
        if (dNorm > 1D-300) then
            dComposition = dComposition / dNorm
            call ProductFractionFromComposition(iKeepPhase, dComposition, dProductFraction)
        else
            dProductFraction = 0D0
        end if

        return
    end subroutine ProductFractionFromBulkComposition


    subroutine ProductFractionFromSite(iSolnPhase, dSiteIn, dProductFraction)
        integer, intent(in) :: iSolnPhase
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(in) :: dSiteIn
        real(8), dimension(:), intent(out) :: dProductFraction

        integer :: iFirst, iLast, iSpecies, iLocal, iSublPhase, iSub, iCon
        real(8) :: dProduct, dNorm

        dProductFraction = 0D0
        iSublPhase = iPhaseSublattice(iSolnPhase)
        if (iSublPhase <= 0) return

        iFirst = nSpeciesPhase(iSolnPhase-1) + 1
        iLast  = nSpeciesPhase(iSolnPhase)
        do iSpecies = iFirst, iLast
            iLocal = iSpecies - iFirst + 1
            if (iLocal > SIZE(dProductFraction)) exit

            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iSublPhase)
                iCon = iConstituentSublattice(iSublPhase,iSub,iLocal)
                if (iCon <= 0) cycle
                dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 0D0)
            end do
            dProductFraction(iLocal) = dProduct
        end do

        dNorm = SUM(dProductFraction)
        if (dNorm > 1D-300) dProductFraction = dProductFraction / dNorm

        return
    end subroutine ProductFractionFromSite


    subroutine CompositionFromSite(iSolnPhase, dSiteIn, dCompositionOut)
        integer, intent(in) :: iSolnPhase
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(in) :: dSiteIn
        real(8), dimension(nElements), intent(out) :: dCompositionOut

        integer :: iFirst, iLast, nLocal, iSpecies, iElement
        real(8) :: dNorm
        real(8), dimension(:), allocatable :: dProductFraction

        dCompositionOut = 0D0
        iFirst = nSpeciesPhase(iSolnPhase-1) + 1
        iLast  = nSpeciesPhase(iSolnPhase)
        nLocal = iLast - iFirst + 1
        if (nLocal <= 0) return

        allocate(dProductFraction(nLocal))
        call ProductFractionFromSite(iSolnPhase, dSiteIn, dProductFraction)

        do iSpecies = iFirst, iLast
            do iElement = 1, nElements
                dCompositionOut(iElement) = dCompositionOut(iElement) + &
                    dProductFraction(iSpecies-iFirst+1) * dStoichSpecies(iSpecies,iElement) / &
                    DFLOAT(iParticlesPerMole(iSpecies))
            end do
        end do

        dNorm = SUM(DABS(dCompositionOut))
        if (dNorm > 1D-300) dCompositionOut = dCompositionOut / dNorm
        deallocate(dProductFraction)

        return
    end subroutine CompositionFromSite


    subroutine ProductFractionFromComposition(iSolnPhase, dCompositionIn, dProductFraction)
        integer, intent(in) :: iSolnPhase
        real(8), dimension(nElements), intent(in) :: dCompositionIn
        real(8), dimension(:), intent(out) :: dProductFraction

        integer :: iFirst, iLast, iSpecies, iLocal, iSublPhase, iSub, iCon, iElement
        real(8) :: dProduct, dNorm

        dProductFraction = 0D0
        iSublPhase = iPhaseSublattice(iSolnPhase)
        if (iSublPhase <= 0) return

        iFirst = nSpeciesPhase(iSolnPhase-1) + 1
        iLast  = nSpeciesPhase(iSolnPhase)
        do iSpecies = iFirst, iLast
            iLocal = iSpecies - iFirst + 1
            if (iLocal > SIZE(dProductFraction)) exit

            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iSublPhase)
                iCon = iConstituentSublattice(iSublPhase,iSub,iLocal)
                if (iCon <= 0) cycle
                if (IsVacancyConstituent(cConstituentNameSUB(iSublPhase,iSub,iCon))) cycle

                iElement = ConstituentElementIndex(cConstituentNameSUB(iSublPhase,iSub,iCon))
                if (iElement <= 0) then
                    dProduct = 0D0
                    exit
                end if
                dProduct = dProduct * DMAX1(dCompositionIn(iElement), 0D0)
            end do
            dProductFraction(iLocal) = dProduct
        end do

        dNorm = SUM(dProductFraction)
        if (dNorm > 1D-300) dProductFraction = dProductFraction / dNorm

        return
    end subroutine ProductFractionFromComposition


    subroutine ClearPhaseSpecies(iSolnPhase)
        integer, intent(in) :: iSolnPhase

        integer :: iFirst, iLast

        if ((iSolnPhase < 1).OR.(iSolnPhase > nSolnPhasesSys)) return
        iFirst = nSpeciesPhase(iSolnPhase-1) + 1
        iLast  = nSpeciesPhase(iSolnPhase)
        if ((iFirst < 1).OR.(iLast > nSpecies)) return
        dMolesSpecies(iFirst:iLast) = 0D0

        return
    end subroutine ClearPhaseSpecies


    subroutine CompactAfterRemovingSlot(iRemoveSlot)
        integer, intent(in) :: iRemoveSlot

        integer :: i, iEntry, iSolnPhase, iFirst, iLast, iDest
        integer :: nConPhasesNew, nSolnPhasesNew
        integer, dimension(nElements) :: iAssemblageNew
        real(8), dimension(nElements) :: dMolesPhaseNew
        real(8), dimension(nSpecies) :: dMolesSpeciesNew

        iAssemblageNew = 0
        dMolesPhaseNew = 0D0
        dMolesSpeciesNew = 0D0
        nConPhasesNew = 0
        nSolnPhasesNew = 0

        do i = 1, nElements
            if (i == iRemoveSlot) cycle
            iEntry = iAssemblage(i)
            if (iEntry <= 0) cycle

            nConPhasesNew = nConPhasesNew + 1
            iAssemblageNew(nConPhasesNew) = iEntry
            dMolesPhaseNew(nConPhasesNew) = dMolesPhase(i)
            dMolesSpeciesNew(iEntry) = dMolesPhase(i)
        end do

        do i = nElements, 1, -1
            if (i == iRemoveSlot) cycle
            iEntry = iAssemblage(i)
            if (iEntry >= 0) cycle

            iSolnPhase = -iEntry
            nSolnPhasesNew = nSolnPhasesNew + 1
            iDest = nElements - nSolnPhasesNew + 1
            iAssemblageNew(iDest) = iEntry
            dMolesPhaseNew(iDest) = dMolesPhase(i)

            iFirst = nSpeciesPhase(iSolnPhase-1) + 1
            iLast  = nSpeciesPhase(iSolnPhase)
            dMolesSpeciesNew(iFirst:iLast) = dMolesPhase(i) * dMolFraction(iFirst:iLast)
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

        return
    end subroutine CompactAfterRemovingSlot


    subroutine RebuildActiveSlotState
        integer :: iSlot, iEntry, iSolnPhase, iFirst, iLast, iPrevSlot, iOrdinal
        integer :: iParentPhase
        real(8) :: dSumLocal

        if (allocated(iActiveSlotThermoPhase)) iActiveSlotThermoPhase = 0
        if (allocated(iActiveSlotDisplayPhase)) iActiveSlotDisplayPhase = 0
        if (allocated(iActiveSlotIdentityOrdinal)) iActiveSlotIdentityOrdinal = 0
        if (allocated(iActiveSlotODClass)) iActiveSlotODClass = OD_CANDIDATE_NOT_EVALUATED
        if (allocated(dActiveSlotMolFraction)) dActiveSlotMolFraction = 0D0
        if (allocated(dActiveSlotSiteFraction)) dActiveSlotSiteFraction = 0D0

        do iSlot = 1, nElements
            iEntry = iAssemblage(iSlot)
            if (iEntry >= 0) cycle

            iSolnPhase = -iEntry
            iParentPhase = ParentPhaseForDisplay(iSolnPhase)
            if (allocated(iActiveSlotThermoPhase)) iActiveSlotThermoPhase(iSlot) = iParentPhase
            if (allocated(iActiveSlotDisplayPhase)) iActiveSlotDisplayPhase(iSlot) = iSolnPhase
            if (allocated(iActiveSlotIdentityOrdinal)) then
                iOrdinal = 1
                do iPrevSlot = 1, iSlot - 1
                    if (allocated(iActiveSlotThermoPhase)) then
                        if (iActiveSlotThermoPhase(iPrevSlot) == iParentPhase) iOrdinal = iOrdinal + 1
                    end if
                end do
                iActiveSlotIdentityOrdinal(iSlot) = iOrdinal
            end if

            iFirst = nSpeciesPhase(iSolnPhase-1) + 1
            iLast  = nSpeciesPhase(iSolnPhase)
            if (allocated(dActiveSlotMolFraction)) then
                dActiveSlotMolFraction(iSlot,iFirst:iLast) = DMAX1(dMolFraction(iFirst:iLast), 0D0)
                dSumLocal = SUM(dActiveSlotMolFraction(iSlot,iFirst:iLast))
                if (dSumLocal > 1D-300) dActiveSlotMolFraction(iSlot,iFirst:iLast) = &
                    dActiveSlotMolFraction(iSlot,iFirst:iLast) / dSumLocal
            end if
            if (allocated(dActiveSlotSiteFraction)) then
                if (TRIM(cSolnPhaseType(iSolnPhase)) == 'SUBOM') then
                    call SiteFractionFromProductFraction(iSolnPhase, dMolFraction(iFirst:iLast), &
                        dActiveSlotSiteFraction(iSlot,:,:))
                end if
            end if
        end do

        return
    end subroutine RebuildActiveSlotState


    integer function ParentPhaseForDisplay(iDisplayPhase)
        integer, intent(in) :: iDisplayPhase

        integer :: iOrderedPhase

        ParentPhaseForDisplay = iDisplayPhase
        if (.NOT.allocated(iDisorderedPhase)) return

        do iOrderedPhase = 1, MIN(nSolnPhasesSys, SIZE(iDisorderedPhase))
            if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') cycle
            if ((CanonicalCompanionPhase(iOrderedPhase) == iDisplayPhase).OR.&
                (lSUBOMTwoSetCandidateEnabled.AND.&
                (iDisorderedPhase(iOrderedPhase) == iDisplayPhase))) then
                ParentPhaseForDisplay = iOrderedPhase
                return
            end if
        end do

        return
    end function ParentPhaseForDisplay


    integer function CanonicalCompanionPhase(iParentPhase)
        integer, intent(in) :: iParentPhase

        CanonicalCompanionPhase = 0
        if (.NOT.lODPartitionUnifiedActive) return
        if (.NOT.allocated(iODCompanionPhase)) return
        if (.NOT.allocated(iODTopologyClass)) return
        if ((iParentPhase < 1).OR.(iParentPhase > SIZE(iODCompanionPhase))) return
        if ((iODTopologyClass(iParentPhase) < OD_TOPOLOGY_HELPER_STANDALONE).OR.&
            (iODTopologyClass(iParentPhase) > OD_TOPOLOGY_HELPER_ONLY)) return
        if (iODCompanionPhase(iParentPhase) == iParentPhase) return
        CanonicalCompanionPhase = iODCompanionPhase(iParentPhase)

        return
    end function CanonicalCompanionPhase


    integer function ConstituentElementIndex(cName)
        character(*), intent(in) :: cName

        integer :: iElement
        character(8) :: cTarget, cElement

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


    logical function IsVacancyConstituent(cName)
        character(*), intent(in) :: cName

        IsVacancyConstituent = TRIM(UpperName(cName)) == 'VA'

        return
    end function IsVacancyConstituent


    character(8) function UpperName(cName)
        character(*), intent(in) :: cName

        integer :: iChar, iCode, nChar

        UpperName = ' '
        nChar = MIN(LEN_TRIM(cName), LEN(UpperName))
        do iChar = 1, nChar
            iCode = IACHAR(cName(iChar:iChar))
            if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                UpperName(iChar:iChar) = ACHAR(iCode - 32)
            else
                UpperName(iChar:iChar) = cName(iChar:iChar)
            end if
        end do

        return
    end function UpperName


    subroutine SiteFractionFromProductFraction(iSolnPhase, dProductFraction, dSiteOut)
        integer, intent(in) :: iSolnPhase
        real(8), dimension(:), intent(in) :: dProductFraction
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(out) :: dSiteOut

        integer :: iSublPhase, iFirst, iLast, iSpecies, iLocal, iSub, iCon
        real(8) :: dSubSum

        dSiteOut = 0D0
        iSublPhase = iPhaseSublattice(iSolnPhase)
        if (iSublPhase <= 0) return

        iFirst = nSpeciesPhase(iSolnPhase-1) + 1
        iLast  = nSpeciesPhase(iSolnPhase)
        do iSpecies = iFirst, iLast
            iLocal = iSpecies - iFirst + 1
            if (iLocal > SIZE(dProductFraction)) exit
            do iSub = 1, nSublatticePhase(iSublPhase)
                iCon = iConstituentSublattice(iSublPhase,iSub,iLocal)
                if (iCon > 0) dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                    DMAX1(dProductFraction(iLocal), 0D0)
            end do
        end do

        do iSub = 1, nSublatticePhase(iSublPhase)
            dSubSum = SUM(dSiteOut(iSub,1:nConstituentSublattice(iSublPhase,iSub)))
            if (dSubSum > 1D-300) then
                dSiteOut(iSub,1:nConstituentSublattice(iSublPhase,iSub)) = &
                    dSiteOut(iSub,1:nConstituentSublattice(iSublPhase,iSub)) / dSubSum
            end if
        end do

        return
    end subroutine SiteFractionFromProductFraction

end subroutine CoalesceDuplicateSUBOMCompositionSets
