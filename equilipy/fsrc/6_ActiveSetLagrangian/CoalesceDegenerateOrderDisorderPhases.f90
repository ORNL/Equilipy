!> \brief Coalesce degenerate random order/disorder active phases.
!!
!! \details If a SUBOM phase has converged to the random-state composition of
!! its DIS_PART companion, both active phases represent the same tangent-plane
!! point.  The helper amount is merged into the ordered phase so the final
!! active set follows the ordered-phase representation used by DIS_PART
!! databases.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    CoalesceDegenerateOrderDisorderPhases.f90
!> \brief   Merge random-state DIS_PART duplicates after CEF Lagrangian convergence.
!> \author  S.Y. Kwon
!> \date    Jun. 27, 2026
!> \sa      RunLagrangianGEM.f90
!> \sa      GEMNewtonCEF.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Coalesced ordered and disordered phases when their equilibrium constitutions become physically degenerate.
    !
    !
! Purpose:
! ========
!
!> \details Order/disorder databases use a DIS_PART companion to define the
!! random state of an ordered SUBOM phase.  At zero order parameter the ordered
!! phase and the DIS_PART companion are a representation duplicate, not a true
!! two-phase assemblage.  This routine detects that degeneracy from active
!! site fractions and bulk compositions, merges the companion amount into the
!! ordered phase, and compacts the active set.  It does not fire for genuinely
!! ordered states or for finite ordered/disordered coexistence with distinct
!! compositions.
!
!
! Required input variables:
! =========================
!
! iAssemblage                 Current converged Lagrangian active set.
! dMolesPhase                 Current active phase amounts.
! dMolFraction                Current phase-local endmember fractions.
! iDisorderedPhase            Ordered phase to DIS_PART companion map.
! cSolnPhaseType              Used to restrict the ordered side to SUBOM.
!
!
! Output/updated variables:
! =========================
!
!> \param[out]  lChanged      True when one order/disorder pair was coalesced.
!
! iAssemblage                 Compacted active set after the companion removal.
! dMolesPhase                 Updated phase amounts.
! dMolesSpecies               Rebuilt species amounts for retained phases.
! dMolFraction                Ordered phase random-state endmember fractions.
! nSolnPhases/lSolnPhases     Updated active solution count and flags.
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
! RunLagrangianGEM            Calls after CEF fixed-active-set convergence.
!
!
! Numerical assumptions:
! ======================
!
! - The routine is representation cleanup after convergence, not phase
!   discovery.  It only merges active order/disorder pairs that already have
!   matching normalized bulk composition.
! - The ordered phase must be at the random state implied by the companion
!   composition on every non-vacancy ordered sublattice.
! - Fixed all-vacancy sublattices are ignored.  Mobile vacancy sublattices are
!   not coalesced by this rule.
!
!-------------------------------------------------------------------------------------------------------------



subroutine CoalesceDegenerateOrderDisorderPhases(lChanged)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lChanged

    integer :: iOrderedPhase, iHelperPhase, iOrderedSlot, iHelperSlot
    real(8) :: dCompositionDiff
    real(8), dimension(nElements) :: dOrderedComposition, dHelperComposition

    lChanged = .FALSE.
    if (.NOT.allocated(iDisorderedPhase)) return
    if (nSolnPhases <= 1) return

    do iOrderedPhase = 1, nSolnPhasesSys
        if (iOrderedPhase > SIZE(iDisorderedPhase)) exit
        if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') cycle

        iHelperPhase = iDisorderedPhase(iOrderedPhase)
        if (iHelperPhase <= 0) cycle
        if (iHelperPhase > nSolnPhasesSys) cycle

        call FindActiveSolutionSlot(iOrderedPhase, iOrderedSlot)
        call FindActiveSolutionSlot(iHelperPhase, iHelperSlot)
        if ((iOrderedSlot == 0).OR.(iHelperSlot == 0)) cycle

        call ComputeNormalizedSolutionComposition(iOrderedPhase, dOrderedComposition)
        call ComputeNormalizedSolutionComposition(iHelperPhase, dHelperComposition)
        dCompositionDiff = SUM(DABS(dOrderedComposition - dHelperComposition)) / DFLOAT(MAX(1,nElements))
        if (dCompositionDiff > 1D-6) cycle
        if (.NOT.OrderedPhaseIsRandomState(iOrderedPhase, dHelperComposition)) cycle

        call MergeHelperIntoOrderedPhase(iOrderedPhase, iOrderedSlot, iHelperPhase, iHelperSlot, &
            dHelperComposition)
        lChanged = .TRUE.
        exit
    end do

    return

contains

    subroutine FindActiveSolutionSlot(iSolnPhaseIndex, iSlotOut)
        integer, intent(in)  :: iSolnPhaseIndex
        integer, intent(out) :: iSlotOut

        integer :: iSoln, iSlot

        iSlotOut = 0
        do iSoln = 1, nSolnPhases
            iSlot = nElements - iSoln + 1
            if (iAssemblage(iSlot) == -iSolnPhaseIndex) then
                iSlotOut = iSlot
                return
            end if
        end do

        return
    end subroutine FindActiveSolutionSlot


    subroutine ComputeNormalizedSolutionComposition(iSolnPhaseIndex, dCompositionOut)
        integer, intent(in) :: iSolnPhaseIndex
        real(8), dimension(nElements), intent(out) :: dCompositionOut

        integer :: iSpecies, iFirst, iLast, iElement
        real(8) :: dNorm

        dCompositionOut = 0D0
        iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
        iLast  = nSpeciesPhase(iSolnPhaseIndex)
        do iSpecies = iFirst, iLast
            do iElement = 1, nElements
                dCompositionOut(iElement) = dCompositionOut(iElement) + &
                    DMAX1(dMolFraction(iSpecies), 0D0) * dStoichSpecies(iSpecies,iElement) / &
                    DFLOAT(iParticlesPerMole(iSpecies))
            end do
        end do

        dNorm = SUM(DABS(dCompositionOut))
        if (dNorm > 1D-300) dCompositionOut = dCompositionOut / dNorm

        return
    end subroutine ComputeNormalizedSolutionComposition


    logical function OrderedPhaseIsRandomState(iOrderedPhase, dHelperComposition)
        integer, intent(in) :: iOrderedPhase
        real(8), dimension(nElements), intent(in) :: dHelperComposition

        integer :: iSublPhase, iSub, iCon, iElement
        real(8) :: dSite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8) :: dSiteTol, dMobileVacancyFraction

        OrderedPhaseIsRandomState = .FALSE.
        iSublPhase = iPhaseSublattice(iOrderedPhase)
        if (iSublPhase <= 0) return

        call BuildSolutionSiteFractions(iOrderedPhase, iSublPhase, dSite)
        dSiteTol = 1D-6

        do iSub = 1, nSublatticePhase(iSublPhase)
            dMobileVacancyFraction = 0D0
            do iCon = 1, nConstituentSublattice(iSublPhase,iSub)
                if (IsVacancyConstituent(cConstituentNameSUB(iSublPhase,iSub,iCon))) then
                    if (nConstituentSublattice(iSublPhase,iSub) > 1) then
                        dMobileVacancyFraction = dMobileVacancyFraction + dSite(iSub,iCon)
                    end if
                    cycle
                end if

                iElement = ConstituentElementIndex(cConstituentNameSUB(iSublPhase,iSub,iCon))
                if (iElement <= 0) return
                if (DABS(dSite(iSub,iCon) - dHelperComposition(iElement)) > dSiteTol) return
            end do
            if (dMobileVacancyFraction > dSiteTol) return
        end do

        OrderedPhaseIsRandomState = .TRUE.

        return
    end function OrderedPhaseIsRandomState


    subroutine BuildSolutionSiteFractions(iSolnPhaseIndex, iSublPhase, dSiteOut)
        integer, intent(in) :: iSolnPhaseIndex, iSublPhase
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(out) :: dSiteOut

        integer :: iSpecies, iFirst, iLast, iLocal, iSub, iCon
        real(8) :: dSubSum

        dSiteOut = 0D0
        iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
        iLast  = nSpeciesPhase(iSolnPhaseIndex)
        do iSpecies = iFirst, iLast
            iLocal = iSpecies - iFirst + 1
            do iSub = 1, nSublatticePhase(iSublPhase)
                iCon = iConstituentSublattice(iSublPhase,iSub,iLocal)
                if (iCon > 0) dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                    DMAX1(dMolFraction(iSpecies), 0D0)
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
    end subroutine BuildSolutionSiteFractions


    subroutine MergeHelperIntoOrderedPhase(iOrderedPhase, iOrderedSlot, iHelperPhase, iHelperSlot, &
        dHelperComposition)
        integer, intent(in) :: iOrderedPhase, iOrderedSlot, iHelperPhase, iHelperSlot
        real(8), dimension(nElements), intent(in) :: dHelperComposition

        integer :: iFirst, iLast, nLocal
        real(8) :: dMergedAmount
        real(8), dimension(:), allocatable :: dProjectedFraction

        iFirst = nSpeciesPhase(iOrderedPhase-1) + 1
        iLast  = nSpeciesPhase(iOrderedPhase)
        nLocal = iLast - iFirst + 1
        if (nLocal <= 0) return

        dMergedAmount = dMolesPhase(iOrderedSlot) + dMolesPhase(iHelperSlot)
        allocate(dProjectedFraction(nLocal))
        call ProjectRandomCompositionToOrderedFractions(iOrderedPhase, dHelperComposition, dProjectedFraction)

        dMolesPhase(iOrderedSlot) = dMergedAmount
        dMolFraction(iFirst:iLast) = dProjectedFraction
        dMolesSpecies(iFirst:iLast) = dMergedAmount * dProjectedFraction

        iFirst = nSpeciesPhase(iHelperPhase-1) + 1
        iLast  = nSpeciesPhase(iHelperPhase)
        dMolesSpecies(iFirst:iLast) = 0D0

        call CompactAfterRemovingSlot(iHelperSlot)

        deallocate(dProjectedFraction)

        return
    end subroutine MergeHelperIntoOrderedPhase


    subroutine ProjectRandomCompositionToOrderedFractions(iOrderedPhase, dCompositionIn, dProjectedFraction)
        integer, intent(in) :: iOrderedPhase
        real(8), dimension(nElements), intent(in) :: dCompositionIn
        real(8), dimension(:), intent(out) :: dProjectedFraction

        integer :: iSpecies, iFirst, iLast, iLocal, iSublPhase, iSub, iCon, iElement
        real(8) :: dProduct, dNorm

        dProjectedFraction = 0D0
        iSublPhase = iPhaseSublattice(iOrderedPhase)
        iFirst = nSpeciesPhase(iOrderedPhase-1) + 1
        iLast  = nSpeciesPhase(iOrderedPhase)

        do iSpecies = iFirst, iLast
            iLocal = iSpecies - iFirst + 1
            if (iLocal > SIZE(dProjectedFraction)) exit

            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iSublPhase)
                iCon = iConstituentSublattice(iSublPhase,iSub,iLocal)
                if (iCon <= 0) cycle

                if (IsVacancyConstituent(cConstituentNameSUB(iSublPhase,iSub,iCon))) then
                    if (nConstituentSublattice(iSublPhase,iSub) == 1) cycle
                    dProduct = 0D0
                    exit
                end if

                iElement = ConstituentElementIndex(cConstituentNameSUB(iSublPhase,iSub,iCon))
                if (iElement <= 0) then
                    dProduct = 0D0
                    exit
                end if
                dProduct = dProduct * DMAX1(dCompositionIn(iElement), 0D0)
            end do
            dProjectedFraction(iLocal) = dProduct
        end do

        dNorm = SUM(dProjectedFraction)
        if (dNorm > 1D-300) dProjectedFraction = dProjectedFraction / dNorm

        return
    end subroutine ProjectRandomCompositionToOrderedFractions


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

end subroutine CoalesceDegenerateOrderDisorderPhases
