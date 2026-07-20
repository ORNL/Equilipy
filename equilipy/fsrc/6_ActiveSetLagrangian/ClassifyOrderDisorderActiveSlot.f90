!> \brief Classify one ordered-parent composition set by exact sublattice symmetry.
!!
!! \details Uses the parser-declared DIS_PART companion and the ordered
!! parent's constituent topology.  No ordering-degree or composition-distance
!! tolerance participates in the ordered/disordered identity decision.
!
subroutine ClassifyOrderDisorderActiveSlot(iParentPhase, dSiteIn, iClassOut, &
    iCompanionPhaseOut)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ClassifyOrderDisorderActiveSlot.f90
    !> \brief   Classify a canonical order/disorder composition set.
    !> \author  S.Y. Kwon
    !> \date    Jul. 16, 2026
    !> \sa      CoalesceDuplicateSUBOMCompositionSets.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Classified active order/disorder slots from exact structural phase identity.
    !
    ! Purpose:
    ! ========
    !
    !> \details A parser-declared ordered/disordered pair is one physical phase.
    !! The ordered parent owns the thermodynamic state.  Its public display is
    !! the disordered companion exactly when every structurally equivalent
    !! ordering sublattice has identical site fractions.  Invalid topology or
    !! non-finite state is returned as a typed ambiguous/failure fact.
    !
    ! Parameters:
    ! ===========
    !
    !> \param[in] iParentPhase       Ordered SUBOM parent phase id.
    !> \param[in] dSiteIn            Parent-model site fractions for one slot.
    !> \param[out] iClassOut         Typed OD_CANDIDATE_* classification.
    !> \param[out] iCompanionPhaseOut Canonical disordered companion phase id.
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Constituent-name equality defines structurally equivalent ordering
    !   sublattices; phase-name substring inference is never used.
    ! - Site-fraction equality is exact.  A nonzero symmetry-breaking component
    !   is ordered; invalid or non-finite data is not coerced onto either branch.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
    USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

    implicit none

    integer, intent(in) :: iParentPhase
    real(8), dimension(*), intent(in) :: dSiteIn
    integer, intent(out) :: iClassOut, iCompanionPhaseOut

    integer :: iPhaseID, iSubA, iSubB, iCon, nConstituent
    logical :: lFoundOrderingGroup

    iClassOut = OD_CANDIDATE_NOT_EVALUATED
    iCompanionPhaseOut = 0
    if (.NOT.lODPartitionUnifiedActive) return
    if (.NOT.allocated(iODCompanionPhase)) return
    if (.NOT.allocated(iODTopologyClass)) return
    if ((iParentPhase < 1).OR.(iParentPhase > nSolnPhasesSys)) return
    if (iParentPhase > SIZE(iODCompanionPhase)) return
    if ((iODTopologyClass(iParentPhase) < OD_TOPOLOGY_HELPER_STANDALONE).OR.&
        (iODTopologyClass(iParentPhase) > OD_TOPOLOGY_HELPER_ONLY)) return
    if (TRIM(cSolnPhaseType(iParentPhase)) /= 'SUBOM') return

    iCompanionPhaseOut = iODCompanionPhase(iParentPhase)
    if ((iCompanionPhaseOut <= 0).OR.(iCompanionPhaseOut > nSolnPhasesSys).OR.&
        (iCompanionPhaseOut == iParentPhase)) then
        iClassOut = OD_CANDIDATE_AMBIGUOUS_COMPANION
        return
    end if

    iPhaseID = iPhaseSublattice(iParentPhase)
    if (iPhaseID <= 0) then
        iClassOut = OD_CANDIDATE_EVALUATION_FAILED
        return
    end if

    lFoundOrderingGroup = .FALSE.
    do iSubA = 1, nSublatticePhase(iPhaseID) - 1
        nConstituent = nConstituentSublattice(iPhaseID,iSubA)
        if (nConstituent <= 1) cycle
        do iSubB = iSubA + 1, nSublatticePhase(iPhaseID)
            if (.NOT.SameConstituentList(iPhaseID, iSubA, iSubB)) cycle
            lFoundOrderingGroup = .TRUE.
            do iCon = 1, nConstituent
                if ((.NOT.ieee_is_finite(SiteValue(iSubA,iCon))).OR.&
                    (.NOT.ieee_is_finite(SiteValue(iSubB,iCon)))) then
                    iClassOut = OD_CANDIDATE_EVALUATION_FAILED
                    return
                end if
                if (SiteValue(iSubA,iCon) /= SiteValue(iSubB,iCon)) then
                    iClassOut = OD_CANDIDATE_ORDERED
                    return
                end if
            end do
        end do
    end do

    if (lFoundOrderingGroup) then
        iClassOut = OD_CANDIDATE_DISORDERED
    else
        iClassOut = OD_CANDIDATE_AMBIGUOUS_COMPANION
    end if

    return

contains

    real(8) function SiteValue(iSublattice, iConstituent)
        integer, intent(in) :: iSublattice, iConstituent

        SiteValue = dSiteIn(iSublattice + (iConstituent - 1) * nMaxSublatticeSys)

        return
    end function SiteValue

    logical function SameConstituentList(iPhaseIDIn, iSubAIn, iSubBIn)
        integer, intent(in) :: iPhaseIDIn, iSubAIn, iSubBIn
        integer :: iConLocal

        SameConstituentList = .FALSE.
        if (nConstituentSublattice(iPhaseIDIn,iSubAIn) /= &
            nConstituentSublattice(iPhaseIDIn,iSubBIn)) return

        do iConLocal = 1, nConstituentSublattice(iPhaseIDIn,iSubAIn)
            if (TRIM(UpperName(cConstituentNameSUB(iPhaseIDIn,iSubAIn,iConLocal))) /= &
                TRIM(UpperName(cConstituentNameSUB(iPhaseIDIn,iSubBIn,iConLocal)))) return
        end do
        SameConstituentList = .TRUE.

        return
    end function SameConstituentList


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

end subroutine ClassifyOrderDisorderActiveSlot
