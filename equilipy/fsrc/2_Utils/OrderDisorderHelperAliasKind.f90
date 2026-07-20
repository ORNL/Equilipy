!> \brief Classify known nonphysical order/disorder helper aliases.
!!
!! \details Returns a small integer identifier for helper phases that exist
!! only to supply a DIS_PART/DISORDER_PART companion.  Physical disordered
!! phases such as BCC_A2 and FCC_A1 return zero.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    OrderDisorderHelperAliasKind.f90
!> \brief   Classify order/disorder helper aliases.
!> \author  S.Y. Kwon
!> \date    Jul. 03, 2026
!> \sa      CheckSystem.f90
!> \sa      RecordGEMIterationDiagnostics.f90
!> \sa      RegisterSUBOMTwoSetCandidateRows.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Classified order/disorder helper aliases from structural partition topology while rejecting ambiguous cases.
    !
    !
! Purpose:
! ========
!
!> \details Centralize the identification of nonphysical order/disorder helper
!! aliases.  These phases should not be treated as ordinary physical
!! disordered phases when the minimizer decides whether same-parent SUBOM
!! composition-set evidence is allowed.
!
!
! Required input variables:
! =========================
!
! iPhaseIndex       System solution-phase index.
! cSolnPhaseName    Current system solution-phase names.
!
!
! Output/updated variables:
! =========================
!
! OrderDisorderHelperAliasKind  Zero for ordinary phases; positive helper kind
!                               for known nonphysical aliases.
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
! CheckSystem
! RecordGEMIterationDiagnostics
! RegisterSUBOMTwoSetCandidateRows
!
!
! Numerical assumptions:
! ======================
!
! - This is an identity helper only; it must not change thermodynamic values.
! - Positive values indicate helper aliases, not physical companion phases.
!
!-------------------------------------------------------------------------------------------------------------



integer function OrderDisorderHelperAliasKind(iPhaseIndex)
    USE ModuleThermo
    USE ModuleGEMSolver, ONLY: OD_TOPOLOGY_HELPER_STANDALONE, &
        OD_TOPOLOGY_HELPER_ONLY, OD_TOPOLOGY_LEGACY_NO_METADATA

    implicit none

    integer, intent(in) :: iPhaseIndex
    integer             :: iOrderedPhase
    character(30)       :: cUpperName

    OrderDisorderHelperAliasKind = 0
    if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) return

    if (allocated(iODTopologyClass).AND.allocated(iDisorderedPhase)) then
        do iOrderedPhase = 1, MIN(nSolnPhasesSys, SIZE(iODTopologyClass), &
            SIZE(iDisorderedPhase))
            if ((iODTopologyClass(iOrderedPhase) /= OD_TOPOLOGY_HELPER_STANDALONE).AND.&
                (iODTopologyClass(iOrderedPhase) /= OD_TOPOLOGY_HELPER_ONLY)) cycle
            if (iDisorderedPhase(iOrderedPhase) /= iPhaseIndex) cycle
            OrderDisorderHelperAliasKind = iOrderedPhase
            return
        end do

        ! A typed non-legacy graph is authoritative.  Do not reinterpret an
        ! ordinary phase through the historical spelling fallback.
        if (ANY((iODTopologyClass > 0).AND.&
            (iODTopologyClass /= OD_TOPOLOGY_LEGACY_NO_METADATA))) return
    end if

    ! Metadata-free DAT calculations retain their historical helper spelling
    ! behavior exactly.  This branch is not used by typed TDB partition paths.
    cUpperName = UpperOrderDisorderAliasName(cSolnPhaseName(iPhaseIndex))
    select case (TRIM(cUpperName))
    case ('A2_BCC')
        OrderDisorderHelperAliasKind = 1
    case ('A1_FCC')
        OrderDisorderHelperAliasKind = 2
    end select

    return

contains

    character(30) function UpperOrderDisorderAliasName(cName)
        character(*), intent(in) :: cName
        integer                  :: iChar, iCode, nChar

        UpperOrderDisorderAliasName = ' '
        nChar = MIN(LEN_TRIM(cName), LEN(UpperOrderDisorderAliasName))
        do iChar = 1, nChar
            iCode = IACHAR(cName(iChar:iChar))
            if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                UpperOrderDisorderAliasName(iChar:iChar) = ACHAR(iCode - 32)
            else
                UpperOrderDisorderAliasName(iChar:iChar) = cName(iChar:iChar)
            end if
        end do

        return
    end function UpperOrderDisorderAliasName

end function OrderDisorderHelperAliasKind


!> \brief Report whether a phase is a nonphysical order/disorder helper alias.
logical function IsOrderDisorderHelperAliasPhase(iPhaseIndex)
    implicit none

    integer, intent(in) :: iPhaseIndex
    integer             :: OrderDisorderHelperAliasKind

    IsOrderDisorderHelperAliasPhase = OrderDisorderHelperAliasKind(iPhaseIndex) > 0

    return

end function IsOrderDisorderHelperAliasPhase
