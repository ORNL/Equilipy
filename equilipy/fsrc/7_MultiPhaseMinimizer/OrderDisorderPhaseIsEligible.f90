!> \brief Return whether an ordered SUBOM phase is eligible as an active PEA witness.
!!
!! \details Ordered/disordered CEF phases need at least two active real
!! constituents to represent an ordering degree of freedom.  If a mapped
!! SUBOM phase contains only one active non-vacancy element, its disordered
!! partner is the physical phase-selection witness and the ordered phase should
!! not be introduced by PEA as a separate pseudo-compound candidate.
!!



logical function OrderDisorderPhaseIsEligible(iSolnPhase)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    OrderDisorderPhaseIsEligible.f90
    !> \brief   Check whether a SUBOM ordered phase has active ordering degrees of freedom.
    !> \author  S.Y. Kwon
    !> \date    Jul. 01, 2026
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/01/2026      S.Y. Kwon           Original helper for filtering one-element ordered phases in PEA.
    !   07/01/2026      S.Y. Kwon           Counted active elements from screened phase stoichiometry instead
    !                                       of CEF constituent names.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details PEA compares solution-phase local minima as pseudo-compounds.
    !! A mapped order/disorder phase with only one active real element has no
    !! ordering variable to minimize and should be represented by its
    !! disordered partner.  This avoids selecting a numerically singular
    !! ordered phase on unary boundaries such as pure Fe in BCC/B2 systems.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! iSolnPhase            Absolute solution phase index in the active system.
    ! cSolnPhaseType        Solution model type.
    ! iDisorderedPhase      Ordered-to-disordered phase mapping.
    ! dStoichSpecies        Screened active-system endmember stoichiometry.
    ! nSpeciesPhase         Active-system endmember boundaries for each solution phase.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! OrderDisorderPhaseIsEligible  FALSE only for mapped SUBOM phases with fewer than
    !                               two active non-vacancy elements.
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
    ! CompInitMinSolnPoint   Filters initial PEA solution candidates.
    ! CompMinSolnPoint       Filters refreshed PEA solution candidates.
    ! CompDrivingForceAll    Filters diagnostic all-phase driving-force checks.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - This is an active-system test; CheckSystem has already screened
    !   inactive elements and endmembers from dStoichSpecies.
    ! - Vacancy-only endmembers have zero active stoichiometry and therefore
    !   do not count as ordering elements.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo

    implicit none

    integer, intent(in) :: iSolnPhase
    integer             :: iFirst, iLast, iSpecies, iElement
    integer             :: nActiveRealElements
    logical, dimension(nElements) :: lRealElementPresent

    OrderDisorderPhaseIsEligible = .TRUE.

    if (iSolnPhase < 1) return
    if (iSolnPhase > nSolnPhasesSys) return
    if (TRIM(cSolnPhaseType(iSolnPhase)) /= 'SUBOM') return
    if (.NOT.allocated(iDisorderedPhase)) return
    if (iSolnPhase > SIZE(iDisorderedPhase)) return
    if (iDisorderedPhase(iSolnPhase) <= 0) return

    lRealElementPresent = .FALSE.
    iFirst = nSpeciesPhase(iSolnPhase-1) + 1
    iLast  = nSpeciesPhase(iSolnPhase)
    do iSpecies = iFirst, iLast
        do iElement = 1, nElements
            if (DABS(dStoichSpecies(iSpecies,iElement)) > 0D0) then
                lRealElementPresent(iElement) = .TRUE.
            end if
        end do
    end do

    nActiveRealElements = COUNT(lRealElementPresent)
    OrderDisorderPhaseIsEligible = nActiveRealElements >= 2

    return

end function OrderDisorderPhaseIsEligible
