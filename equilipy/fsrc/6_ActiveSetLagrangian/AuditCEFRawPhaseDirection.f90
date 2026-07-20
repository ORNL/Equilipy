!> \brief Classify raw CEF phase-removal Newton directions.
!!
!! \details Checks the pre-line-search CEF Newton phase-amount direction and
!! records whether a full Newton step would make an active phase amount
!! negative.  The event is diagnostic at this stage; RunLagrangianGEM may use
!! it only after the positivity-preserving CEF line search fails to find descent.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    AuditCEFRawPhaseDirection.f90
!> \brief   Audit raw CEF phase-removal Newton directions.
!> \author  S.Y. Kwon
!> \date    Jun. 25, 2026
!> \sa      RunLagrangianGEM.f90
!> \sa      GEMNewtonCEF.f90
!> \sa      GEMLineSearchCEF.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Classified raw CEF phase directions against phase-boundary complementarity during active-set repair.
    !
    !
! Purpose:
! ========
!
!> \details The purpose of this routine is to preserve visibility of raw
!! phase-amount boundary directions without deleting phases from an undamped
!! Newton step.  Lambda correction and line-search clipping are useful for
!! keeping site fractions positive, but an undamped raw direction alone is not
!! sufficient evidence that a phase should leave the active set.
!
!
! Required input variables:
! =========================
!
! dMolesPhase                 Current active-set phase amounts.
! dGEMCEFPhaseDirection       Phase-amount Newton direction from GEMNewtonCEF.
! iGEMCEFPhaseSlot            Map from CEF phase variable to active assemblage slot.
! nGEMCEFPhaseVariables       Number of active phase-amount variables.
! iterGlobal                  Current Lagrangian iteration index.
!
!
! Output/updated variables:
! =========================
!
! dGEMLSRawPhaseMoles         Raw full-step phase amounts for diagnostics.
! iGEMLineSearchNegativePhaseCount
!                             Number of raw full-step negative phase amounts.
! iGEMRawNegativePhase*       Most negative raw phase-removal target for delayed repair.
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
! RunLagrangianGEM            Calls after GEMNewtonCEF and before GEMLineSearchCEF.
!
!
! Numerical assumptions:
! ======================
!
! - Postprocess heat-capacity solves must remain fixed assemblage and therefore
!   must not call phase repair.
! - A raw full-step negative phase amount is not an active-set event by itself.
!   The damped CEF line search first tries a feasible step.  Only if that line
!   search fails to find descent can RunLagrangianGEM use this audit to remove
!   the selected active phase during normal Lagrangian solves.  During
!   PEA-polish solves, the same evidence rejects the polish and returns
!   active-set control to PEA.
!
!-------------------------------------------------------------------------------------------------------------



subroutine AuditCEFRawPhaseDirection

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: iPhaseVar, iSlot
    real(8) :: dRawPhaseAmount
    real(8) :: dMostNegativeRawPhaseAmount
    real(8) :: dComplementarityTolerance

    if (lPostProcess) return
    if (.NOT.lGEMCEFSiteDirectionActive) return
    if (nGEMCEFPhaseVariables <= 0) return

    iGEMLineSearchNegativePhaseCount = 0
    iGEMRawNegativePhaseSlot = 0
    iGEMRawNegativePhaseSoln = 0
    iGEMRawNegativePhaseSpecies = 0
    iGEMRawNegComp = 0
    dGEMRawNegativePhaseAmount = 0D0
    dGEMRawNegativePhaseDirection = 0D0
    dGEMRawNegPhaseResidual = 0D0
    dMostNegativeRawPhaseAmount = 0D0
    dComplementarityTolerance = 1D-10
    dGEMLineSearchMinRawPhaseMoles = MINVAL(dMolesPhase)
    if (allocated(dGEMLSRawPhaseMoles)) dGEMLSRawPhaseMoles = dMolesPhase

    do iPhaseVar = 1, nGEMCEFPhaseVariables
        iSlot = iGEMCEFPhaseSlot(iPhaseVar)
        if ((iSlot < 1).OR.(iSlot > nElements)) cycle

        dRawPhaseAmount = dMolesPhase(iSlot) + dGEMCEFPhaseDirection(iPhaseVar)
        if (allocated(dGEMLSRawPhaseMoles)) dGEMLSRawPhaseMoles(iSlot) = dRawPhaseAmount
        dGEMLineSearchMinRawPhaseMoles = DMIN1(dGEMLineSearchMinRawPhaseMoles, dRawPhaseAmount)
        if (dRawPhaseAmount < -dTolerance(8)) then
            iGEMLineSearchNegativePhaseCount = iGEMLineSearchNegativePhaseCount + 1
            if ((iGEMRawNegativePhaseSlot == 0).OR.(dRawPhaseAmount < dMostNegativeRawPhaseAmount)) then
                iGEMRawNegativePhaseSlot = iSlot
                iGEMRawNegativePhaseSoln = iGEMCEFPhaseSoln(iPhaseVar)
                iGEMRawNegativePhaseSpecies = iGEMCEFPhaseSpecies(iPhaseVar)
                dGEMRawNegativePhaseAmount = dRawPhaseAmount
                dGEMRawNegativePhaseDirection = dGEMCEFPhaseDirection(iPhaseVar)
                if (allocated(dGEMCEFPhaseResidual)) then
                    if (iPhaseVar <= SIZE(dGEMCEFPhaseResidual)) then
                        dGEMRawNegPhaseResidual = dGEMCEFPhaseResidual(iPhaseVar)
                    end if
                end if
                if (dGEMRawNegPhaseResidual >= -dComplementarityTolerance) then
                    iGEMRawNegComp = 1
                else
                    iGEMRawNegComp = -1
                end if
                dMostNegativeRawPhaseAmount = dRawPhaseAmount
            end if
        end if
    end do

    return

end subroutine AuditCEFRawPhaseDirection
