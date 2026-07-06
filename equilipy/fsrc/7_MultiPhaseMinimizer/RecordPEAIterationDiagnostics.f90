!> \brief Record one CheckPhaseAssemblage/PEA iteration for active-set audits.
!!
!! \details Stores the Leveling assemblage, embedded Lagrangian polish status,
!! elemental potentials, and convergence metrics for the current PEA iteration.
!!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    RecordPEAIterationDiagnostics.f90
!> \brief   Record per-iteration PEA diagnostics.
!> \author  S.Y. Kwon
!> \date    Jun. 26, 2026
!> \sa      CheckPhaseAssemblage.f90
!> \sa      PEALagrangianPolish.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   06/26/2026      S.Y. Kwon           Original per-PEA-iteration diagnostic recorder.
!   07/04/2026      S.Y. Kwon           Recorded passive C4-a PEA gate, direct-witness, and
!                                       repeated-handoff exit counters.
!
!
! Purpose:
! ========
!
!> \details The purpose of this routine is to make PEA active-set evolution
!! inspectable without changing the minimization path.  Each recorded row
!! captures the Leveling assemblage, whether embedded Lagrangian polishing was
!! attempted and accepted, the polish rejection reason, and the current
!! driving-force and elemental-potential convergence measures.
!
!
! Required input variables:
! =========================
!
!> \param[in]   iterLevel          Number of Leveling iterations used in the
!!                                 current PEA outer iteration.
!> \param[in]   iPolishAttempted   One when PEALagrangianPolish actually ran a
!!                                 Lagrangian attempt in this PEA iteration.
!
! iterPEA                       Current PEA outer iteration.
! iAssemblage                   Current Leveling assemblage after the Leveling pass.
! dElementPotential             Elemental potentials used for the next driving-force pass.
! dMinPhasePotential            Minimum phase potential after solution minima are refreshed.
! dMaxElementPotential          Elemental-potential change from the previous PEA iteration.
! dGEMFunctionNorm              Current GEM residual norm after the refreshed driving-force pass.
!                                Compare with dPEALagrangianPolishNormAfter for the embedded polish norm.
! lPEALagrangianPolishAccepted  True when the embedded Lagrangian polish was accepted.
! iPEALagrangianPolishReason    Reason code from the latest embedded Lagrangian polish.
! dPEALagrangianPolish*         Norm and potential-change diagnostics from the polish.
!
!
! Output/updated variables:
! =========================
!
! nPEARecorded                         Highest PEA iteration recorded in the current run.
! iPEAAssemblageHist                   PEA assemblage history after Leveling.
! iPEALevelIterHist                    Leveling iteration count history.
! iPEAPolishAttemptHist                Embedded Lagrangian attempt flags.
! iPEAPolishAcceptedHist               Embedded Lagrangian acceptance flags.
! iPEAPolishReasonHist                 Embedded Lagrangian rejection/exit reason codes.
! iPEAPolishIterGlobalHist             Lagrangian iteration count used by each polish.
! iPEAGate*Hist                        C4-a passive histories for the polish-entry gate.
! iPEADirect*Hist                      C4-a passive histories for direct inactive-witness handoffs.
! iPEARepeatExit*Hist                  C4-a passive histories for repeated-handoff exits.
! dPEAMinPhasePotentialHist            Minimum refreshed phase potential per PEA iteration.
! dPEAMaxElementPotentialHist          Elemental-potential plateau metric per PEA iteration.
! dPEAGEMFunctionNormHist              GEM norm after the refreshed driving-force pass.
! dPEAPolishNormHist                   Embedded Lagrangian polish norm per PEA iteration.
! dPEAPolishPotentialChangeHist        Embedded polish elemental-potential change.
! dPEAElementPotentialHist             Elemental-potential vector per PEA iteration.
! iPEALagrangianHandoffHist            Translated Lagrangian handoff assemblage, recorded by
!                                      PEALagrangianPolish after Level2Lagrange.
! iPEALagrangianHandoffRepeatHist      One when the translated handoff repeats an earlier PEA handoff.
! iPEALagrangianHandoffFirstIterHist   First earlier PEA iteration with the same translated handoff.
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
! CheckPhaseAssemblage       Calls once per completed PEA outer iteration.
!
!
! Numerical assumptions:
! ======================
!
! - This routine is diagnostic only.  It must not change active-set selection,
!   convergence criteria, phase amounts, or elemental potentials.
! - Histories are capped by iterPEAMax.  Iterations beyond that cap are not
!   recorded because CheckPhaseAssemblage should not exceed the same cap.
!
!-------------------------------------------------------------------------------------------------------------
subroutine RecordPEAIterationDiagnostics(iterLevel, iPolishAttempted)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in) :: iterLevel, iPolishAttempted
    integer             :: iIter

    iIter = iterPEA
    if ((iIter < 1).OR.(iIter > iterPEAMax)) return

    nPEARecorded = MAX(nPEARecorded, iIter)

    if (allocated(iPEAAssemblageHist)) iPEAAssemblageHist(:,iIter) = iAssemblage
    if (allocated(iPEALevelIterHist)) iPEALevelIterHist(iIter) = iterLevel
    if (allocated(iPEAPolishAttemptHist)) iPEAPolishAttemptHist(iIter) = iPolishAttempted
    if (allocated(iPEAPolishAcceptedHist)) then
        if ((iPolishAttempted > 0).AND.lPEALagrangianPolishAccepted) then
            iPEAPolishAcceptedHist(iIter) = 1
        else
            iPEAPolishAcceptedHist(iIter) = 0
        end if
    end if
    if (allocated(iPEAPolishReasonHist)) iPEAPolishReasonHist(iIter) = iPEALagrangianPolishReason
    if (allocated(iPEAPolishIterGlobalHist)) iPEAPolishIterGlobalHist(iIter) = iPEALagrangianPolishIterGlobal
    if (allocated(iPEAGatePassHist)) iPEAGatePassHist(iIter) = iPEAGatePass
    if (allocated(iPEAGateBlockHist)) iPEAGateBlockHist(iIter) = iPEAGateBlock
    if (allocated(iPEAGatePotBlockHist)) then
        iPEAGatePotBlockHist(iIter) = iPEAGatePotBlock
    end if
    if (allocated(iPEAGateDFBlockHist)) then
        iPEAGateDFBlockHist(iIter) = iPEAGateDFBlock
    end if
    if (allocated(iPEADirectHandoffHist)) then
        iPEADirectHandoffHist(iIter) = iPEADirectHandoff
    end if
    if (allocated(iPEADirectPhaseHist)) iPEADirectPhaseHist(iIter) = iPEADirectPhase
    if (allocated(iPEARepeatExitHist)) then
        iPEARepeatExitHist(iIter) = iPEARepeatExit
    end if
    if (allocated(iPEARepeatExitReasonHist)) then
        iPEARepeatExitReasonHist(iIter) = iPEARepeatExitReason
    end if

    if (allocated(dPEAMinPhasePotentialHist)) dPEAMinPhasePotentialHist(iIter) = dMinPhasePotential
    if (allocated(dPEAMaxElementPotentialHist)) dPEAMaxElementPotentialHist(iIter) = dMaxElementPotential
    if (allocated(dPEAGEMFunctionNormHist)) dPEAGEMFunctionNormHist(iIter) = dGEMFunctionNorm
    if (allocated(dPEAPolishNormHist)) dPEAPolishNormHist(iIter) = dPEALagrangianPolishNormAfter
    if (allocated(dPEAPolishPotentialChangeHist)) then
        dPEAPolishPotentialChangeHist(iIter) = dPEALagrangianPolishPotentialChange
    end if
    if (allocated(dPEAElementPotentialHist)) dPEAElementPotentialHist(:,iIter) = dElementPotential

    return

end subroutine RecordPEAIterationDiagnostics
