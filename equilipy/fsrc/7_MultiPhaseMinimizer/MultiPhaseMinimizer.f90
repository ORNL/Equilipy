!> \brief Orchestrate Leveling, fixed-active-set Lagrangian GEM, and PEA repair.
!!
!! \details Runs the high-level minimization sequence used by GEMSolver.  The
!! routine starts from classical Leveling, checks the active set with
!! CheckPhaseAssemblage/PEA, then uses Lagrangian GEM and event-triggered PEA
!! repair to finish the active-set solve.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    MultiPhaseMinimizer.f90
    !> \brief   Orchestrate the high-level Gibbs energy minimization workflow.
    !> \author  S.Y. Kwon
    !> \date    Jun. 24, 2026
    !> \sa      RunLeveling.f90
    !> \sa      RunLagrangianGEM.f90
    !> \sa      CheckPhaseAssemblage.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
!   06/24/2026      S.Y. Kwon           Documented the staged minimizer workflow.
!   06/24/2026      S.Y. Kwon           Stopped repeated PEA repair when PEA returns the same assemblage.
!   06/24/2026      S.Y. Kwon           Restored initial Leveling-to-PEA assemblage refinement before Lagrangian GEM.
!   06/28/2026      S.Y. Kwon           Skipped PEA for unary active-component systems and converted
!                                        Leveling directly to the Lagrangian state.
!   06/29/2026      S.Y. Kwon           Added initial Leveling-to-Lagrangian polish before PEA so
!                                        order-disorder branch information can refine the first PEA plane.
!   06/29/2026      S.Y. Kwon           Restricted initial Leveling-to-Lagrangian polish to
!                                        order/disorder seeds to avoid damaging ordinary PEA phase discovery.
!   07/01/2026      S.Y. Kwon           Removed the multi-component pre-PEA Lagrangian polish so PEA sees
!                                        the original classical Leveling plane.
!   07/02/2026      S.Y. Kwon           Removed stale initial-polish helper routines that were no longer
!                                        called by the production minimizer flow.
!   07/05/2026      S.Y. Kwon           Added passive substage timing buckets for Scheil warm-start
!                                        expected-value census.
!   07/05/2026      S.Y. Kwon           Reset passive inactive-candidate certification timing for WS-a2.
!   07/05/2026      S.Y. Kwon           Exposed final unconverged GEM exits through iGEMExitStatus even
!                                        when the last internal reason was stagnation or repeated repair.
    !
    ! Purpose:
    ! ========
    !
    !> \details This routine is the production minimizer tree:
    !! Leveling -> CheckPhaseAssemblage/PEA
    !! -> Lagrangian GEM -> event-triggered PEA repair -> Lagrangian GEM retry.  PEA is expensive because it subminimizes
    !! candidate solution phases.  When PEA returns the same active assemblage
    !! and the following Lagrangian retry still requests repair, repeating PEA
    !! does not add new thermodynamic information; the current state is left for
    !! the final norm/error handling instead of looping to the hard repair cap.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! dMolesElement      Normalized bulk composition prepared by InitThermo.
    ! dLeveling*         Leveling arrays prepared by RunLeveling/InitGEMSolver.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! iAssemblage        Active phase assemblage after Leveling, PEA, and Lagrangian retries.
    ! dMolesPhase        Active phase amounts after the final retry.
    ! lConverged         True when Lagrangian GEM accepts the active set.
    ! lPhaseChange       True when the final retry still requests active-set repair.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! RunLeveling                    Builds the initial active-set candidate.
    ! RunLagrangianGEM               Refines the current fixed active set.
    ! CheckPhaseAssemblage           Repairs active-set candidates with PEA.
    !
    ! Primary callers:
    ! ================
    !
    ! GEMSolver                      Delegates the complete GEM workflow here.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Classical Leveling only chooses among endmembers and compounds.
    !   Multi-component systems pass that original Leveling plane directly to
    !   PEA so missing phases such as SIGMA can still be discovered.
    ! - Unary active-component systems have no PEA active-set search.  They are
    !   converted from Leveling directly to Lagrangian GEM.
    ! - CheckPhaseAssemblage is useful only when it changes the active set or
    !   produces a state that lets the next Lagrangian solve converge.
    ! - If PEA returns the same assemblage on a repeated repair pass, another
    !   identical repair pass is treated as no progress rather than a new phase
    !   selection event.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine MultiPhaseMinimizer

    USE ModuleThermoIO
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none
    integer :: iterAssemblageRepair
    integer, parameter :: iterAssemblageRepairMax = 100
    integer, dimension(nElements) :: iAssemblageBeforePEA, iAssemblageAfterPEA
    real(8) :: dTimerStart, dTimerStop
    logical :: lRepeatedPEAAssemblage, lPEAAllowed
    
    iterGlobal = 0
    dGEMTimingLeveling = 0D0
    dGEMTimingHandoff = 0D0
    dGEMTimingPEA = 0D0
    dGEMTimingLagrangian = 0D0
    dGEMTimingCertification = 0D0
    nGEMCertificationSweep = 0

    ! Stage 1: initialize GEM and run coarse leveling.
    call cpu_time(dTimerStart)
    call RunLeveling
    call cpu_time(dTimerStop)
    dGEMTimingLeveling = dGEMTimingLeveling + MAX(0D0, dTimerStop - dTimerStart)
    if (lCompbdOnly) return
    lPEAAllowed = (nElements > 1)

    ! Stage 2: run PEA from the original Leveling plane for multi-component
    ! active-set discovery.  Unary active-component systems have no PEA search
    ! and are translated directly to Lagrangian form.
    if (lPhaseChange .AND. (iPhaseChangeReason == PHASE_CHANGE_REASON_LEVELING_INITIAL)) then
        if (.NOT.lPEAAllowed) then
            call cpu_time(dTimerStart)
            call Level2Lagrange
            call cpu_time(dTimerStop)
            dGEMTimingHandoff = dGEMTimingHandoff + MAX(0D0, dTimerStop - dTimerStart)
            call cpu_time(dTimerStart)
            call RunLagrangianGEM
            call cpu_time(dTimerStop)
            dGEMTimingLagrangian = dGEMTimingLagrangian + MAX(0D0, dTimerStop - dTimerStart)
            if (lCompbdOnly) return
        end if

        call cpu_time(dTimerStart)
        call CheckPhaseAssemblage
        call cpu_time(dTimerStop)
        dGEMTimingPEA = dGEMTimingPEA + MAX(0D0, dTimerStop - dTimerStart)
        if ((nSolnPhases==0 .and. dGEMFunctionNorm<1D-7 ).or.lCompbdOnly) return
    end if

    ! Stage 3: try the PEA-refined active set directly with Lagrangian GEM.
    call cpu_time(dTimerStart)
    call RunLagrangianGEM
    call cpu_time(dTimerStop)
    dGEMTimingLagrangian = dGEMTimingLagrangian + MAX(0D0, dTimerStop - dTimerStart)
    if (lConverged .or. lCompbdOnly) return
    if (.NOT.lPEAAllowed) return

    ! If the fixed active set failed without a specific trigger, run one PEA
    ! repair pass so the failure is explicit instead of silently maxing out.
    if (.NOT.lPhaseChange) then
        lPhaseChange = .TRUE.
        iPhaseChangeReason = PHASE_CHANGE_REASON_LAGRANGIAN_UNCONVERGED
        iGEMExitStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
    end if

    ! Stage 4: event-triggered active-set repair followed by Lagrangian retry.
    LOOP_AssemblageRepair: do iterAssemblageRepair = 1, iterAssemblageRepairMax
        if (.NOT.lPhaseChange) exit LOOP_AssemblageRepair

        iAssemblageBeforePEA = iAssemblage
        call cpu_time(dTimerStart)
        call CheckPhaseAssemblage
        call cpu_time(dTimerStop)
        dGEMTimingPEA = dGEMTimingPEA + MAX(0D0, dTimerStop - dTimerStart)
        iAssemblageAfterPEA = iAssemblage
        if ((nSolnPhases==0 .and. dGEMFunctionNorm<1D-7 ).or.lCompbdOnly) return

        call cpu_time(dTimerStart)
        call RunLagrangianGEM
        call cpu_time(dTimerStop)
        dGEMTimingLagrangian = dGEMTimingLagrangian + MAX(0D0, dTimerStop - dTimerStart)
        if (lConverged .or. lCompbdOnly) exit LOOP_AssemblageRepair

        lRepeatedPEAAssemblage = ALL(iAssemblageAfterPEA == iAssemblageBeforePEA)
        if (lPhaseChange .AND. lRepeatedPEAAssemblage) exit LOOP_AssemblageRepair
    end do LOOP_AssemblageRepair

    if ((.NOT.lConverged).AND.(dGEMFunctionNorm > 1D-5)) then
        iGEMExitStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
    end if

    return

end subroutine MultiPhaseMinimizer
