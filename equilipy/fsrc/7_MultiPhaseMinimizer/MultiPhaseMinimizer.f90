!> \brief Orchestrate Leveling, fixed-active-set Lagrangian GEM, and PEA repair.
!!
!! \details Runs the high-level minimization sequence used by GEMSolver.  The
!! routine starts from classical Leveling, checks the active set with
!! CheckPhaseAssemblage/PEA, then uses Lagrangian GEM and event-triggered PEA
!! repair to finish the active-set solve.  Optional grid rows augment the
!! Leveling and PEA candidate pool; they never replace this workflow.
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
    !   07/20/2026      S.Y. Kwon           Integrated optional static-grid discovery into the standard certified PEA pipeline with honest terminal status.
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
    ! CertifyLevelingCompoundsOnly   Certifies a provisional compounds-only
    !                                Leveling plane against all solution minima.
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
    ! - A compounds-only Leveling result is terminal without a sweep only when
    !   the system contains no solution phases.  Otherwise one initial PEA
    !   solution-minimum sweep must certify the Leveling plane.
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
    integer :: iterAssemblageRepair, iSlot, iParentPhase
    integer, parameter :: iterAssemblageRepairMax = 100
    integer, dimension(nElements) :: iAssemblageBeforePEA, iAssemblageAfterPEA
    real(8) :: dTimerStart, dTimerStop
    logical :: lCanonicalPartitionWitnessPending
    logical :: lRepeatedPEAAssemblage, lPEAAllowed
    
    iterGlobal = 0
    dGEMTimingLeveling = 0D0
    dGEMTimingHandoff = 0D0
    dGEMTimingPEA = 0D0
    dGEMTimingLagrangian = 0D0
    dGEMTimingCertification = 0D0
    dGEMTimingGridGeneration = 0D0
    dGEMTimingGridRefinement = 0D0
    nSUBLHessianEndmemberProjectionCall = 0
    nGEMCertificationSweep = 0
    iPEAUncertifiedHandoffSeen = 0

    ! Step 1. Initialize GEM and run coarse or static-grid Leveling.
    call cpu_time(dTimerStart)
    call RunLeveling
    call cpu_time(dTimerStop)
    dGEMTimingLeveling = dGEMTimingLeveling + MAX(0D0, dTimerStop - dTimerStart)
    if (lCompbdOnly) then
        if (nSolnPhasesSys == 0) return
        call CertifyLevelingCompoundsOnly
        if (lCompbdOnly) return
    end if
    lPEAAllowed = (nElements > 1)

    ! Step 2. Run PEA from the original Leveling plane for multi-component
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
        if (iPEAExitStatus /= PEA_EXIT_STATUS_CERTIFIED_SETTLED) then
            iPEAUncertifiedHandoffSeen = 1
        end if
        if ((nSolnPhases==0 .and. dGEMFunctionNorm<1D-7 ).or.lCompbdOnly) return
    end if

    ! Step 3. Try the PEA-refined active set directly with Lagrangian GEM.
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

    ! Step 4. Repair classified active-set events and retry Lagrangian GEM.
    LOOP_AssemblageRepair: do iterAssemblageRepair = 1, iterAssemblageRepairMax
        if (.NOT.lPhaseChange) exit LOOP_AssemblageRepair

        iAssemblageBeforePEA = iAssemblage
        call cpu_time(dTimerStart)
        call CheckPhaseAssemblage
        call cpu_time(dTimerStop)
        dGEMTimingPEA = dGEMTimingPEA + MAX(0D0, dTimerStop - dTimerStart)
        if (iPEAExitStatus /= PEA_EXIT_STATUS_CERTIFIED_SETTLED) then
            iPEAUncertifiedHandoffSeen = 1
        end if
        iAssemblageAfterPEA = iAssemblage
        if ((nSolnPhases==0 .and. dGEMFunctionNorm<1D-7 ).or.lCompbdOnly) return

        call cpu_time(dTimerStart)
        call RunLagrangianGEM
        call cpu_time(dTimerStop)
        dGEMTimingLagrangian = dGEMTimingLagrangian + MAX(0D0, dTimerStop - dTimerStart)
        if (lConverged .or. lCompbdOnly) exit LOOP_AssemblageRepair

        lRepeatedPEAAssemblage = ALL(iAssemblageAfterPEA == iAssemblageBeforePEA)
        lCanonicalPartitionWitnessPending = .FALSE.
        if (lODPartitionUnifiedActive.AND.&
            (iPhaseChangeReason == PHASE_CHANGE_REASON_PHASE_POTENTIAL).AND.&
            allocated(iActiveSlotThermoPhase).AND.&
            allocated(iActiveSlotODClass).AND.&
            allocated(iODCompanionPhase)) then
            do iSlot = 1, MIN(nElements, SIZE(iActiveSlotThermoPhase))
                iParentPhase = iActiveSlotThermoPhase(iSlot)
                if ((iParentPhase <= 0).OR.(iParentPhase > SIZE(iODCompanionPhase))) cycle
                if (iODCompanionPhase(iParentPhase) <= 0) cycle
                if (iActiveSlotODClass(iSlot) /= OD_CANDIDATE_ORDERED) cycle
                lCanonicalPartitionWitnessPending = .TRUE.
                exit
            end do
        end if
        ! A repeated phase id is not a repeated physical set when the ordered
        ! parent has just exposed a favorable second composition set.  Let PEA
        ! price that same-parent witness once more; the bounded repair loop and
        ! honest final status remain authoritative if no new set is admitted.
        if (lPhaseChange .AND. lRepeatedPEAAssemblage.AND.&
            (.NOT.lCanonicalPartitionWitnessPending)) exit LOOP_AssemblageRepair
    end do LOOP_AssemblageRepair

    if (.NOT.lConverged) then
        if (iPEAUncertifiedHandoffSeen == 1) then
            iGEMExitStatus = GEM_EXIT_STATUS_PEA_UNCERTIFIED_HANDOFF
        else
            iGEMExitStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
        end if
    end if

    return

end subroutine MultiPhaseMinimizer
