!> \brief Run fixed-active-set Lagrangian Gibbs energy minimization.
!!
!! \details Solves the current active phase assemblage with Newton updates,
!! line search, trace-species handling, and convergence checks.  Phase
!! assemblage repair is reported through lPhaseChange and iPhaseChangeReason;
!! this routine does not call CheckPhaseAssemblage internally.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    RunLagrangianGEM.f90
    !> \brief   Run the fixed-active-set Lagrangian GEM iteration.
    !> \author  S.Y. Kwon
    !> \date    Jun. 24, 2026
    !> \sa      MultiPhaseMinimizer.f90
    !> \sa      GEMNewton.f90
    !> \sa      GEMLineSearch.f90
    !> \sa      CheckConvergence.f90
    !> \sa      RecordGEMIterationDiagnostics.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/24/2026      S.Y. Kwon           Documented fixed-assemblage GEM stage
    !   06/24/2026      S.Y. Kwon           Recorded per-iteration Lagrangian diagnostics
    !   06/25/2026      S.Y. Kwon           Dispatched CEF fixed assemblages to site-fraction Newton/line search
    !   06/25/2026      S.Y. Kwon           Allow compound-only postprocess solves to update the elemental
    !                                        potential plane at perturbed temperature.
!   06/25/2026      S.Y. Kwon           Delayed raw CEF phase removal until the positivity-preserving line
!                                        search fails to find descent.
!   06/25/2026      S.Y. Kwon           Removed overfull CEF active phases pinned at the phase-amount
!                                        boundary after stalled line search.
!   06/25/2026      S.Y. Kwon           Marked when the CEF site-fraction KKT path is active so
!                                        CompFunctionNorm uses matching residual coordinates.
!   06/26/2026      S.Y. Kwon           Added a PEA-polish mode for refining elemental potentials
!                                        inside CheckPhaseAssemblage.
!   06/26/2026      S.Y. Kwon           Rejected CEF PEA-polish attempts on boundary/stagnation events
!                                        instead of removing phases inside the embedded Lagrangian solve.
!   06/27/2026      S.Y. Kwon           Coalesced degenerate random order/disorder active sets before
!                                        accepting CEF Lagrangian convergence.
!   06/27/2026      S.Y. Kwon           Let PEA-polish Lagrangian use the same coalescing, trace handling,
!                                        and CEF boundary-removal rules as the final Lagrangian solve.
!   06/29/2026      S.Y. Kwon           Added fixed-assemblage postprocess stagnation exit without
!                                        requesting phase-assemblage repair.
!   06/29/2026      S.Y. Kwon           Added a one-shot ideal-log retry when a SUBG/SUBQ analytical
!                                        model-Hessian direction gives no line-search descent.
!   06/30/2026      S.Y. Kwon           Removed tiny boundary solution phases before declaring stagnation.
!   07/01/2026      S.Y. Kwon           Try the CEF residual-LM direction before line search when the
!                                        analytical CEF Newton step would make a finite phase negative,
!                                        including fixed-assemblage postprocess solves.
!   06/30/2026      S.Y. Kwon           Activated CEF site-fraction retry only for stalled underfull sets.
!   06/30/2026      S.Y. Kwon           Added debug-only max-iteration warnings for Lagrangian GEM.
!   07/01/2026      S.Y. Kwon           Audited converged SUBOM phases for lower ordered active-phase
!                                        minima before accepting a random-state stationary point.
!   07/01/2026      S.Y. Kwon           Added a residual Levenberg fallback for stalled CEF site-fraction
!                                        Newton directions while preserving trace-site complementarity.
!   07/01/2026      S.Y. Kwon           Allowed the residual Levenberg fallback for single active CEF
!                                        phases when the analytical site-fraction direction gives no descent.
!   07/02/2026      S.Y. Kwon           Split Phase 0 scaffold counters for trace handling, CEF retry,
!                                        raw-negative removal, and boundary-pinned removal.
!   07/02/2026      S.Y. Kwon           Limited CEF residual-LM activation to SUBOM active sets so
!                                        ordinary SUBL/SUBLM phases can use the positivity-preserving
!                                        analytical KKT line search first.
!   07/03/2026      S.Y. Kwon           Split residual-LM scaffold diagnostics by raw-negative and
!                                        no-descent activation reason.
!   07/03/2026      S.Y. Kwon           Tried a complementarity-valid bound-active CEF phase KKT solve
!                                        before raw-negative residual-LM fallback.
!   07/03/2026      S.Y. Kwon           Protected switch-gated two-composition-set SUBOM slots from the
!                                        legacy order/disorder coalescer and endpoint stability restart.
!   07/04/2026      S.Y. Kwon           Counted tiny-boundary removals separately from boundary-pinned
!                                        removals for C1-a passive diagnostics.
!   07/03/2026      S.Y. Kwon           Routed switch-gated two-set SUBOM raw-negative overshoots through
!                                        the bound-active KKT retry before residual-LM fallback.
!   07/03/2026      S.Y. Kwon           Let converged two-set SUBOM class-2 no-descent states reach
!                                        CheckConvergence instead of residual-LM.
!   07/03/2026      S.Y. Kwon           Coalesced duplicate same-parent SUBOM composition sets after
!                                        convergence when their slot-local site fractions collapse.
!   07/03/2026      S.Y. Kwon           Added per-minimization trial-7 and residual-LM activation totals
!                                        so scaffold dependence is visible from Python tests.
!   07/03/2026      S.Y. Kwon           Split passive totals for duplicate SUBOM coalescence, legacy
!                                        order/disorder coalescence, and stability restarts.
!   07/03/2026      S.Y. Kwon           Removed the SUBG/SUBQ ideal-log no-descent retry after
!                                        validating the analytical pair derivative basis.
!   07/04/2026      S.Y. Kwon           Routed invalid-complementarity raw-negative CEF phase events through
!                                        a certified bound-active KKT multiplier before residual-LM fallback.
!   07/04/2026      S.Y. Kwon           Snapshotted analytical-direction diagnostics before residual-LM
!                                        fallback for C3-a2 event-buffer reporting.
!   07/04/2026      S.Y. Kwon           Tried C3-c1 inertia-regularized CEF KKT directions before
!                                        no-descent residual-LM fallback.
!   07/05/2026      S.Y. Kwon           Set explicit GEM exit status when a nonconverged fixed-active-set
!                                        solve exits through stagnation.
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to minimize Gibbs energy for
    !! the phase assemblage already selected by Leveling or CheckPhaseAssemblage.
    !! It performs repeated Newton, line-search, trace-control, and convergence
    !! checks until the fixed active set converges or a classified active-set
    !! repair trigger is raised.
    !
    ! Required input variables:
    ! =========================
    !
    ! iAssemblage          Current active phase assemblage.
    ! dMolesSpecies        Current active-set species amounts.
    ! dMolesPhase          Current active-set phase amounts.
    ! dElementPotential    Current elemental potentials.
    ! lPostProcess         True for fixed-assemblage Cp perturbation solves.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! lConverged           True when fixed-active-set Lagrangian GEM converges.
    ! lPhaseChange         True when the active set needs CheckPhaseAssemblage repair.
    ! iPhaseChangeReason   Classified active-set repair reason.
    ! dGEM*History         Per-iteration residual and line-search diagnostics.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! UseCEFLagrangian                  Selects the CEF site-fraction GEM path when safe.
    ! GEMNewton/GEMNewtonCEF            Builds and solves the fixed-active-set Newton system.
    ! AuditCEFRawPhaseDirection         Records raw CEF phase-removal directions for diagnostics.
    ! RemoveRawNegativeCEFPhase         Removes a raw-negative CEF phase only after line-search failure.
    ! RemoveCEFBoundaryPhaseAfterLineSearch
    !                                   Removes a CEF phase pinned at the phase-amount boundary.
    ! CoalesceDegenerateOrderDisorderPhases
    !                                   Merges random-state DIS_PART duplicates into the ordered phase.
    ! CoalesceDuplicateSUBOMCompositionSets
    !                                   Merges duplicate same-parent SUBOM composition-set slots.
    ! StabilizeActiveOrderDisorderPhases
    !                                   Restarts an active SUBOM phase when endpoint starts find a lower
    !                                   ordered minimum under the current elemental-potential plane.
    ! GEMLineSearch/GEMLineSearchCEF    Applies the damped Newton step.
    ! RemoveTraceSpeciesBelowThreshold  Temporarily removes tiny active solution endmembers.
    ! CheckConvergence                  Classifies convergence or active-set repair triggers.
    ! ShouldReinjectTraceSpeciesForSlowProgress
    !                                   Requests reinjection after slow reduced progress.
    ! ReinjectTraceSpecies              Reintroduces temporarily removed trace endmembers.
    ! RecordGEMIterationDiagnostics     Stores the iteration state for debugging.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! MultiPhaseMinimizer       Production minimizer driver.
    ! PostProcess               Fixed-assemblage Cp perturbation path.
    ! Python run_lagrangian_gem Staged debugging entry point.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - CheckPhaseAssemblage/PEA is owned by the caller for phase discovery,
    !   removal, and replacement.  Raw CEF phase amount directions are recorded
    !   as diagnostics.  A raw undamped Newton overshoot can remove a phase only
    !   after the positivity-preserving CEF line search fails to find descent.
    ! - In overfull CEF active sets, a phase whose raw and accepted line-search
    !   amounts are both below the phase-mole boundary tolerance can be removed
    !   after line-search stagnation; otherwise the positivity line search can
    !   waste iterations shrinking an inactive phase geometrically.
    ! - During PEA-polish solves, Lagrangian uses the same numerical solve rules
    !   as the final fixed-active-set solve.  PEALagrangianPolish owns what is
    !   persisted back to PEA after this routine returns.
    ! - A SUBOM phase at its random state is thermodynamically degenerate with
    !   its DIS_PART companion.  Normal Lagrangian solves coalesce that
    !   duplicate representation into the ordered phase and then resolve the
    !   reduced active set once more.
    ! - A first-order CEF KKT residual can be small at a random ordered-phase
    !   saddle.  Before accepting convergence, active SUBOM phases are checked
    !   against endpoint-start phase-local minimization under the same
    !   elemental-potential plane.  A lower ordered point is used only as a
    !   new seed for another fixed-active-set Lagrangian solve.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine RunLagrangianGEM

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer :: INFO, iterReason, iterCurrentMax
    integer :: iSolnSlot, iSolnPhase
    integer :: iCEFStabilityRestartCount
    logical :: lTraceChanged, lTraceSlowReinject, lUseCEFSiteGEM
    logical :: lCEFPhaseRemoved, lCompEverything
    logical :: lOrderDisorderCoalesced
    logical :: lDuplicateSUBOMCoalesced
    logical :: lOrderDisorderStabilized
    logical :: lBoundaryPhaseRemoved, lCEFRetryActivated
    logical :: lActiveSetContainsSUBOM
    logical :: lUseResidualLMForRawNegative, lUseResidualLMForNoDescent
    logical :: lUseBoundPhaseRetry
    logical :: lInvalidCompBoundRetry
    logical :: lDistinctSUBOMCompositionSets
    integer :: iInvalidCompInfo
    integer, parameter :: iCEFStabilityRestartMax = 4
    real(8), parameter :: dPostProcessStagnationRelativeTolerance = 1D-12
    real(8), parameter :: dPostProcessStagnationNormTolerance = 1D-10
    real(8) :: dNormBeforeIteration, dNormScale, dRelativeProgress
    real(8) :: dInvalidCompPhi
    logical :: lStagnationCandidate, lHitLagrangianMax, lResidualLMUsedThisIteration

    if (lCompbdOnly .AND. (.NOT.lPostProcess)) return

    ! Trace masks are initialized by InitGEMSolver and intentionally persist
    ! across PEA/Lagrangian retries so removed trace species can be reinjected
    ! by the slow-progress logic instead of being forced back immediately.
    iterGlobal          = 0
    lConverged          = .False.
    lPhaseChange        = .False.
    iPhaseChangeReason  = PHASE_CHANGE_REASON_NONE
    iGEMExitStatus      = GEM_EXIT_STATUS_OK
    iTraceSpeciesSlowProgressCount = 0
    iGEMStagnationCount = 0
    iGEMCEFResidualLMFallbackCount = 0
    iterTraceSpeciesLastRemoval = 0
    iterTraceSpeciesLastReinject = 0
    dTraceSpeciesReducedNormLast = HUGE(1D0)
    lGEMCEFBndPhaseActive = .FALSE.
    iGEMCEFBndPhaseSlot = 0
    dGEMCEFBndPhaseStep = 0D0
    lCompEverything = .FALSE.
    iCEFStabilityRestartCount = 0
    iterCurrentMax = iterGlobalMax
    if (lPEALagrangianPolishActive) iterCurrentMax = iterPEALagrangianPolishMax
    call UseCEFLagrangian(lUseCEFSiteGEM)
    lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM

    LOOP_GEMSolver: do iterGlobal = 1, iterCurrentMax

        dNormBeforeIteration = dGEMFunctionNorm

        call ResetGEMIterationDiagnostics
        call UseCEFLagrangian(lUseCEFSiteGEM)
        lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM

        ! Construct the Hessian matrix and compute the direction vector.
        lResidualLMUsedThisIteration = .FALSE.
        if (lUseCEFSiteGEM) then
            call GEMNewtonCEF(INFO)
        else
            call GEMNewton(INFO)
        end if
        if (INFO /= 0) then
            lPhaseChange = .True.
            iPhaseChangeReason = PHASE_CHANGE_REASON_NEWTON_FAILURE
            lPhaseChangeHistory(iterGlobal) = .TRUE.
            iPhaseChangeReasonHistory(iterGlobal) = iPhaseChangeReason
            call RecordGEMIterationDiagnostics
            exit LOOP_GEMSolver
        end if

        if (lUseCEFSiteGEM) then
            call AuditCEFRawPhaseDirection
            lActiveSetContainsSUBOM = .FALSE.
            LOOP_RESIDUAL_LM_SUBOM_CHECK: do iSolnSlot = 1, nSolnPhases
                iSolnPhase = -iAssemblage(nElements - iSolnSlot + 1)
                if ((iSolnPhase <= 0).OR.(iSolnPhase > nSolnPhasesSys)) cycle LOOP_RESIDUAL_LM_SUBOM_CHECK
                if (TRIM(cSolnPhaseType(iSolnPhase)) == 'SUBOM') then
                    lActiveSetContainsSUBOM = .TRUE.
                    exit LOOP_RESIDUAL_LM_SUBOM_CHECK
                end if
            end do LOOP_RESIDUAL_LM_SUBOM_CHECK
            lUseBoundPhaseRetry = .FALSE.
            lInvalidCompBoundRetry = .FALSE.
            lDistinctSUBOMCompositionSets = DistinctRepeatedSUBOMCompositionSetsActive()
            ! Trial 7 is the first recovery path for complementarity-valid
            ! raw-negative CEF phase directions.  Residual-LM remains a final
            ! fallback only when this bound-active KKT retry is not available.
            if ((iGEMRawNegativePhaseSlot > 0).AND.(nSolnPhases > 1).AND.&
                lActiveSetContainsSUBOM.AND.&
                ((iGEMRawNegComp /= 0).OR.lDistinctSUBOMCompositionSets)) then
                lInvalidCompBoundRetry = (iGEMRawNegComp < 0)
                iGEMTrial7BoundRetryAttemptedTotal = iGEMTrial7BoundRetryAttemptedTotal + 1
                if (lInvalidCompBoundRetry) then
                    iGEMInvalidCompBoundAttemptedUsed = iGEMInvalidCompBoundAttemptedUsed + 1
                    iGEMInvalidCompBoundAttemptedTotal = iGEMInvalidCompBoundAttemptedTotal + 1
                end if
                lGEMCEFBndPhaseActive = .TRUE.
                iGEMCEFBndPhaseSlot = iGEMRawNegativePhaseSlot
                dGEMCEFBndPhaseStep = DMAX1(10D0 * dTolerance(8), dTraceSpeciesRemoveFraction) - &
                    dMolesPhase(iGEMRawNegativePhaseSlot)
                call GEMNewtonCEF(INFO)
                if (INFO == 0) then
                    if (lInvalidCompBoundRetry) then
                        call ComputeInvalidCompBoundPhi(iGEMRawNegativePhaseSlot, &
                            dInvalidCompPhi, iInvalidCompInfo)
                        dGEMInvalidCompBoundPhi = dInvalidCompPhi
                        if ((iInvalidCompInfo == 0).AND.(dInvalidCompPhi >= -1D-10)) then
                            lUseBoundPhaseRetry = .TRUE.
                            iGEMInvalidCompBoundAcceptedUsed = iGEMInvalidCompBoundAcceptedUsed + 1
                            iGEMInvalidCompBoundAcceptedTotal = iGEMInvalidCompBoundAcceptedTotal + 1
                            iGEMInvalidCompBoundVerdict = 1
                        else
                            iGEMInvalidCompBoundRejectedUsed = iGEMInvalidCompBoundRejectedUsed + 1
                            iGEMInvalidCompBoundRejectedTotal = iGEMInvalidCompBoundRejectedTotal + 1
                            iGEMInvalidCompBoundVerdict = -1
                        end if
                    else
                        lUseBoundPhaseRetry = .TRUE.
                    end if

                    if (lUseBoundPhaseRetry) then
                        iGEMTrial7BoundRetryAcceptedTotal = iGEMTrial7BoundRetryAcceptedTotal + 1
                    end if
                end if

                if (.NOT.lUseBoundPhaseRetry) then
                    if (lInvalidCompBoundRetry.AND.(INFO /= 0)) then
                        iGEMInvalidCompBoundRejectedUsed = iGEMInvalidCompBoundRejectedUsed + 1
                        iGEMInvalidCompBoundRejectedTotal = iGEMInvalidCompBoundRejectedTotal + 1
                        iGEMInvalidCompBoundVerdict = -2
                    end if
                    lGEMCEFBndPhaseActive = .FALSE.
                    iGEMCEFBndPhaseSlot = 0
                    dGEMCEFBndPhaseStep = 0D0
                end if
            end if
            lUseResidualLMForRawNegative = .FALSE.
            if ((.NOT.lUseBoundPhaseRetry).AND.(iGEMRawNegativePhaseSlot > 0).AND.(nSolnPhases > 1)) then
                lUseResidualLMForRawNegative = lActiveSetContainsSUBOM
            end if
            if (lUseResidualLMForRawNegative) then
                iGEMPreLMNoDescent = 0
                iGEMPreLMNoDescentClass = 0
                call SnapshotGEMPreLMDiagnostics
                lGEMCEFResidualLMDirection = .TRUE.
                call GEMNewtonCEF(INFO)
                lGEMCEFResidualLMDirection = .FALSE.
                if (INFO == 0) then
                    iGEMCEFResidualLMFallbackCount = iGEMCEFResidualLMFallbackCount + 1
                    iGEMResidualLMUsed = iGEMResidualLMUsed + 1
                    iGEMResidualLMRawNegativeUsed = iGEMResidualLMRawNegativeUsed + 1
                    iGEMResidualLMTotalUsed = iGEMResidualLMTotalUsed + 1
                    iGEMResidualLMRawNegativeTotalUsed = iGEMResidualLMRawNegativeTotalUsed + 1
                    call RecordGEMResidualLMEvent(1)
                    lResidualLMUsedThisIteration = .TRUE.
                end if
            end if
        end if

        ! Perform a line search using the direction vector.
        if (lUseCEFSiteGEM) then
            call GEMLineSearchCEF
            if (lUseBoundPhaseRetry) then
                lGEMCEFBndPhaseActive = .FALSE.
                iGEMCEFBndPhaseSlot = 0
                dGEMCEFBndPhaseStep = 0D0
            end if
            lUseResidualLMForNoDescent = .FALSE.
            if ((iGEMLineSearchNoDescent > 0).AND.&
                (iGEMLineSearchNegativePhaseCount == 0).AND.&
                (.NOT.lResidualLMUsedThisIteration)) then
                lUseResidualLMForNoDescent = lActiveSetContainsSUBOM
                if (lSUBOMTwoSetCandidateEnabled.AND.lDistinctSUBOMCompositionSets.AND.&
                    (iGEMLineSearchNoDescentClass == 2).AND.&
                    (dGEMFunctionNorm < 1D-5)) then
                    lUseResidualLMForNoDescent = .FALSE.
                end if
            end if
            if (lUseResidualLMForNoDescent) then
                iGEMPreLMNoDescent = iGEMLineSearchNoDescent
                iGEMPreLMNoDescentClass = iGEMLineSearchNoDescentClass
                call SnapshotGEMPreLMDiagnostics
                if (lGEMCEFInertiaRegularizationEnabled.AND.CurrentGEMCEFInertiaWrong()) then
                    lGEMCEFInertiaRegularizationActive = .TRUE.
                    call GEMNewtonCEF(INFO)
                    lGEMCEFInertiaRegularizationActive = .FALSE.
                    if ((INFO == 0).AND.(iGEMInertiaRegularizationAcceptedUsed > 0)) then
                        call GEMLineSearchCEF
                        if (iGEMLineSearchNoDescent == 0) then
                            lUseResidualLMForNoDescent = .FALSE.
                        else
                            call SnapshotGEMPreLMDiagnostics
                        end if
                    end if
                end if
            end if
            if (lUseResidualLMForNoDescent) then
                lGEMCEFResidualLMDirection = .TRUE.
                call GEMNewtonCEF(INFO)
                lGEMCEFResidualLMDirection = .FALSE.
                if (INFO == 0) then
                    iGEMCEFResidualLMFallbackCount = iGEMCEFResidualLMFallbackCount + 1
                    iGEMResidualLMUsed = iGEMResidualLMUsed + 1
                    iGEMResidualLMNoDescentUsed = iGEMResidualLMNoDescentUsed + 1
                    iGEMResidualLMNoDescentClass = iGEMLineSearchNoDescentClass
                    iGEMResidualLMTotalUsed = iGEMResidualLMTotalUsed + 1
                    iGEMResidualLMNoDescentTotalUsed = iGEMResidualLMNoDescentTotalUsed + 1
                    call RecordGEMResidualLMEvent(2)
                    call GEMLineSearchCEF
                end if
            end if
            if (.NOT.lPostProcess) then
                lCEFPhaseRemoved = .FALSE.
                if ((iGEMLineSearchNegativePhaseCount > 0).AND.(iGEMLineSearchNoDescent > 0)) then
                    call RemoveRawNegativeCEFPhase(lCEFPhaseRemoved)
                    if (lCEFPhaseRemoved) iGEMRawNegativeRemovalUsed = iGEMRawNegativeRemovalUsed + 1
                end if
                if (.NOT.lCEFPhaseRemoved) then
                    call RemoveCEFBoundaryPhaseAfterLineSearch(lCEFPhaseRemoved)
                    if (lCEFPhaseRemoved) then
                        iGEMBoundaryPinnedRemovalUsed = iGEMBoundaryPinnedRemovalUsed + 1
                    end if
                end if
                if (lCEFPhaseRemoved) then
                    iGEMBoundaryRemovalUsed = iGEMBoundaryRemovalUsed + 1
                    lPhaseChange = .FALSE.
                    iPhaseChangeReason = PHASE_CHANGE_REASON_NONE
                    call CompChemicalPotential(lCompEverything)
                    call CompFunctionNorm
                    call UseCEFLagrangian(lUseCEFSiteGEM)
                    lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM
                    call RecordGEMIterationDiagnostics
                    cycle LOOP_GEMSolver
                end if
            end if
            lTraceChanged = .FALSE.
        else
            call GEMLineSearch
            call RemoveTraceSpeciesBelowThreshold(lTraceChanged)
            if (lTraceChanged) iGEMTraceRemoveUsed = iGEMTraceRemoveUsed + 1
        end if

        call CheckConvergence

        if (lUseCEFSiteGEM.AND.lConverged) then
            call CoalesceDuplicateSUBOMCompositionSets(lDuplicateSUBOMCoalesced)
            if (lDuplicateSUBOMCoalesced) then
                iGEMCoalesceUsed = iGEMCoalesceUsed + 1
                iGEMDuplicateSUBOMCoalesceTotalUsed = iGEMDuplicateSUBOMCoalesceTotalUsed + 1
                lConverged = .FALSE.
                call CompChemicalPotential(lCompEverything)
                call CompFunctionNorm
                call UseCEFLagrangian(lUseCEFSiteGEM)
                lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM
                call RecordGEMIterationDiagnostics
                cycle LOOP_GEMSolver
            end if
        end if

        lDistinctSUBOMCompositionSets = DistinctRepeatedSUBOMCompositionSetsActive()

        if (lUseCEFSiteGEM.AND.lConverged.AND.(.NOT.lDistinctSUBOMCompositionSets)) then
            call CoalesceDegenerateOrderDisorderPhases(lOrderDisorderCoalesced)
            if (lOrderDisorderCoalesced) then
                iGEMCoalesceUsed = iGEMCoalesceUsed + 1
                iGEMDegenerateODCoalesceTotalUsed = iGEMDegenerateODCoalesceTotalUsed + 1
                lConverged = .FALSE.
                call CompChemicalPotential(lCompEverything)
                call CompFunctionNorm
                call UseCEFLagrangian(lUseCEFSiteGEM)
                lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM
                call RecordGEMIterationDiagnostics
                cycle LOOP_GEMSolver
            end if
        end if

        if (lUseCEFSiteGEM.AND.lConverged.AND.(.NOT.lPostProcess).AND.&
            (.NOT.lDistinctSUBOMCompositionSets)) then
            lOrderDisorderStabilized = .FALSE.
            if (iCEFStabilityRestartCount < iCEFStabilityRestartMax) then
                call StabilizeActiveOrderDisorderPhases(lOrderDisorderStabilized)
            end if
            if (lOrderDisorderStabilized) then
                iGEMStabilizeUsed = iGEMStabilizeUsed + 1
                iGEMStabilizeTotalUsed = iGEMStabilizeTotalUsed + 1
                iCEFStabilityRestartCount = iCEFStabilityRestartCount + 1
                lConverged = .FALSE.
                call CompChemicalPotential(lCompEverything)
                call CompFunctionNorm
                call UseCEFLagrangian(lUseCEFSiteGEM)
                lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM
                call RecordGEMIterationDiagnostics
                cycle LOOP_GEMSolver
            end if
        end if

        lTraceSlowReinject = .FALSE.
        if ((.NOT.lUseCEFSiteGEM).AND.(.NOT.lConverged)) then
            call ShouldReinjectTraceSpeciesForSlowProgress(lTraceSlowReinject)
        end if

        if ((.NOT.lUseCEFSiteGEM).AND.(lConverged .OR. lTraceSlowReinject)) then
            call ReinjectTraceSpecies(lTraceChanged)
            if (lTraceChanged) then
                iGEMTraceReinjectUsed = iGEMTraceReinjectUsed + 1
                lConverged = .FALSE.
                call RecordGEMIterationDiagnostics
                cycle LOOP_GEMSolver
            end if
            if (lConverged) then
                call RecordGEMIterationDiagnostics
                exit LOOP_GEMSolver
            end if
        else if (lConverged) then
            call RecordGEMIterationDiagnostics
            exit LOOP_GEMSolver
        end if

        ! Surface a bad fixed assemblage to Python instead of invoking PEA here.
        if (lPhaseChange) then
            call RecordGEMIterationDiagnostics
            exit LOOP_GEMSolver
        end if

        dNormScale = MAX(DABS(dNormBeforeIteration), 1D-300)
        dRelativeProgress = (dNormBeforeIteration - dGEMFunctionNorm) / dNormScale
        lStagnationCandidate = .FALSE.

        if (lPostProcess) then
            lStagnationCandidate = (dGEMFunctionNorm > dPostProcessStagnationNormTolerance).AND.&
                ((iGEMLineSearchNoDescent > 0).OR.&
                (DABS(dRelativeProgress) < dPostProcessStagnationRelativeTolerance))
        else if (dGEMFunctionNorm > dGEMStagnationNormTolerance) then
            lStagnationCandidate = dRelativeProgress < dGEMStagnationRelativeTolerance
        end if

        if (lStagnationCandidate) then
            iGEMStagnationCount = iGEMStagnationCount + 1
        else
            iGEMStagnationCount = 0
        end if

        if (iGEMStagnationCount >= iGEMStagnationWindow) then
            call RemoveTinyBoundarySolutionPhase(lBoundaryPhaseRemoved)
            if (lBoundaryPhaseRemoved) then
                iGEMBoundaryRemovalUsed = iGEMBoundaryRemovalUsed + 1
                iGEMBoundaryPinnedRemovalUsed = iGEMBoundaryPinnedRemovalUsed + 1
                iGEMTinyBoundaryRemovalUsed = iGEMTinyBoundaryRemovalUsed + 1
                lPhaseChange = .FALSE.
                iPhaseChangeReason = PHASE_CHANGE_REASON_NONE
                iGEMStagnationCount = 0
                call CompChemicalPotential(lCompEverything)
                call CompFunctionNorm
                call UseCEFLagrangian(lUseCEFSiteGEM)
                lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM
                call RecordGEMIterationDiagnostics
                cycle LOOP_GEMSolver
            end if

            lCEFRetryActivated = .FALSE.
            if ((.NOT.lPostProcess).AND.(.NOT.lUseCEFSiteGEM).AND.&
                (nSolnPhases > 0).AND.(nConPhases + nSolnPhases < nElements)) then
                lGEMCEFSiteLagrangianEnabled = .TRUE.
                call UseCEFLagrangian(lCEFRetryActivated)
                if (lCEFRetryActivated) then
                    iGEMCEFRetryActivationUsed = iGEMCEFRetryActivationUsed + 1
                    iGEMCEFRetryActivationTotalUsed = iGEMCEFRetryActivationTotalUsed + 1
                    lPhaseChange = .FALSE.
                    iPhaseChangeReason = PHASE_CHANGE_REASON_NONE
                    iGEMStagnationCount = 0
                    lUseCEFSiteGEM = .TRUE.
                    lGEMCEFSiteLagrangianActive = .TRUE.
                    call RecordGEMIterationDiagnostics
                    cycle LOOP_GEMSolver
                end if
                lGEMCEFSiteLagrangianEnabled = .FALSE.
            end if

            iPhaseChangeReason = PHASE_CHANGE_REASON_GEM_STAGNATION
            iPhaseChangeReasonHistory(iterGlobal) = iPhaseChangeReason
            if (lPostProcess) then
                lPhaseChange = .FALSE.
                lPhaseChangeHistory(iterGlobal) = .FALSE.
            else
                lPhaseChange = .TRUE.
                lPhaseChangeHistory(iterGlobal) = .TRUE.
            end if
            call RecordGEMIterationDiagnostics
            exit LOOP_GEMSolver
        end if

        call RecordGEMIterationDiagnostics

    end do LOOP_GEMSolver

    lHitLagrangianMax = (.NOT.lConverged).AND.(.NOT.lPhaseChange).AND.&
        (INFOThermo == 0).AND.(iterGlobal > iterCurrentMax)

    if (lDebugMode.AND.lHitLagrangianMax) then
        print *, 'WARNING: Lagrangian GEM hit max iterations ', &
            ' iterPEA=', iterPEA, &
            ' iterMax=', iterCurrentMax, &
            ' polishActive=', lPEALagrangianPolishActive, &
            ' useCEFSite=', lUseCEFSiteGEM, &
            ' T=', dTemperature, &
            ' P=', dPressure, &
            ' dGEMFunctionNorm=', dGEMFunctionNorm, &
            ' massNorm=', dGEMMassBalanceNorm, &
            ' chemNorm=', dGEMChemicalPotentialNorm, &
            ' solnChemNorm=', dGEMSolutionChemicalPotentialNorm, &
            ' condChemNorm=', dGEMCondensedChemicalPotentialNorm, &
            ' sublatticeExchangeNorm=', dSublatticeExchangeNorm
        print *, 'WARNING: Lagrangian GEM line search ', &
            ' iterations=', iGEMLineSearchIterationCount, &
            ' negativeFactors=', iGEMLineSearchNegativeFactorCount, &
            ' negativePhases=', iGEMLineSearchNegativePhaseCount, &
            ' floorCount=', iGEMLineSearchFloorCount, &
            ' noDescent=', iGEMLineSearchNoDescent, &
            ' noDescentClass=', iGEMLineSearchNoDescentClass, &
            ' initialNorm=', dGEMLineSearchInitialNorm, &
            ' bestNorm=', dGEMLineSearchBestNorm, &
            ' finalNorm=', dGEMLineSearchFinalNorm, &
            ' bestStep=', dGEMLineSearchBestStep, &
            ' initialG=', dGEMLineSearchInitialGibbs, &
            ' bestG=', dGEMLineSearchBestGibbs, &
            ' dirNormSlope=', dGEMNewtonDirNormSlope, &
            ' dirGSlope=', dGEMNewtonDirGibbsSlope
        print *, 'WARNING: Lagrangian GEM Newton diagnostics ', &
            ' kktSize=', iGEMNewtonKKTSize, &
            ' dsysvInfo=', iGEMNewtonDSYSVInfo, &
            ' pivot1x1=', iGEMNewtonPivot1x1Count, &
            ' pivot2x2=', iGEMNewtonPivot2x2Count, &
            ' pivotPositive=', iGEMNewtonPivotPositiveCount, &
            ' pivotNegative=', iGEMNewtonPivotNegativeCount, &
            ' pivotZero=', iGEMNewtonPivotZeroCount, &
            ' minPivot=', dGEMNewtonMinPivotScale, &
            ' directionNorm=', dGEMNewtonDirectionNorm
        print *, 'WARNING: Lagrangian GEM scaffold diagnostics ', &
            ' residualLM=', iGEMResidualLMUsed, &
            ' coalesce=', iGEMCoalesceUsed, &
            ' stabilize=', iGEMStabilizeUsed, &
            ' boundaryRemoval=', iGEMBoundaryRemovalUsed, &
            ' rawNegativeRemoval=', iGEMRawNegativeRemovalUsed, &
            ' boundaryPinnedRemoval=', iGEMBoundaryPinnedRemovalUsed, &
            ' traceRemove=', iGEMTraceRemoveUsed, &
            ' traceReinject=', iGEMTraceReinjectUsed, &
            ' cefRetry=', iGEMCEFRetryActivationUsed, &
            ' subgqRideAlong=', iGEMSUBGQRideAlongUsed
    end if

    if ((.NOT.lConverged).AND.(.NOT.lPhaseChange).AND.&
        (iPhaseChangeReason == PHASE_CHANGE_REASON_NONE)) then
        iPhaseChangeReason = PHASE_CHANGE_REASON_LAGRANGIAN_UNCONVERGED
        iGEMExitStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
        iterReason = MIN(MAX(iterGlobal,1),iterGlobalMax)
        iPhaseChangeReasonHistory(iterReason) = iPhaseChangeReason
    end if
    if ((.NOT.lConverged).AND.(dGEMFunctionNorm > 1D-5)) then
        iGEMExitStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
    end if

    lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM
    lGEMCEFBndPhaseActive = .FALSE.
    iGEMCEFBndPhaseSlot = 0
    dGEMCEFBndPhaseStep = 0D0

    if (allocated(dElementPotentialLast)) deallocate(dElementPotentialLast)
    if (allocated(dMolesSpeciesLast)) deallocate(dMolesSpeciesLast)
    if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)

    return

contains

    logical function DistinctRepeatedSUBOMCompositionSetsActive()

        implicit none

        integer :: iSlotA, iSlotB, iPhaseA, iPhaseB
        real(8) :: dSiteDifference
        real(8), parameter :: dSameSetTolerance = 1D-8

        DistinctRepeatedSUBOMCompositionSetsActive = .FALSE.
        if (.NOT.lSUBOMTwoSetCandidateEnabled) return
        if (.NOT.allocated(iActiveSlotThermoPhase)) return
        if (.NOT.allocated(dActiveSlotSiteFraction)) return

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
                if (dSiteDifference > dSameSetTolerance) then
                    DistinctRepeatedSUBOMCompositionSetsActive = .TRUE.
                    return
                end if
            end do
        end do

        return

    end function DistinctRepeatedSUBOMCompositionSetsActive

    logical function CurrentGEMCEFInertiaWrong()

        implicit none

        integer :: nPrimalLocal, nConstraintLocal

        nPrimalLocal = nGEMCEFPhaseVariables + nGEMCEFSiteVariables
        nConstraintLocal = nElements
        CurrentGEMCEFInertiaWrong = (iGEMNewtonDSYSVInfo /= 0).OR.&
            (iGEMNewtonPivotZeroCount /= 0).OR.&
            (iGEMNewtonPivotPositiveCount /= nPrimalLocal).OR.&
            (iGEMNewtonPivotNegativeCount /= nConstraintLocal)

        return

    end function CurrentGEMCEFInertiaWrong

    subroutine ComputeInvalidCompBoundPhi(iSlotIn, dPhiOut, iInfoOut)

        implicit none

        integer, intent(in) :: iSlotIn
        real(8), intent(out) :: dPhiOut
        integer, intent(out) :: iInfoOut

        real(8) :: dPhaseComposition(nElements)

        dPhiOut = 0D0
        iInfoOut = 0
        call ComputeCertificatePhaseComposition(iSlotIn, dPhaseComposition, iInfoOut)
        if (iInfoOut /= 0) return

        dPhiOut = dGEMRawNegPhaseResidual - &
            SUM(dGEMCEFElementDirection(1:nElements) * dPhaseComposition(1:nElements))

        return

    end subroutine ComputeInvalidCompBoundPhi

    subroutine ComputeCertificatePhaseComposition(iSlotIn, dCompositionOut, iInfoOut)

        implicit none

        integer, intent(in) :: iSlotIn
        real(8), intent(out) :: dCompositionOut(nElements)
        integer, intent(out) :: iInfoOut

        integer :: iSolnPhaseLocal, iThermoPhaseLocal, iPhaseIDLocal
        real(8) :: dSiteLocal(nMaxSublatticeSys,nMaxConstituentSys)

        dCompositionOut = 0D0
        iInfoOut = 0
        if ((iSlotIn <= 0).OR.(iSlotIn > nElements)) then
            iInfoOut = 1
            return
        end if

        if (iAssemblage(iSlotIn) > 0) then
            dCompositionOut = dStoichSpecies(iAssemblage(iSlotIn),1:nElements)
            return
        end if

        iSolnPhaseLocal = -iAssemblage(iSlotIn)
        iThermoPhaseLocal = iSolnPhaseLocal
        if (allocated(iActiveSlotThermoPhase)) then
            if ((iSlotIn <= SIZE(iActiveSlotThermoPhase)).AND.&
                (iActiveSlotThermoPhase(iSlotIn) > 0)) then
                iThermoPhaseLocal = iActiveSlotThermoPhase(iSlotIn)
            end if
        end if
        if ((iThermoPhaseLocal <= 0).OR.(iThermoPhaseLocal > nSolnPhasesSys)) then
            iInfoOut = 1
            return
        end if

        if (IsCertificateCEFPhase(iThermoPhaseLocal)) then
            iPhaseIDLocal = iPhaseSublattice(iThermoPhaseLocal)
            call CertificateSlotSiteFractions(iSlotIn, dSiteLocal, iInfoOut)
            if (iInfoOut /= 0) return
            call ComputeCertificateCEFComposition(iThermoPhaseLocal, iPhaseIDLocal, &
                dSiteLocal, dCompositionOut)
        else
            call CompStoichSolnPhase(iThermoPhaseLocal)
            dCompositionOut = dEffStoichSolnPhase(iThermoPhaseLocal,1:nElements)
        end if

        return

    end subroutine ComputeCertificatePhaseComposition

    subroutine CertificateSlotSiteFractions(iSlotIn, dSiteOut, iInfoOut)

        implicit none

        integer, intent(in) :: iSlotIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(out) :: iInfoOut

        integer :: iPhaseVarLocal

        dSiteOut = 0D0
        iInfoOut = 0
        if (allocated(dActiveSlotSiteFraction)) then
            if ((iSlotIn <= SIZE(dActiveSlotSiteFraction,1)).AND.&
                (SUM(dActiveSlotSiteFraction(iSlotIn,:,:)) > 0D0)) then
                dSiteOut = dActiveSlotSiteFraction(iSlotIn,:,:)
                return
            end if
        end if

        do iPhaseVarLocal = 1, nGEMCEFPhaseVariables
            if (iGEMCEFPhaseSlot(iPhaseVarLocal) == iSlotIn) then
                dSiteOut = dGEMCEFPhaseSiteLast(iPhaseVarLocal,:,:)
                return
            end if
        end do

        iInfoOut = 1
        return

    end subroutine CertificateSlotSiteFractions

    subroutine ComputeCertificateCEFComposition(iSolnPhaseIn, iPhaseIDIn, dSiteIn, dCompositionOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dCompositionOut(nElements)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dProductSum
        real(8), allocatable :: dProductSpecies(:)

        dCompositionOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIn)
        allocate(dProductSpecies(iLastLocal-iFirstLocal+1))
        dProductSpecies = 0D0
        dProductSum = 0D0

        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 1D-75)
            end do
            dProductSpecies(iLocal) = dProduct
            dProductSum = dProductSum + dProduct
        end do

        if (dProductSum > 0D0) then
            do iLocalSpecies = iFirstLocal, iLastLocal
                iLocal = iLocalSpecies - iFirstLocal + 1
                dCompositionOut = dCompositionOut + &
                    dProductSpecies(iLocal) * dStoichSpecies(iLocalSpecies,1:nElements) / dProductSum
            end do
        end if

        if (allocated(dProductSpecies)) deallocate(dProductSpecies)

        return

    end subroutine ComputeCertificateCEFComposition

    logical function IsCertificateCEFPhase(iSolnPhaseIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIn

        IsCertificateCEFPhase = .FALSE.
        if ((iSolnPhaseIn <= 0).OR.(iSolnPhaseIn > nSolnPhasesSys)) return
        if (.NOT.allocated(iPhaseSublattice)) return
        if (iSolnPhaseIn > SIZE(iPhaseSublattice)) return
        if (iPhaseSublattice(iSolnPhaseIn) <= 0) return

        IsCertificateCEFPhase = (TRIM(cSolnPhaseType(iSolnPhaseIn)) == 'SUBL').OR.&
            (TRIM(cSolnPhaseType(iSolnPhaseIn)) == 'SUBLM').OR.&
            (TRIM(cSolnPhaseType(iSolnPhaseIn)) == 'SUBOM')

        return

    end function IsCertificateCEFPhase

end subroutine RunLagrangianGEM
