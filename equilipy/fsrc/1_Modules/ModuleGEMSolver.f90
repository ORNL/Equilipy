!-------------------------------------------------------------------------------------------------------------
!
!> \file    ModuleGEMSolver.f90
!> \brief   Shared state for Leveling, Lagrangian GEM, and CheckPhaseAssemblage.
!> \author  M.H.A. Piro
!> \date    Apr. 25, 2012
!> \sa      MultiPhaseMinimizer.f90
!> \sa      RunLagrangianGEM.f90
!> \sa      CheckPhaseAssemblage.f90
!> \sa      TraceSpeciesControl.f90
!> \sa      CompFunctionNorm.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   04/25/2012      M.H.A. Piro         Original GEM solver state module
!   06/23/2026      S.Y. Kwon           Added phase-change reasons and trace-species controls
!   06/23/2026      S.Y. Kwon           Collapsed immutable trace/stagnation defaults into parameters
!   06/23/2026      S.Y. Kwon           Organized declarations by rank and role
!   06/23/2026      S.Y. Kwon           Consolidated GEM active-set and PEA candidate state
!   06/24/2026      S.Y. Kwon           Added line-search diagnostics for Lagrangian stagnation analysis
!   06/24/2026      S.Y. Kwon           Added per-iteration Lagrangian residual and line-search histories
!   06/25/2026      S.Y. Kwon           Removed legacy SUBOM residual state after moving order/disorder
!                                        refinement to site-fraction minimization
!   06/25/2026      S.Y. Kwon           Renamed trace-species enable flag for production use
!   06/25/2026      S.Y. Kwon           Added CEF site-fraction Lagrangian direction state
!   06/25/2026      S.Y. Kwon           Kept CEF site-fraction Lagrangian disabled by default pending
!                                        multiphase KKT validation.
!   06/25/2026      S.Y. Kwon           Added raw CEF phase-removal target state for Lagrangian active-set repair.
!   06/25/2026      S.Y. Kwon           Added an active CEF KKT flag so residual norms match the solver
!                                        coordinates used for the current fixed assemblage.
!   06/26/2026      S.Y. Kwon           Added Level2Lagrange handoff diagnostics for the embedded PEA-Lagrangian
!                                        minimizer implementation path.
!   06/26/2026      S.Y. Kwon           Added PEA-internal Lagrangian polish controls and diagnostics.
!   06/26/2026      S.Y. Kwon           Clarified that PEA-polish accepts converged Lagrangian planes
!                                        without bounding them to the PEA Leveling plane.
!   06/26/2026      S.Y. Kwon           Removed inactive sampled/grid repair controls from the production
!                                        minimizer state.
!   06/26/2026      S.Y. Kwon           Added per-iteration PEA diagnostics for Lagrangian-polished
!                                        active-set development.
!   06/26/2026      S.Y. Kwon           Added explicit Subminimization candidate status so PEA Leveling can
!                                        reject unknown pseudo-compound rows.
!   06/27/2026      S.Y. Kwon           Generalized negative driving-force witnesses for PEA candidate rows.
!   06/28/2026      S.Y. Kwon           Added an elemental-potential plateau gate before PEA-Lagrangian polish.
!   06/30/2026      S.Y. Kwon           Tightened GEM stagnation detection to avoid premature active-set repair.
!   07/01/2026      S.Y. Kwon           Documented that Leveling/PEA endmember vectors are handoff state,
!                                        while active CEF directions are site-fraction KKT state.
!   07/02/2026      S.Y. Kwon           Added passive KKT, pivot, and scaffold-activation diagnostics
!                                        for the minimizer Phase 0 baseline.
!   07/02/2026      S.Y. Kwon           Added passive active order/disorder pair diagnostics for
!                                        A2/B2 representation-null debugging.
!   07/02/2026      S.Y. Kwon           Added passive active-slot model/display/identity arrays for
!                                        SUBOM helper-safe representation work.
!   07/02/2026      S.Y. Kwon           Added Phase 0 directional-derivative and line-search energy
!                                        diagnostics, and split scaffold counters by activation type.
!   07/02/2026      S.Y. Kwon           Renamed the active-slot composition-set placeholder to an
!                                        identity ordinal until real parent composition sets are built.
!   07/02/2026      S.Y. Kwon           Added passive per-active-slot constitution storage for the
!                                        SUBOM helper-safe identity scaffold.
!   07/02/2026      S.Y. Kwon           Split passive active-slot thermodynamic parent identity from
!                                        display phase identity for order/disorder scaffolding.
!   07/02/2026      S.Y. Kwon           Added passive per-active-slot parent-site-fraction storage.
!   07/03/2026      S.Y. Kwon           Added pre-residual-LM CEF line-search diagnostics for
!                                        class-1 no-descent trials.
!   07/03/2026      S.Y. Kwon           Added a default-off SUBOM two-composition-set handoff
!                                        switch and diagnostics for Phase 2e-c slice c2.
!   07/04/2026      S.Y. Kwon           Split tiny-boundary phase-removal diagnostics from
!                                        boundary-pinned removals for C1-a passive census.
!   07/04/2026      S.Y. Kwon           Added default-off standalone SUBG/SUBQ CEF-routing switch
!                                        for C1-b0 plumbing.
!   07/04/2026      S.Y. Kwon           Promoted standalone SUBG/SUBQ mixed-coordinate routing to
!                                        the production default after C1 validation.
!   07/04/2026      S.Y. Kwon           Added call-site-stable residual-LM event storage for C3-a2
!                                        globalization diagnostics.
!   07/04/2026      S.Y. Kwon           Added C3-c1 primal-inertia regularization counters.
!   07/04/2026      S.Y. Kwon           Added default-off C3-c2 funnel line-search switch and counters.
!   07/05/2026      S.Y. Kwon           Added passive MultiPhaseMinimizer timing buckets for Scheil
!                                        warm-start expected-value census.
!   07/05/2026      S.Y. Kwon           Added passive inactive-candidate certification timing and
!                                        sweep count for WS-a2.
!   07/05/2026      S.Y. Kwon           Added explicit GEM exit status and max CEF exchange residual
!                                        diagnostics for trace-phase Scheil failures.
!
!
! Purpose:
! ========
!
!> \details This module stores global minimizer state shared by initial
!! Leveling, Lagrangian GEM, CheckPhaseAssemblage/PEA, trace-species control,
!! and postprocess fixed-assemblage solves.  It is runtime solver state, not
!! parser/database state.
!
!
! Required input variables:
! =========================
!
! lCompbdOnly                    True when only compound thermodynamic properties are requested.
! lPostProcess                   True during fixed-assemblage Cp perturbation solves.
! lSolnPhases                    Solution phase active flags for the current system.
! lMiscibility                   Miscibility-gap flags for solution phases.
! iAssemblageGEMinit             Initial assemblage snapshot for Lagrangian translation.
! dMolesPhaseLast                Previous phase amounts used by line search/update diagnostics.
! dMolesSpeciesLast              Previous species amounts used by line search/update diagnostics.
! dElementPotentialLast          Previous elemental potentials used by line search/update diagnostics.
!
!
! Output/updated variables:
! =========================
!
! iterGlobal                     Current Lagrangian GEM iteration.
! iterPEA                        Current CheckPhaseAssemblage/PEA iteration.
! iterHistory                    Phase assemblage history.
! lConverged                     True when Lagrangian GEM satisfies convergence criteria.
! lPhaseChange                   True when the current active set must be repaired by PEA.
! iPhaseChangeReason             Reason code for the latest active-set repair request.
! lPhaseChangeHistory            Per-iteration phase-change flag history.
! iPhaseChangeReasonHistory      Per-iteration phase-change reason history.
! dGEMFunctionNorm               Combined GEM residual norm.
! dGEMFunctionNormLast           Previous GEM residual norm.
! dGEMMassBalanceNorm            Mass-balance contribution to the GEM residual norm.
! dGEMChemicalPotentialNorm      Chemical-potential contribution to the GEM residual norm.
! dGEMSolutionChemicalPotentialNorm
!                                 Solution-species stationarity contribution to the chemical norm.
! dGEMCondensedChemicalPotentialNorm
!                                 Pure condensed phase stationarity contribution to the chemical norm.
! dSublatticeExchangeNorm        Order/disorder sublattice exchange residual norm.
! dUpdateVar                     Lagrangian Newton update vector.
! dChemicalPotentialGEM          GEM active-set chemical potentials translated from Leveling/PEA.
! dStoichSpeciesGEM              GEM active-set stoichiometry translated from Leveling/PEA.
! dAtomFractionSpeciesGEM        GEM active-set atom fractions translated from Leveling/PEA.
! dMolFractionGEM                GEM active-set solution compositions translated from Leveling/PEA.
!                                 For CEF phases this is only a handoff/candidate representation.
!                                 GEMNewtonCEF converts it to site fractions before building directions.
! dLevelCandidateMolFraction     Optional PEA candidate compositions retained for Leveling rows.
! iLevelCandidateParentPhase     Thermodynamic parent phase id for each candidate-pool row.
! iLevelCandidateDisplayPhase    Display/reporting phase id for each candidate-pool row.
! iLevelCandidateIdentityOrdinal Composition-set ordinal under the parent phase for each candidate row.
! iGEMLineSearchIterationCount   Number of damping trials in the latest Lagrangian line search.
! iGEMLineSearchNegativeFactorCount
!                                 Number of raw species update factors below zero before lambda damping.
! iGEMLineSearchFloorCount       Number of raw full-step species amounts below the numerical mole floor.
! iGEMLineSearchNoDescent        One when the latest line search could not reduce the starting norm.
! iGEMLineSearchNoDescentClass   1 when residual-only rejection still lowered G; 2 when G also failed to descend.
! iGEMPreLMNoDescent            Preserves the CEF line-search no-descent flag before residual-LM retry.
! iGEMPreLMNoDescentClass       Preserves the CEF no-descent class before residual-LM retry.
! iGEMNewtonKKTSize              Dimension of the latest Lagrangian Newton linear system.
! iGEMNewtonPivot*Count          Bunch-Kaufman pivot/inertia summary from the latest DSYSV factorization.
! iGEMInertiaRegularization*    C3-c1 counters for primal-block Levenberg retries when KKT inertia is wrong.
! iGEM*Used                     Per-iteration scaffold activation counters recorded for Phase 0 auditing.
! iGEMResidualLM*Used           Split per-iteration residual-LM activation counters by raw CEF phase-boundary
!                                 evidence and post-line-search no-descent evidence.
! iGEMResidualLM*TotalUsed      Per-minimization residual-LM activation totals exposed to Python tests.
! iGEMResidualLMEvent*          Call-site-stable residual-LM event buffer for C3-a2 passive diagnostics.
! iGEMTrial7BoundRetry*Total    Per-minimization count of bound-active KKT retries attempted before
!                                 residual-LM fallback.
! iGEMInvalidCompBound*Total    Per-minimization count of invalid-complementarity raw-negative CEF
!                                 phase events handled by a recomputed bound-active KKT multiplier.
! iGEM*CoalesceTotalUsed       Per-minimization duplicate/degenerate active-set coalescer totals.
! iGEMStabilizeTotalUsed       Per-minimization SUBOM endpoint-stability restart total.
! iODPairCount
!                                 Number of active ordered/DIS_PART helper pairs in the current GEM set.
! dGEMLineSearchInitialNorm      GEM norm before the latest line search.
! dGEMLineSearchBestNorm         Best GEM norm found during the latest line search.
! dGEMLineSearchFinalNorm        GEM norm after restoring the best line-search state.
! dGEMLineSearchInitialStep      Initial step length after lambda correction.
! dGEMLineSearchBestStep         Step length that produced dGEMLineSearchBestNorm.
! dGEMLineSearchFinalStep        Last step length evaluated by the latest line search.
! dGEMLineSearch*Gibbs           Internal active-set Gibbs energy diagnostics for line-search classification.
! dGEMLineSearch*Merit           Gibbs plus mass-balance penalty diagnostics for line-search classification.
! dGEMNewtonDir*Slope            First feasible line-search slope along the latest Newton direction.
! dGEMNormHist        Per-Lagrangian-iteration combined GEM norm.
! dGEMMassNormHist     Per-Lagrangian-iteration mass-balance norm.
! dGEMChemNormHist
!                                 Per-Lagrangian-iteration chemical-potential norm.
! dGEMLSRawPhaseMolesHist
!                                 Raw full-step phase amounts before lambda damping.
! dGEMLSFinalPhaseMolesHist
!                                 Phase amounts after accepted line search and trace handling.
! iGEMRawNegativePhaseSlotHist
!                                 Per-iteration CEF raw-negative slot selected before residual-LM fallback.
! iGEMInvalidCompBoundVerdictHist
!                                 Per-iteration invalid-complementarity lower-bound certificate:
!                                 1 accepted, -1 rejected by phi, -2 rejected by KKT failure.
! dGEMRawNegativePhase*Hist
!                                 Per-iteration raw-negative phase amount and direction selected by the
!                                 CEF phase-boundary audit before any fallback direction overwrites
!                                 line-search raw-mole diagnostics.
! dGEMInvalidCompBoundPhiHist
!                                 Per-iteration recomputed phase residual under the bound-active
!                                 KKT element-potential update.
! dODCompDist
!                                 Maximum normalized composition distance across active ordered/helper pairs.
! dODOrderNorm  Maximum ordered-phase site-fraction departure from its random manifold.
! iODOrderingModeCount
!                                 Number of active SUBOM phases with a passive projected ordering-mode Hessian.
! dODOrderingEigenMin
!                                 Smallest passive fixed-average ordering-mode Hessian eigenvalue across
!                                 active SUBOM phases.
! lGEMCEFSiteLagrangianActive     True when the current fixed assemblage uses the CEF site KKT path.
! lGEMCEFSiteDirectionActive      True when the latest Lagrangian Newton step uses CEF site variables.
! lGEMCEFInertiaRegularizationActive
!                                 True only while C3-c1 tests a primal-regularized KKT direction.
! lGEMCEFInertiaRegularizationEnabled
!                                 Internal default-off switch for the C3-c1 inertia-regularization
!                                 prototype.  The first behavior attempt was rejected because it
!                                 moved the row-458 thermodynamic state.
! dGEMCEF*Direction               Mixed-coordinate CEF Newton directions for phase amount, site fraction,
!                                 and elemental potentials.  These arrays, not dMolFractionGEM, own
!                                 the active CEF minimizer direction.
! dGEMAnalyticalSpeciesDirection  Additive species-mole direction from analytical non-CEF solution Hessians.
! iGEMRawNegativePhase*           Raw CEF phase-removal target used by Lagrangian GEM active-set repair.
! nLevel2Lagrange*                Counters describing the latest Leveling-to-Lagrangian candidate handoff.
! iLevel2Lagrange*                Phase, candidate, source, and assemblage identities from the latest handoff.
! dLevel2Lagrange*                Phase amounts and elemental potentials from the latest handoff.
! lSUBOMTwoSetCandidateEnabled    Default-off Phase 2e-c/c2 switch that permits a DIS_PART helper
!                                 candidate to become a second composition-set slot of its ordered
!                                 SUBOM parent when the c1 ordering-mode diagnostic is unstable.
! lSUBQStandaloneEnabled          Production default permitting standalone SUBG/SUBQ active sets
!                                 to use the mixed-coordinate Lagrangian path after the SUBQ/SUBG
!                                 derivative-basis validation and C1 routing gates.
! iSUBOMOrderingGate*             Per-PEA-iteration cache for the expensive projected ordering-mode
!                                 eigenvalue used only by the two-set candidate gate.
! nLevel2LagrangeTwoSetCreated    Number of second SUBOM composition-set handoff slots created
!                                 during the latest Level2Lagrange call.
! iActiveSlotThermoPhase          Thermodynamic model phase id for each active Lagrangian slot.
! iActiveSlotDisplayPhase         Display/reporting phase id for each active Lagrangian slot.
! iActiveSlotIdentityOrdinal      Passive per-thermodynamic-parent composition-set ordinal for each active slot.
! dActiveSlotMolFraction          Passive per-active-slot solution endmember/product-fraction snapshot.
! dActiveSlotSiteFraction         Passive per-active-slot parent-model site-fraction snapshot.
! lPEALagrangianPolish*           Flags for embedded Lagrangian polishing inside PEA.
! nPEALagrangianPolish*           Counters and iteration diagnostics for PEA-internal polishing.
! dPEALagrangianPolish*           Norm and elemental-potential diagnostics for PEA-internal polishing.
! iSubMinCandidateStatusSoln      Per-solution validity/status of the latest subminimized PEA candidate row.
! iPEA*Hist                       Per-PEA-iteration assemblage, Leveling, and polish status diagnostics.
! dPEA*Hist                       Per-PEA-iteration potential, driving-force, and polish norm diagnostics.
!
!
! Numerical parameters:
! =====================
!
! iterGlobalMax                       Maximum Lagrangian GEM iterations per fixed active set.
! PHASE_CHANGE_REASON_*               Named active-set repair trigger codes.
! dTraceSpeciesRemoveFraction         Trace endmember removal threshold.
! dTraceSpeciesReinjectFraction       Reintroduced trace endmember fraction.
! dTraceSpeciesResidualTolerance      Trace residual tolerance for convergence checks.
! dTraceSpeciesSlowProgressTolerance  Slow-progress threshold that triggers reinjection.
! dTraceSpeciesReducedSolveTolerance  Reduced active-set convergence threshold.
! iTraceSpeciesSlowProgressWindow     Consecutive slow-progress count before reinjection.
! iTraceSpeciesMinIterBeforeReinject  Minimum wait after removal/reinjection.
! dGEMStagnationRelativeTolerance     Relative GEM norm progress threshold for stagnation.
! dGEMStagnationNormTolerance         GEM norm floor below which stagnation is ignored.
! iGEMStagnationWindow                Consecutive stagnation count before PEA repair.
! dPEALagrangianPolishMaxElementPotential
!                                    Loose elemental-potential plateau threshold before embedded PEA polish.
! dPEALagrangianPolishMaxDrivingForce
!                                    Driving-force threshold before embedded PEA polish is attempted.
! nPEAGate* / iPEAGate*
!                                    C4-a passive counters for the PEA-polish elemental-potential gate.
! nPEADirectHandoff*          C4-a passive counters for direct inactive-witness PEA handoffs.
! nPEARepeatExit*           C4-a passive counters for repeated-handoff PEA exits.
! nPEARepeatGuard        C4-a context counter for the KEEP repeated-Leveling-assemblage guard.
! dGEMTiming*             Per-minimization passive CPU-time buckets for Leveling, direct
!                         Level2Lagrange handoff, CheckPhaseAssemblage/PEA, and
!                         fixed-active-set Lagrangian GEM.
! nGEMCertificationSweep  Per-minimization count of inactive-candidate driving-force
!                         refresh/certification sweeps.
!
!
! Primary callers/users:
! ======================
!
! MultiPhaseMinimizer       Owns the top-level minimization sequence.
! RunLagrangianGEM          Uses Newton/line-search state and phase-change diagnostics.
! CheckPhaseAssemblage      Uses phase-change reason state and PEA iteration state.
! TraceSpeciesControl       Removes/reinjects trace endmembers during Lagrangian GEM.
! SubMinTraceSpeciesControl Reuses trace thresholds inside PEA subminimization.
! CompFunctionNorm          Updates dGEMFunctionNorm and trace-aware residuals.
! ResetThermo               Resets mutable solver state between calculations.
!
!
! Numerical assumptions:
! ======================
!
! - Trace-species thresholds are fixed numerical parameters.  They are not
!   runtime tunables and should not be reset through separate default values.
! - Lagrangian GEM normally requests PEA repair by setting lPhaseChange and
!   iPhaseChangeReason; PEA owns phase discovery and replacement.
! - CEF Lagrangian may directly remove an active phase during normal final GEM
!   solves when the raw Newton phase-amount direction crosses the zero boundary
!   before line-search damping.  Embedded PEA-polish solves reject that polish
!   instead and leave phase evolution to PEA.
! - Postprocess Cp perturbation runs with lPostProcess true and should not
!   trigger trace handling or phase assemblage changes.
!
!
!-------------------------------------------------------------------------------------------------------------
module ModuleGEMSolver

    implicit none

    SAVE

    ! Integer parameters:
    integer, parameter                         :: PHASE_CHANGE_REASON_NONE = 0
    integer, parameter                         :: PHASE_CHANGE_REASON_LEVELING_INITIAL = 1
    integer, parameter                         :: PHASE_CHANGE_REASON_NEGATIVE_PHASE_AMOUNT = 2
    integer, parameter                         :: PHASE_CHANGE_REASON_GEM_STAGNATION = 3
    integer, parameter                         :: PHASE_CHANGE_REASON_PHASE_POTENTIAL = 4
    integer, parameter                         :: PHASE_CHANGE_REASON_DUMMY_PHASE = 5
    integer, parameter                         :: PHASE_CHANGE_REASON_NEWTON_FAILURE = 6
    integer, parameter                         :: PHASE_CHANGE_REASON_LAGRANGIAN_UNCONVERGED = 7
    integer, parameter                         :: GEM_EXIT_STATUS_OK = 0
    integer, parameter                         :: GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED = 1
    integer, parameter                         :: SUBMIN_CANDIDATE_UNKNOWN = 0
    integer, parameter                         :: SUBMIN_CANDIDATE_CONVERGED = 1
    integer, parameter                         :: SUBMIN_CANDIDATE_NEGATIVE_WITNESS = 2
    integer, parameter                         :: SUBMIN_CANDIDATE_REJECTED = 3
    integer, parameter                         :: SUBMIN_CANDIDATE_DUPLICATE = 4
    integer, parameter                         :: SUBMIN_CANDIDATE_MAX_ITER = 5
    integer, parameter                         :: SUBOM_TRACE_REGISTER = 1
    integer, parameter                         :: SUBOM_TRACE_REFRESH = 2
    integer, parameter                         :: SUBOM_TRACE_HANDOFF = 3
    integer, parameter                         :: SUBOM_TRACE_SOLVER_SLOT = 4
    integer, parameter                         :: iterGlobalMax = 1000
    integer, parameter                         :: iTraceSpeciesSlowProgressWindow = 3
    integer, parameter                         :: iTraceSpeciesMinIterBeforeReinject = 3
    integer, parameter                         :: iGEMStagnationWindow = 3
    integer, parameter                         :: iterPEALagrangianPolishMax = 25
    integer, parameter                         :: iterPEAMax = 100

    ! Real parameters:
    real(8), parameter                         :: dTraceSpeciesRemoveFraction = 1D-30
    real(8), parameter                         :: dTraceSpeciesReinjectFraction = 1D-10
    real(8), parameter                         :: dTraceSpeciesResidualTolerance = 1D-3
    real(8), parameter                         :: dTraceSpeciesSlowProgressTolerance = 5D-3
    real(8), parameter                         :: dTraceSpeciesReducedSolveTolerance = 1D-5
    real(8), parameter                         :: dGEMStagnationRelativeTolerance = 1D-15
    real(8), parameter                         :: dGEMStagnationNormTolerance = 1D-5
    real(8), parameter                         :: dPEALagrangianPolishMaxElementPotential = 5D-1
    real(8), parameter                         :: dPEALagrangianPolishMaxDrivingForce = 1D-4


    ! Integer scalars:
    integer                                    :: iterLast, iterStep, iterRevert, iterGlobal, iterPEA, iterLG
    integer                                    :: iterLastCon, iterLastSoln, iterSwap, iterLastMiscGapCheck
    integer                                    :: iterTraceSpeciesLastRemoval, iterTraceSpeciesLastReinject
    integer                                    :: iConPhaseLast, iSolnPhaseLast, iSolnSwap, iPureConSwap
    integer                                    :: iMinDrivingForceStoich, iMinDrivingForceSoln, iSpeciesRemove
    integer                                    :: iPhaseChangeReason
    integer                                    :: iGEMExitStatus
    integer                                    :: iGEMNewtonSolver, iGEMNewtonDSYSVInfo
    integer                                    :: iGEMNewtonKKTSize
    integer                                    :: iGEMNewtonPivot1x1Count, iGEMNewtonPivot2x2Count
    integer                                    :: iGEMNewtonPivotPositiveCount, iGEMNewtonPivotNegativeCount
    integer                                    :: iGEMNewtonPivotZeroCount
    integer                                    :: iGEMMaxSublExchangeSlot
    integer                                    :: iGEMMaxSublExchangePhase
    integer                                    :: iGEMMaxSublExchangeSite
    integer                                    :: iGEMMaxSublExchangeConstituent
    integer                                    :: iGEMInertiaRegularizationAttemptedUsed
    integer                                    :: iGEMInertiaRegularizationAcceptedUsed
    integer                                    :: iGEMInertiaRegularizationFailedUsed
    integer                                    :: iGEMInertiaRegularizationAttemptedTotal
    integer                                    :: iGEMInertiaRegularizationAcceptedTotal
    integer                                    :: iGEMInertiaRegularizationFailedTotal
    integer                                    :: iGEMInertiaRegularizationStepTotal
    integer                                    :: iGEMInertiaRegularizationStepLast
    integer                                    :: iTraceSpeciesSlowProgressCount
    integer                                    :: iGEMStagnationCount
    integer                                    :: iGEMLineSearchIterationCount
    integer                                    :: iGEMLineSearchNegativeFactorCount
    integer                                    :: iGEMLineSearchNegativePhaseCount
    integer                                    :: iGEMLineSearchFloorCount
    integer                                    :: iGEMLineSearchNoDescent
    integer                                    :: iGEMLineSearchNoDescentClass
    integer                                    :: iGEMPreLMNoDescent
    integer                                    :: iGEMPreLMNoDescentClass
    integer                                    :: iGEMPreLMNewtonInfo
    integer                                    :: iGEMPreLMKKTSize
    integer                                    :: iGEMPreLMPivot1x1
    integer                                    :: iGEMPreLMPivot2x2
    integer                                    :: iGEMPreLMPivotPositive
    integer                                    :: iGEMPreLMPivotNegative
    integer                                    :: iGEMPreLMPivotZero
    integer                                    :: iGEMPreLMODOrdPhase
    integer                                    :: iGEMPreLMODCompPhase
    integer                                    :: iGEMAnalyticalHessianFallbackCount
    integer                                    :: iGEMCEFResidualLMFallbackCount
    integer                                    :: iGEMResidualLMUsed
    integer                                    :: iGEMResidualLMRawNegativeUsed
    integer                                    :: iGEMResidualLMNoDescentUsed
    integer                                    :: iGEMResidualLMNoDescentClass
    integer                                    :: iGEMResidualLMTotalUsed
    integer                                    :: iGEMResidualLMRawNegativeTotalUsed
    integer                                    :: iGEMResidualLMNoDescentTotalUsed
    integer                                    :: iGEMLagrangianCallSite
    integer                                    :: nGEMResidualLMEvent = 0
    integer                                    :: nGEMResidualLMEventCapacity = 0
    integer                                    :: iGEMTrial7BoundRetryAttemptedTotal
    integer                                    :: iGEMTrial7BoundRetryAcceptedTotal
    integer                                    :: iGEMInvalidCompBoundAttemptedUsed
    integer                                    :: iGEMInvalidCompBoundAcceptedUsed
    integer                                    :: iGEMInvalidCompBoundRejectedUsed
    integer                                    :: iGEMInvalidCompBoundVerdict
    integer                                    :: iGEMInvalidCompBoundAttemptedTotal
    integer                                    :: iGEMInvalidCompBoundAcceptedTotal
    integer                                    :: iGEMInvalidCompBoundRejectedTotal
    integer                                    :: iGEMCoalesceUsed
    integer                                    :: iGEMDuplicateSUBOMCoalesceTotalUsed
    integer                                    :: iGEMDegenerateODCoalesceTotalUsed
    integer                                    :: iGEMStabilizeUsed
    integer                                    :: iGEMStabilizeTotalUsed
    integer                                    :: iGEMBoundaryRemovalUsed
    integer                                    :: iGEMRawNegativeRemovalUsed
    integer                                    :: iGEMBoundaryPinnedRemovalUsed
    integer                                    :: iGEMTinyBoundaryRemovalUsed
    integer                                    :: iGEMTraceRemoveUsed
    integer                                    :: iGEMTraceReinjectUsed
    integer                                    :: iGEMCEFRetryActivationUsed
    integer                                    :: iGEMCEFRetryActivationTotalUsed
    integer                                    :: iGEMSUBGQRideAlongUsed
    integer                                    :: iODPairCount
    integer                                    :: iODOrderingModeCount
    integer                                    :: nLevelCandidate, nLevelCandidateCapacity
    integer                                    :: nLevelCandidateRowOffset = 0
    integer                                    :: nGEMCEFPhaseVariables = 0
    integer                                    :: nGEMCEFSiteVariables = 0
    integer                                    :: iGEMPreLMNPrimal = 0
    integer                                    :: iGEMPreLMNConstraint = 0
    integer                                    :: iGEMRawNegativePhaseSlot = 0
    integer                                    :: iGEMRawNegativePhaseSoln = 0
    integer                                    :: iGEMRawNegativePhaseSpecies = 0
    integer                                    :: iGEMRawNegComp = 0
    integer                                    :: iGEMCEFBndPhaseSlot = 0
    integer                                    :: iGEMBoundaryRemovalSlot = 0
    integer                                    :: iGEMBoundaryRemovalPhase = 0
    integer                                    :: iGEMTinyBoundaryRemovalSlot = 0
    integer                                    :: iGEMTinyBoundaryRemovalPhase = 0
    integer                                    :: iGEMBoundaryRankGuardUsed = 0
    integer                                    :: iGEMBoundaryRankGuardSlot = 0
    integer                                    :: iGEMBoundaryRankGuardPhase = 0
    integer                                    :: iGEMTraceRemoveSpecies = 0
    integer                                    :: iGEMTraceRemovePhase = 0
    integer                                    :: iGEMTraceRemoveCount = 0
    integer                                    :: iGEMTraceReinjectSpecies = 0
    integer                                    :: iGEMTraceReinjectPhase = 0
    integer                                    :: iGEMTraceReinjectCount = 0
    integer                                    :: nLevel2LagrangeInput = 0
    integer                                    :: nLevel2LagrangePruned = 0
    integer                                    :: nLevel2LagrangeStoichSelected = 0
    integer                                    :: nLevel2LagrangeSolnSelected = 0
    integer                                    :: nLevel2LagrangeSolnMerged = 0
    integer                                    :: nLevel2LagrangeSolnAdded = 0
    integer                                    :: nLevel2LagrangeOrderProjected = 0
    integer                                    :: nLevel2LagrangeTwoSetCreated = 0
    integer                                    :: nPEALagrangianPolishAttempt = 0
    integer                                    :: nPEALagrangianPolishAccepted = 0
    integer                                    :: nPEALagrangianPolishRejected = 0
    integer                                    :: iPEALagrangianPolishIterGlobal = 0
    integer                                    :: iPEALagrangianPolishReason = PHASE_CHANGE_REASON_NONE
    integer                                    :: iPEALagrangianHandoffRepeat = 0
    integer                                    :: iPEALagrangianHandoffFirstIter = 0
    integer                                    :: nPEALagrangianHandoffRepeated = 0
    integer                                    :: iPEAGatePass = 0
    integer                                    :: iPEAGateBlock = 0
    integer                                    :: iPEAGatePotBlock = 0
    integer                                    :: iPEAGateDFBlock = 0
    integer                                    :: nPEAGatePassed = 0
    integer                                    :: nPEAGateBlocked = 0
    integer                                    :: nPEAGatePotBlocked = 0
    integer                                    :: nPEAGateDFBlocked = 0
    integer                                    :: iPEADirectHandoff = 0
    integer                                    :: iPEADirectPhase = 0
    integer                                    :: nPEADirectHandoff = 0
    integer                                    :: iPEARepeatExit = 0
    integer                                    :: iPEARepeatExitReason = 0
    integer                                    :: nPEARepeatExit = 0
    integer                                    :: nPEARepeatGuard = 0
    integer                                    :: nPEARecorded = 0
    integer                                    :: nSUBOMTwoSetStored = 0
    integer                                    :: nSUBOMTwoSetTrace = 0
    integer                                    :: nSUBOMTwoSetTraceCapacity = 0
    integer                                    :: nSUBOMOrderingGateEvaluated = 0
    integer                                    :: nGEMCertificationSweep = 0

    ! Real scalars:
    real(8)                                    :: dGEMFunctionNorm, dGEMFunctionNormLast
    real(8)                                    :: dGEMMassBalanceNorm, dGEMChemicalPotentialNorm
    real(8)                                    :: dGEMSolutionChemicalPotentialNorm, dGEMCondensedChemicalPotentialNorm
    real(8)                                    :: dMaxSpeciesChange, dMinGibbs
    real(8)                                    :: dSublatticeExchangeNorm
    real(8)                                    :: dGEMMaxSublExchangeResidual
    real(8)                                    :: dGEMMaxSublExchangeWeightedResidual
    real(8)                                    :: dGEMNewtonSymmetryResidual
    real(8)                                    :: dGEMNewtonMinPivotScale, dGEMNewtonMaxPivotScale
    real(8)                                    :: dGEMNewtonDirectionNorm
    real(8)                                    :: dMinDrivingForceStoich, dMinDrivingForceSoln
    real(8)                                    :: dSpeciesRemove, dPlateau
    real(8)                                    :: dTraceSpeciesReducedNormLast
    real(8)                                    :: dMinPhasePotential, dToleranceLevel, dMaxElementPotential
    real(8)                                    :: dPEATol, dMinMolesPhase, dMaxPotentialTol
    real(8)                                    :: dGEMLineSearchInitialNorm, dGEMLineSearchBestNorm
    real(8)                                    :: dGEMLineSearchFinalNorm
    real(8)                                    :: dGEMLineSearchInitialStep, dGEMLineSearchBestStep
    real(8)                                    :: dGEMLineSearchFinalStep
    real(8)                                    :: dGEMLineSearchMinRawPhaseMoles
    real(8)                                    :: dGEMLineSearchMinFinalPhaseMoles
    real(8)                                    :: dGEMLineSearchInitialGibbs
    real(8)                                    :: dGEMLineSearchBestGibbs
    real(8)                                    :: dGEMLineSearchFinalGibbs
    real(8)                                    :: dGEMLineSearchInitialMerit
    real(8)                                    :: dGEMLineSearchBestMerit
    real(8)                                    :: dGEMLineSearchFinalMerit
    real(8)                                    :: dGEMLineSearchMeritCandNorm
    real(8)                                    :: dGEMLineSearchMeritCandMass
    real(8)                                    :: dGEMLineSearchMeritCandStep
    real(8)                                    :: dGEMLineSearchMeritCandGibbs
    real(8)                                    :: dGEMLineSearchMeritCandMerit
    real(8)                                    :: dGEMPreLMInitNorm
    real(8)                                    :: dGEMPreLMBestNorm
    real(8)                                    :: dGEMPreLMInitGibbs
    real(8)                                    :: dGEMPreLMBestGibbs
    real(8)                                    :: dGEMPreLMInitMerit
    real(8)                                    :: dGEMPreLMBestMerit
    real(8)                                    :: dGEMPreLMMeritCandNorm
    real(8)                                    :: dGEMPreLMMeritCandMass
    real(8)                                    :: dGEMPreLMMeritCandStep
    real(8)                                    :: dGEMPreLMMeritCandGibbs
    real(8)                                    :: dGEMPreLMMeritCandMerit
    real(8)                                    :: dGEMPreLMFinalNorm
    real(8)                                    :: dGEMPreLMFinalGibbs
    real(8)                                    :: dGEMPreLMFinalMerit
    real(8)                                    :: dGEMPreLMDirNormSlope
    real(8)                                    :: dGEMPreLMDirGibbsSlope
    real(8)                                    :: dGEMPreLMDirMeritSlope
    real(8)                                    :: dGEMPreLMDirectionNorm
    real(8)                                    :: dGEMPreLMMinPivotScale
    real(8)                                    :: dGEMPreLMMaxPivotScale
    real(8)                                    :: dGEMPreLMODAlign
    real(8)                                    :: dGEMPreLMODEigen
    real(8)                                    :: dGEMNewtonDirNormSlope
    real(8)                                    :: dGEMNewtonDirGibbsSlope
    real(8)                                    :: dGEMNewtonDirMeritSlope
    real(8)                                    :: dGEMRawNegativePhaseAmount
    real(8)                                    :: dGEMRawNegativePhaseDirection
    real(8)                                    :: dGEMRawNegPhaseResidual
    real(8)                                    :: dGEMInvalidCompBoundPhi
    real(8)                                    :: dGEMCEFBndPhaseStep
    real(8)                                    :: dODCompDist
    real(8)                                    :: dODOrderNorm
    real(8)                                    :: dODOrderingEigenMin
    real(8)                                    :: dPEALagrangianPolishNormBefore = 0D0
    real(8)                                    :: dPEALagrangianPolishNormAfter = 0D0
    real(8)                                    :: dPEALagrangianPolishPotentialChange = 0D0
    real(8)                                    :: dGEMTimingLeveling = 0D0
    real(8)                                    :: dGEMTimingHandoff = 0D0
    real(8)                                    :: dGEMTimingPEA = 0D0
    real(8)                                    :: dGEMTimingLagrangian = 0D0
    real(8)                                    :: dGEMTimingCertification = 0D0

    ! Logical scalars:
    logical                                    :: lDebugMode, lRevertSystem, lConverged, lSubConverged
    logical                                    :: lNegativeMolesPhase, lGibbsMinCheck, lPhaseChange
    logical                                    :: lCompbdOnly, lPostProcess, lTraceSpeciesControlEnabled
    logical                                    :: lSampledLevelingThermoExtended = .FALSE.
    logical                                    :: lGEMCEFSiteLagrangianEnabled = .FALSE.
    logical                                    :: lGEMCEFSiteLagrangianActive = .FALSE.
    logical                                    :: lGEMCEFSiteDirectionActive = .FALSE.
    logical                                    :: lGEMCEFResidualLMDirection = .FALSE.
    logical                                    :: lGEMCEFInertiaRegularizationEnabled = .FALSE.
    logical                                    :: lGEMCEFInertiaRegularizationActive = .FALSE.
    logical                                    :: lGEMCEFBndPhaseActive = .FALSE.
    logical                                    :: lPEALagrangianPolishEnabled = .TRUE.
    logical                                    :: lPEALagrangianPolishActive = .FALSE.
    logical                                    :: lPEALagrangianPolishAccepted = .FALSE.
    logical                                    :: lSUBOMTwoSetCandidateEnabled = .FALSE.
    logical                                    :: lSUBQStandaloneEnabled = .TRUE.


    ! Integer 1D arrays:
    integer, dimension(:), allocatable         :: iAssemblageGEM, iAssemblageGEMinit
    integer, dimension(:), allocatable         :: iPhaseChangeReasonHistory
    integer, dimension(:), allocatable         :: iGEMNewtonSolverHist, iGEMNewtonInfoHist
    integer, dimension(:), allocatable         :: iGEMNewtonKKTSizeHist
    integer, dimension(:), allocatable         :: iGEMNewtonPivot1x1Hist, iGEMNewtonPivot2x2Hist
    integer, dimension(:), allocatable         :: iGEMNewtonPivotPositiveHist, iGEMNewtonPivotNegativeHist
    integer, dimension(:), allocatable         :: iGEMNewtonPivotZeroHist
    integer, dimension(:), allocatable         :: iGEMLSIterHist
    integer, dimension(:), allocatable         :: iGEMLSNegFactorHist
    integer, dimension(:), allocatable         :: iGEMLSNegPhaseHist
    integer, dimension(:), allocatable         :: iGEMLSFloorHist
    integer, dimension(:), allocatable         :: iGEMLSNoDescentHist
    integer, dimension(:), allocatable         :: iGEMLSNoDescentClassHist
    integer, dimension(:), allocatable         :: iGEMPreLMNoDescentHist
    integer, dimension(:), allocatable         :: iGEMPreLMNoDescentClassHist
    integer, dimension(:), allocatable         :: iGEMTraceInactiveHist
    integer, dimension(:), allocatable         :: iGEMTraceReinjectedHist
    integer, dimension(:), allocatable         :: iGEMTraceRemoveHist
    integer, dimension(:), allocatable         :: iGEMTraceReinjectHist
    integer, dimension(:), allocatable         :: iGEMResidualLMHist
    integer, dimension(:), allocatable         :: iGEMResidualLMRawNegativeHist
    integer, dimension(:), allocatable         :: iGEMResidualLMNoDescentHist
    integer, dimension(:), allocatable         :: iGEMResidualLMNoDescentClassHist
    integer, dimension(:), allocatable         :: iGEMResidualLMEventCallSite
    integer, dimension(:), allocatable         :: iGEMResidualLMEventIterPEA
    integer, dimension(:), allocatable         :: iGEMResidualLMEventIterGlobal
    integer, dimension(:), allocatable         :: iGEMResidualLMEventReason
    integer, dimension(:), allocatable         :: iGEMResidualLMEventNoDescentClass
    integer, dimension(:), allocatable         :: iGEMResidualLMEventNewtonInfo
    integer, dimension(:), allocatable         :: iGEMResidualLMEventKKTSize
    integer, dimension(:), allocatable         :: iGEMResidualLMEventNPrimal
    integer, dimension(:), allocatable         :: iGEMResidualLMEventNConstraint
    integer, dimension(:), allocatable         :: iGEMResidualLMEventPivot1x1
    integer, dimension(:), allocatable         :: iGEMResidualLMEventPivot2x2
    integer, dimension(:), allocatable         :: iGEMResidualLMEventPivotPositive
    integer, dimension(:), allocatable         :: iGEMResidualLMEventPivotNegative
    integer, dimension(:), allocatable         :: iGEMResidualLMEventPivotZero
    integer, dimension(:), allocatable         :: iGEMLMEventODOrdPhase
    integer, dimension(:), allocatable         :: iGEMLMEventODCompPhase
    integer, dimension(:), allocatable         :: iGEMCoalesceHist
    integer, dimension(:), allocatable         :: iGEMStabilizeHist
    integer, dimension(:), allocatable         :: iGEMBoundaryRemovalHist
    integer, dimension(:), allocatable         :: iGEMRawNegativeRemovalHist
    integer, dimension(:), allocatable         :: iGEMRawNegativePhaseSlotHist
    integer, dimension(:), allocatable         :: iGEMRawNegativePhaseSolnHist
    integer, dimension(:), allocatable         :: iGEMRawNegativePhaseSpeciesHist
    integer, dimension(:), allocatable         :: iGEMRawNegCompHist
    integer, dimension(:), allocatable         :: iGEMInvalidCompBoundAttemptedHist
    integer, dimension(:), allocatable         :: iGEMInvalidCompBoundAcceptedHist
    integer, dimension(:), allocatable         :: iGEMInvalidCompBoundRejectedHist
    integer, dimension(:), allocatable         :: iGEMInvalidCompBoundVerdictHist
    integer, dimension(:), allocatable         :: iGEMBoundaryPinnedRemovalHist
    integer, dimension(:), allocatable         :: iGEMBoundaryRemovalSlotHist
    integer, dimension(:), allocatable         :: iGEMBoundaryRemovalPhaseHist
    integer, dimension(:), allocatable         :: iGEMTinyBoundaryRemovalHist
    integer, dimension(:), allocatable         :: iGEMTinyBoundaryRemovalSlotHist
    integer, dimension(:), allocatable         :: iGEMTinyBoundaryRemovalPhaseHist
    integer, dimension(:), allocatable         :: iGEMBoundaryRankGuardHist
    integer, dimension(:), allocatable         :: iGEMBoundaryRankGuardSlotHist
    integer, dimension(:), allocatable         :: iGEMBoundaryRankGuardPhaseHist
    integer, dimension(:), allocatable         :: iGEMTraceRemoveSpeciesHist
    integer, dimension(:), allocatable         :: iGEMTraceRemovePhaseHist
    integer, dimension(:), allocatable         :: iGEMTraceRemoveCountHist
    integer, dimension(:), allocatable         :: iGEMTraceReinjectSpeciesHist
    integer, dimension(:), allocatable         :: iGEMTraceReinjectPhaseHist
    integer, dimension(:), allocatable         :: iGEMTraceReinjectCountHist
    integer, dimension(:), allocatable         :: iGEMCEFRetryActivationHist
    integer, dimension(:), allocatable         :: iGEMSUBGQRideAlongHist
    integer, dimension(:), allocatable         :: iODPairHist
    integer, dimension(:), allocatable         :: iODOrderingModeHist
    integer, dimension(:,:), allocatable       :: iGEMResidualLMEventAssemblage
    integer, dimension(:), allocatable         :: iPhaseGEM, iPhaseGEMinit
    integer, dimension(:), allocatable         :: iLevelCandidatePhase, iLevelCandidateSource
    integer, dimension(:), allocatable         :: iLevelCandidateParentPhase, iLevelCandidateDisplayPhase
    integer, dimension(:), allocatable         :: iLevelCandidateIdentityOrdinal
    integer, dimension(:), allocatable         :: iLevelCandidateFromLevel
    integer, dimension(:), allocatable         :: iSubMinCandidateStatusSoln
    integer, dimension(:), allocatable         :: iGEMCEFPhaseSlot, iGEMCEFPhaseSoln, iGEMCEFPhaseSpecies
    integer, dimension(:), allocatable         :: iGEMCEFVarPhaseVar, iGEMCEFVarSolnPhase, iGEMCEFVarPhaseID
    integer, dimension(:), allocatable         :: iGEMCEFVarSub, iGEMCEFVarCon, iGEMCEFVarRef
    integer, dimension(:), allocatable         :: iLevel2LagrangeInputAssemblage
    integer, dimension(:), allocatable         :: iLevel2LagrangeInputPhase
    integer, dimension(:), allocatable         :: iLevel2LagrangeInputCandidate
    integer, dimension(:), allocatable         :: iLevel2LagrangeInputSource
    integer, dimension(:), allocatable         :: iLevel2LagrangeOutputAssemblage
    integer, dimension(:), allocatable         :: iActiveSlotThermoPhase
    integer, dimension(:), allocatable         :: iActiveSlotDisplayPhase
    integer, dimension(:), allocatable         :: iActiveSlotIdentityOrdinal
    integer, dimension(:), allocatable         :: iPEALevelIterHist
    integer, dimension(:), allocatable         :: iPEAPolishAttemptHist
    integer, dimension(:), allocatable         :: iPEAPolishAcceptedHist
    integer, dimension(:), allocatable         :: iPEAPolishReasonHist
    integer, dimension(:), allocatable         :: iPEAPolishIterGlobalHist
    integer, dimension(:), allocatable         :: iPEALagrangianHandoffRepeatHist
    integer, dimension(:), allocatable         :: iPEALagrangianHandoffFirstIterHist
    integer, dimension(:), allocatable         :: iPEAGatePassHist
    integer, dimension(:), allocatable         :: iPEAGateBlockHist
    integer, dimension(:), allocatable         :: iPEAGatePotBlockHist
    integer, dimension(:), allocatable         :: iPEAGateDFBlockHist
    integer, dimension(:), allocatable         :: iPEADirectHandoffHist
    integer, dimension(:), allocatable         :: iPEADirectPhaseHist
    integer, dimension(:), allocatable         :: iPEARepeatExitHist
    integer, dimension(:), allocatable         :: iPEARepeatExitReasonHist
    integer, dimension(:), allocatable         :: iSUBOMTwoSetTraceStage
    integer, dimension(:), allocatable         :: iSUBOMTwoSetTraceIterPEA
    integer, dimension(:), allocatable         :: iSUBOMTwoSetTraceIterGlobal
    integer, dimension(:), allocatable         :: iSUBOMTwoSetTracePhase
    integer, dimension(:), allocatable         :: iSUBOMTwoSetTraceOrdinal
    integer, dimension(:), allocatable         :: iSUBOMTwoSetTraceSlot
    integer, dimension(:), allocatable         :: iSUBOMTwoSetStoredPhase
    integer, dimension(:), allocatable         :: iSUBOMTwoSetStoredOrdinal
    integer, dimension(:), allocatable         :: iSUBOMOrderingGateIterPEA
    integer, dimension(:), allocatable         :: iSUBOMOrderingGateModeCount
    integer, dimension(:), allocatable         :: iSUBOMOrderingGateInfo

    ! Real 1D arrays:
    real(8), dimension(:), allocatable         :: dSumMolFractionSoln, dMolesPhaseLast
    real(8), dimension(:), allocatable         :: dUpdateVar, dDrivingForceSoln
    real(8), dimension(:), allocatable         :: dUpdateVarLast
    real(8), dimension(:), allocatable         :: dMolesSpeciesLast, dElementPotentialLast
    real(8), dimension(:), allocatable         :: dPhasePotential, dChemicalPotentialGEM
    real(8), dimension(:), allocatable         :: dChemicalPotentialGEMinit
    real(8), dimension(:), allocatable         :: dGEMNormHist
    real(8), dimension(:), allocatable         :: dGEMMassNormHist
    real(8), dimension(:), allocatable         :: dGEMChemNormHist
    real(8), dimension(:), allocatable         :: dGEMSolnChemNormHist
    real(8), dimension(:), allocatable         :: dGEMCondChemNormHist
    real(8), dimension(:), allocatable         :: dGEMSublExchangeNormHist
    real(8), dimension(:), allocatable         :: dGEMNewtonSymResidualHist
    real(8), dimension(:), allocatable         :: dGEMNewtonMinPivotScaleHist
    real(8), dimension(:), allocatable         :: dGEMNewtonMaxPivotScaleHist
    real(8), dimension(:), allocatable         :: dGEMNewtonDirectionNormHist
    real(8), dimension(:), allocatable         :: dGEMLSInitialNormHist
    real(8), dimension(:), allocatable         :: dGEMLSBestNormHist
    real(8), dimension(:), allocatable         :: dGEMLSFinalNormHist
    real(8), dimension(:), allocatable         :: dGEMLSInitialStepHist
    real(8), dimension(:), allocatable         :: dGEMLSBestStepHist
    real(8), dimension(:), allocatable         :: dGEMLSFinalStepHist
    real(8), dimension(:), allocatable         :: dGEMLSMinRawPhaseMolesHist
    real(8), dimension(:), allocatable         :: dGEMLSMinFinalPhaseMolesHist
    real(8), dimension(:), allocatable         :: dGEMLSInitialGibbsHist
    real(8), dimension(:), allocatable         :: dGEMLSBestGibbsHist
    real(8), dimension(:), allocatable         :: dGEMLSFinalGibbsHist
    real(8), dimension(:), allocatable         :: dGEMLSInitialMeritHist
    real(8), dimension(:), allocatable         :: dGEMLSBestMeritHist
    real(8), dimension(:), allocatable         :: dGEMLSFinalMeritHist
    real(8), dimension(:), allocatable         :: dGEMPreLMInitNormHist
    real(8), dimension(:), allocatable         :: dGEMPreLMBestNormHist
    real(8), dimension(:), allocatable         :: dGEMPreLMInitGibbsHist
    real(8), dimension(:), allocatable         :: dGEMPreLMBestGibbsHist
    real(8), dimension(:), allocatable         :: dGEMPreLMInitMeritHist
    real(8), dimension(:), allocatable         :: dGEMPreLMBestMeritHist
    real(8), dimension(:), allocatable         :: dGEMPreLMMeritCandNormHist
    real(8), dimension(:), allocatable         :: dGEMPreLMMeritCandMassHist
    real(8), dimension(:), allocatable         :: dGEMPreLMMeritCandStepHist
    real(8), dimension(:), allocatable         :: dGEMPreLMMeritCandGibbsHist
    real(8), dimension(:), allocatable         :: dGEMPreLMMeritCandMeritHist
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventNormSlope
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventGibbsSlope
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventMeritSlope
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventInitialNorm
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventBestNorm
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventFinalNorm
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventInitialGibbs
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventBestGibbs
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventFinalGibbs
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventInitialMerit
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventBestMerit
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventFinalMerit
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventMeritCandNorm
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventMeritCandMass
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventMeritCandStep
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventMeritCandGibbs
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventMeritCandMerit
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventDirectionNorm
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventMinPivot
    real(8), dimension(:), allocatable         :: dGEMResidualLMEventMaxPivot
    real(8), dimension(:), allocatable         :: dGEMLMEventODAlign
    real(8), dimension(:), allocatable         :: dGEMLMEventODEigen
    real(8), dimension(:), allocatable         :: dGEMNewtonDirNormSlopeHist
    real(8), dimension(:), allocatable         :: dGEMNewtonDirGibbsSlopeHist
    real(8), dimension(:), allocatable         :: dGEMNewtonDirMeritSlopeHist
    real(8), dimension(:), allocatable         :: dGEMRawNegativePhaseAmountHist
    real(8), dimension(:), allocatable         :: dGEMRawNegativePhaseDirectionHist
    real(8), dimension(:), allocatable         :: dGEMRawNegPhaseResidualHist
    real(8), dimension(:), allocatable         :: dGEMInvalidCompBoundPhiHist
    real(8), dimension(:), allocatable         :: dODCompDistHist
    real(8), dimension(:), allocatable         :: dODOrderNormHist
    real(8), dimension(:), allocatable         :: dODOrderingEigenMinHist
    real(8), dimension(:), allocatable         :: dGEMLSRawPhaseMoles
    real(8), dimension(:), allocatable         :: dGEMLSFinalPhaseMoles
    real(8), dimension(:), allocatable         :: dGEMCEFPhaseDirection, dGEMCEFSiteDirection
    real(8), dimension(:), allocatable         :: dGEMCEFElementDirection
    real(8), dimension(:), allocatable         :: dGEMCEFPhaseResidual
    real(8), dimension(:), allocatable         :: dGEMAnalyticalSpeciesDirection
    real(8), dimension(:), allocatable         :: dLevel2LagrangeInputMoles
    real(8), dimension(:), allocatable         :: dLevel2LagrangeOutputMoles
    real(8), dimension(:), allocatable         :: dLevel2LagrangeElementPotentialIn
    real(8), dimension(:), allocatable         :: dLevel2LagrangeElementPotentialOut
    real(8), dimension(:), allocatable         :: dPEAMinPhasePotentialHist
    real(8), dimension(:), allocatable         :: dPEAMaxElementPotentialHist
    real(8), dimension(:), allocatable         :: dPEAGEMFunctionNormHist
    real(8), dimension(:), allocatable         :: dPEAPolishNormHist
    real(8), dimension(:), allocatable         :: dPEAPolishPotentialChangeHist
    real(8), dimension(:), allocatable         :: dSUBOMTwoSetTraceAmount
    real(8), dimension(:), allocatable         :: dSUBOMOrderingGateEigenMin

    ! Logical 1D arrays:
    logical, dimension(:), allocatable         :: lSolnPhases, lMiscibility, lPhaseChangeHistory
    logical, dimension(:), allocatable         :: lTraceSpeciesInactive, lTraceSpeciesReinjected
    logical, dimension(:), allocatable         :: lGEMAnalyticalSpeciesDirection
    logical, dimension(:), allocatable         :: lSUBOMOrderingGateUnstable


    ! Integer 2D arrays:
    integer, dimension(:,:), allocatable       :: iterHistory
    integer, dimension(:,:), allocatable       :: iterHistoryLevel
    integer, dimension(:,:), allocatable       :: iPEAAssemblageHist
    integer, dimension(:,:), allocatable       :: iPEALagrangianHandoffHist

    ! Real 2D arrays:
    real(8), dimension(:,:), allocatable       :: dStoichSpeciesGEM, dAtomFractionSpeciesGEM
    real(8), dimension(:,:), allocatable       :: dMolesPhaseHistory, dEffStoichSolnPhase
    real(8), dimension(:,:), allocatable       :: dMolFractionGEM
    real(8), dimension(:,:), allocatable       :: dActiveSlotMolFraction
    real(8), dimension(:,:,:), allocatable     :: dActiveSlotSiteFraction
    real(8), dimension(:,:), allocatable       :: dStoichSpeciesGEMinit, dAtomFractionSpeciesGEMinit
    real(8), dimension(:,:), allocatable       :: dMolFractionGEMinit
    real(8), dimension(:,:), allocatable       :: dLevelCandidateMolFraction
    real(8), dimension(:,:), allocatable       :: dSUBOMTwoSetStoredMol
    real(8), dimension(:,:), allocatable       :: dGEMLSRawPhaseMolesHist
    real(8), dimension(:,:), allocatable       :: dGEMLSFinalPhaseMolesHist
    real(8), dimension(:,:), allocatable       :: dPEAElementPotentialHist
    real(8), dimension(:,:), allocatable       :: dSUBOMTwoSetTraceMol

    ! Real 3D arrays:
    real(8), dimension(:,:,:), allocatable     :: dGEMCEFSiteLast
    real(8), dimension(:,:,:), allocatable     :: dGEMCEFPhaseSiteLast
    real(8), dimension(:,:,:), allocatable     :: dSUBOMTwoSetTraceSite

end module ModuleGEMSolver
