!-------------------------------------------------------------------------------------------------------------
!
!> \file    ModuleSubMin.f90
!> \brief   Shared state for solution-phase Subminimization.
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   08/21/2012      M.H.A. Piro         Original subminimization state module
    !   07/20/2026      S.Y. Kwon           Added subminimization state and diagnostics for CEF site fractions, trace estimates, and analytical solution curvature.
    !
    !
! Purpose:
! ========
!
!> \details This module stores process-global state for the solution-phase
!! Subminimization stage used by CheckPhaseAssemblage/PEA.  Subminimization
!! minimizes the phase driving force for one candidate solution phase at the
!! current elemental potentials.  The active path is Subminimization ->
!! SubMinInit, then either SubMinSiteFractionCEF for CEF SUBL/SUBLM/SUBOM
!! phases or SubMinNewton -> SubMinLineSearch -> SubMinFunctionNorm for the
!! non-CEF path, with optional phase-local trace-species removal/reinjection.
!! SubMinNewton may assemble either the classic ideal-mixing endmember Hessian
!! or a solution-model analytical Hessian when one exists.  CEF phases never
!! enter the classic path; a feasible negative-driving-force CEF point is
!! treated as a phase-entry witness.
!
!
! Required setup variables:
! =========================
!
! nVar                 Number of endmembers/species in the phase being subminimized.
! iFirstSUB            First absolute species index for the active solution phase.
! iLastSUB             Last absolute species index for the active solution phase.
! iterSub              Current Subminimization iteration.  Used by Newton damping and trace reinjection.
! iSolnPhaseIndexOther Matching miscibility-gap phase index used by SubMinCheckDuplicate.
!
!
! Output/updated variables:
! =========================
!
! dDrivingForce        Current subminimized driving force for the active phase.
! dDrivingForceLast    Previous driving force used by SubMinLineSearch acceptance.
! dSubMinFunctionNorm  Mole-fraction, mass-balance, and driving-force residual norm.
! dSubminGibbsEst      Mole-fraction-weighted residual estimate used for potential-vector updates.
! dMaxPotentialVector  Maximum active endmember chemical-potential residual.
! lSubMinConverged     True when the active subminimization solve satisfies the current criteria.
! lSubMinMaxHit        True after any Subminimization call reaches the max-iteration guard.
! iSubMinMaxPhaseIndex Solution-phase index from the latest max-iteration event.
! iSubMinMaxIter       Iteration counter from the latest max-iteration event.
! dSubMinMaxFunctionNorm
!                      Function norm from the latest max-iteration event.
! dSubMinMaxPotentialVector
!                      Potential-vector residual from the latest max-iteration event.
! dSubMinMaxMinFraction
!                      Minimum endmember fraction from the latest max-iteration event.
! dSubMinMaxMaxFraction
!                      Maximum endmember fraction from the latest max-iteration event.
! dSubMinMaxElementPotential
!                      Element-potential vector captured at the latest max-iteration event.
! dSubMinInitialMolFraction
!                      Normalized starting endmember fractions for the active submin call.
! dSubMinMaxInitialMolFraction
!                      Starting endmember fractions from the latest max-iteration event.
! dSubMinMaxFinalMolFraction
!                      Final endmember fractions from the latest max-iteration event.
! lSubMinCEFAttempted True when the current submin call entered the CEF site-fraction path.
! lSubMinCEFHandled   True when the current CEF path converged tightly enough to pass the result.
! iSubMinCEFIter      Final CEF site-fraction iteration for the current submin call.
! dSubMinCEFGradientNorm
!                      Final active site-gradient norm from the current CEF path.
! lSubMinNewtonAnalyticalHessianAttempted
!                      True after this submin call enters a solution-model analytical Hessian branch.
! lSubMinNewtonAnalyticalHessianAccepted
!                      True when the latest Newton direction keeps the analytical Hessian branch after guards.
!
!
! Working arrays:
! ===============
!
! dChemicalPotentialStar  Element-potential projection for each active phase endmember.
! dRHS                    Newton right-hand side on entry; Newton update direction on return.
! dHessian                Symmetric Newton/KKT system or projected analytical Hessian.
! iHessian                DSYSV pivot/work index array for the symmetric indefinite solve.
! dPotentialVector        Endmember chemical-potential residual vector.
!
!
! Trace-species state:
! ====================
!
! lSubMinTraceInactive          Marks endmembers temporarily removed from the active Newton solve.
! lSubMinTraceReinjected        Prevents immediate remove/reinject cycling for a trace endmember.
! dSubMinTraceReducedNormLast   Previous reduced active-set potential norm for slow-progress checks.
! iSubMinTraceSlowProgressCount Number of consecutive low-progress reduced solves.
! iterSubMinTraceLastRemoval    Iteration of the last trace endmember removal.
! iterSubMinTraceLastReinject   Iteration of the last trace endmember reinjection.
! iSubMinHenrianTraceEstimateCount
!                                Number of accepted Henrian trace estimates in the current submin call.
!
!
! Newton diagnostics:
! ===================
!
! dSubMinNewtonSymmetryResidual Largest symmetry mismatch detected before DSYSV.
! iSubMinNewtonDSYSVInfo        LAPACK DSYSV status from the latest Newton solve.
! lSubMinNewtonAnalyticalHessianAttempted
!                                Analytical solution-Hessian branch was attempted during this submin call.
! lSubMinNewtonAnalyticalHessianAccepted
!                                Latest Newton direction accepted the analytical Hessian after descent checks.
!
!
! Primary callers/users:
! ======================
!
! Subminimization              Owns the phase-level subminimization workflow.
! SubMinInit                   Allocates and initializes the module state for one phase.
! SubMinSiteFractionCEF        Handles CEF SUBL/SUBLM/SUBOM phases in site-fraction variables.
! SubMinNewton                 Builds and solves the classic or analytical-Hessian Newton system.
! SubMinLineSearch             Accepts, damps, or multiplicatively repairs classic Newton updates.
! SubMinFunctionNorm           Updates the subminimization residual norm.
! SubMinTraceSpeciesControl    Removes/reinjects trace endmembers during slow active-set progress.
! SubMinCheckDuplicate         Uses miscibility-gap partner state to reject duplicate minima.
! ResetThermo                  Deallocates remaining Subminimization allocatable arrays.
!
!
! Numerical assumptions:
! ======================
!
! - Active endmember mole fractions are normalized to unity.  Trace-inactive
!   endmembers are set to zero and skipped by Newton, line-search, norm, and
!   driving-force updates.
! - CEF SUBL/SUBLM/SUBOM phases use site-fraction variables.  Positive
!   unconverged CEF candidates are rejected without classic fallback.  A
!   feasible negative-driving-force CEF point is a phase-entry witness even
!   when the site-gradient has not fully converged.  Magnetic RKMPM/QKTOM and
!   unsupported ordinary non-CEF solution models use the classic endmember path
!   directly.
! - SUBG/SUBQ/RKMP/QKTO phases use the analytical solution-Hessian branch
!   inside SubMinNewton when all variables are active and no ionic charge row is
!   required.  Generic trace controls remain disabled for SUBG/SUBQ because
!   their variables are MQM pair/quadruplet populations.  A pair-aware trace
!   update must include the coordination-number factors in the pair
!   chemical-potential relation.
! - dMinMoleFraction is 1D-30, a physical trace floor used to avoid exactly
!   zero mole fractions during active subminimization updates.
! - dSubMinHenrianTraceThreshold is the local trace-region limit used by the
!   classic Henrian estimate.  It should stay far above dMinMoleFraction and
!   below ordinary constituent fractions.
! - The Newton matrix is expected to be symmetric.  Ionic phases add a
!   charge-neutrality equation/column before the DSYSV solve.
! - dMaxPotentialVector is the primary chemical-potential stationarity
!   residual for active endmembers; dSubMinFunctionNorm remains the
!   thermodynamic/mass-balance residual guard.
! - Active non-CEF subminimization uses Newton plus line search.
!
!-------------------------------------------------------------------------------------------------------------
module ModuleSubMin
    implicit none

    SAVE

    ! Real parameters:
    real(8), parameter                  :: dSubMinTolerance = 1D-8, dMinMoleFraction = 1D-30
    real(8), parameter                  :: dSubMinHenrianTraceThreshold = 1D-5
    real(8), parameter                  :: dTolEuclideanNorm = 1D-4


    ! Integer scalars:
    integer                             :: nVar, iFirstSUB, iLastSUB, iSolnPhaseIndexOther, iterSub
    integer                             :: iSubMinMaxPhaseIndex = 0, iSubMinMaxIter = 0
    integer                             :: iSubMinTraceSlowProgressCount
    integer                             :: iSubMinHenrianTraceEstimateCount
    integer                             :: iSubMinNewtonDSYSVInfo
    integer                             :: iterSubMinTraceLastRemoval, iterSubMinTraceLastReinject
    integer                             :: iSubMinCEFPhaseIndex = 0, iSubMinCEFIter = 0
    integer                             :: nSubMinCEFIndependent = 0, nSubMinCEFActiveIndependent = 0

    ! Real scalars:
    real(8)                             :: dDrivingForce, dDrivingForceLast, dSubMinFunctionNorm
    real(8)                             :: dSubminGibbsEst, dMaxPotentialVector
    real(8)                             :: dSubMinMaxFunctionNorm = 0D0
    real(8)                             :: dSubMinMaxPotentialVector = 0D0
    real(8)                             :: dSubMinMaxMinFraction = 0D0
    real(8)                             :: dSubMinMaxMaxFraction = 0D0
    real(8)                             :: dSubMinTraceReducedNormLast
    real(8)                             :: dSubMinNewtonSymmetryResidual
    real(8)                             :: dSubMinCEFObjective = 0D0
    real(8)                             :: dSubMinCEFGradientNorm = 0D0
    real(8)                             :: dSubMinCEFStepNorm = 0D0
    real(8)                             :: dSubMinCEFChargeResidual = 0D0

    ! Logical scalars:
    logical                             :: lSubMinConverged
    logical                             :: lSubMinMaxHit = .FALSE.
    logical                             :: lSubMinCEFAttempted = .FALSE.
    logical                             :: lSubMinCEFHandled = .FALSE.
    logical                             :: lSubMinCEFLineSearchAccepted = .FALSE.
    logical                             :: lSubMinNewtonAnalyticalHessianAttempted = .FALSE.
    logical                             :: lSubMinNewtonAnalyticalHessianAccepted = .FALSE.


    ! Integer 1D arrays:
    integer, dimension(:), allocatable  :: iHessian

    ! Real 1D arrays:
    real(8), dimension(:), allocatable  :: dChemicalPotentialStar, dRHS
    real(8), dimension(:), allocatable  :: dPotentialVector
    real(8), dimension(:), allocatable  :: dSubMinInitialMolFraction
    real(8), dimension(:), allocatable  :: dSubMinMaxElementPotential
    real(8), dimension(:), allocatable  :: dSubMinMaxInitialMolFraction, dSubMinMaxFinalMolFraction

    ! Logical 1D arrays:
    logical, dimension(:), allocatable  :: lSubMinTraceInactive, lSubMinTraceReinjected


    ! Real 2D arrays:
    real(8), dimension(:,:), allocatable :: dHessian

end module ModuleSubMin
