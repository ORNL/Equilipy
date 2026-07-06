!> \brief Polish PEA elemental potentials with an embedded Lagrangian solve.
!!
!! \details Runs Level2Lagrange and RunLagrangianGEM on the current PEA
!! Leveling active set, then restores PEA-owned arrays and keeps only the
!! polished elemental-potential plane when the Lagrangian solve gives a clean
!! converged result.
!!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    PEALagrangianPolish.f90
!> \brief   Polish PEA elemental potentials without constraining the Lagrangian plane.
!> \author  S.Y. Kwon
!> \date    Jun. 26, 2026
!> \sa      CheckPhaseAssemblage.f90
!> \sa      Level2Lagrange.f90
!> \sa      RunLagrangianGEM.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   06/26/2026      S.Y. Kwon           Original PEA-internal Lagrangian polish routine.
!   06/26/2026      S.Y. Kwon           Rejected PEA-polish attempts on Lagrangian active-set
!                                        invalidation instead of accepting reduced active sets.
!   06/26/2026      S.Y. Kwon           Reset per-iteration polish state before guard returns.
!   06/27/2026      S.Y. Kwon           Skipped embedded PEA polish for CEF active sets until the
!                                        CEF polish path preserves metallic SUBL/SUBOM regressions.
!   06/27/2026      S.Y. Kwon           Enabled embedded PEA polish for CEF active sets through the
!                                        site-fraction Lagrangian path.
!   06/27/2026      S.Y. Kwon           Matched PEA-polish Lagrangian solve rules to the final Lagrangian
!                                        path while preserving PEA-owned trace masks on restore.
!   06/27/2026      S.Y. Kwon           Recorded translated PEA-to-Lagrangian handoff assemblages and
!                                        flagged repeated handoffs for oscillation diagnostics.
!   07/02/2026      S.Y. Kwon           Snapshotted active-slot CEF composition-set state so rejected
!                                        embedded polish attempts leave no slot-local state behind.
!   07/03/2026      S.Y. Kwon           Persisted accepted SUBOM two-set slot constitutions as candidate-pool
!                                        evidence before restoring the PEA state.
!
!
! Purpose:
! ========
!
!> \details PEA is responsible for active-set discovery, but its Leveling plane
!! can be less precise than a Lagrangian solve.  This routine
!! temporarily translates the current PEA Leveling assemblage to Lagrangian
!! variables, runs a fixed-active-set Lagrangian polish, and feeds back only
!! the refined elemental potentials.  The polished plane is not bounded by the
!! PEA Leveling plane, but PEA remains the active-set owner.  If the embedded
!! Lagrangian solve detects a phase-boundary, Newton, stagnation, or
!! convergence failure, the polish is rejected and PEA continues from its own
!! candidate rows.  Ordinary Lagrangian-refined constitutions are not stored as
!! persistent PEA candidates; accepted SUBOM two-set constitutions may be stored
!! as switch-gated candidate-pool evidence.
!
!
! Required input variables:
! =========================
!
! iAssemblage                Current PEA Leveling candidate assemblage.
! dMolesPhase                PEA phase amounts for the candidate assemblage.
! dElementPotential          Current PEA elemental-potential plane.
! dStoichSpeciesGEM          PEA candidate stoichiometry rows.
! dChemicalPotentialGEM      PEA candidate Gibbs rows.
! dMolFractionGEM            PEA solution candidate constitutions.
!
!
! Output/updated variables:
! =========================
!
! dElementPotential          Replaced by the Lagrangian-polished plane only if the polish is accepted.
! dPhasePotential            Recomputed from the restored PEA rows and current elemental potentials.
! nPEALagrangianPolish*      Attempt/accept/reject counters and iteration diagnostics.
! iPEALagrangianHandoff*     Latest and per-iteration repeated-handoff diagnostics.
! dPEALagrangianPolish*      Norm and elemental-potential change diagnostics.
!
!
! Called subroutines/functions:
! =============================
!
! Level2Lagrange             Converts the selected PEA rows to Lagrangian active-set variables.
! RunLagrangianGEM           Polishes the PEA active-set estimate in PEA-polish mode.
!
!
! Primary callers:
! ================
!
! CheckPhaseAssemblage       Calls after each PEA Leveling pass and before CompMinSolnPoint.
!
!
! Numerical assumptions:
! ======================
!
! - PEA remains the owner of candidate rows, candidate pool, and active-set discovery.
! - A polish is attempted only for a full non-empty PEA assemblage with nonnegative phase amounts.
! - The polished plane is accepted when the Lagrangian solve converges without
!   requesting active-set repair.  The plane is allowed to move far from the
!   PEA Leveling plane, because the Lagrangian solve is the more accurate
!   cotangent-plane estimate.
! - RunLagrangianGEM uses the same coalescing, trace handling, and CEF
!   boundary-removal rules here as in the final Lagrangian solve.  This routine
!   still restores PEA-owned arrays afterward and persists accepted elemental
!   potentials plus switch-gated two-set candidate-pool evidence back to PEA.
! - A phase-potential, Newton-failure, stagnation, or unconverged trigger is
!   rejected because it means Lagrangian did not produce a complete equilibrium
!   plane for the tested active set.
! - A non-CEF trace-reduced polish is accepted only after trace species are
!   fully reinjected; reduced trace masks are local to the polish attempt.
!
!-------------------------------------------------------------------------------------------------------------



subroutine PEALagrangianPolish

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer  :: nSpeciesLevelSave, nConPhasesSave, nSolnPhasesSave
    integer  :: iterGlobalSave, iPhaseChangeReasonSave
    integer  :: iTraceSlowSave, iStagnationSave
    integer  :: iterTraceRemovalSave, iterTraceReinjectSave
    integer  :: nGEMCEFPhaseVariablesSave, nGEMCEFSiteVariablesSave
    integer  :: iPolishReason, iPolishIter
    logical  :: lConvergedSave, lPhaseChangeSave, lSkipLagrangeSave
    logical  :: lSampledThermoSave, lGEMCEFActiveSave, lGEMCEFDirectionSave
    logical  :: lPolishConverged, lPolishAccepted
    logical  :: lHadAssemblageGEM, lHadMolesPhaseLast, lHadMolFractionOld
    logical  :: lHadTraceInactive, lHadTraceReinjected, lTraceInactiveAfter
    logical  :: lHadActiveSlotThermo, lHadActiveSlotDisplay, lHadActiveSlotIdentity
    logical  :: lHadActiveSlotMol, lHadActiveSlotSite
    logical  :: lHadGEMCEFSiteLast, lHadGEMCEFPhaseSiteLast
    real(8)  :: dNormBeforeSave, dNormLastSave, dTraceReducedNormSave
    real(8)  :: dMassNormSave, dChemNormSave, dSolnChemNormSave, dCondChemNormSave
    real(8)  :: dSublatticeNormSave, dPolishNorm
    real(8), dimension(nElements) :: dElementPotentialBefore, dElementPotentialAfter

    integer, dimension(:), allocatable   :: iAssemblageSave, iPhaseGEMSave
    integer, dimension(:), allocatable   :: iPhaseLevelSave, iAssemblageGEMSave
    integer, dimension(:), allocatable   :: iActiveSlotThermoSave, iActiveSlotDisplaySave
    integer, dimension(:), allocatable   :: iActiveSlotIdentitySave
    real(8), dimension(:), allocatable   :: dMolesPhaseSave, dMolesSpeciesSave
    real(8), dimension(:), allocatable   :: dMolFractionSave, dElementPotentialSave
    real(8), dimension(:), allocatable   :: dElementPotentialLastSave, dPhasePotentialSave
    real(8), dimension(:), allocatable   :: dChemicalPotentialSave, dLevelingChemicalPotentialSave
    real(8), dimension(:), allocatable   :: dMolesPhaseLastSave, dMolFractionOldSave
    real(8), dimension(:), allocatable   :: dUpdateVarSave, dUpdateVarLastSave
    real(8), dimension(:), allocatable   :: dChemicalPotentialGEMSave
    real(8), dimension(:), allocatable   :: dGibbsSolnPhaseSave, dSumMolFractionSolnSave
    real(8), dimension(:), allocatable   :: dDrivingForceSolnSave
    real(8), dimension(:), allocatable   :: dPartialExcessGibbsSave
    real(8), dimension(:), allocatable   :: dPartialEnthalpyXSSave
    real(8), dimension(:), allocatable   :: dPartialEntropyXSSave
    real(8), dimension(:), allocatable   :: dPartialHeatCapacityXSSave
    real(8), dimension(:,:), allocatable :: dAtomFractionSpeciesSave
    real(8), dimension(:,:), allocatable :: dLevelingCompositionSpeciesSave
    real(8), dimension(:,:), allocatable :: dStoichSpeciesLevelSave
    real(8), dimension(:,:), allocatable :: dStoichSpeciesGEMSave
    real(8), dimension(:,:), allocatable :: dAtomFractionSpeciesGEMSave
    real(8), dimension(:,:), allocatable :: dMolFractionGEMSave
    real(8), dimension(:,:), allocatable :: dActiveSlotMolSave
    real(8), dimension(:,:), allocatable :: dMolesPhaseHistorySave
    real(8), dimension(:,:), allocatable :: dEffStoichSolnPhaseSave
    real(8), dimension(:,:,:), allocatable :: dActiveSlotSiteSave
    real(8), dimension(:,:,:), allocatable :: dGEMCEFSiteLastSave
    real(8), dimension(:,:,:), allocatable :: dGEMCEFPhaseSiteLastSave
    integer, dimension(:,:), allocatable :: iterHistoryLevelSave
    logical, dimension(:), allocatable   :: lSolnPhasesSave
    logical, dimension(:), allocatable   :: lTraceSpeciesInactiveSave
    logical, dimension(:), allocatable   :: lTraceSpeciesReinjectedSave

    lPEALagrangianPolishAccepted = .FALSE.
    iPEALagrangianPolishReason = PHASE_CHANGE_REASON_NONE
    iPEALagrangianPolishIterGlobal = 0
    dPEALagrangianPolishNormBefore = 0D0
    dPEALagrangianPolishNormAfter = 0D0
    dPEALagrangianPolishPotentialChange = 0D0
    iPEALagrangianHandoffRepeat = 0
    iPEALagrangianHandoffFirstIter = 0

    if (.NOT.lPEALagrangianPolishEnabled) return
    if (nSolnPhasesSys <= 0) return
    if (.NOT.allocated(iAssemblage)) return
    if (.NOT.allocated(dMolesPhase)) return
    if (ANY(iAssemblage == 0)) return
    if (MINVAL(dMolesPhase) < -1D-10) return
    if (INFOThermo /= 0) return

    call SnapshotPEAState

    nPEALagrangianPolishAttempt = nPEALagrangianPolishAttempt + 1
    lPEALagrangianPolishAccepted = .FALSE.
    dPEALagrangianPolishNormBefore = dGEMFunctionNorm
    dElementPotentialBefore = dElementPotential

    lPEALagrangianPolishActive = .TRUE.
    call Level2Lagrange
    call RecordPEALagrangianHandoff
    call RunLagrangianGEM
    lPEALagrangianPolishActive = .FALSE.

    dElementPotentialAfter = dElementPotential
    dPolishNorm = dGEMFunctionNorm
    lPolishConverged = lConverged
    iPolishReason = iPhaseChangeReason
    iPolishIter = iterGlobal
    lTraceInactiveAfter = .FALSE.
    if (allocated(lTraceSpeciesInactive)) lTraceInactiveAfter = ANY(lTraceSpeciesInactive)
    lPolishAccepted = lFiniteVector(dElementPotentialAfter).AND.lPolishConverged.AND.&
        (iPolishReason == PHASE_CHANGE_REASON_NONE).AND.(.NOT.lTraceInactiveAfter)

    if (lPolishAccepted) call StoreAcceptedSUBOMTwoSetSlots
    call RestorePEAState

    iPEALagrangianPolishIterGlobal = iPolishIter
    iPEALagrangianPolishReason = iPolishReason
    dPEALagrangianPolishNormAfter = dPolishNorm
    dPEALagrangianPolishPotentialChange = MAXVAL(DABS(dElementPotentialAfter - dElementPotentialBefore))

    if (lPolishAccepted) then
        dElementPotential = dElementPotentialAfter
        lPEALagrangianPolishAccepted = .TRUE.
        nPEALagrangianPolishAccepted = nPEALagrangianPolishAccepted + 1
    else
        dElementPotential = dElementPotentialBefore
        lPEALagrangianPolishAccepted = .FALSE.
        nPEALagrangianPolishRejected = nPEALagrangianPolishRejected + 1
    end if

    dPhasePotential = dLevelingChemicalPotential - MATMUL(dLevelingCompositionSpecies,dElementPotential)

    return

contains

    subroutine RecordPEALagrangianHandoff
        integer :: iIter, iPrevious

        iPEALagrangianHandoffRepeat = 0
        iPEALagrangianHandoffFirstIter = 0

        iIter = iterPEA
        if ((iIter < 1).OR.(iIter > iterPEAMax)) return
        if (.NOT.allocated(iLevel2LagrangeOutputAssemblage)) return
        if (.NOT.allocated(iPEALagrangianHandoffHist)) return

        iPEALagrangianHandoffHist(:,iIter) = iLevel2LagrangeOutputAssemblage

        do iPrevious = 1, iIter - 1
            if (ALL(iPEALagrangianHandoffHist(:,iPrevious) == iLevel2LagrangeOutputAssemblage)) then
                iPEALagrangianHandoffRepeat = 1
                iPEALagrangianHandoffFirstIter = iPrevious
                nPEALagrangianHandoffRepeated = nPEALagrangianHandoffRepeated + 1
                exit
            end if
        end do

        if (allocated(iPEALagrangianHandoffRepeatHist)) then
            iPEALagrangianHandoffRepeatHist(iIter) = iPEALagrangianHandoffRepeat
        end if
        if (allocated(iPEALagrangianHandoffFirstIterHist)) then
            iPEALagrangianHandoffFirstIterHist(iIter) = iPEALagrangianHandoffFirstIter
        end if

        return
    end subroutine RecordPEALagrangianHandoff


    subroutine StoreAcceptedSUBOMTwoSetSlots
        integer :: iSlotLocal, iSolnPhaseLocal, iFirstLocal, iLastLocal
        real(8) :: dNormLocal

        if (.NOT.lSUBOMTwoSetCandidateEnabled) return
        if (.NOT.allocated(iActiveSlotThermoPhase)) return
        if (.NOT.allocated(iActiveSlotIdentityOrdinal)) return
        if (.NOT.allocated(dActiveSlotMolFraction)) return

        if (.NOT.allocated(iSUBOMTwoSetStoredPhase)) then
            allocate(iSUBOMTwoSetStoredPhase(MAX(1,nSolnPhasesSys)))
        end if
        if (.NOT.allocated(iSUBOMTwoSetStoredOrdinal)) then
            allocate(iSUBOMTwoSetStoredOrdinal(MAX(1,nSolnPhasesSys)))
        end if
        if (.NOT.allocated(dSUBOMTwoSetStoredMol)) then
            allocate(dSUBOMTwoSetStoredMol(MAX(1,nSolnPhasesSys),nSpecies))
        end if

        nSUBOMTwoSetStored = 0
        iSUBOMTwoSetStoredPhase = 0
        iSUBOMTwoSetStoredOrdinal = 0
        dSUBOMTwoSetStoredMol = 0D0

        do iSlotLocal = 1, MIN(nElements, SIZE(iActiveSlotThermoPhase))
            iSolnPhaseLocal = iActiveSlotThermoPhase(iSlotLocal)
            if ((iSolnPhaseLocal < 1).OR.(iSolnPhaseLocal > nSolnPhasesSys)) cycle
            if (TRIM(cSolnPhaseType(iSolnPhaseLocal)) /= 'SUBOM') cycle
            if (iActiveSlotIdentityOrdinal(iSlotLocal) <= 1) cycle
            if (iSolnPhaseLocal > SIZE(iSUBOMTwoSetStoredPhase)) cycle

            iFirstLocal = nSpeciesPhase(iSolnPhaseLocal-1) + 1
            iLastLocal  = nSpeciesPhase(iSolnPhaseLocal)
            if ((iFirstLocal < 1).OR.(iLastLocal > SIZE(dActiveSlotMolFraction,2))) cycle
            if (iLastLocal > SIZE(dSUBOMTwoSetStoredMol,2)) cycle

            dNormLocal = SUM(DMAX1(dActiveSlotMolFraction(iSlotLocal,iFirstLocal:iLastLocal), 0D0))
            if (dNormLocal <= 1D-300) cycle

            iSUBOMTwoSetStoredPhase(iSolnPhaseLocal) = iSolnPhaseLocal
            iSUBOMTwoSetStoredOrdinal(iSolnPhaseLocal) = iActiveSlotIdentityOrdinal(iSlotLocal)
            dSUBOMTwoSetStoredMol(iSolnPhaseLocal,iFirstLocal:iLastLocal) = &
                DMAX1(dActiveSlotMolFraction(iSlotLocal,iFirstLocal:iLastLocal), 0D0) / dNormLocal
            nSUBOMTwoSetStored = nSUBOMTwoSetStored + 1
        end do

        return
    end subroutine StoreAcceptedSUBOMTwoSetSlots

    subroutine SnapshotPEAState

        nSpeciesLevelSave = nSpeciesLevel
        nConPhasesSave = nConPhases
        nSolnPhasesSave = nSolnPhases
        iterGlobalSave = iterGlobal
        iPhaseChangeReasonSave = iPhaseChangeReason
        iTraceSlowSave = iTraceSpeciesSlowProgressCount
        iStagnationSave = iGEMStagnationCount
        iterTraceRemovalSave = iterTraceSpeciesLastRemoval
        iterTraceReinjectSave = iterTraceSpeciesLastReinject
        nGEMCEFPhaseVariablesSave = nGEMCEFPhaseVariables
        nGEMCEFSiteVariablesSave = nGEMCEFSiteVariables
        lConvergedSave = lConverged
        lPhaseChangeSave = lPhaseChange
        lSkipLagrangeSave = lSkipLagrange
        lSampledThermoSave = lSampledLevelingThermoExtended
        lGEMCEFActiveSave = lGEMCEFSiteLagrangianActive
        lGEMCEFDirectionSave = lGEMCEFSiteDirectionActive
        dNormBeforeSave = dGEMFunctionNorm
        dNormLastSave = dGEMFunctionNormLast
        dTraceReducedNormSave = dTraceSpeciesReducedNormLast
        dMassNormSave = dGEMMassBalanceNorm
        dChemNormSave = dGEMChemicalPotentialNorm
        dSolnChemNormSave = dGEMSolutionChemicalPotentialNorm
        dCondChemNormSave = dGEMCondensedChemicalPotentialNorm
        dSublatticeNormSave = dSublatticeExchangeNorm

        lHadAssemblageGEM = allocated(iAssemblageGEM)
        lHadMolesPhaseLast = allocated(dMolesPhaseLast)
        lHadMolFractionOld = allocated(dMolFractionOld)
        lHadTraceInactive = allocated(lTraceSpeciesInactive)
        lHadTraceReinjected = allocated(lTraceSpeciesReinjected)
        lHadActiveSlotThermo = allocated(iActiveSlotThermoPhase)
        lHadActiveSlotDisplay = allocated(iActiveSlotDisplayPhase)
        lHadActiveSlotIdentity = allocated(iActiveSlotIdentityOrdinal)
        lHadActiveSlotMol = allocated(dActiveSlotMolFraction)
        lHadActiveSlotSite = allocated(dActiveSlotSiteFraction)
        lHadGEMCEFSiteLast = allocated(dGEMCEFSiteLast)
        lHadGEMCEFPhaseSiteLast = allocated(dGEMCEFPhaseSiteLast)

        iAssemblageSave = iAssemblage
        iPhaseGEMSave = iPhaseGEM
        iPhaseLevelSave = iPhaseLevel
        dMolesPhaseSave = dMolesPhase
        dMolesSpeciesSave = dMolesSpecies
        dMolFractionSave = dMolFraction
        dElementPotentialSave = dElementPotential
        dElementPotentialLastSave = dElementPotentialLast
        dPhasePotentialSave = dPhasePotential
        dChemicalPotentialSave = dChemicalPotential
        dLevelingChemicalPotentialSave = dLevelingChemicalPotential
        dUpdateVarSave = dUpdateVar
        dUpdateVarLastSave = dUpdateVarLast
        dGibbsSolnPhaseSave = dGibbsSolnPhase
        dSumMolFractionSolnSave = dSumMolFractionSoln
        dDrivingForceSolnSave = dDrivingForceSoln
        dPartialExcessGibbsSave = dPartialExcessGibbs
        dPartialEnthalpyXSSave = dPartialEnthalpyXS
        dPartialEntropyXSSave = dPartialEntropyXS
        dPartialHeatCapacityXSSave = dPartialHeatCapacityXS
        dAtomFractionSpeciesSave = dAtomFractionSpecies
        dLevelingCompositionSpeciesSave = dLevelingCompositionSpecies
        dStoichSpeciesLevelSave = dStoichSpeciesLevel
        dChemicalPotentialGEMSave = dChemicalPotentialGEM
        dStoichSpeciesGEMSave = dStoichSpeciesGEM
        dAtomFractionSpeciesGEMSave = dAtomFractionSpeciesGEM
        dMolFractionGEMSave = dMolFractionGEM
        dMolesPhaseHistorySave = dMolesPhaseHistory
        dEffStoichSolnPhaseSave = dEffStoichSolnPhase
        iterHistoryLevelSave = iterHistoryLevel
        lSolnPhasesSave = lSolnPhases

        if (lHadAssemblageGEM) iAssemblageGEMSave = iAssemblageGEM
        if (lHadMolesPhaseLast) dMolesPhaseLastSave = dMolesPhaseLast
        if (lHadMolFractionOld) dMolFractionOldSave = dMolFractionOld
        if (lHadTraceInactive) lTraceSpeciesInactiveSave = lTraceSpeciesInactive
        if (lHadTraceReinjected) lTraceSpeciesReinjectedSave = lTraceSpeciesReinjected
        if (lHadActiveSlotThermo) iActiveSlotThermoSave = iActiveSlotThermoPhase
        if (lHadActiveSlotDisplay) iActiveSlotDisplaySave = iActiveSlotDisplayPhase
        if (lHadActiveSlotIdentity) iActiveSlotIdentitySave = iActiveSlotIdentityOrdinal
        if (lHadActiveSlotMol) dActiveSlotMolSave = dActiveSlotMolFraction
        if (lHadActiveSlotSite) dActiveSlotSiteSave = dActiveSlotSiteFraction
        if (lHadGEMCEFSiteLast) dGEMCEFSiteLastSave = dGEMCEFSiteLast
        if (lHadGEMCEFPhaseSiteLast) dGEMCEFPhaseSiteLastSave = dGEMCEFPhaseSiteLast

        return
    end subroutine SnapshotPEAState

    subroutine RestorePEAState

        nSpeciesLevel = nSpeciesLevelSave
        nConPhases = nConPhasesSave
        nSolnPhases = nSolnPhasesSave
        iterGlobal = iterGlobalSave
        iPhaseChangeReason = iPhaseChangeReasonSave
        iTraceSpeciesSlowProgressCount = iTraceSlowSave
        iGEMStagnationCount = iStagnationSave
        iterTraceSpeciesLastRemoval = iterTraceRemovalSave
        iterTraceSpeciesLastReinject = iterTraceReinjectSave
        nGEMCEFPhaseVariables = nGEMCEFPhaseVariablesSave
        nGEMCEFSiteVariables = nGEMCEFSiteVariablesSave
        lConverged = lConvergedSave
        lPhaseChange = lPhaseChangeSave
        lSkipLagrange = lSkipLagrangeSave
        lSampledLevelingThermoExtended = lSampledThermoSave
        lGEMCEFSiteLagrangianActive = lGEMCEFActiveSave
        lGEMCEFSiteDirectionActive = lGEMCEFDirectionSave
        lPEALagrangianPolishActive = .FALSE.
        dGEMFunctionNorm = dNormBeforeSave
        dGEMFunctionNormLast = dNormLastSave
        dTraceSpeciesReducedNormLast = dTraceReducedNormSave
        dGEMMassBalanceNorm = dMassNormSave
        dGEMChemicalPotentialNorm = dChemNormSave
        dGEMSolutionChemicalPotentialNorm = dSolnChemNormSave
        dGEMCondensedChemicalPotentialNorm = dCondChemNormSave
        dSublatticeExchangeNorm = dSublatticeNormSave

        call RestoreInteger1D(iAssemblage, iAssemblageSave)
        call RestoreInteger1D(iPhaseGEM, iPhaseGEMSave)
        call RestoreInteger1D(iPhaseLevel, iPhaseLevelSave)
        call RestoreReal1D(dMolesPhase, dMolesPhaseSave)
        call RestoreReal1D(dMolesSpecies, dMolesSpeciesSave)
        call RestoreReal1D(dMolFraction, dMolFractionSave)
        call RestoreReal1D(dElementPotential, dElementPotentialSave)
        call RestoreReal1D(dElementPotentialLast, dElementPotentialLastSave)
        call RestoreReal1D(dPhasePotential, dPhasePotentialSave)
        call RestoreReal1D(dChemicalPotential, dChemicalPotentialSave)
        call RestoreReal1D(dLevelingChemicalPotential, dLevelingChemicalPotentialSave)
        call RestoreReal1D(dUpdateVar, dUpdateVarSave)
        call RestoreReal1D(dUpdateVarLast, dUpdateVarLastSave)
        call RestoreReal1D(dGibbsSolnPhase, dGibbsSolnPhaseSave)
        call RestoreReal1D(dSumMolFractionSoln, dSumMolFractionSolnSave)
        call RestoreReal1D(dDrivingForceSoln, dDrivingForceSolnSave)
        call RestoreReal1D(dPartialExcessGibbs, dPartialExcessGibbsSave)
        call RestoreReal1D(dPartialEnthalpyXS, dPartialEnthalpyXSSave)
        call RestoreReal1D(dPartialEntropyXS, dPartialEntropyXSSave)
        call RestoreReal1D(dPartialHeatCapacityXS, dPartialHeatCapacityXSSave)
        call RestoreReal2D(dAtomFractionSpecies, dAtomFractionSpeciesSave)
        call RestoreReal2D(dLevelingCompositionSpecies, dLevelingCompositionSpeciesSave)
        call RestoreReal2D(dStoichSpeciesLevel, dStoichSpeciesLevelSave)
        call RestoreReal1D(dChemicalPotentialGEM, dChemicalPotentialGEMSave)
        call RestoreReal2D(dStoichSpeciesGEM, dStoichSpeciesGEMSave)
        call RestoreReal2D(dAtomFractionSpeciesGEM, dAtomFractionSpeciesGEMSave)
        call RestoreReal2D(dMolFractionGEM, dMolFractionGEMSave)
        call RestoreReal2D(dMolesPhaseHistory, dMolesPhaseHistorySave)
        call RestoreReal2D(dEffStoichSolnPhase, dEffStoichSolnPhaseSave)
        call RestoreInteger2D(iterHistoryLevel, iterHistoryLevelSave)
        call RestoreLogical1D(lSolnPhases, lSolnPhasesSave)

        if (lHadAssemblageGEM) then
            call RestoreInteger1D(iAssemblageGEM, iAssemblageGEMSave)
        else
            if (allocated(iAssemblageGEM)) deallocate(iAssemblageGEM)
        end if

        if (lHadMolesPhaseLast) then
            call RestoreReal1D(dMolesPhaseLast, dMolesPhaseLastSave)
        else
            if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)
        end if

        if (lHadMolFractionOld) then
            call RestoreReal1D(dMolFractionOld, dMolFractionOldSave)
        else
            if (allocated(dMolFractionOld)) deallocate(dMolFractionOld)
        end if

        if (lHadTraceInactive) then
            call RestoreLogical1D(lTraceSpeciesInactive, lTraceSpeciesInactiveSave)
        else
            if (allocated(lTraceSpeciesInactive)) deallocate(lTraceSpeciesInactive)
        end if

        if (lHadTraceReinjected) then
            call RestoreLogical1D(lTraceSpeciesReinjected, lTraceSpeciesReinjectedSave)
        else
            if (allocated(lTraceSpeciesReinjected)) deallocate(lTraceSpeciesReinjected)
        end if

        if (lHadActiveSlotThermo) then
            call RestoreInteger1D(iActiveSlotThermoPhase, iActiveSlotThermoSave)
        else
            if (allocated(iActiveSlotThermoPhase)) deallocate(iActiveSlotThermoPhase)
        end if

        if (lHadActiveSlotDisplay) then
            call RestoreInteger1D(iActiveSlotDisplayPhase, iActiveSlotDisplaySave)
        else
            if (allocated(iActiveSlotDisplayPhase)) deallocate(iActiveSlotDisplayPhase)
        end if

        if (lHadActiveSlotIdentity) then
            call RestoreInteger1D(iActiveSlotIdentityOrdinal, iActiveSlotIdentitySave)
        else
            if (allocated(iActiveSlotIdentityOrdinal)) deallocate(iActiveSlotIdentityOrdinal)
        end if

        if (lHadActiveSlotMol) then
            call RestoreReal2D(dActiveSlotMolFraction, dActiveSlotMolSave)
        else
            if (allocated(dActiveSlotMolFraction)) deallocate(dActiveSlotMolFraction)
        end if

        if (lHadActiveSlotSite) then
            call RestoreReal3D(dActiveSlotSiteFraction, dActiveSlotSiteSave)
        else
            if (allocated(dActiveSlotSiteFraction)) deallocate(dActiveSlotSiteFraction)
        end if

        if (lHadGEMCEFSiteLast) then
            call RestoreReal3D(dGEMCEFSiteLast, dGEMCEFSiteLastSave)
        else
            if (allocated(dGEMCEFSiteLast)) deallocate(dGEMCEFSiteLast)
        end if

        if (lHadGEMCEFPhaseSiteLast) then
            call RestoreReal3D(dGEMCEFPhaseSiteLast, dGEMCEFPhaseSiteLastSave)
        else
            if (allocated(dGEMCEFPhaseSiteLast)) deallocate(dGEMCEFPhaseSiteLast)
        end if

        return
    end subroutine RestorePEAState

    subroutine RestoreInteger1D(iTarget, iSource)
        integer, dimension(:), allocatable, intent(inout) :: iTarget
        integer, dimension(:), intent(in)                 :: iSource

        if (allocated(iTarget)) deallocate(iTarget)
        allocate(iTarget(SIZE(iSource)))
        iTarget = iSource

        return
    end subroutine RestoreInteger1D

    subroutine RestoreInteger2D(iTarget, iSource)
        integer, dimension(:,:), allocatable, intent(inout) :: iTarget
        integer, dimension(:,:), intent(in)                 :: iSource

        if (allocated(iTarget)) deallocate(iTarget)
        allocate(iTarget(SIZE(iSource,1),SIZE(iSource,2)))
        iTarget = iSource

        return
    end subroutine RestoreInteger2D

    subroutine RestoreReal1D(dTarget, dSource)
        real(8), dimension(:), allocatable, intent(inout) :: dTarget
        real(8), dimension(:), intent(in)                 :: dSource

        if (allocated(dTarget)) deallocate(dTarget)
        allocate(dTarget(SIZE(dSource)))
        dTarget = dSource

        return
    end subroutine RestoreReal1D

    subroutine RestoreReal2D(dTarget, dSource)
        real(8), dimension(:,:), allocatable, intent(inout) :: dTarget
        real(8), dimension(:,:), intent(in)                 :: dSource

        if (allocated(dTarget)) deallocate(dTarget)
        allocate(dTarget(SIZE(dSource,1),SIZE(dSource,2)))
        dTarget = dSource

        return
    end subroutine RestoreReal2D

    subroutine RestoreReal3D(dTarget, dSource)
        real(8), dimension(:,:,:), allocatable, intent(inout) :: dTarget
        real(8), dimension(:,:,:), intent(in)                 :: dSource

        if (allocated(dTarget)) deallocate(dTarget)
        allocate(dTarget(SIZE(dSource,1),SIZE(dSource,2),SIZE(dSource,3)))
        dTarget = dSource

        return
    end subroutine RestoreReal3D

    subroutine RestoreLogical1D(lTarget, lSource)
        logical, dimension(:), allocatable, intent(inout) :: lTarget
        logical, dimension(:), intent(in)                 :: lSource

        if (allocated(lTarget)) deallocate(lTarget)
        allocate(lTarget(SIZE(lSource)))
        lTarget = lSource

        return
    end subroutine RestoreLogical1D

    logical function lFiniteVector(dVector)
        real(8), dimension(:), intent(in) :: dVector
        integer                           :: i

        lFiniteVector = .TRUE.
        do i = 1, SIZE(dVector)
            if (dVector(i) /= dVector(i)) then
                lFiniteVector = .FALSE.
                return
            end if
            if (DABS(dVector(i)) > HUGE(1D0) / 100D0) then
                lFiniteVector = .FALSE.
                return
            end if
        end do

        return
    end function lFiniteVector

end subroutine PEALagrangianPolish
