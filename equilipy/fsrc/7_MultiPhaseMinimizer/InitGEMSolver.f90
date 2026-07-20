!> \brief Initialize Leveling and Lagrangian GEM solver state.
!!
!! \details Allocates and resets the global minimizer arrays used by Leveling,
!! Leveling-to-Lagrangian translation, Lagrangian GEM, trace-species handling,
!! and per-iteration diagnostics.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    InitGEMSolver.f90
    !> \brief   Initialize shared GEM solver state before Leveling.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      RunLeveling.f90
    !> \sa      GetFirstAssemblage.f90
    !> \sa      LevelingSolver.f90
    !> \sa      RunLagrangianGEM.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original GEM solver initialization
    !   07/20/2026      S.Y. Kwon           Allocated typed candidate, order/disorder, reporting, grid-certificate, and atomic PEA generation state.
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to prepare the minimizer for
    !! a fresh equilibrium calculation.  It creates the initial Leveling
    !! candidate arrays, records the initial active assemblage, allocates
    !! Lagrangian GEM work arrays, and resets diagnostic histories.
    !
    ! Required input variables:
    ! =========================
    !
    ! nSpecies             Number of species in the current system.
    ! nElements            Number of conserved system components.
    ! nSolnPhasesSys       Number of solution phases in the current system.
    ! dStoichSpecies       Species stoichiometry from InitThermo/CheckSystem.
    ! dLevelingChemicalPotential
    !                      Leveling chemical potentials after order/disorder corrections.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! iAssemblage          Initial Leveling assemblage.
    ! dMolesSpecies        Solver species moles initialized for Leveling/Lagrangian.
    ! dElementPotential    Element potentials initialized for Leveling/Lagrangian.
    ! dGEM*History         Per-iteration Lagrangian diagnostic histories.
    ! dPEATol              Named PEA residual tolerance available to RunLeveling.
    ! dToleranceLevel      Signed phase-potential certificate tolerance.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! ApplyOrderDisorderLevelingPotentials   Updates Leveling rows for SUBOM phases.
    ! InitLevelCandidatePool                 Resets optional Leveling candidate handoff storage.
    ! GetFirstAssemblage                     Selects the first Leveling assemblage.
    ! InitTraceSpeciesMask                   Resets Lagrangian trace-species state.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! RunLeveling       Starts the production Leveling stage.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Sampled Leveling is disabled here; current production minimization
    !   starts from classical Leveling.
    ! - Diagnostic histories must be reset here so repeated staged probes do
    !   not carry stale Lagrangian information.
    ! - RunLeveling executes before InitCheckPhaseAssemblage, so both named
    !   certificate tolerances must already hold their standard PEA values.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine InitGEMSolver
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    implicit none
    integer::   i, j, k, l
    ! Initialize variables:

    lPostProcess = .False.
    nGridRecoveryAttempt = 0
    lGridRecoveryPassActive = .FALSE.
    if (allocated(lGridRecoveryPhase)) deallocate(lGridRecoveryPhase)
!
    nSpeciesLevel = nSpecies

    ! Check to see if allocatable arrays have already been allocated:
    if (allocated(dMolesPhase))             deallocate(dMolesPhase)
    if (allocated(dMolesSpecies))           deallocate(dMolesSpecies)
    if (allocated(dElementPotential))       deallocate(dElementPotential)
    if (allocated(iAssemblage))             deallocate(iAssemblage)
    if (allocated(dPhasePotential))         deallocate(dPhasePotential)
    if (allocated(iShuffled))               deallocate(iShuffled)
    allocate(dMolesPhase(nElements),dMolesSpecies(nSpecies),dElementPotential(nElements),&
        iAssemblage(nElements),iShuffled(nElements),dPhasePotential(nSpeciesLevel))
    
    ! Initialize allocatable variables:
    dElementPotential = 0D0
    iAssemblage       = 0
    iShuffled         = (/(i, i = 1,nElements, 1)/)
    dMolesPhase       = 0D0
    dMolesSpecies     = 0D0


    ! Initialize variables for GetFirstAssemblage and GetNewAssemblage
    if (allocated(iPhaseLevel))             deallocate(iPhaseLevel)
    if (allocated(dStoichSpeciesLevel))     deallocate(dStoichSpeciesLevel)
    if (allocated(iterHistoryLevel))        deallocate(iterHistoryLevel)
    if (allocated(dMolesPhaseHistory))      deallocate(dMolesPhaseHistory)
    

    if (allocated(dChemicalPotentialGEM))   deallocate(dChemicalPotentialGEM)
    if (allocated(dAtomFractionSpeciesGEM)) deallocate(dAtomFractionSpeciesGEM)
    if (allocated(dStoichSpeciesGEM))       deallocate(dStoichSpeciesGEM)
    if (allocated(dMolFractionGEM))         deallocate(dMolFractionGEM)
    if (allocated(dActiveSlotMolFraction))  deallocate(dActiveSlotMolFraction)
    if (allocated(dActiveSlotSiteFraction)) deallocate(dActiveSlotSiteFraction)
    if (allocated(dActiveSlotChemPot))       deallocate(dActiveSlotChemPot)
    if (allocated(dActiveSlotPartialH))      deallocate(dActiveSlotPartialH)
    if (allocated(dActiveSlotPartialS))      deallocate(dActiveSlotPartialS)
    if (allocated(dActiveSlotPartialCp))     deallocate(dActiveSlotPartialCp)
    if (allocated(lActiveSlotPropValid))     deallocate(lActiveSlotPropValid)
    if (allocated(dMolFraction))        deallocate(dMolFraction)
    if (allocated(iPhaseGEM))               deallocate(iPhaseGEM)
    if (allocated(iActiveSlotThermoPhase))  deallocate(iActiveSlotThermoPhase)
    if (allocated(iActiveSlotDisplayPhase)) deallocate(iActiveSlotDisplayPhase)
    if (allocated(iActiveSlotIdentityOrdinal)) deallocate(iActiveSlotIdentityOrdinal)
    if (allocated(iActiveSlotODClass)) deallocate(iActiveSlotODClass)

    ! Allocate variables for GetFirstAssemblage and GetNewAssemblage
    allocate(iPhaseLevel(nSpeciesLevel),dStoichSpeciesLevel(nSpeciesLevel,nElements),&
    iterHistoryLevel(nElements,500),dMolesPhaseHistory(nElements,500),dChemicalPotentialGEM(nElements),&
    dAtomFractionSpeciesGEM(nElements,nElements),dStoichSpeciesGEM(nElements,nElements),&
    dMolFractionGEM(nElements,nSpecies),dActiveSlotMolFraction(nElements,nSpecies),&
    dActiveSlotSiteFraction(nElements,nMaxSublatticeSys,nMaxConstituentSys),&
    dActiveSlotChemPot(nElements,nSpecies),dActiveSlotPartialH(nElements,nSpecies),&
    dActiveSlotPartialS(nElements,nSpecies),dActiveSlotPartialCp(nElements,nSpecies),&
    lActiveSlotPropValid(nElements),&
    iPhaseGEM(nElements),dMolFraction(nSpecies),&
    iActiveSlotThermoPhase(nElements),iActiveSlotDisplayPhase(nElements),&
    iActiveSlotIdentityOrdinal(nElements),iActiveSlotODClass(nElements))
    
    iPhaseLevel          = 0
    dStoichSpeciesLevel  = 0D0
    iPhaseLevel(:nSpecies)          = iPhase
    dStoichSpeciesLevel(:nSpecies,:)  = dStoichSpecies
    iterHistoryLevel     = 0
    dMolFractionGEM      = 1D0
    dActiveSlotMolFraction = 0D0
    dActiveSlotSiteFraction = 0D0
    dActiveSlotChemPot = 0D0
    dActiveSlotPartialH = 0D0
    dActiveSlotPartialS = 0D0
    dActiveSlotPartialCp = 0D0
    lActiveSlotPropValid = .FALSE.
    dMolFraction = 1D0
    iActiveSlotThermoPhase = 0
    iActiveSlotDisplayPhase = 0
    iActiveSlotIdentityOrdinal = 0
    iActiveSlotODClass = OD_CANDIDATE_NOT_EVALUATED

    call ApplyOrderDisorderLevelingPotentials

    call InitLevelCandidatePool(0)

    ! Establish the very first phase assemblage:
    call GetFirstAssemblage

    do i =1,nElements
        iPhaseGEM(i)= iPhaseLevel(iAssemblage(i))
        if (iPhaseGEM(i) > 0) then
            iActiveSlotThermoPhase(i) = iPhaseGEM(i)
            if (allocated(iODCompanionPhase)) then
                do k = 1, MIN(nSolnPhasesSys, SIZE(iODCompanionPhase))
                    if (TRIM(cSolnPhaseType(k)) /= 'SUBOM') cycle
                    if (iODCompanionPhase(k) == iPhaseGEM(i)) then
                        if (allocated(iODStandalonePhase)) then
                            if (k <= SIZE(iODStandalonePhase)) then
                                if (iODStandalonePhase(k) == iPhaseGEM(i)) cycle
                            end if
                        end if
                        iActiveSlotThermoPhase(i) = k
                        exit
                    end if
                end do
            end if
            iActiveSlotDisplayPhase(i) = iPhaseGEM(i)
            iActiveSlotIdentityOrdinal(i) = 1
            do j = 1, i - 1
                if (iActiveSlotThermoPhase(j) == iActiveSlotThermoPhase(i)) then
                    iActiveSlotIdentityOrdinal(i) = iActiveSlotIdentityOrdinal(i) + 1
                end if
            end do
        end if
    end do

    do j = 1,nElements
        do i = 1,nElements
            dStoichSpeciesGEM(j,i)=dStoichSpeciesLevel(iAssemblage(j),i)
            dAtomFractionSpeciesGEM(j,i)=dLevelingCompositionSpecies(iAssemblage(j),i)
        end do
        dChemicalPotentialGEM(j) = dLevelingChemicalPotential(iAssemblage(j))
    end do

    ! Store the conditions for initial phase assemblage
    if(allocated(iPhaseGEMinit)) deallocate(iPhaseGEMinit)
    if(allocated(dStoichSpeciesGEMinit)) deallocate(dStoichSpeciesGEMinit)
    if(allocated(dAtomFractionSpeciesGEMinit)) deallocate(dAtomFractionSpeciesGEMinit)
    if(allocated(dChemicalPotentialGEMinit)) deallocate(dChemicalPotentialGEMinit)
    if(allocated(dMolFractionGEMinit)) deallocate(dMolFractionGEMinit)
    if(allocated(iAssemblageGEMinit)) deallocate(iAssemblageGEMinit)

    allocate(iPhaseGEMinit(nElements),dStoichSpeciesGEMinit(nElements,nElements),&
    dAtomFractionSpeciesGEMinit(nElements,nElements),dChemicalPotentialGEMinit(nElements),&
    dMolFractionGEMinit(nElements,nSpecies),iAssemblageGEMinit(nElements))
    iPhaseGEMinit = iPhaseGEM
    dStoichSpeciesGEMinit = dStoichSpeciesGEM
    dAtomFractionSpeciesGEMinit = dAtomFractionSpeciesGEM
    dChemicalPotentialGEMinit = dChemicalPotentialGEM
    dMolFractionGEMinit = dMolFractionGEM
    iAssemblageGEMinit = iAssemblage
    

    ! This part is for Lagrangian Multipler method
        ! Allocate memory:
    l = MAX(1,nSolnPhasesSys)
    if(allocated(dUpdateVar)) deallocate(dUpdateVar)
    if(allocated(iterHistory)) deallocate(iterHistory)
    if(allocated(dUpdateVarLast)) deallocate(dUpdateVarLast)
    if(allocated(lPhaseChangeHistory)) deallocate(lPhaseChangeHistory)
    if(allocated(iPhaseChangeReasonHistory)) deallocate(iPhaseChangeReasonHistory)
    if(allocated(iGEMNewtonSolverHist)) deallocate(iGEMNewtonSolverHist)
    if(allocated(iGEMNewtonInfoHist)) deallocate(iGEMNewtonInfoHist)
    if(allocated(iGEMNewtonKKTSizeHist)) deallocate(iGEMNewtonKKTSizeHist)
    if(allocated(iGEMNewtonPivot1x1Hist)) deallocate(iGEMNewtonPivot1x1Hist)
    if(allocated(iGEMNewtonPivot2x2Hist)) deallocate(iGEMNewtonPivot2x2Hist)
    if(allocated(iGEMNewtonPivotPositiveHist)) deallocate(iGEMNewtonPivotPositiveHist)
    if(allocated(iGEMNewtonPivotNegativeHist)) deallocate(iGEMNewtonPivotNegativeHist)
    if(allocated(iGEMNewtonPivotZeroHist)) deallocate(iGEMNewtonPivotZeroHist)
    if(allocated(iGEMLSIterHist)) deallocate(iGEMLSIterHist)
    if(allocated(iGEMLSNegFactorHist)) deallocate(iGEMLSNegFactorHist)
    if(allocated(iGEMLSNegPhaseHist)) deallocate(iGEMLSNegPhaseHist)
    if(allocated(iGEMLSFloorHist)) deallocate(iGEMLSFloorHist)
    if(allocated(iGEMLSNoDescentHist)) deallocate(iGEMLSNoDescentHist)
    if(allocated(iGEMLSNoDescentClassHist)) deallocate(iGEMLSNoDescentClassHist)
    if(allocated(iGEMPreLMNoDescentHist)) deallocate(iGEMPreLMNoDescentHist)
    if(allocated(iGEMPreLMNoDescentClassHist)) deallocate(iGEMPreLMNoDescentClassHist)
    if(allocated(iGEMTraceInactiveHist)) deallocate(iGEMTraceInactiveHist)
    if(allocated(iGEMTraceReinjectedHist)) deallocate(iGEMTraceReinjectedHist)
    if(allocated(iGEMTraceRemoveHist)) deallocate(iGEMTraceRemoveHist)
    if(allocated(iGEMTraceReinjectHist)) deallocate(iGEMTraceReinjectHist)
    if(allocated(iGEMResidualLMHist)) deallocate(iGEMResidualLMHist)
    if(allocated(iGEMResidualLMRawNegativeHist)) deallocate(iGEMResidualLMRawNegativeHist)
    if(allocated(iGEMResidualLMNoDescentHist)) deallocate(iGEMResidualLMNoDescentHist)
    if(allocated(iGEMResidualLMNoDescentClassHist)) deallocate(iGEMResidualLMNoDescentClassHist)
    if(allocated(iGEMResidualLMEventCallSite)) deallocate(iGEMResidualLMEventCallSite)
    if(allocated(iGEMResidualLMEventIterPEA)) deallocate(iGEMResidualLMEventIterPEA)
    if(allocated(iGEMResidualLMEventIterGlobal)) deallocate(iGEMResidualLMEventIterGlobal)
    if(allocated(iGEMResidualLMEventReason)) deallocate(iGEMResidualLMEventReason)
    if(allocated(iGEMResidualLMEventNoDescentClass)) deallocate(iGEMResidualLMEventNoDescentClass)
    if(allocated(iGEMResidualLMEventNewtonInfo)) deallocate(iGEMResidualLMEventNewtonInfo)
    if(allocated(iGEMResidualLMEventKKTSize)) deallocate(iGEMResidualLMEventKKTSize)
    if(allocated(iGEMResidualLMEventNPrimal)) deallocate(iGEMResidualLMEventNPrimal)
    if(allocated(iGEMResidualLMEventNConstraint)) deallocate(iGEMResidualLMEventNConstraint)
    if(allocated(iGEMResidualLMEventPivot1x1)) deallocate(iGEMResidualLMEventPivot1x1)
    if(allocated(iGEMResidualLMEventPivot2x2)) deallocate(iGEMResidualLMEventPivot2x2)
    if(allocated(iGEMResidualLMEventPivotPositive)) deallocate(iGEMResidualLMEventPivotPositive)
    if(allocated(iGEMResidualLMEventPivotNegative)) deallocate(iGEMResidualLMEventPivotNegative)
    if(allocated(iGEMResidualLMEventPivotZero)) deallocate(iGEMResidualLMEventPivotZero)
    if(allocated(iGEMLMEventODOrdPhase)) then
        deallocate(iGEMLMEventODOrdPhase)
    end if
    if(allocated(iGEMLMEventODCompPhase)) then
        deallocate(iGEMLMEventODCompPhase)
    end if
    if(allocated(iGEMResidualLMEventAssemblage)) deallocate(iGEMResidualLMEventAssemblage)
    if(allocated(iGEMCoalesceHist)) deallocate(iGEMCoalesceHist)
    if(allocated(iGEMStabilizeHist)) deallocate(iGEMStabilizeHist)
    if(allocated(iGEMBoundaryRemovalHist)) deallocate(iGEMBoundaryRemovalHist)
    if(allocated(iGEMRawNegativeRemovalHist)) deallocate(iGEMRawNegativeRemovalHist)
    if(allocated(iGEMRawNegativePhaseSlotHist)) deallocate(iGEMRawNegativePhaseSlotHist)
    if(allocated(iGEMRawNegativePhaseSolnHist)) deallocate(iGEMRawNegativePhaseSolnHist)
    if(allocated(iGEMRawNegativePhaseSpeciesHist)) deallocate(iGEMRawNegativePhaseSpeciesHist)
    if(allocated(iGEMRawNegCompHist)) then
        deallocate(iGEMRawNegCompHist)
    end if
    if(allocated(iGEMInvalidCompBoundAttemptedHist)) deallocate(iGEMInvalidCompBoundAttemptedHist)
    if(allocated(iGEMInvalidCompBoundAcceptedHist)) deallocate(iGEMInvalidCompBoundAcceptedHist)
    if(allocated(iGEMInvalidCompBoundRejectedHist)) deallocate(iGEMInvalidCompBoundRejectedHist)
    if(allocated(iGEMInvalidCompBoundVerdictHist)) deallocate(iGEMInvalidCompBoundVerdictHist)
    if(allocated(iGEMBoundaryPinnedRemovalHist)) deallocate(iGEMBoundaryPinnedRemovalHist)
    if(allocated(iGEMBoundaryRemovalSlotHist)) deallocate(iGEMBoundaryRemovalSlotHist)
    if(allocated(iGEMBoundaryRemovalPhaseHist)) deallocate(iGEMBoundaryRemovalPhaseHist)
    if(allocated(iGEMTinyBoundaryRemovalHist)) deallocate(iGEMTinyBoundaryRemovalHist)
    if(allocated(iGEMTinyBoundaryRemovalSlotHist)) deallocate(iGEMTinyBoundaryRemovalSlotHist)
    if(allocated(iGEMTinyBoundaryRemovalPhaseHist)) deallocate(iGEMTinyBoundaryRemovalPhaseHist)
    if(allocated(iGEMBoundaryRankGuardHist)) deallocate(iGEMBoundaryRankGuardHist)
    if(allocated(iGEMBoundaryRankGuardSlotHist)) deallocate(iGEMBoundaryRankGuardSlotHist)
    if(allocated(iGEMBoundaryRankGuardPhaseHist)) deallocate(iGEMBoundaryRankGuardPhaseHist)
    if(allocated(iGEMTraceRemoveSpeciesHist)) deallocate(iGEMTraceRemoveSpeciesHist)
    if(allocated(iGEMTraceRemovePhaseHist)) deallocate(iGEMTraceRemovePhaseHist)
    if(allocated(iGEMTraceRemoveCountHist)) deallocate(iGEMTraceRemoveCountHist)
    if(allocated(iGEMTraceReinjectSpeciesHist)) deallocate(iGEMTraceReinjectSpeciesHist)
    if(allocated(iGEMTraceReinjectPhaseHist)) deallocate(iGEMTraceReinjectPhaseHist)
    if(allocated(iGEMTraceReinjectCountHist)) deallocate(iGEMTraceReinjectCountHist)
    if(allocated(iGEMCEFRetryActivationHist)) deallocate(iGEMCEFRetryActivationHist)
    if(allocated(iGEMSUBGQRideAlongHist)) deallocate(iGEMSUBGQRideAlongHist)
    if(allocated(iODPairHist)) deallocate(iODPairHist)
    if(allocated(iODOrderingModeHist)) deallocate(iODOrderingModeHist)
    if(allocated(dGEMNormHist)) deallocate(dGEMNormHist)
    if(allocated(dGEMMassNormHist)) deallocate(dGEMMassNormHist)
    if(allocated(dGEMChemNormHist)) deallocate(dGEMChemNormHist)
    if(allocated(dGEMSolnChemNormHist)) deallocate(dGEMSolnChemNormHist)
    if(allocated(dGEMCondChemNormHist)) deallocate(dGEMCondChemNormHist)
    if(allocated(dGEMSublExchangeNormHist)) deallocate(dGEMSublExchangeNormHist)
    if(allocated(dGEMNewtonSymResidualHist)) deallocate(dGEMNewtonSymResidualHist)
    if(allocated(dGEMNewtonMinPivotScaleHist)) deallocate(dGEMNewtonMinPivotScaleHist)
    if(allocated(dGEMNewtonMaxPivotScaleHist)) deallocate(dGEMNewtonMaxPivotScaleHist)
    if(allocated(dGEMNewtonDirectionNormHist)) deallocate(dGEMNewtonDirectionNormHist)
    if(allocated(dGEMLSInitialNormHist)) deallocate(dGEMLSInitialNormHist)
    if(allocated(dGEMLSBestNormHist)) deallocate(dGEMLSBestNormHist)
    if(allocated(dGEMLSFinalNormHist)) deallocate(dGEMLSFinalNormHist)
    if(allocated(dGEMLSInitialStepHist)) deallocate(dGEMLSInitialStepHist)
    if(allocated(dGEMLSBestStepHist)) deallocate(dGEMLSBestStepHist)
    if(allocated(dGEMLSFinalStepHist)) deallocate(dGEMLSFinalStepHist)
    if(allocated(dGEMLSMinRawPhaseMolesHist)) deallocate(dGEMLSMinRawPhaseMolesHist)
    if(allocated(dGEMLSMinFinalPhaseMolesHist)) deallocate(dGEMLSMinFinalPhaseMolesHist)
    if(allocated(dGEMLSInitialGibbsHist)) deallocate(dGEMLSInitialGibbsHist)
    if(allocated(dGEMLSBestGibbsHist)) deallocate(dGEMLSBestGibbsHist)
    if(allocated(dGEMLSFinalGibbsHist)) deallocate(dGEMLSFinalGibbsHist)
    if(allocated(dGEMLSInitialMeritHist)) deallocate(dGEMLSInitialMeritHist)
    if(allocated(dGEMLSBestMeritHist)) deallocate(dGEMLSBestMeritHist)
    if(allocated(dGEMLSFinalMeritHist)) deallocate(dGEMLSFinalMeritHist)
    if(allocated(dGEMPreLMInitNormHist)) deallocate(dGEMPreLMInitNormHist)
    if(allocated(dGEMPreLMBestNormHist)) deallocate(dGEMPreLMBestNormHist)
    if(allocated(dGEMPreLMInitGibbsHist)) deallocate(dGEMPreLMInitGibbsHist)
    if(allocated(dGEMPreLMBestGibbsHist)) deallocate(dGEMPreLMBestGibbsHist)
    if(allocated(dGEMPreLMInitMeritHist)) deallocate(dGEMPreLMInitMeritHist)
    if(allocated(dGEMPreLMBestMeritHist)) deallocate(dGEMPreLMBestMeritHist)
    if(allocated(dGEMPreLMMeritCandNormHist)) deallocate(dGEMPreLMMeritCandNormHist)
    if(allocated(dGEMPreLMMeritCandMassHist)) deallocate(dGEMPreLMMeritCandMassHist)
    if(allocated(dGEMPreLMMeritCandStepHist)) deallocate(dGEMPreLMMeritCandStepHist)
    if(allocated(dGEMPreLMMeritCandGibbsHist)) deallocate(dGEMPreLMMeritCandGibbsHist)
    if(allocated(dGEMPreLMMeritCandMeritHist)) deallocate(dGEMPreLMMeritCandMeritHist)
    if(allocated(dGEMResidualLMEventNormSlope)) deallocate(dGEMResidualLMEventNormSlope)
    if(allocated(dGEMResidualLMEventGibbsSlope)) deallocate(dGEMResidualLMEventGibbsSlope)
    if(allocated(dGEMResidualLMEventMeritSlope)) deallocate(dGEMResidualLMEventMeritSlope)
    if(allocated(dGEMResidualLMEventInitialNorm)) deallocate(dGEMResidualLMEventInitialNorm)
    if(allocated(dGEMResidualLMEventBestNorm)) deallocate(dGEMResidualLMEventBestNorm)
    if(allocated(dGEMResidualLMEventFinalNorm)) deallocate(dGEMResidualLMEventFinalNorm)
    if(allocated(dGEMResidualLMEventInitialGibbs)) deallocate(dGEMResidualLMEventInitialGibbs)
    if(allocated(dGEMResidualLMEventBestGibbs)) deallocate(dGEMResidualLMEventBestGibbs)
    if(allocated(dGEMResidualLMEventFinalGibbs)) deallocate(dGEMResidualLMEventFinalGibbs)
    if(allocated(dGEMResidualLMEventInitialMerit)) deallocate(dGEMResidualLMEventInitialMerit)
    if(allocated(dGEMResidualLMEventBestMerit)) deallocate(dGEMResidualLMEventBestMerit)
    if(allocated(dGEMResidualLMEventFinalMerit)) deallocate(dGEMResidualLMEventFinalMerit)
    if(allocated(dGEMResidualLMEventMeritCandNorm)) deallocate(dGEMResidualLMEventMeritCandNorm)
    if(allocated(dGEMResidualLMEventMeritCandMass)) deallocate(dGEMResidualLMEventMeritCandMass)
    if(allocated(dGEMResidualLMEventMeritCandStep)) deallocate(dGEMResidualLMEventMeritCandStep)
    if(allocated(dGEMResidualLMEventMeritCandGibbs)) deallocate(dGEMResidualLMEventMeritCandGibbs)
    if(allocated(dGEMResidualLMEventMeritCandMerit)) deallocate(dGEMResidualLMEventMeritCandMerit)
    if(allocated(dGEMResidualLMEventDirectionNorm)) deallocate(dGEMResidualLMEventDirectionNorm)
    if(allocated(dGEMResidualLMEventMinPivot)) deallocate(dGEMResidualLMEventMinPivot)
    if(allocated(dGEMResidualLMEventMaxPivot)) deallocate(dGEMResidualLMEventMaxPivot)
    if(allocated(dGEMLMEventODAlign)) then
        deallocate(dGEMLMEventODAlign)
    end if
    if(allocated(dGEMLMEventODEigen)) then
        deallocate(dGEMLMEventODEigen)
    end if
    if(allocated(dGEMNewtonDirNormSlopeHist)) deallocate(dGEMNewtonDirNormSlopeHist)
    if(allocated(dGEMNewtonDirGibbsSlopeHist)) deallocate(dGEMNewtonDirGibbsSlopeHist)
    if(allocated(dGEMNewtonDirMeritSlopeHist)) deallocate(dGEMNewtonDirMeritSlopeHist)
    if(allocated(dGEMRawNegativePhaseAmountHist)) deallocate(dGEMRawNegativePhaseAmountHist)
    if(allocated(dGEMRawNegativePhaseDirectionHist)) deallocate(dGEMRawNegativePhaseDirectionHist)
    if(allocated(dGEMRawNegPhaseResidualHist)) deallocate(dGEMRawNegPhaseResidualHist)
    if(allocated(dGEMInvalidCompBoundPhiHist)) deallocate(dGEMInvalidCompBoundPhiHist)
    if(allocated(dODCompDistHist)) deallocate(dODCompDistHist)
    if(allocated(dODOrderNormHist)) deallocate(dODOrderNormHist)
    if(allocated(dODOrderingEigenMinHist)) deallocate(dODOrderingEigenMinHist)
    if(allocated(dGEMLSRawPhaseMoles)) deallocate(dGEMLSRawPhaseMoles)
    if(allocated(dGEMLSFinalPhaseMoles)) deallocate(dGEMLSFinalPhaseMoles)
    if(allocated(dGEMLSRawPhaseMolesHist)) deallocate(dGEMLSRawPhaseMolesHist)
    if(allocated(dGEMLSFinalPhaseMolesHist)) deallocate(dGEMLSFinalPhaseMolesHist)
    if(allocated(iGEMCEFPhaseSlot)) deallocate(iGEMCEFPhaseSlot)
    if(allocated(iGEMCEFPhaseSoln)) deallocate(iGEMCEFPhaseSoln)
    if(allocated(iGEMCEFPhaseSpecies)) deallocate(iGEMCEFPhaseSpecies)
    if(allocated(iGEMCEFVarPhaseVar)) deallocate(iGEMCEFVarPhaseVar)
    if(allocated(iGEMCEFVarSolnPhase)) deallocate(iGEMCEFVarSolnPhase)
    if(allocated(iGEMCEFVarPhaseID)) deallocate(iGEMCEFVarPhaseID)
    if(allocated(iGEMCEFVarSub)) deallocate(iGEMCEFVarSub)
    if(allocated(iGEMCEFVarCon)) deallocate(iGEMCEFVarCon)
    if(allocated(iGEMCEFVarRef)) deallocate(iGEMCEFVarRef)
    if(allocated(nGEMCEFVarTie)) deallocate(nGEMCEFVarTie)
    if(allocated(iGEMCEFVarTieSub)) deallocate(iGEMCEFVarTieSub)
    if(allocated(dGEMCEFPhaseDirection)) deallocate(dGEMCEFPhaseDirection)
    if(allocated(dGEMCEFSiteDirection)) deallocate(dGEMCEFSiteDirection)
    if(allocated(dGEMCEFElementDirection)) deallocate(dGEMCEFElementDirection)
    if(allocated(dGEMCEFPhaseResidual)) deallocate(dGEMCEFPhaseResidual)
    if(allocated(dGEMCEFSiteLast)) deallocate(dGEMCEFSiteLast)
    if(allocated(dGEMCEFPhaseSiteLast)) deallocate(dGEMCEFPhaseSiteLast)
    if(allocated(dGEMAnalyticalSpeciesDirection)) deallocate(dGEMAnalyticalSpeciesDirection)
    if(allocated(lGEMAnalyticalSpeciesDirection)) deallocate(lGEMAnalyticalSpeciesDirection)
    if(allocated(iPEALevelIterHist)) deallocate(iPEALevelIterHist)
    if(allocated(iPEAPlaneGenerationHist)) deallocate(iPEAPlaneGenerationHist)
    if(allocated(iPEADFSweepGenerationHist)) deallocate(iPEADFSweepGenerationHist)
    if(allocated(iPEADFSweepOutcomeHist)) deallocate(iPEADFSweepOutcomeHist)
    if(allocated(iPEADFPendingWitnessHist)) deallocate(iPEADFPendingWitnessHist)
    if(allocated(iPEAPolishAttemptHist)) deallocate(iPEAPolishAttemptHist)
    if(allocated(iPEAPolishAcceptedHist)) deallocate(iPEAPolishAcceptedHist)
    if(allocated(iPEAPolishReasonHist)) deallocate(iPEAPolishReasonHist)
    if(allocated(iPEAPolishIterGlobalHist)) deallocate(iPEAPolishIterGlobalHist)
    if(allocated(iPEAGatePassHist)) deallocate(iPEAGatePassHist)
    if(allocated(iPEAGateBlockHist)) deallocate(iPEAGateBlockHist)
    if(allocated(iPEAGatePotBlockHist)) deallocate(iPEAGatePotBlockHist)
    if(allocated(iPEAGateDFBlockHist)) deallocate(iPEAGateDFBlockHist)
    if(allocated(iPEADirectHandoffHist)) deallocate(iPEADirectHandoffHist)
    if(allocated(iPEADirectPhaseHist)) deallocate(iPEADirectPhaseHist)
    if(allocated(iPEARepeatExitHist)) deallocate(iPEARepeatExitHist)
    if(allocated(iPEARepeatExitReasonHist)) deallocate(iPEARepeatExitReasonHist)
    if(allocated(iPEALagrangianHandoffRepeatHist)) deallocate(iPEALagrangianHandoffRepeatHist)
    if(allocated(iPEALagrangianHandoffFirstIterHist)) deallocate(iPEALagrangianHandoffFirstIterHist)
    if(allocated(iSUBOMTwoSetTraceStage)) deallocate(iSUBOMTwoSetTraceStage)
    if(allocated(iSUBOMTwoSetTraceIterPEA)) deallocate(iSUBOMTwoSetTraceIterPEA)
    if(allocated(iSUBOMTwoSetTraceIterGlobal)) deallocate(iSUBOMTwoSetTraceIterGlobal)
    if(allocated(iSUBOMTwoSetTracePhase)) deallocate(iSUBOMTwoSetTracePhase)
    if(allocated(iSUBOMTwoSetTraceOrdinal)) deallocate(iSUBOMTwoSetTraceOrdinal)
    if(allocated(iSUBOMTwoSetTraceSlot)) deallocate(iSUBOMTwoSetTraceSlot)
    if(allocated(iSUBOMTwoSetStoredPhase)) deallocate(iSUBOMTwoSetStoredPhase)
    if(allocated(iSUBOMTwoSetStoredOrdinal)) deallocate(iSUBOMTwoSetStoredOrdinal)
    if(allocated(iSUBOMOrderingGateIterPEA)) deallocate(iSUBOMOrderingGateIterPEA)
    if(allocated(iSUBOMOrderingGateModeCount)) deallocate(iSUBOMOrderingGateModeCount)
    if(allocated(iSUBOMOrderingGateInfo)) deallocate(iSUBOMOrderingGateInfo)
    if(allocated(iODCandidateClass)) deallocate(iODCandidateClass)
    if(allocated(iODCandidateCompanionPhase)) deallocate(iODCandidateCompanionPhase)
    if(allocated(iPEAAssemblageHist)) deallocate(iPEAAssemblageHist)
    if(allocated(iPEALagrangianHandoffHist)) deallocate(iPEALagrangianHandoffHist)
    if(allocated(dPEAMinPhasePotentialHist)) deallocate(dPEAMinPhasePotentialHist)
    if(allocated(dPEAMaxElementPotentialHist)) deallocate(dPEAMaxElementPotentialHist)
    if(allocated(dPEAGEMFunctionNormHist)) deallocate(dPEAGEMFunctionNormHist)
    if(allocated(dPEAPolishNormHist)) deallocate(dPEAPolishNormHist)
    if(allocated(dPEAPolishPotentialChangeHist)) deallocate(dPEAPolishPotentialChangeHist)
    if(allocated(dPEAElementPotentialHist)) deallocate(dPEAElementPotentialHist)
    if(allocated(dSUBOMTwoSetTraceAmount)) deallocate(dSUBOMTwoSetTraceAmount)
    if(allocated(dSUBOMTwoSetStoredMol)) deallocate(dSUBOMTwoSetStoredMol)
    if(allocated(dSUBOMTwoSetTraceMol)) deallocate(dSUBOMTwoSetTraceMol)
    if(allocated(dSUBOMTwoSetTraceSite)) deallocate(dSUBOMTwoSetTraceSite)
    if(allocated(dSUBOMOrderingGateEigenMin)) deallocate(dSUBOMOrderingGateEigenMin)
    if(allocated(dODCandidateCurrentGibbs)) deallocate(dODCandidateCurrentGibbs)
    if(allocated(dODCandidateDisorderedGibbs)) deallocate(dODCandidateDisorderedGibbs)
    if(allocated(dODCandidateOrderingEigenMin)) deallocate(dODCandidateOrderingEigenMin)
    if(allocated(dODCompanionEigenMin)) then
        deallocate(dODCompanionEigenMin)
    end if
    if(allocated(lSUBOMOrderingGateUnstable)) deallocate(lSUBOMOrderingGateUnstable)
    allocate(lPhaseChangeHistory(iterGlobalMax))
    allocate(iPhaseChangeReasonHistory(iterGlobalMax))
    allocate(iGEMNewtonSolverHist(iterGlobalMax))
    allocate(iGEMNewtonInfoHist(iterGlobalMax))
    allocate(iGEMNewtonKKTSizeHist(iterGlobalMax))
    allocate(iGEMNewtonPivot1x1Hist(iterGlobalMax))
    allocate(iGEMNewtonPivot2x2Hist(iterGlobalMax))
    allocate(iGEMNewtonPivotPositiveHist(iterGlobalMax))
    allocate(iGEMNewtonPivotNegativeHist(iterGlobalMax))
    allocate(iGEMNewtonPivotZeroHist(iterGlobalMax))
    allocate(iGEMLSIterHist(iterGlobalMax))
    allocate(iGEMLSNegFactorHist(iterGlobalMax))
    allocate(iGEMLSNegPhaseHist(iterGlobalMax))
    allocate(iGEMLSFloorHist(iterGlobalMax))
    allocate(iGEMLSNoDescentHist(iterGlobalMax))
    allocate(iGEMLSNoDescentClassHist(iterGlobalMax))
    allocate(iGEMPreLMNoDescentHist(iterGlobalMax))
    allocate(iGEMPreLMNoDescentClassHist(iterGlobalMax))
    allocate(iGEMTraceInactiveHist(iterGlobalMax))
    allocate(iGEMTraceReinjectedHist(iterGlobalMax))
    allocate(iGEMTraceRemoveHist(iterGlobalMax))
    allocate(iGEMTraceReinjectHist(iterGlobalMax))
    allocate(iGEMResidualLMHist(iterGlobalMax))
    allocate(iGEMResidualLMRawNegativeHist(iterGlobalMax))
    allocate(iGEMResidualLMNoDescentHist(iterGlobalMax))
    allocate(iGEMResidualLMNoDescentClassHist(iterGlobalMax))
    allocate(iGEMCoalesceHist(iterGlobalMax))
    allocate(iGEMStabilizeHist(iterGlobalMax))
    allocate(iGEMBoundaryRemovalHist(iterGlobalMax))
    allocate(iGEMRawNegativeRemovalHist(iterGlobalMax))
    allocate(iGEMRawNegativePhaseSlotHist(iterGlobalMax))
    allocate(iGEMRawNegativePhaseSolnHist(iterGlobalMax))
    allocate(iGEMRawNegativePhaseSpeciesHist(iterGlobalMax))
    allocate(iGEMRawNegCompHist(iterGlobalMax))
    allocate(iGEMInvalidCompBoundAttemptedHist(iterGlobalMax))
    allocate(iGEMInvalidCompBoundAcceptedHist(iterGlobalMax))
    allocate(iGEMInvalidCompBoundRejectedHist(iterGlobalMax))
    allocate(iGEMInvalidCompBoundVerdictHist(iterGlobalMax))
    allocate(iGEMBoundaryPinnedRemovalHist(iterGlobalMax))
    allocate(iGEMBoundaryRemovalSlotHist(iterGlobalMax))
    allocate(iGEMBoundaryRemovalPhaseHist(iterGlobalMax))
    allocate(iGEMTinyBoundaryRemovalHist(iterGlobalMax))
    allocate(iGEMTinyBoundaryRemovalSlotHist(iterGlobalMax))
    allocate(iGEMTinyBoundaryRemovalPhaseHist(iterGlobalMax))
    allocate(iGEMBoundaryRankGuardHist(iterGlobalMax))
    allocate(iGEMBoundaryRankGuardSlotHist(iterGlobalMax))
    allocate(iGEMBoundaryRankGuardPhaseHist(iterGlobalMax))
    allocate(iGEMTraceRemoveSpeciesHist(iterGlobalMax))
    allocate(iGEMTraceRemovePhaseHist(iterGlobalMax))
    allocate(iGEMTraceRemoveCountHist(iterGlobalMax))
    allocate(iGEMTraceReinjectSpeciesHist(iterGlobalMax))
    allocate(iGEMTraceReinjectPhaseHist(iterGlobalMax))
    allocate(iGEMTraceReinjectCountHist(iterGlobalMax))
    allocate(iGEMCEFRetryActivationHist(iterGlobalMax))
    allocate(iGEMSUBGQRideAlongHist(iterGlobalMax))
    allocate(iODPairHist(iterGlobalMax))
    allocate(iODOrderingModeHist(iterGlobalMax))
    allocate(dGEMNormHist(iterGlobalMax))
    allocate(dGEMMassNormHist(iterGlobalMax))
    allocate(dGEMChemNormHist(iterGlobalMax))
    allocate(dGEMSolnChemNormHist(iterGlobalMax))
    allocate(dGEMCondChemNormHist(iterGlobalMax))
    allocate(dGEMSublExchangeNormHist(iterGlobalMax))
    allocate(dGEMNewtonSymResidualHist(iterGlobalMax))
    allocate(dGEMNewtonMinPivotScaleHist(iterGlobalMax))
    allocate(dGEMNewtonMaxPivotScaleHist(iterGlobalMax))
    allocate(dGEMNewtonDirectionNormHist(iterGlobalMax))
    allocate(dGEMLSInitialNormHist(iterGlobalMax))
    allocate(dGEMLSBestNormHist(iterGlobalMax))
    allocate(dGEMLSFinalNormHist(iterGlobalMax))
    allocate(dGEMLSInitialStepHist(iterGlobalMax))
    allocate(dGEMLSBestStepHist(iterGlobalMax))
    allocate(dGEMLSFinalStepHist(iterGlobalMax))
    allocate(dGEMLSMinRawPhaseMolesHist(iterGlobalMax))
    allocate(dGEMLSMinFinalPhaseMolesHist(iterGlobalMax))
    allocate(dGEMLSInitialGibbsHist(iterGlobalMax))
    allocate(dGEMLSBestGibbsHist(iterGlobalMax))
    allocate(dGEMLSFinalGibbsHist(iterGlobalMax))
    allocate(dGEMLSInitialMeritHist(iterGlobalMax))
    allocate(dGEMLSBestMeritHist(iterGlobalMax))
    allocate(dGEMLSFinalMeritHist(iterGlobalMax))
    allocate(dGEMPreLMInitNormHist(iterGlobalMax))
    allocate(dGEMPreLMBestNormHist(iterGlobalMax))
    allocate(dGEMPreLMInitGibbsHist(iterGlobalMax))
    allocate(dGEMPreLMBestGibbsHist(iterGlobalMax))
    allocate(dGEMPreLMInitMeritHist(iterGlobalMax))
    allocate(dGEMPreLMBestMeritHist(iterGlobalMax))
    allocate(dGEMPreLMMeritCandNormHist(iterGlobalMax))
    allocate(dGEMPreLMMeritCandMassHist(iterGlobalMax))
    allocate(dGEMPreLMMeritCandStepHist(iterGlobalMax))
    allocate(dGEMPreLMMeritCandGibbsHist(iterGlobalMax))
    allocate(dGEMPreLMMeritCandMeritHist(iterGlobalMax))
    allocate(dGEMNewtonDirNormSlopeHist(iterGlobalMax))
    allocate(dGEMNewtonDirGibbsSlopeHist(iterGlobalMax))
    allocate(dGEMNewtonDirMeritSlopeHist(iterGlobalMax))
    allocate(dGEMRawNegativePhaseAmountHist(iterGlobalMax))
    allocate(dGEMRawNegativePhaseDirectionHist(iterGlobalMax))
    allocate(dGEMRawNegPhaseResidualHist(iterGlobalMax))
    allocate(dGEMInvalidCompBoundPhiHist(iterGlobalMax))
    allocate(dODCompDistHist(iterGlobalMax))
    allocate(dODOrderNormHist(iterGlobalMax))
    allocate(dODOrderingEigenMinHist(iterGlobalMax))
    allocate(dGEMLSRawPhaseMoles(nElements))
    allocate(dGEMLSFinalPhaseMoles(nElements))
    allocate(dGEMLSRawPhaseMolesHist(nElements,iterGlobalMax))
    allocate(dGEMLSFinalPhaseMolesHist(nElements,iterGlobalMax))
    allocate(iPEALevelIterHist(iterPEAMax))
    allocate(iPEAPlaneGenerationHist(iterPEAMax))
    allocate(iPEADFSweepGenerationHist(iterPEAMax))
    allocate(iPEADFSweepOutcomeHist(iterPEAMax))
    allocate(iPEADFPendingWitnessHist(iterPEAMax))
    allocate(iPEAPolishAttemptHist(iterPEAMax))
    allocate(iPEAPolishAcceptedHist(iterPEAMax))
    allocate(iPEAPolishReasonHist(iterPEAMax))
    allocate(iPEAPolishIterGlobalHist(iterPEAMax))
    allocate(iPEAGatePassHist(iterPEAMax))
    allocate(iPEAGateBlockHist(iterPEAMax))
    allocate(iPEAGatePotBlockHist(iterPEAMax))
    allocate(iPEAGateDFBlockHist(iterPEAMax))
    allocate(iPEADirectHandoffHist(iterPEAMax))
    allocate(iPEADirectPhaseHist(iterPEAMax))
    allocate(iPEARepeatExitHist(iterPEAMax))
    allocate(iPEARepeatExitReasonHist(iterPEAMax))
    allocate(iPEALagrangianHandoffRepeatHist(iterPEAMax))
    allocate(iPEALagrangianHandoffFirstIterHist(iterPEAMax))
    allocate(iPEAAssemblageHist(nElements,iterPEAMax))
    allocate(iPEALagrangianHandoffHist(nElements,iterPEAMax))
    allocate(dPEAMinPhasePotentialHist(iterPEAMax))
    allocate(dPEAMaxElementPotentialHist(iterPEAMax))
    allocate(dPEAGEMFunctionNormHist(iterPEAMax))
    allocate(dPEAPolishNormHist(iterPEAMax))
    allocate(dPEAPolishPotentialChangeHist(iterPEAMax))
    allocate(dPEAElementPotentialHist(nElements,iterPEAMax))
    nSUBOMTwoSetTraceCapacity = MAX(256, 2*iterGlobalMax + 12*iterPEAMax)
    allocate(iSUBOMTwoSetTraceStage(nSUBOMTwoSetTraceCapacity))
    allocate(iSUBOMTwoSetTraceIterPEA(nSUBOMTwoSetTraceCapacity))
    allocate(iSUBOMTwoSetTraceIterGlobal(nSUBOMTwoSetTraceCapacity))
    allocate(iSUBOMTwoSetTracePhase(nSUBOMTwoSetTraceCapacity))
    allocate(iSUBOMTwoSetTraceOrdinal(nSUBOMTwoSetTraceCapacity))
    allocate(iSUBOMTwoSetTraceSlot(nSUBOMTwoSetTraceCapacity))
    allocate(iSUBOMTwoSetStoredPhase(MAX(1,nSolnPhasesSys)))
    allocate(iSUBOMTwoSetStoredOrdinal(MAX(1,nSolnPhasesSys)))
    allocate(iSUBOMOrderingGateIterPEA(MAX(1,nSolnPhasesSys)))
    allocate(iSUBOMOrderingGateModeCount(MAX(1,nSolnPhasesSys)))
    allocate(iSUBOMOrderingGateInfo(MAX(1,nSolnPhasesSys)))
    allocate(iODCandidateClass(MAX(1,nSolnPhasesSys)))
    allocate(iODCandidateCompanionPhase(MAX(1,nSolnPhasesSys)))
    allocate(dSUBOMTwoSetTraceAmount(nSUBOMTwoSetTraceCapacity))
    allocate(dSUBOMTwoSetStoredMol(MAX(1,nSolnPhasesSys),nSpecies))
    allocate(dSUBOMTwoSetTraceMol(nSUBOMTwoSetTraceCapacity,nSpecies))
    allocate(dSUBOMTwoSetTraceSite(nSUBOMTwoSetTraceCapacity,nMaxSublatticeSys,nMaxConstituentSys))
    allocate(dSUBOMOrderingGateEigenMin(MAX(1,nSolnPhasesSys)))
    allocate(dODCandidateCurrentGibbs(MAX(1,nSolnPhasesSys)))
    allocate(dODCandidateDisorderedGibbs(MAX(1,nSolnPhasesSys)))
    allocate(dODCandidateOrderingEigenMin(MAX(1,nSolnPhasesSys)))
    allocate(dODCompanionEigenMin(MAX(1,nSolnPhasesSys)))
    allocate(lSUBOMOrderingGateUnstable(MAX(1,nSolnPhasesSys)))
    lPhaseChangeHistory = .FALSE.
    iPhaseChangeReasonHistory = PHASE_CHANGE_REASON_NONE
    iGEMNewtonSolverHist = 0
    iGEMNewtonInfoHist = 0
    iGEMNewtonKKTSizeHist = 0
    iGEMNewtonPivot1x1Hist = 0
    iGEMNewtonPivot2x2Hist = 0
    iGEMNewtonPivotPositiveHist = 0
    iGEMNewtonPivotNegativeHist = 0
    iGEMNewtonPivotZeroHist = 0
    iGEMLSIterHist = 0
    iGEMLSNegFactorHist = 0
    iGEMLSNegPhaseHist = 0
    iGEMLSFloorHist = 0
    iGEMLSNoDescentHist = 0
    iGEMLSNoDescentClassHist = 0
    iGEMPreLMNoDescentHist = 0
    iGEMPreLMNoDescentClassHist = 0
    iGEMTraceInactiveHist = 0
    iGEMTraceReinjectedHist = 0
    iGEMTraceRemoveHist = 0
    iGEMTraceReinjectHist = 0
    iGEMResidualLMHist = 0
    iGEMResidualLMRawNegativeHist = 0
    iGEMResidualLMNoDescentHist = 0
    iGEMResidualLMNoDescentClassHist = 0
    iGEMCoalesceHist = 0
    iGEMStabilizeHist = 0
    iGEMBoundaryRemovalHist = 0
    iGEMRawNegativeRemovalHist = 0
    iGEMRawNegativePhaseSlotHist = 0
    iGEMRawNegativePhaseSolnHist = 0
    iGEMRawNegativePhaseSpeciesHist = 0
    iGEMRawNegCompHist = 0
    iGEMInvalidCompBoundAttemptedHist = 0
    iGEMInvalidCompBoundAcceptedHist = 0
    iGEMInvalidCompBoundRejectedHist = 0
    iGEMInvalidCompBoundVerdictHist = 0
    iGEMBoundaryPinnedRemovalHist = 0
    iGEMBoundaryRemovalSlotHist = 0
    iGEMBoundaryRemovalPhaseHist = 0
    iGEMTinyBoundaryRemovalHist = 0
    iGEMTinyBoundaryRemovalSlotHist = 0
    iGEMTinyBoundaryRemovalPhaseHist = 0
    iGEMBoundaryRankGuardHist = 0
    iGEMBoundaryRankGuardSlotHist = 0
    iGEMBoundaryRankGuardPhaseHist = 0
    iGEMTraceRemoveSpeciesHist = 0
    iGEMTraceRemovePhaseHist = 0
    iGEMTraceRemoveCountHist = 0
    iGEMTraceReinjectSpeciesHist = 0
    iGEMTraceReinjectPhaseHist = 0
    iGEMTraceReinjectCountHist = 0
    iGEMCEFRetryActivationHist = 0
    iGEMSUBGQRideAlongHist = 0
    iODPairHist = 0
    iODOrderingModeHist = 0
    dGEMNormHist = 0D0
    dGEMMassNormHist = 0D0
    dGEMChemNormHist = 0D0
    dGEMSolnChemNormHist = 0D0
    dGEMCondChemNormHist = 0D0
    dGEMSublExchangeNormHist = 0D0
    dGEMNewtonSymResidualHist = 0D0
    dGEMNewtonMinPivotScaleHist = 0D0
    dGEMNewtonMaxPivotScaleHist = 0D0
    dGEMNewtonDirectionNormHist = 0D0
    dGEMLSInitialNormHist = 0D0
    dGEMLSBestNormHist = 0D0
    dGEMLSFinalNormHist = 0D0
    dGEMLSInitialStepHist = 0D0
    dGEMLSBestStepHist = 0D0
    dGEMLSFinalStepHist = 0D0
    dGEMLSMinRawPhaseMolesHist = 0D0
    dGEMLSMinFinalPhaseMolesHist = 0D0
    dGEMLSInitialGibbsHist = 0D0
    dGEMLSBestGibbsHist = 0D0
    dGEMLSFinalGibbsHist = 0D0
    dGEMLSInitialMeritHist = 0D0
    dGEMLSBestMeritHist = 0D0
    dGEMLSFinalMeritHist = 0D0
    dGEMPreLMInitNormHist = 0D0
    dGEMPreLMBestNormHist = 0D0
    dGEMPreLMInitGibbsHist = 0D0
    dGEMPreLMBestGibbsHist = 0D0
    dGEMPreLMInitMeritHist = 0D0
    dGEMPreLMBestMeritHist = 0D0
    dGEMPreLMMeritCandNormHist = 0D0
    dGEMPreLMMeritCandMassHist = 0D0
    dGEMPreLMMeritCandStepHist = 0D0
    dGEMPreLMMeritCandGibbsHist = 0D0
    dGEMPreLMMeritCandMeritHist = 0D0
    dGEMNewtonDirNormSlopeHist = 0D0
    dGEMNewtonDirGibbsSlopeHist = 0D0
    dGEMNewtonDirMeritSlopeHist = 0D0
    dGEMRawNegativePhaseAmountHist = 0D0
    dGEMRawNegativePhaseDirectionHist = 0D0
    dGEMRawNegPhaseResidualHist = 0D0
    dGEMInvalidCompBoundPhiHist = 0D0
    dODCompDistHist = 0D0
    dODOrderNormHist = 0D0
    dODOrderingEigenMinHist = 0D0
    dGEMLSRawPhaseMoles = 0D0
    dGEMLSFinalPhaseMoles = 0D0
    dGEMLSRawPhaseMolesHist = 0D0
    dGEMLSFinalPhaseMolesHist = 0D0
    iPEALevelIterHist = 0
    iPEAPlaneGenerationHist = 0
    iPEADFSweepGenerationHist = 0
    iPEADFSweepOutcomeHist = PEA_DF_SWEEP_NOT_RUN
    iPEADFPendingWitnessHist = 0
    iPEAPolishAttemptHist = 0
    iPEAPolishAcceptedHist = 0
    iPEAPolishReasonHist = PHASE_CHANGE_REASON_NONE
    iPEAPolishIterGlobalHist = 0
    iPEAGatePassHist = 0
    iPEAGateBlockHist = 0
    iPEAGatePotBlockHist = 0
    iPEAGateDFBlockHist = 0
    iPEADirectHandoffHist = 0
    iPEADirectPhaseHist = 0
    iPEARepeatExitHist = 0
    iPEARepeatExitReasonHist = 0
    iPEALagrangianHandoffRepeatHist = 0
    iPEALagrangianHandoffFirstIterHist = 0
    iPEAAssemblageHist = 0
    iPEALagrangianHandoffHist = 0
    dPEAMinPhasePotentialHist = 0D0
    dPEAMaxElementPotentialHist = 0D0
    dPEAGEMFunctionNormHist = 0D0
    dPEAPolishNormHist = 0D0
    dPEAPolishPotentialChangeHist = 0D0
    dPEAElementPotentialHist = 0D0
    nSUBOMTwoSetTrace = 0
    nSUBOMTwoSetStored = 0
    nSUBOMOrderingGateEvaluated = 0
    iSUBOMTwoSetTraceStage = 0
    iSUBOMTwoSetTraceIterPEA = 0
    iSUBOMTwoSetTraceIterGlobal = 0
    iSUBOMTwoSetTracePhase = 0
    iSUBOMTwoSetTraceOrdinal = 0
    iSUBOMTwoSetTraceSlot = 0
    iSUBOMTwoSetStoredPhase = 0
    iSUBOMTwoSetStoredOrdinal = 0
    iSUBOMOrderingGateIterPEA = -1
    iSUBOMOrderingGateModeCount = 0
    iSUBOMOrderingGateInfo = 0
    iODCandidateClass = OD_CANDIDATE_NOT_EVALUATED
    lODCommittedOwnershipApplied = .FALSE.
    iODCandidateCompanionPhase = 0
    dSUBOMTwoSetTraceAmount = 0D0
    dSUBOMTwoSetStoredMol = 0D0
    dSUBOMTwoSetTraceMol = 0D0
    dSUBOMTwoSetTraceSite = 0D0
    dSUBOMOrderingGateEigenMin = 0D0
    dODCandidateCurrentGibbs = 0D0
    dODCandidateDisorderedGibbs = 0D0
    dODCandidateOrderingEigenMin = 0D0
    dODCompanionEigenMin = 0D0
    lSUBOMOrderingGateUnstable = .FALSE.

    !Solution related variables
    if (allocated(dPartialExcessGibbs)) then
        deallocate(dPartialExcessGibbs,dPartialEnthalpyXS,dPartialEntropyXS,dPartialHeatCapacityXS)
    end if
    if (allocated(dEffStoichSolnPhase)) deallocate(dEffStoichSolnPhase)
    if (allocated(dDrivingForceSoln)) deallocate(dDrivingForceSoln)
    if (allocated(iSubMinCandidateStatusSoln)) deallocate(iSubMinCandidateStatusSoln)
    if (allocated(iGEMCertPhase)) deallocate(iGEMCertPhase)
    if (allocated(iGEMCertCopy)) deallocate(iGEMCertCopy)
    if (allocated(iGEMCertLevelRow)) deallocate(iGEMCertLevelRow)
    if (allocated(iGEMCertBasis)) deallocate(iGEMCertBasis)
    if (allocated(iGEMCertNorm)) deallocate(iGEMCertNorm)
    if (allocated(iGEMCertStatus)) deallocate(iGEMCertStatus)
    if (allocated(iGEMCertProof)) deallocate(iGEMCertProof)
    if (allocated(iGEMCertSubMinStatus)) deallocate(iGEMCertSubMinStatus)
    if (allocated(iGEMCertRank)) deallocate(iGEMCertRank)
    if (allocated(iGEMCertFirstSpecies)) deallocate(iGEMCertFirstSpecies)
    if (allocated(iGEMCertLastSpecies)) deallocate(iGEMCertLastSpecies)
    if (allocated(dGEMCertDrivingForce)) deallocate(dGEMCertDrivingForce)
    if (allocated(dGEMCertTiming)) deallocate(dGEMCertTiming)
    if (allocated(dGEMCertPlane)) deallocate(dGEMCertPlane)
    if (allocated(dSumMolFractionSoln)) deallocate(dSumMolFractionSoln)

    nGEMCertCapacity = MAX(1, l*MAX(1,nElements))
    allocate(dUpdateVar(nElements*2),dUpdateVarLast(nElements*2))
    allocate(dGEMAnalyticalSpeciesDirection(nSpecies),lGEMAnalyticalSpeciesDirection(l))
    allocate(iterHistory(nElements,iterGlobalMax))
    allocate(&
        dEffStoichSolnPhase(l,nElements),&
        dSumMolFractionSoln(l),&
        dDrivingForceSoln(l),&
        iSubMinCandidateStatusSoln(l),&
        dPartialExcessGibbs(nSpecies),&
        dPartialEnthalpyXS(nSpecies),&
        dPartialEntropyXS(nSpecies),&
        dPartialHeatCapacityXS(nSpecies))
    allocate(&
        iGEMCertPhase(nGEMCertCapacity),&
        iGEMCertCopy(nGEMCertCapacity),&
        iGEMCertLevelRow(nGEMCertCapacity),&
        iGEMCertBasis(nGEMCertCapacity),&
        iGEMCertNorm(nGEMCertCapacity),&
        iGEMCertStatus(nGEMCertCapacity),&
        iGEMCertProof(nGEMCertCapacity),&
        iGEMCertSubMinStatus(nGEMCertCapacity),&
        iGEMCertRank(nGEMCertCapacity),&
        iGEMCertFirstSpecies(nGEMCertCapacity),&
        iGEMCertLastSpecies(nGEMCertCapacity),&
        dGEMCertDrivingForce(nGEMCertCapacity),&
        dGEMCertTiming(nGEMCertCapacity),&
        dGEMCertPlane(nGEMCertCapacity,nElements))
    lConverged              = .FALSE.
    lRevertSystem           = .FALSE.
    iPhaseChangeReason      = PHASE_CHANGE_REASON_NONE
    iGEMExitStatus          = GEM_EXIT_STATUS_OK
    iGEMNewtonSolver        = 0
    iGEMNewtonDSYSVInfo     = 0
    iGEMNewtonKKTSize       = 0
    iGEMMaxSublExchangeSlot = 0
    iGEMMaxSublExchangePhase = 0
    iGEMMaxSublExchangeSite = 0
    iGEMMaxSublExchangeConstituent = 0
    dGEMMaxSublExchangeResidual = 0D0
    dGEMMaxSublExchangeWeightedResidual = 0D0
    iGEMNewtonPivot1x1Count = 0
    iGEMNewtonPivot2x2Count = 0
    iGEMNewtonPivotPositiveCount = 0
    iGEMNewtonPivotNegativeCount = 0
    iGEMNewtonPivotZeroCount = 0
    iGEMInertiaRegularizationAttemptedUsed = 0
    iGEMInertiaRegularizationAcceptedUsed = 0
    iGEMInertiaRegularizationFailedUsed = 0
    iGEMInertiaRegularizationAttemptedTotal = 0
    iGEMInertiaRegularizationAcceptedTotal = 0
    iGEMInertiaRegularizationFailedTotal = 0
    iGEMInertiaRegularizationStepTotal = 0
    iGEMInertiaRegularizationStepLast = 0
    iODPairCount = 0
    iODOrderingModeCount = 0
    nLevel2LagrangeTwoSetCreated = 0
    dGEMNewtonSymmetryResidual = 0D0
    dGEMNewtonMinPivotScale = 0D0
    dGEMNewtonMaxPivotScale = 0D0
    dGEMNewtonDirectionNorm = 0D0
    dODCompDist = 0D0
    dODOrderNorm = 0D0
    dODOrderingEigenMin = 0D0
    iGEMLineSearchIterationCount = 0
    iGEMLineSearchNegativeFactorCount = 0
    iGEMLineSearchFloorCount = 0
    iGEMLineSearchNoDescent = 0
    iGEMLineSearchNoDescentClass = 0
    iGEMPreLMNoDescent = 0
    iGEMPreLMNoDescentClass = 0
    iGEMPreLMNewtonInfo = 0
    iGEMPreLMKKTSize = 0
    iGEMPreLMPivot1x1 = 0
    iGEMPreLMPivot2x2 = 0
    iGEMPreLMPivotPositive = 0
    iGEMPreLMPivotNegative = 0
    iGEMPreLMPivotZero = 0
    iGEMPreLMODOrdPhase = 0
    iGEMPreLMODCompPhase = 0
    dGEMLineSearchInitialNorm = 0D0
    dGEMLineSearchBestNorm = 0D0
    dGEMLineSearchFinalNorm = 0D0
    dGEMLineSearchInitialStep = 0D0
    dGEMLineSearchBestStep = 0D0
    dGEMLineSearchFinalStep = 0D0
    dGEMLineSearchInitialGibbs = 0D0
    dGEMLineSearchBestGibbs = 0D0
    dGEMLineSearchFinalGibbs = 0D0
    dGEMLineSearchInitialMerit = 0D0
    dGEMLineSearchBestMerit = 0D0
    dGEMLineSearchFinalMerit = 0D0
    dGEMLineSearchMeritCandNorm = 0D0
    dGEMLineSearchMeritCandMass = 0D0
    dGEMLineSearchMeritCandStep = 0D0
    dGEMLineSearchMeritCandGibbs = 0D0
    dGEMLineSearchMeritCandMerit = 0D0
    dGEMPreLMInitNorm = 0D0
    dGEMPreLMBestNorm = 0D0
    dGEMPreLMInitGibbs = 0D0
    dGEMPreLMBestGibbs = 0D0
    dGEMPreLMInitMerit = 0D0
    dGEMPreLMBestMerit = 0D0
    dGEMPreLMMeritCandNorm = 0D0
    dGEMPreLMMeritCandMass = 0D0
    dGEMPreLMMeritCandStep = 0D0
    dGEMPreLMMeritCandGibbs = 0D0
    dGEMPreLMMeritCandMerit = 0D0
    dGEMPreLMFinalNorm = 0D0
    dGEMPreLMFinalGibbs = 0D0
    dGEMPreLMFinalMerit = 0D0
    dGEMPreLMDirNormSlope = 0D0
    dGEMPreLMDirGibbsSlope = 0D0
    dGEMPreLMDirMeritSlope = 0D0
    dGEMPreLMDirectionNorm = 0D0
    dGEMPreLMMinPivotScale = 0D0
    dGEMPreLMMaxPivotScale = 0D0
    dGEMPreLMODAlign = 0D0
    dGEMPreLMODEigen = 0D0
    dGEMNewtonDirNormSlope = 0D0
    dGEMNewtonDirGibbsSlope = 0D0
    dGEMNewtonDirMeritSlope = 0D0
    dMaxPotentialTol        = 1D-5
    dPEATol                 = 1D-8
    dToleranceLevel         = -dPEATol
    dUpdateVarLast       = 0D0
    dGEMAnalyticalSpeciesDirection = 0D0
    lGEMAnalyticalSpeciesDirection = .FALSE.
    dPartialExcessGibbs=0D0
    dPartialEnthalpyXS=0D0
    dPartialEntropyXS=0D0
    dPartialHeatCapacityXS=0D0
    dGEMFunctionNorm     = 10D0
    dGEMMassBalanceNorm = 0D0
    dGEMChemicalPotentialNorm = 0D0
    dGEMSolutionChemicalPotentialNorm = 0D0
    dGEMCondensedChemicalPotentialNorm = 0D0
    dSublatticeExchangeNorm = 0D0
    lGEMCEFSiteLagrangianEnabled = .TRUE.
    lGEMCEFSiteLagrangianActive = .FALSE.
    lGEMCEFSiteDirectionActive = .FALSE.
    nGEMCEFPhaseVariables = 0
    nGEMCEFSiteVariables = 0
    iGEMRawNegativePhaseSlot = 0
    iGEMRawNegativePhaseSoln = 0
    iGEMRawNegativePhaseSpecies = 0
    iGEMRawNegComp = 0
    iGEMCEFBndPhaseSlot = 0
    iGEMBoundaryRemovalSlot = 0
    iGEMBoundaryRemovalPhase = 0
    iGEMTinyBoundaryRemovalSlot = 0
    iGEMTinyBoundaryRemovalPhase = 0
    iGEMBoundaryRankGuardUsed = 0
    iGEMBoundaryRankGuardSlot = 0
    iGEMBoundaryRankGuardPhase = 0
    iGEMTraceRemoveSpecies = 0
    iGEMTraceRemovePhase = 0
    iGEMTraceRemoveCount = 0
    iGEMTraceReinjectSpecies = 0
    iGEMTraceReinjectPhase = 0
    iGEMTraceReinjectCount = 0
    dGEMRawNegativePhaseAmount = 0D0
    dGEMRawNegativePhaseDirection = 0D0
    dGEMRawNegPhaseResidual = 0D0
    dGEMInvalidCompBoundPhi = 0D0
    dGEMCEFBndPhaseStep = 0D0
    nPEALagrangianPolishAttempt = 0
    nPEALagrangianPolishAccepted = 0
    nPEALagrangianPolishRejected = 0
    iPEALagrangianPolishIterGlobal = 0
    iPEALagrangianPolishReason = PHASE_CHANGE_REASON_NONE
    iPEALagrangianHandoffRepeat = 0
    iPEALagrangianHandoffFirstIter = 0
    nPEALagrangianHandoffRepeated = 0
    iPEAGatePass = 0
    iPEAGateBlock = 0
    iPEAGatePotBlock = 0
    iPEAGateDFBlock = 0
    nPEAGatePassed = 0
    nPEAGateBlocked = 0
    nPEAGatePotBlocked = 0
    nPEAGateDFBlocked = 0
    iPEADirectHandoff = 0
    iPEADirectPhase = 0
    nPEADirectHandoff = 0
    iPEARepeatExit = 0
    iPEARepeatExitReason = 0
    nPEARepeatExit = 0
    nPEARepeatGuard = 0
    iPEAExitStatus = PEA_EXIT_STATUS_UNKNOWN
    iPEAExitReason = PEA_EXIT_REASON_NONE
    iPEAExitFreshMinPointSweep = 0
    iPEAUncertifiedHandoffSeen = 0
    dPEAExitMinPhasePotential = 0D0
    dPEAExitTolerance = 0D0
    nPEARecorded = 0
    nSUBOMTwoSetTrace = 0
    nSUBOMTwoSetStored = 0
    iGEMAnalyticalHessianFallbackCount = 0
    iGEMCEFResidualLMFallbackCount = 0
    iGEMResidualLMUsed = 0
    iGEMResidualLMRawNegativeUsed = 0
    iGEMResidualLMNoDescentUsed = 0
    iGEMResidualLMNoDescentClass = 0
    iGEMResidualLMTotalUsed = 0
    iGEMResidualLMRawNegativeTotalUsed = 0
    iGEMResidualLMNoDescentTotalUsed = 0
    iGEMLagrangianCallSite = 1
    nGEMResidualLMEvent = 0
    nGEMResidualLMEventCapacity = 0
    iGEMTrial7BoundRetryAttemptedTotal = 0
    iGEMTrial7BoundRetryAcceptedTotal = 0
    iGEMInvalidCompBoundAttemptedUsed = 0
    iGEMInvalidCompBoundAcceptedUsed = 0
    iGEMInvalidCompBoundRejectedUsed = 0
    iGEMInvalidCompBoundVerdict = 0
    iGEMInvalidCompBoundAttemptedTotal = 0
    iGEMInvalidCompBoundAcceptedTotal = 0
    iGEMInvalidCompBoundRejectedTotal = 0
    iGEMCoalesceUsed = 0
    iGEMDuplicateSUBOMCoalesceTotalUsed = 0
    iGEMDegenerateODCoalesceTotalUsed = 0
    iGEMStabilizeUsed = 0
    iGEMStabilizeTotalUsed = 0
    iGEMBoundaryRemovalUsed = 0
    iGEMRawNegativeRemovalUsed = 0
    iGEMBoundaryPinnedRemovalUsed = 0
    iGEMTinyBoundaryRemovalUsed = 0
    iGEMTraceRemoveUsed = 0
    iGEMTraceReinjectUsed = 0
    iGEMCEFRetryActivationUsed = 0
    iGEMCEFRetryActivationTotalUsed = 0
    iGEMSUBGQRideAlongUsed = 0
    dPEALagrangianPolishNormBefore = 0D0
    dPEALagrangianPolishNormAfter = 0D0
    dPEALagrangianPolishPotentialChange = 0D0
    lPEALagrangianPolishEnabled = .TRUE.
    lPEALagrangianPolishActive = .FALSE.
    lPEALagrangianPolishAccepted = .FALSE.
    lGEMCEFBndPhaseActive = .FALSE.
    lGEMCEFInertiaRegularizationEnabled = .FALSE.
    lGEMCEFInertiaRegularizationActive = .FALSE.

    dDrivingForceSoln =0d0
    iSubMinCandidateStatusSoln = SUBMIN_CANDIDATE_UNKNOWN
    nGEMCertCount = 0
    nGEMCertEmissionCount = 0
    nGEMCertDropped = 0
    iGEMCertPhase = 0
    iGEMCertCopy = 0
    iGEMCertLevelRow = 0
    iGEMCertBasis = 0
    iGEMCertNorm = GEM_CERT_NORMALIZATION_PER_MOLE_ATOMS
    iGEMCertStatus = GEM_CERT_STATUS_UNKNOWN
    iGEMCertProof = GEM_CERT_PROOF_UNKNOWN
    iGEMCertSubMinStatus = SUBMIN_CANDIDATE_UNKNOWN
    iGEMCertRank = 0
    iGEMCertFirstSpecies = 0
    iGEMCertLastSpecies = 0
    dGEMCertDrivingForce = 9D5
    dGEMCertTiming = 0D0
    dGEMCertPlane = 0D0

    call InitTraceSpeciesMask
    call InitPEADFSweepStorage

    return
end subroutine InitGEMSolver
!
!
