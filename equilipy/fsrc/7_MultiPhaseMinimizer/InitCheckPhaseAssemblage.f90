!> Initialize PEA work arrays from Leveling or Lagrangian active-set state.
subroutine InitCheckPhaseAssemblage
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    InitCheckPhaseAssemblage.f90
    !> \brief   Initialize CheckPhaseAssemblage/PEA arrays from the current active set.
    !> \author  S.Y. Kwon
    !> \date    Jun. 25, 2026
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Integrated immutable grid rows into PEA initialization and used selected refined starts without grid-mode initialization multistart.
    !
    ! Purpose:
    ! ========
    !
    !> \details This routine prepares the PEA candidate arrays used by
    !! CheckPhaseAssemblage.  A Lagrangian solution assemblage may contain fewer
    !! active phases than system components; in that case empty rows are left for
    !! GetNewAssemblage to fill and the existing Lagrangian elemental potentials
    !! are retained.
    !
    !-------------------------------------------------------------------------------------------------------------
!
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE GridDiscovery, ONLY: RestoreStaticGridLevelingRows, RegisterSelectedGridSeeds, RefineSelectedGridSeeds
    implicit none
    integer::   i, j, k, l, m, n, iter, INFO
    integer,dimension(nElements)           :: IPIV
    integer,dimension(nElements)           :: iAssemblageLeveling
    real(8),dimension(nElements)           :: dMolesPhaseInput
    real(8),dimension(nElements,nElements) :: A
    real(8)                                :: dLevelingDenom, dLevelingFormulaAtomDenom
    logical                                :: lHasEmptyAssemblageSlot
    ! Step 1. Size the dynamic candidate rows after any immutable grid rows.
    IPIV   = 0
    A      = 0D0
    nConPhases        = nElements
    nSpeciesLevel     = nSpecies + nSolnPhasesSys
    if (lSUBOMTwoSetCandidateEnabled.OR.lODPartitionUnifiedActive) then
        nSpeciesLevel = nSpeciesLevel + nSolnPhasesSys
    end if
    if (lGridFrontEndActive) then
        nSpeciesLevel = nSpecies + nGridPoint + 2*nSolnPhasesSys + nElements
    end if
    l = MAX(1,nSolnPhasesSys)

    ! Step 2. Preserve the current Leveling assemblage before reallocating
    ! the PEA work arrays.  On the first PEA call, iAssemblageGEM has not been
    ! created by Level2Lagrange and must not be used as an implicit dependency.
    iAssemblageLeveling = 0
    if (allocated(iAssemblage)) iAssemblageLeveling = iAssemblage
    dMolesPhaseInput = 0D0
    if (allocated(dMolesPhase)) dMolesPhaseInput = dMolesPhase
!
    ! Step 3. Reallocate the PEA and candidate arrays for the current row table.
    if (allocated(dMolesPhase))             deallocate(dMolesPhase)
    if (allocated(dElementPotentialLast))       deallocate(dElementPotentialLast)
    if (allocated(iAssemblage))             deallocate(iAssemblage)
    if (allocated(iterHistoryLevel))        deallocate(iterHistoryLevel)
    if (allocated(dPhasePotential))         deallocate(dPhasePotential)
    if (allocated(dAtomFractionSpeciesOld)) deallocate(dAtomFractionSpeciesOld)
    if (allocated(dLevelingCompositionSpeciesOld)) deallocate(dLevelingCompositionSpeciesOld)
    if (allocated(dChemicalPotentialOld))   deallocate(dChemicalPotentialOld)
    if (allocated(dLevelingChemicalPotentialOld))   deallocate(dLevelingChemicalPotentialOld)
    if (allocated(dStoichSpeciesLevel))     deallocate(dStoichSpeciesLevel)
    if (allocated(iPhaseLevel))             deallocate(iPhaseLevel)
    if (allocated(iShuffled))               deallocate(iShuffled)
    if (allocated(dMolesPhaseHistory))      deallocate(dMolesPhaseHistory)
    if (allocated(dEffStoichSolnPhase))     deallocate(dEffStoichSolnPhase)
    if (allocated(iCandidate))              deallocate(iCandidate)
    allocate(&
        dMolesPhase(nElements),&
        dElementPotentialLast(nElements),&
        iAssemblage(nElements),&
        dAtomFractionSpeciesOld(nSpecies,nElements),&
        dLevelingCompositionSpeciesOld(nSpecies,nElements),&
        dChemicalPotentialOld(nSpecies),&
        dLevelingChemicalPotentialOld(nSpecies),&
        iterHistoryLevel(nElements,500),&
        dStoichSpeciesLevel(nSpeciesLevel,nElements),&
        iPhaseLevel(nSpeciesLevel),&
        iShuffled(nElements),&
        dMolesPhaseHistory(nElements,500),&
        iCandidate(nSpecies))
    allocate(dEffStoichSolnPhase(l,nElements))
    
    ! Step 4. Restore the incoming active set and initialize fresh PEA diagnostics.
    dElementPotentialLast=dElementPotential
    if (allocated(iAssemblageGEM)) then
        iAssemblage = iAssemblageGEM
    else
        iAssemblage = iAssemblageLeveling
    end if
    lHasEmptyAssemblageSlot = ANY(iAssemblage == 0)
    iShuffled         = (/(i, i = 1,nElements, 1)/)
    iCandidate        = 0
    dAtomFractionSpeciesOld = dAtomFractionSpecies(:nSpecies,:)
    dLevelingCompositionSpeciesOld = dLevelingCompositionSpecies(:nSpecies,:)
    dChemicalPotentialOld   =(dStdGibbsEnergy) &
    *DFLOAT(iParticlesPerMole)/ dSpeciesTotalAtoms
    dLevelingChemicalPotentialOld = dLevelingChemicalPotential(:nSpecies)
    dMolesPhase       = 0D0
    dEffStoichSolnPhase=0D0
    dPEATol = 1D-8
    call InitLevelCandidatePool(0)
    if (lGridFrontEndActive) nLevelCandidateRowOffset = nGridPoint
    
    iterPEA = 0
    nPEARecorded = 0
    iPEALagrangianHandoffRepeat = 0
    iPEALagrangianHandoffFirstIter = 0
    nPEALagrangianHandoffRepeated = 0
    iPEAExitStatus = PEA_EXIT_STATUS_UNKNOWN
    iPEAExitReason = PEA_EXIT_REASON_NONE
    iPEAExitFreshMinPointSweep = 0
    dPEAExitMinPhasePotential = 0D0
    dPEAExitTolerance = 0D0
    if (allocated(iPEALevelIterHist)) iPEALevelIterHist = 0
    if (allocated(iPEAPlaneGenerationHist)) iPEAPlaneGenerationHist = 0
    if (allocated(iPEADFSweepGenerationHist)) iPEADFSweepGenerationHist = 0
    if (allocated(iPEADFSweepOutcomeHist)) iPEADFSweepOutcomeHist = PEA_DF_SWEEP_NOT_RUN
    if (allocated(iPEADFPendingWitnessHist)) iPEADFPendingWitnessHist = 0
    if (allocated(iPEAPolishAttemptHist)) iPEAPolishAttemptHist = 0
    if (allocated(iPEAPolishAcceptedHist)) iPEAPolishAcceptedHist = 0
    if (allocated(iPEAPolishReasonHist)) iPEAPolishReasonHist = PHASE_CHANGE_REASON_NONE
    if (allocated(iPEAPolishIterGlobalHist)) iPEAPolishIterGlobalHist = 0
    if (allocated(iPEALagrangianHandoffRepeatHist)) iPEALagrangianHandoffRepeatHist = 0
    if (allocated(iPEALagrangianHandoffFirstIterHist)) iPEALagrangianHandoffFirstIterHist = 0
    if (allocated(iPEAAssemblageHist)) iPEAAssemblageHist = 0
    if (allocated(iPEALagrangianHandoffHist)) iPEALagrangianHandoffHist = 0
    if (allocated(dPEAMinPhasePotentialHist)) dPEAMinPhasePotentialHist = 0D0
    if (allocated(dPEAMaxElementPotentialHist)) dPEAMaxElementPotentialHist = 0D0
    if (allocated(dPEAGEMFunctionNormHist)) dPEAGEMFunctionNormHist = 0D0
    if (allocated(dPEAPolishNormHist)) dPEAPolishNormHist = 0D0
    if (allocated(dPEAPolishPotentialChangeHist)) dPEAPolishPotentialChangeHist = 0D0
    if (allocated(dPEAElementPotentialHist)) dPEAElementPotentialHist = 0D0
!
    
    ! Step 5. Rebuild the physical-species rows, followed by immutable grid rows.
    if (allocated(dAtomFractionSpecies))     deallocate(dAtomFractionSpecies)
    if (allocated(dLevelingCompositionSpecies)) deallocate(dLevelingCompositionSpecies)
    if (allocated(dChemicalPotential))       deallocate(dChemicalPotential)
    if (allocated(dLevelingChemicalPotential)) deallocate(dLevelingChemicalPotential)
    allocate(dAtomFractionSpecies(nSpeciesLevel,nElements),&
    dLevelingCompositionSpecies(nSpeciesLevel,nElements),&
    dChemicalPotential(nSpeciesLevel),dLevelingChemicalPotential(nSpeciesLevel),&
    dPhasePotential(nSpeciesLevel))
!

    dAtomFractionSpecies      = 0D0
    dLevelingCompositionSpecies = 0D0
    dStoichSpeciesLevel       = 0D0
    dChemicalPotential        = 5D9
    dLevelingChemicalPotential = 5D9
    dPhasePotential           = 5D9
    iPhaseLevel               = 0

    dChemicalPotential(:nSpecies)     = dChemicalPotentialOld
    dLevelingChemicalPotential(:nSpecies) = dLevelingChemicalPotentialOld(:nSpecies)

    dAtomFractionSpecies(:nSpecies,:) = dAtomFractionSpeciesOld(:nSpecies,:)
    dLevelingCompositionSpecies(:nSpecies,:) = &
        dLevelingCompositionSpeciesOld(:nSpecies,:)
    dStoichSpeciesLevel(:nSpecies,:)  = dStoichSpecies
    iPhaseLevel(:nSpecies)            = iPhase
    if (lGridFrontEndActive) then
        call RestoreStaticGridLevelingRows
        if (INFOThermo /= 0) return
        call RegisterSelectedGridSeeds
    end if

    ! Step 6. Reconstruct the current active thermodynamic rows and refine all
    ! selected solution candidates from their own constitutions.
    if((iterGlobal>1).AND.(.NOT.lGridRecoveryPassActive)) then
        ! Case 1: No rigorous phase assemblage check & Strict PEA tolerance
        if(lPostProcess) lPostProcess =.False.
        dPEATol = 1D-8

        ! Rebuild the PEA candidate rows from the current Lagrangian elemental
        ! potentials before restarting the Leveling search from the stored
        ! initial assemblage.
        if (lGridFrontEndActive) then
            call RefineSelectedGridSeeds
        else
            call CompInitMinSolnPoint
        end if

        dStoichSpeciesGEM       = dStoichSpeciesGEMinit
        dChemicalPotentialGEM   = dChemicalPotentialGEMinit
        iPhaseGEM               = iPhaseGEMinit
        iAssemblage             = iAssemblageGEMinit

        ! Update dMolesPhase based on the first assemblage.
        do i = 1, nElements
            do j = 1, nElements
                A(j,i) = dStoichSpeciesGEM(i,j)
            end do
            dMolesPhase(i) = dMolesElement(i)
        end do
        ! Call the linear equation solver to compute molar quantities of the phase assemblage:
        INFO = 0
        IPIV = 0
        call DGESV( nElements, 1, A, nElements, IPIV, dMolesPhase, nElements, INFO )
        if (INFO /= 0) INFOThermo = 10
        
        
        ! Update elemental potentials based on the first assemblage.
        dAtomFractionSpeciesGEM = dAtomFractionSpeciesGEMinit
        dMolFractionGEM = dMolFractionGEMinit
        
    
        ! Re-Calculate elemental potentiall based on new assemblage
        do j = 1,nElements
            do i = 1,nElements
                A(i,j) = dAtomFractionSpeciesGEM(i,j)
            end do
            dElementPotential(j) = dChemicalPotentialGEM(j)
        end do
    
        ! Reinitialize variables:
        INFO = 0
        IPIV = 0
    
        ! Call linear equation solver to solve the adjustments applied to the Gibbs Plane:
        call DGESV( nElements, 1, A, nElements, IPIV, dElementPotential, nElements, INFO )
        if (INFO /= 0) INFOThermo = 10


        
    else
        ! This part updates dMolesPhase based on GEM assemblage:
        ! We use elemental potentials from Leveling (or previous calculation) to rigorously testing phase assemblage

        ! Define PEA tolerance
        if (lPostProcess) then
            ! Case 2: Rigorous phase assemblage check & Strict PEA tolerance
            dPEATol = 1D-8
            lPostProcess =.False.
        else
            ! Case 3: Rigorous phase assemblage check & Loose PEA tolerance
            dPEATol = 1D-8
        end if

        if (lHasEmptyAssemblageSlot) then
            dMolesPhase = dMolesPhaseInput
        else
            ! Calculate dMolesPhase (phase amount)
            do i = 1, nElements
                do j = 1, nElements
                    A(j,i) = dStoichSpeciesGEM(i,j)
                end do
                dMolesPhase(i) = dMolesElement(i)
            end do

            INFO = 0
            IPIV = 0
            call DGESV( nElements, 1, A, nElements, IPIV, dMolesPhase, nElements, INFO )
            if (INFO /= 0) INFOThermo = 10
        end if

        ! Update ChemicalPotentialGEM
        do i = 1, nElements
            j=iAssemblage(i)
            if (j == 0) then
                dChemicalPotentialGEM(i) = 0D0
                dAtomFractionSpeciesGEM(i,:) = 0D0
            else if ((j<=nSpecies).AND.&
                ((.NOT.lGridFrontEndActive).OR.(iPhaseGEM(i) <= 0))) then
                !Update Sthoichometric compounds
                dChemicalPotentialGEM(i)=dLevelingChemicalPotentialOld(j)
            else
                l = iPhaseGEM(i)
                if (l <= 0) l = j - nSpecies - nLevelCandidateRowOffset
                if ((l < 1).OR.(l > nSolnPhasesSys)) then
                    INFOThermo = 42
                    return
                end if
                m = nSpeciesPhase(l-1) + 1      ! First constituent in phase.
                n = nSpeciesPhase(l)            ! Last  constituent in phase.
                dMolFraction(m:n) = dMolFractionGEM(i,m:n)
                call CompExcessGibbsEnergy(l)
                dLevelingDenom = SUM(dMolFraction(m:n) * &
                    dLevelingSpeciesTotalAtoms(m:n) / DFLOAT(iParticlesPerMole(m:n)))
                dLevelingFormulaAtomDenom = SUM(dMolFraction(m:n) * &
                    dLevelingSpeciesFormulaAtoms(m:n) / DFLOAT(iParticlesPerMole(m:n)))
                if (dLevelingDenom > 1D-300) then
                    dAtomFractionSpeciesGEM(i,:) = &
                        matmul(dMolFraction(m:n)/DFLOAT(iParticlesPerMole(m:n)), dStoichSpecies(m:n,:)) / &
                        dLevelingDenom
                    if (dLevelingFormulaAtomDenom > 1D-300) then
                        dChemicalPotentialGEM(i) = dot_product(dMolFraction(m:n),dChemicalPotential(m:n)) / &
                            dLevelingFormulaAtomDenom
                    else
                        dChemicalPotentialGEM(i) = 5D9
                    end if
                else
                    dAtomFractionSpeciesGEM(i,:) = 0D0
                    dChemicalPotentialGEM(i) = 5D9
                end if
            end if
        end do

        if (.NOT.lHasEmptyAssemblageSlot) then
            ! Refine elemental potential based on assemblage (This is necessary particularly for PostProcess)
            do j = 1,nElements
                do i = 1,nElements
                    A(i,j) = dAtomFractionSpeciesGEM(i,j)
                end do
                dElementPotential(j) = dChemicalPotentialGEM(j)
            end do

            ! Reinitialize variables:
            INFO = 0
            IPIV = 0

            ! Call linear equation solver to solve the adjustments applied to the Gibbs Plane:
            call DGESV( nElements, 1, A, nElements, IPIV, dElementPotential, nElements, INFO )

            ! Calculate phase potential of each phases
            if (INFO /= 0) INFOThermo = 10
        end if
        if (lGridFrontEndActive) then
            call RefineSelectedGridSeeds
        else
            call CompInitMinSolnPoint
        end if

    end if

    ! Initialization classification precedes the first committed PEA plane.
    lODCommittedOwnershipApplied = .FALSE.

    return
end subroutine InitCheckPhaseAssemblage
!
!
