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
    !   06/25/2026      S.Y. Kwon           Preserved under-filled Lagrangian active sets for PEA repair without
    !                                       solving singular zero-slot phase-amount and potential matrices.
    !   06/26/2026      S.Y. Kwon           Reset per-iteration PEA diagnostics for each PEA entry.
!   06/27/2026      S.Y. Kwon           Reset Lagrangian-handoff repeat diagnostics for each PEA entry.
!   07/03/2026      S.Y. Kwon           Reserved switch-gated extra Leveling rows for same-parent SUBOM
!                                       composition-set candidates.
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
    implicit none
    integer::   i, j, k, l, m, n, iter, INFO
    integer,dimension(nElements)           :: IPIV
    integer,dimension(nElements)           :: iAssemblageLeveling
    real(8),dimension(nElements)           :: dMolesPhaseInput
    real(8),dimension(nElements,nElements) :: A
    real(8)                                :: dLevelingDenom, dLevelingFormulaAtomDenom
    logical                                :: lHasEmptyAssemblageSlot
    ! Initialize variables:
    IPIV   = 0
    A      = 0D0
    nConPhases        = nElements
    nSpeciesLevel     = nSpecies + nSolnPhasesSys
    if (lSUBOMTwoSetCandidateEnabled) nSpeciesLevel = nSpeciesLevel + nSolnPhasesSys
    l = MAX(1,nSolnPhasesSys)

    ! Preserve the current Leveling assemblage before this routine reallocates
    ! the PEA work arrays.  On the first PEA call, iAssemblageGEM has not been
    ! created by Level2Lagrange and must not be used as an implicit dependency.
    iAssemblageLeveling = 0
    if (allocated(iAssemblage)) iAssemblageLeveling = iAssemblage
    dMolesPhaseInput = 0D0
    if (allocated(dMolesPhase)) dMolesPhaseInput = dMolesPhase
!
    ! Check to see if allocatable arrays have already been allocated:
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
    
    ! Initialize allocatable variables:
    dElementPotentialLast=dElementPotential
    if (allocated(iAssemblageGEM)) then
        iAssemblage = iAssemblageGEM
    else
        iAssemblage = iAssemblageLeveling
    end if
    lHasEmptyAssemblageSlot = ANY(iAssemblage == 0)
    iShuffled         = (/(i, i = 1,nElements, 1)/)
    iCandidate        = 0
    dAtomFractionSpeciesOld = dAtomFractionSpecies
    dLevelingCompositionSpeciesOld = dLevelingCompositionSpecies
    dChemicalPotentialOld   =(dStdGibbsEnergy) &
    *DFLOAT(iParticlesPerMole)/ dSpeciesTotalAtoms
    dLevelingChemicalPotentialOld = dLevelingChemicalPotential
    dMolesPhase       = 0D0
    dEffStoichSolnPhase=0D0
    dPEATol = 1D-8
    call InitLevelCandidatePool(0)
    
    iterPEA = 0
    nPEARecorded = 0
    iPEALagrangianHandoffRepeat = 0
    iPEALagrangianHandoffFirstIter = 0
    nPEALagrangianHandoffRepeated = 0
    if (allocated(iPEALevelIterHist)) iPEALevelIterHist = 0
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
    
    !SY: Modifying dAtomFractionSpecies, dChemicalPotential, dPhasePotential
    !These needs to be modified back to the size of nSpecies
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
    dLevelingChemicalPotential(:nSpecies) = dLevelingChemicalPotentialOld

    dAtomFractionSpecies(:nSpecies,:) = dAtomFractionSpeciesOld
    dLevelingCompositionSpecies(:nSpecies,:) = dLevelingCompositionSpeciesOld
    dStoichSpeciesLevel(:nSpecies,:)  = dStoichSpecies
    iPhaseLevel(:nSpecies)            = iPhase

    ! !========================================================================================================================

    if(iterGlobal>1) then
        ! Case 1: No rigorous phase assemblage check & Strict PEA tolerance
        if(lPostProcess) lPostProcess =.False.
        dPEATol = 1D-8

        ! Rebuild the PEA candidate rows from the current Lagrangian elemental
        ! potentials before restarting the Leveling search from the stored
        ! initial assemblage.
        call CompInitMinSolnPoint

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
            else if (j<=nSpecies) then
                !Update Sthoichometric compounds
                dChemicalPotentialGEM(i)=dLevelingChemicalPotentialOld(j)
            else
                l = j - nSpecies
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
        call CompInitMinSolnPoint

    end if

    !========================================================================================================================
    return
end subroutine InitCheckPhaseAssemblage
!
!
