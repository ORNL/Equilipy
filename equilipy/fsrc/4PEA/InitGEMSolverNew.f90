!
subroutine InitGEMSolverNew
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    implicit none
    integer::   i, j, k, l, m, n, iter, INFO
    integer,dimension(nElements)           :: IPIV
    real(8),dimension(nElements,nElements) :: A
    ! Initialize variables:
    IPIV   = 0
    A      = 0D0
    nConPhases        = nElements
    nSpeciesLevel     = nSpecies + nSolnPhasesSys
    l = MAX(1,nSolnPhasesSys)

    ! print*,' '
    ! print*,'initUBC',iterGlobal,dGEMFunctionNorm,dElementPotential
    ! print*,'dMolesPhase',dMolesPhase
!
    ! Check to see if allocatable arrays have already been allocated:
    if (allocated(dMolesPhase))             deallocate(dMolesPhase)
    if (allocated(dElementPotentialLast))       deallocate(dElementPotentialLast)
    if (allocated(iAssemblage))             deallocate(iAssemblage)
    if (allocated(iterHistoryLevel))        deallocate(iterHistoryLevel)
    if (allocated(dPhasePotential))         deallocate(dPhasePotential)
    if (allocated(dAtomFractionSpeciesOld)) deallocate(dAtomFractionSpeciesOld)
    if (allocated(dChemicalPotentialOld))   deallocate(dChemicalPotentialOld)
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
        dChemicalPotentialOld(nSpecies),&
        iterHistoryLevel(nElements,500),&
        dStoichSpeciesLevel(nSpeciesLevel,nElements),&
        iPhaseLevel(nSpeciesLevel),&
        iShuffled(nElements),&
        dMolesPhaseHistory(nElements,500),&
        iCandidate(nSpecies))
    allocate(dEffStoichSolnPhase(l,nElements))
    
    ! Initialize allocatable variables:
    dElementPotentialLast=dElementPotential
    iAssemblage       = iAssemblageGEM
    iShuffled         = (/(i, i = 1,nElements, 1)/)
    iCandidate        = 0
    dAtomFractionSpeciesOld = dAtomFractionSpecies
    dChemicalPotentialOld   =(dStdGibbsEnergy) &
    *DFLOAT(iParticlesPerMole)/ dSpeciesTotalAtoms
    dMolesPhase       = 0D0
    dEffStoichSolnPhase=0D0
    dPEATol = 1D-7
    
    iterPEA = 0
!
    
    !SY: Modifying dAtomFractionSpecies, dChemicalPotential, dPhasePotential
    !These needs to be modified back to the size of nSpecies
    if (allocated(dAtomFractionSpecies))     deallocate(dAtomFractionSpecies)
    if (allocated(dChemicalPotential))       deallocate(dChemicalPotential)
    allocate(dAtomFractionSpecies(nSpeciesLevel,nElements),&
    dChemicalPotential(nSpeciesLevel),dPhasePotential(nSpeciesLevel))
!

    dChemicalPotential(:nSpecies)     = dChemicalPotentialOld

    dAtomFractionSpecies(:nSpecies,:) = dAtomFractionSpeciesOld
    dStoichSpeciesLevel(:nSpecies,:)  = dStoichSpecies
    iPhaseLevel(:nSpecies)            = iPhase

    ! !========================================================================================================================

    if(iterGlobal>1) then
        ! Case 1: No rigorous phase assemblage check & Strict PEA tolerance
        if(lPostProcess) lPostProcess =.False.
        dPEATol = 1D-10

        ! Goal
        ! Debug: Consider changing this part from the scratch. Take first assemblage and level through them based on the LG element potential.
        ! 1. Update dStoichSpeciesGEM, dMolFractionGEM, dAtomFractionGEM, iAssemblageGEM from LG assemblages and mole fraction
        ! 2. Update dMolesPhase based on the first assemblage
        ! 3. Update ChemicalPotentialGEM based on the first assemblage
        ! 4. Calculate CompMinSolnPoint based on the LG element potential
        ! 5. Update element potential based on the first assemblage
        ! 5. Move on to leveling

        ! 1. Calculate Gibbs energy of solution as stoichimatric phases using LG element potential and driving force
        call CompMinSolnPoint

        ! 2. Initialize leveling condition to the first assemblage
        dStoichSpeciesGEM       = dStoichSpeciesGEMinit
        dChemicalPotentialGEM   = dChemicalPotentialGEMinit
        iPhaseGEM               = iPhaseGEMinit
        iAssemblage             = iAssemblageGEMinit

        ! 3. Update dMolesPhase based on the first assemblage
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
        
        
        ! 4. Update elemental potential based on the first assemblage
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
            dPEATol = 1D-10
            lPostProcess =.False.
        else
            ! Case 3: Rigorous phase assemblage check & Loose PEA tolerance
            dPEATol = 1D-7
        end if

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

        ! Update ChemicalPotentialGEM
        do i = 1, nElements
            j=iAssemblage(i)
            if (j<=nSpecies) then
                !Update Sthoichometric compounds
                dChemicalPotentialGEM(i)=dChemicalPotentialOld(j)
            else
                l = j - nSpecies
                m = nSpeciesPhase(l-1) + 1      ! First constituent in phase.
                n = nSpeciesPhase(l)            ! Last  constituent in phase.
                dMolFraction(m:n) = dMolFractionGEM(i,m:n)
                call CompExcessGibbsEnergy(l)
                dChemicalPotentialGEM(i) = dot_product(dMolFraction(m:n),dChemicalPotential(m:n))/&
                sum(matmul(dMolFraction(m:n)/DFLOAT(iParticlesPerMole(m:n)), dStoichSpecies(m:n,:)))
            end if
        end do

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
        call CompMinSolnPoint

    end if

    !========================================================================================================================
    return
end subroutine InitGEMSolverNew
!
!
