!
subroutine InitGEMSolver
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    implicit none
    integer::   i, j, k, l
    ! Initialize variables:

    lPostProcess = .False.
!
    ! Check to see if allocatable arrays have already been allocated:
    if (allocated(dMolesPhase))             deallocate(dMolesPhase)
    if (allocated(dMolesSpecies))           deallocate(dMolesSpecies)
    if (allocated(dElementPotential))       deallocate(dElementPotential)
    if (allocated(iAssemblage))             deallocate(iAssemblage)
    if (allocated(dPhasePotential))         deallocate(dPhasePotential)
    if (allocated(iShuffled))               deallocate(iShuffled)
    allocate(dMolesPhase(nElements),dMolesSpecies(nSpecies),dElementPotential(nElements),&
        iAssemblage(nElements),iShuffled(nElements),dPhasePotential(nSpecies))
    
    ! Initialize allocatable variables:
    dElementPotential = 0D0
    iAssemblage       = 0
    iShuffled         = (/(i, i = 1,nElements, 1)/)
    dMolesPhase       = 0D0
    dMolesSpecies     = 0D0


    ! Initialize variables for GetFirstAssemblage and GetNewAssemblage
    nSpeciesLevel     = nSpecies !For GetFirstAssemblage
    if (allocated(iPhaseLevel))             deallocate(iPhaseLevel)
    if (allocated(dStoichSpeciesLevel))     deallocate(dStoichSpeciesLevel)
    if (allocated(iterHistoryLevel))        deallocate(iterHistoryLevel)
    if (allocated(dMolesPhaseHistory))      deallocate(dMolesPhaseHistory)
    

    if (allocated(dChemicalPotentialGEM))   deallocate(dChemicalPotentialGEM)
    if (allocated(dAtomFractionSpeciesGEM)) deallocate(dAtomFractionSpeciesGEM)
    if (allocated(dStoichSpeciesGEM))       deallocate(dStoichSpeciesGEM)
    if (allocated(dMolFractionGEM))         deallocate(dMolFractionGEM)
    if (allocated(dMolFraction))        deallocate(dMolFraction)
    if (allocated(iPhaseGEM))               deallocate(iPhaseGEM)

    ! Allocate variables for GetFirstAssemblage and GetNewAssemblage
    allocate(iPhaseLevel(nSpeciesLevel),dStoichSpeciesLevel(nSpeciesLevel,nElements),&
    iterHistoryLevel(nElements,500),dMolesPhaseHistory(nElements,500),dChemicalPotentialGEM(nElements),&
    dAtomFractionSpeciesGEM(nElements,nElements),dStoichSpeciesGEM(nElements,nElements),&
    dMolFractionGEM(nElements,nSpecies),iPhaseGEM(nElements),dMolFraction(nSpecies))
    
    iPhaseLevel          = iPhase
    dStoichSpeciesLevel  = dStoichSpecies
    iterHistoryLevel     = 0
    dMolFractionGEM      = 1D0
    dMolFraction = 1D0

    ! Establish the very first phase assemblage:
    call GetFirstAssemblage

    do i =1,nElements
        iPhaseGEM(i)= iPhaseLevel(iAssemblage(i))
    end do

    do j = 1,nElements
        do i = 1,nElements
            dStoichSpeciesGEM(j,i)=dStoichSpeciesLevel(iAssemblage(j),i)
            dAtomFractionSpeciesGEM(j,i)=dAtomFractionSpecies(iAssemblage(j),i)
        end do
        dChemicalPotentialGEM(j) = dChemicalPotential(iAssemblage(j))
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
    allocate(lPhaseChangeHistory(iterGlobalMax))
    lPhaseChangeHistory = .FALSE.

    !Solution related variables
    if (allocated(dPartialExcessGibbs)) then
        deallocate(dPartialExcessGibbs,dPartialEnthalpyXS,dPartialEntropyXS,dPartialHeatCapacityXS)
    end if
    if (allocated(dEffStoichSolnPhase)) deallocate(dEffStoichSolnPhase)
    if (allocated(dDrivingForceSoln)) deallocate(dDrivingForceSoln)
    if (allocated(dSumMolFractionSoln)) deallocate(dSumMolFractionSoln)

    allocate(dUpdateVar(nElements*2),dUpdateVarLast(nElements*2))
    allocate(iterHistory(nElements,iterGlobalMax))
    allocate(&
        dEffStoichSolnPhase(l,nElements),&
        dSumMolFractionSoln(l),&
        dDrivingForceSoln(l),&
        dPartialExcessGibbs(nSpecies),&
        dPartialEnthalpyXS(nSpecies),&
        dPartialEntropyXS(nSpecies),&
        dPartialHeatCapacityXS(nSpecies))
    lConverged              = .FALSE.
    lRevertSystem           = .FALSE.
    dMaxPotentialTol        = 1D-5
    dUpdateVarLast       = 0D0
    dPartialExcessGibbs=0D0
    dPartialEnthalpyXS=0D0
    dPartialEntropyXS=0D0
    dPartialHeatCapacityXS=0D0
    dGEMFunctionNorm     = 10D0

    dDrivingForceSoln =0d0

    return
end subroutine InitGEMSolver
!
!
