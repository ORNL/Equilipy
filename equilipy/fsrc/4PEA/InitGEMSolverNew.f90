!
subroutine InitGEMSolverNew
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    implicit none
    integer::   i, j, k, l, m, n
    ! Initialize variables:
    nConPhases        = nElements
    nSpeciesLevel     = nSpecies + nSolnPhasesSys
    l = MAX(1,nSolnPhasesSys)
!
    ! Check to see if allocatable arrays have already been allocated:
    if (allocated(dMolesPhase))             deallocate(dMolesPhase)
    if (allocated(dElementPotential))       deallocate(dElementPotential)
    if (allocated(dElementPotentialLast))       deallocate(dElementPotentialLast)
    if (allocated(iAssemblage))             deallocate(iAssemblage)
    if (allocated(iterHistoryLevel))        deallocate(iterHistoryLevel)
    if (allocated(dPhasePotential))         deallocate(dPhasePotential)
    if (allocated(dAtomFractionSpeciesOld)) deallocate(dAtomFractionSpeciesOld)
    if (allocated(dChemicalPotentialOld))   deallocate(dChemicalPotentialOld)
    if (allocated(dStoichSpeciesLevel))     deallocate(dStoichSpeciesLevel)
    if (allocated(iPhaseLevel))             deallocate(iPhaseLevel)
    if (allocated(dChemicalPotentialGEM))   deallocate(dChemicalPotentialGEM)
    if (allocated(dAtomFractionSpeciesGEM)) deallocate(dAtomFractionSpeciesGEM)
    if (allocated(dStoichSpeciesGEM))       deallocate(dStoichSpeciesGEM)
    if (allocated(dMolFractionGEM))         deallocate(dMolFractionGEM)
    if (allocated(iShuffled))               deallocate(iShuffled)
    if (allocated(dMolesPhaseHistory))      deallocate(dMolesPhaseHistory)
    if (allocated(dEffStoichSolnPhase))     deallocate(dEffStoichSolnPhase)
    if (allocated(iPhaseGEM))               deallocate(iPhaseGEM)
    if (allocated(iCandidate))              deallocate(iCandidate)
    allocate(dMolesPhase(nElements),dElementPotential(nElements),dElementPotentialLast(nElements),&
        iAssemblage(nElements),dAtomFractionSpeciesOld(nSpecies,nElements),&
        dChemicalPotentialOld(nSpecies),iterHistoryLevel(nElements,500),&
        dStoichSpeciesLevel(nSpeciesLevel,nElements),iPhaseLevel(nSpeciesLevel),&
        dAtomFractionSpeciesGEM(nElements,nElements),dChemicalPotentialGEM(nElements),&
        dStoichSpeciesGEM(nElements,nElements),dMolFractionGEM(nElements,nSpecies),&
        iShuffled(nElements),dMolesPhaseHistory(nElements,500),iCandidate(nSpecies),iPhaseGEM(nElements))
    allocate(dEffStoichSolnPhase(l,nElements))
    
    ! Initialize allocatable variables:
    dElementPotential = 0D0
    dElementPotentialLast=0D0
    iAssemblage       = 0
    iShuffled         = (/(i, i = 1,nElements, 1)/)
    iCandidate        = 0
    iPhaseGEM         = 0
    dAtomFractionSpeciesOld = dAtomFractionSpecies
    dChemicalPotentialOld   = dChemicalPotential
    dMolFractionGEM   = 0D0
    dMolesPhase       = 0D0
    dEffStoichSolnPhase=0D0
!
    !SY: Modifying dAtomFractionSpecies, dChemicalPotential, dPhasePotential
    !These needs to be modified back to the size of nSpecies
    if (allocated(dAtomFractionSpecies))     deallocate(dAtomFractionSpecies)
    if (allocated(dChemicalPotential))       deallocate(dChemicalPotential)
    ! deallocate(dAtomFractionSpecies,dChemicalPotential)
    allocate(dAtomFractionSpecies(nSpeciesLevel,nElements),&
    dChemicalPotential(nSpeciesLevel),dPhasePotential(nSpeciesLevel))
!
!
    dAtomFractionSpecies = 0D0
    dChemicalPotential   = 0D0
    dPhasePotential      = 0D0
    dStoichSpeciesLevel  = 0D0
    dStoichSpeciesGEM    = 0D0
    dAtomFractionSpeciesGEM= 0D0
    dChemicalPotentialGEM = 0D0
    dGEMFunctionNorm     = 10D0
    
!
!
    dAtomFractionSpecies(:nSpecies,:) = dAtomFractionSpeciesOld
    dChemicalPotential(:nSpecies)     = dChemicalPotentialOld
    dStoichSpeciesLevel(:nSpecies,:)  = dStoichSpecies
    iPhaseLevel(:nSpecies)            = iPhase
!
    ! Establish the very first phase assemblage:
    call GetFirstAssemblage
!
    do i =1,nElements
        iPhaseGEM(i)= iPhaseLevel(iAssemblage(i))
    end do
    ! Calculate the minimum point of each solution phase assuming zero ElementPotential
    call CompInitMinSolnPoint
!
!

    do j = 1,nElements
        do i = 1,nElements
            dStoichSpeciesGEM(j,i)=dStoichSpeciesLevel(iAssemblage(j),i)
            dAtomFractionSpeciesGEM(j,i)=dAtomFractionSpecies(iAssemblage(j),i)
        end do
        dChemicalPotentialGEM(j) = dChemicalPotential(iAssemblage(j))
        k = iAssemblage(j)
        l = iPhaseLevel(k)
        if(l>0) then
            m = nSpeciesPhase(l-1) + 1      ! First constituent in phase.
            n = nSpeciesPhase(l)            ! Last  constituent in phase.
            dMolFractionGEM(j,m:n) = dMolFraction(m:n)
!
        end if
    end do

!
    return
end subroutine InitGEMSolverNew
!
!
