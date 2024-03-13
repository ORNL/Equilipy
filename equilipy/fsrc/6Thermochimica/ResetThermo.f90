


subroutine ResetThermo

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ResetThermo.f90
    !> \brief   Deallocate allocatable variables used by the ModuleThermo.f90, ModulePGESolver.f90 modules.
    !> \author  M.H.A. Piro
    !> \date    Apr. 26, 2012
    !> \sa      ModuleThermo.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer         Description of change
    !    ----          ----------         ---------------------
    !    11/04/2011    M.H.A. Piro        Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to attempt to gracefully exit Thermochimica.  Allocatable
    !! arrays are deallocated and memory is stored for output to external packages.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO                  An error is returned if deallocation is unsuccessful.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                       it encounters an error.  A description for each error is given in ThermoDebug.f90.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo, dElementMass
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer::   i, INFO


    ! Initialize variables:
    i = 0
    if(allocated(dMagEnthalpy)) deallocate(dMagEnthalpy)
    if(allocated(dMagHeatCapacity)) deallocate(dMagHeatCapacity)
    if(allocated(dMagEntropy)) deallocate(dMagEntropy)
    if(allocated(iRegularParam)) deallocate(iRegularParam)
    if(allocated(dCoeffGibbsMagnetic)) deallocate(dCoeffGibbsMagnetic)
    if(allocated(cRegularParam)) deallocate(cRegularParam)
    if(allocated(iPhaseLevel)) deallocate(iPhaseLevel)
    if(allocated(dMagGibbsEnergy)) deallocate(dMagGibbsEnergy)
    if(allocated(iParticlesPerMole)) deallocate(iParticlesPerMole)
    if(allocated(dSpeciesTotalAtoms)) deallocate(dSpeciesTotalAtoms)
    if(allocated(iSpeciesPass)) deallocate(iSpeciesPass)
    if(allocated(iElementSystem)) deallocate(iElementSystem)
    if(allocated(nParamPhase)) deallocate(nParamPhase)
    if(allocated(nSpeciesPhase)) deallocate(nSpeciesPhase)
    if(allocated(iPhase)) deallocate(iPhase)
    if(allocated(dStoichSpeciesLevel)) deallocate(dStoichSpeciesLevel)
    if(allocated(dStoichSpecies)) deallocate(dStoichSpecies)
    if (allocated(dChemicalPotentialGEM))        deallocate(dChemicalPotentialGEM)
    if (allocated(dStoichSpeciesGEM))            deallocate(dStoichSpeciesGEM)
    if (allocated(dAtomFractionSpeciesGEM))      deallocate(dAtomFractionSpeciesGEM)
    if (allocated(iShuffled))                    deallocate(iShuffled)
    if (allocated(iCandidate))                   deallocate(iCandidate)
    if (allocated(iSolnPS))                      deallocate(iSolnPS)
    if (allocated(iAssemblage)) deallocate (iAssemblage)
    if(allocated(dChemicalPotential)) deallocate(dChemicalPotential)
    if(allocated(dMolesElement)) deallocate(dMolesElement)
    if(allocated(dAtomFractionSpecies)) deallocate(dAtomFractionSpecies)
    if(allocated(nMagParamPhase)) deallocate(nMagParamPhase)
    if(allocated(iMagneticParam)) deallocate(iMagneticParam)
    if(allocated(dMagneticParam)) deallocate(dMagneticParam)
    if(allocated(iSUBLParamData)) deallocate(iSUBLParamData)
    if(allocated(dPartialEnthalpy))  deallocate(dPartialEnthalpy)
    if(allocated(dPartialEntropy))  deallocate(dPartialEntropy)
    if(allocated(dPartialHeatCapacity))  deallocate(dPartialHeatCapacity)
    if(allocated(dChemicalPotentialOld)) deallocate(dChemicalPotentialOld)
    if(allocated(dAtomFractionSpeciesOld)) deallocate(dAtomFractionSpeciesOld)
    if(allocated(dExcessGibbsParam)) deallocate(dExcessGibbsParam)
    if(allocated(dExcessHParam)) deallocate(dExcessHParam)
    if(allocated(dExcessSParam)) deallocate(dExcessSParam)
    if(allocated(dExcessCpParam)) deallocate(dExcessCpParam)
    if(allocated(dStdGibbsEnergy)) deallocate(dStdGibbsEnergy)
    if(allocated(dStdEnthalpy)) deallocate(dStdEnthalpy)
    if(allocated(dStdEntropy)) deallocate(dStdEntropy)
    if(allocated(dStdHeatCapacity)) deallocate(dStdHeatCapacity)
    if(allocated(dMolesSpecies)) deallocate(dMolesSpecies)
    if(allocated(dElementPotential)) deallocate(dElementPotential)
    if(allocated(dMolesPhase)) deallocate(dMolesPhase)
    if (allocated(dStoichDependent)) deallocate (dStoichDependent)
    dElementMass=0.0d0
    if(allocated(cElementName)) deallocate(cElementName)
    if(allocated(cSpeciesName)) deallocate(cSpeciesName)
    if(allocated(cSolnPhaseType)) deallocate(cSolnPhaseType)
    if(allocated(cSolnPhaseName)) deallocate(cSolnPhaseName)
    if(allocated(dAtomicMass)) deallocate(dAtomicMass)
    if (allocated(iterHistoryLevel)) deallocate (iterHistoryLevel)
    if (allocated(dPhasePotential)) deallocate (dPhasePotential)
    if (allocated(iterHistory)) deallocate (iterHistory)
    if(allocated(dSumMolFractionSoln)) deallocate(dSumMolFractionSoln)
    if(allocated(dDrivingForceSoln)) deallocate(dDrivingForceSoln)
    if(allocated(dEffStoichSolnPhase)) deallocate(dEffStoichSolnPhase)
    if(allocated(dMolesSpeciesPEA)) deallocate(dMolesSpeciesPEA)
    if(allocated(dChemicalPotentialPEA)) deallocate(dChemicalPotentialPEA)
    if(allocated(dMolFractionPEA)) deallocate(dMolFractionPEA)
    if(allocated(dMolesPhasePEA)) deallocate(dMolesPhasePEA)
    if(allocated(dEffStoichSolnPhasePEA)) deallocate(dEffStoichSolnPhasePEA)
    if(allocated(dGibbsEnergySysHist)) deallocate(dGibbsEnergySysHist)
    if(allocated(dPartialExcessGibbsLast)) deallocate(dPartialExcessGibbsLast)
    if(allocated(dPartialExcessGibbs)) deallocate(dPartialExcessGibbs)
    if(allocated(dPartialEnthalpyXS)) deallocate(dPartialEnthalpyXS)
    if(allocated(dPartialEntropyXS)) deallocate(dPartialEntropyXS)
    if(allocated(dPartialHeatCapacityXS)) deallocate(dPartialHeatCapacityXS)
    if(allocated(dMolFraction)) deallocate(dMolFraction)
    if(allocated(dUpdateVar)) deallocate(dUpdateVar)
    if(allocated(lSolnPhases)) deallocate(lSolnPhases)
    if(allocated(dGibbsSolnPhase)) deallocate(dGibbsSolnPhase)
    if(allocated(lMiscibility)) deallocate(lMiscibility)
    if(allocated(iPhaseSublattice)) deallocate(iPhaseSublattice)
    if(allocated(dStoichSublattice)) deallocate(dStoichSublattice)
    if(allocated(cConstituentNameSUB)) deallocate(cConstituentNameSUB)
    if(allocated(iConstituentSublattice)) deallocate(iConstituentSublattice)
    if(allocated(nSublatticePhase)) deallocate(nSublatticePhase)
    if(allocated(nConstituentSublattice)) deallocate(nConstituentSublattice)
    if(allocated(nSublatticeElements)) deallocate(nSublatticeElements)
    if(allocated(dSublatticeCharge)) deallocate(dSublatticeCharge)
    if(allocated(iConstituentPass)) deallocate(iConstituentPass)
    if(allocated(dSiteFraction)) deallocate(dSiteFraction)
    if(allocated(iPhaseElectronID)) deallocate(iPhaseElectronID)
    if(allocated(dZetaSpecies)) deallocate(dZetaSpecies)
    if(allocated(dConstituentCoefficients)) deallocate(dConstituentCoefficients)
    if(allocated(iChemicalGroup)) deallocate(iChemicalGroup)
    if(allocated(cPairName)) deallocate(cPairName)
    if(allocated(dStoichPairs)) deallocate(dStoichPairs)
    if(allocated(iHessian)) deallocate(iHessian)
    if(allocated(dChemicalPotentialStar)) deallocate(dChemicalPotentialStar)
    if(allocated(dChemicalPotentialDiff)) deallocate(dChemicalPotentialDiff)
    if(allocated(dRHS)) deallocate(dRHS)
    if(allocated(dHessian)) deallocate(dHessian)
    if (allocated(iPairID)) deallocate(iPairID)
    if (allocated(dCoordinationNumber)) deallocate(dCoordinationNumber)
    if (allocated(nPairsSRO)) deallocate(nPairsSRO)
    if (i > 0) then
        INFOThermo = 15
    end if

    return

end subroutine ResetThermo
