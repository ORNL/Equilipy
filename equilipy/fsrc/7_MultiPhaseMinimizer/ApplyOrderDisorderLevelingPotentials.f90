subroutine ApplyOrderDisorderLevelingPotentials

    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer :: i, k, m, n, nConstituents
    real(8) :: dLevelingFormulaAtomDenom
    real(8), dimension(:), allocatable :: dMolFractionSaved, dMolesSpeciesSaved
    real(8), dimension(:), allocatable :: dChemicalPotentialSaved, dGibbsSolnPhaseSaved
    real(8), dimension(:), allocatable :: dPartialEnthalpySaved, dPartialEntropySaved
    real(8), dimension(:), allocatable :: dPartialHeatCapacitySaved
    real(8), dimension(:), allocatable :: dPartialExcessGibbsSaved, dPartialEnthalpyXSSaved
    real(8), dimension(:), allocatable :: dPartialEntropyXSSaved, dPartialHeatCapacityXSSaved

    if (nSolnPhasesSys <= 0) return
    if (.NOT.allocated(cSolnPhaseType)) return
    if (.NOT.allocated(dMolFraction)) return
    if (.NOT.allocated(dMolesSpecies)) return
    if (.NOT.allocated(dPartialExcessGibbs)) then
        allocate(dPartialExcessGibbs(nSpecies), dPartialEnthalpyXS(nSpecies), &
            dPartialEntropyXS(nSpecies), dPartialHeatCapacityXS(nSpecies))
        dPartialExcessGibbs = 0D0
        dPartialEnthalpyXS = 0D0
        dPartialEntropyXS = 0D0
        dPartialHeatCapacityXS = 0D0
    end if

    allocate(dMolFractionSaved(nSpecies), dMolesSpeciesSaved(nSpecies))
    allocate(dChemicalPotentialSaved(nSpecies), dGibbsSolnPhaseSaved(nSolnPhasesSys))
    allocate(dPartialEnthalpySaved(nSpecies), dPartialEntropySaved(nSpecies))
    allocate(dPartialHeatCapacitySaved(nSpecies))
    allocate(dPartialExcessGibbsSaved(nSpecies), dPartialEnthalpyXSSaved(nSpecies))
    allocate(dPartialEntropyXSSaved(nSpecies), dPartialHeatCapacityXSSaved(nSpecies))

    dMolFractionSaved = dMolFraction
    dMolesSpeciesSaved = dMolesSpecies
    dChemicalPotentialSaved = dChemicalPotential
    dGibbsSolnPhaseSaved = dGibbsSolnPhase
    dPartialEnthalpySaved = dPartialEnthalpy
    dPartialEntropySaved = dPartialEntropy
    dPartialHeatCapacitySaved = dPartialHeatCapacity
    dPartialExcessGibbsSaved = dPartialExcessGibbs
    dPartialEnthalpyXSSaved = dPartialEnthalpyXS
    dPartialEntropyXSSaved = dPartialEntropyXS
    dPartialHeatCapacityXSSaved = dPartialHeatCapacityXS

    do k = 1, nSolnPhasesSys
        if (TRIM(cSolnPhaseType(k)) /= 'SUBOM') cycle

        m = nSpeciesPhase(k-1) + 1
        n = nSpeciesPhase(k)
        nConstituents = n - m + 1
        if (nConstituents <= 0) cycle

        do i = m, n
            dMolFraction(m:n) = 0D0
            dMolesSpecies(m:n) = 0D0
            dMolFraction(i) = 1D0
            dMolesSpecies(i) = 1D0
            dGibbsSolnPhase(k) = 0D0

            call CompExcessGibbsEnergy(k)
            if (INFOThermo /= 0) exit

            dLevelingFormulaAtomDenom = dLevelingSpeciesFormulaAtoms(i) / DBLE(iParticlesPerMole(i))
            if (dLevelingFormulaAtomDenom > 1D-300) then
                dLevelingChemicalPotential(i) = dChemicalPotential(i) / dLevelingFormulaAtomDenom
            else
                dLevelingChemicalPotential(i) = 5D9
            end if
        end do
        if (INFOThermo /= 0) exit
    end do

    dMolFraction = dMolFractionSaved
    dMolesSpecies = dMolesSpeciesSaved
    dChemicalPotential = dChemicalPotentialSaved
    dGibbsSolnPhase = dGibbsSolnPhaseSaved
    dPartialEnthalpy = dPartialEnthalpySaved
    dPartialEntropy = dPartialEntropySaved
    dPartialHeatCapacity = dPartialHeatCapacitySaved
    dPartialExcessGibbs = dPartialExcessGibbsSaved
    dPartialEnthalpyXS = dPartialEnthalpyXSSaved
    dPartialEntropyXS = dPartialEntropyXSSaved
    dPartialHeatCapacityXS = dPartialHeatCapacityXSSaved

    deallocate(dMolFractionSaved, dMolesSpeciesSaved, dChemicalPotentialSaved)
    deallocate(dGibbsSolnPhaseSaved, dPartialEnthalpySaved, dPartialEntropySaved)
    deallocate(dPartialHeatCapacitySaved, dPartialExcessGibbsSaved)
    deallocate(dPartialEnthalpyXSSaved, dPartialEntropyXSSaved, dPartialHeatCapacityXSSaved)

    return

end subroutine ApplyOrderDisorderLevelingPotentials
