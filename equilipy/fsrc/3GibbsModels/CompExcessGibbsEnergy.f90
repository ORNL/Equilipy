!
!
subroutine CompExcessGibbsEnergy(iSolnIndex)
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompExcessGibbsEnergy.f90
    !> \brief   Compute the partial molar excess Gibbs energy of mixing of solution phase constituents by calling
    !!           model specific subroutines.
    !> \author  M.H.A. Piro
    !> \date    April 1, 2018
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   01/14/2013      M.H.A. Piro         Original code.
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to call a specific subroutine to compute the partial molar
    !! excess Gibbs energy of mixing for each constituent in a solution phase based on the model type.
    !
    ! Note: the chemical potential term is computed in this subroutine as opposed to CompChemicalPotential.f90
    ! because the chemical potential of a phase component is model dependent.  See SUBL for an example below.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in] iSolnIndex    An integer scalar representing the absolute solution phase index.
    !
    ! iFirst                    An integer scalar representing the first species in a solution phase.
    ! iLast                     An integer scalar representing the last species in a solution phase.
    ! dChemicalPotential        A double real vector represening the chemical potential for every species in the
    !                            system.
    ! dStdGibbsEnergy           A double real vector represending the standard molar Gibbs energy of every pure
    !                            species in the system.
    ! dMolFraction              A double real vector representing the mole fraction for every species in the
    !                            system.
    ! dPartialExcessGibbs       A double real vector representing the partial molar excess Gibbs energy of mixing
    !                            for every species in the system.
    ! dMolesSpecies             A double real vector representing the number of moles for every species in the
    !                            system.
    ! dGibbsSolnPhase          A double real vector represending the molar Gibbs energy for every solution phase
    !                            in the system.
    !
    !-------------------------------------------------------------------------------------------------------------
!
    USE ModuleThermo
    USE ModuleThermoIO!, ONLY: INFOThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer :: i, j, iSolnIndex, iFirst, iLast
    real(8) :: dSpeciesTotalMole
    
!
!
    ! Temporary variables used for convenience:
    iFirst = nSpeciesPhase(iSolnIndex - 1) + 1
    iLast  = nSpeciesPhase(iSolnIndex)
!
    ! Do not compute when none of the species are considered.
    if (iLast-iFirst<1) then
        dChemicalPotential(iFirst:iLast)    = dStdGibbsEnergy(iFirst:iLast)
        dPartialEnthalpy(iFirst:iLast)      = dStdEnthalpy(iFirst:iLast)
        dPartialEntropy(iFirst:iLast)       = dStdEntropy(iFirst:iLast)
        dPartialHeatCapacity(iFirst:iLast)  = dStdHeatCapacity(iFirst:iLast)
        return
    end if

    dChemicalPotential(iFirst:iLast)        = 0D0
    dPartialEnthalpy(iFirst:iLast)          = 0D0
    dPartialEntropy(iFirst:iLast)           = 0D0
    dPartialHeatCapacity(iFirst:iLast)      = 0D0

    dSpeciesTotalMole = 0D0
    
    dPartialExcessGibbs(iFirst:iLast) = 0D0 !Added by SY
    dPartialEnthalpyXS(iFirst:iLast)  = 0D0 !Added by SY
    dPartialEntropyXS(iFirst:iLast)   = 0D0 !Added by SY
    dPartialHeatCapacityXS(iFirst:iLast) = 0D0 !Added by SY
!
!
    ! Compute excess terms based on solution phase type:
    select case (cSolnPhaseType(iSolnIndex))
        case ('IDMX')
!
            ! Compute the chemical potentials of each species and the molar Gibbs energy of the phase:
            do i = iFirst, iLast
                dChemicalPotential(i)       = dStdGibbsEnergy(i) + DLOG(DMAX1(dMolFraction(i), 1D-75))
                dPartialEnthalpy(i)         = dStdEnthalpy(i) 
                dPartialEntropy(i)          = dStdEntropy(i)  - DLOG(DMAX1(dMolFraction(i), 1D-75))
                dPartialHeatCapacity(i)     = dStdHeatCapacity(i) 
                dGibbsSolnPhase(iSolnIndex) = dGibbsSolnPhase(iSolnIndex) + dChemicalPotential(i) * dMolesSpecies(i)
            end do
!
        case ('QKTO')
!
            ! Compute the excess terms for a Quasichemical Kohler-TOop (QKTO) model:
            call CompExcessGibbsEnergyQKTO(iSolnIndex)
!
            ! Compute the chemical potentials of each species and the molar Gibbs energy of the phase:
            do i = iFirst, iLast
                dChemicalPotential(i)       = dStdGibbsEnergy(i) + DLOG(DMAX1(dMolFraction(i), 1D-75)) + dPartialExcessGibbs(i)
                dPartialEnthalpy(i)         = dStdEnthalpy(i) + dPartialEnthalpyXS(i)
                dPartialEntropy(i)          = dStdEntropy(i)  - DLOG(DMAX1(dMolFraction(i), 1D-75)) + dPartialEntropyXS(i)
                dPartialHeatCapacity(i)     = dStdHeatCapacity(i) + dPartialHeatCapacityXS(i)

                dGibbsSolnPhase(iSolnIndex) = dGibbsSolnPhase(iSolnIndex) + dChemicalPotential(i) * dMolesSpecies(i)
            end do
!
        case ('RKMP','RKMPM')
!
            ! Compute magnetic terms if this is a magnetic phase:
            if (cSolnPhaseType(iSolnIndex) == 'RKMPM') call CompGibbsMagneticSoln(iSolnIndex)
!
            ! Compute the excess terms for a Redlich-Kister-Muggiano-Polynomial (RKMP) model:
            call CompExcessGibbsEnergyRKMP(iSolnIndex)
!
            ! Compute the chemical potentials of each species and the molar Gibbs energy of the phase:
            do i = iFirst, iLast
                dChemicalPotential(i)  = dStdGibbsEnergy(i) + DLOG(DMAX1(dMolFraction(i), 1D-75)) + dPartialExcessGibbs(i)&
                                        + dMagGibbsEnergy(i)
                dPartialEnthalpy(i)    = dStdEnthalpy(i) + dPartialEnthalpyXS(i) + dMagEnthalpy(i)
                
                dPartialEntropy(i)     = dStdEntropy(i)  - DLOG(DMAX1(dMolFraction(i), 1D-75)) + dPartialEntropyXS(i)&
                                        + dMagEntropy(i)
                dPartialHeatCapacity(i)     = dStdHeatCapacity(i) + dPartialHeatCapacityXS(i) + dMagHeatCapacity(i)
                
                dGibbsSolnPhase(iSolnIndex) = dGibbsSolnPhase(iSolnIndex) + dChemicalPotential(i) * dMolesSpecies(i)
                
            end do
!
        case ('SUBL','SUBLM')
!
            ! Compute the excess terms for a Compound Energy Formalism model:
            call CompExcessGibbsEnergySUBL(iSolnIndex)
!
            ! Compute the chemical potentials of each species and the molar Gibbs energy of the phase:
            do i = iFirst, iLast
                ! print*,dSpeciesTotalMole
                dChemicalPotential(i)       = dChemicalPotential(i) + dPartialExcessGibbs(i) + dMagGibbsEnergy(i)
                dPartialEnthalpy(i)         = dPartialEnthalpy(i) + dPartialEnthalpyXS(i) + dMagEnthalpy(i)
                dPartialEntropy(i)          = dPartialEntropy(i)  + dPartialEntropyXS(i) + dMagEntropy(i)
                dPartialHeatCapacity(i)     = dPartialHeatCapacity(i) + dPartialHeatCapacityXS(i) + dMagHeatCapacity(i)
                dGibbsSolnPhase(iSolnIndex) = dGibbsSolnPhase(iSolnIndex) + dChemicalPotential(i) * dMolesSpecies(i)

                dSpeciesTotalMole = MAX1(dSpeciesTotalMole,sum(dStoichSpecies(i,:)))

                
            end do
!
        case ('SUBG','SUBQ')
!
            ! Compute the excess terms for a phase represented by the Modified Quasichemical Model:
            call CompExcessGibbsEnergySUBG(iSolnIndex)
!
            ! Compute the chemical potentials of each species and the molar Gibbs energy of the phase:
            do i = iFirst, iLast
                dChemicalPotential(i)       = dChemicalPotential(i) + dPartialExcessGibbs(i)
                dPartialEnthalpy(i)         = dPartialEnthalpy(i) + dPartialEnthalpyXS(i)
                dPartialEntropy(i)          = dPartialEntropy(i)  + dPartialEntropyXS(i)
                dPartialHeatCapacity(i)     = dPartialHeatCapacity(i) + dPartialHeatCapacityXS(i)
                dGibbsSolnPhase(iSolnIndex) = dGibbsSolnPhase(iSolnIndex) + dChemicalPotential(i) * dMolesSpecies(i)
            end do
        case default
            ! Report an error (solution phase type unsupported).  Note that there is an earlier check when parsing
            ! the data-file; although this is redundant, it is conservative.
!
            INFOThermo = 17
            return
!
    end select
!
    return
!
end subroutine CompExcessGibbsEnergy
!
!
