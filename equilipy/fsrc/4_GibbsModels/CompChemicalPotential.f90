!
!
subroutine CompChemicalPotential(lCompEverything)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompChemicalPotential.f90
    !> \brief   Compute the chemical potentials of all solution phase constituents.
    !> \author  M.H.A. Piro
    !> \date    Apr. 25, 2012
    !> \sa      CompMolFraction.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   04/25/2012      M.H.A. Piro         Original code
    !   06/13/2012      M.H.A. Piro         Include excess parameters
    !   01/18/2013      M.H.A. Piro         The chemical potential terms for each phase component are now
    !                                        computed in the CompExcessGibbsEnergy.f90 subroutine because
    !                                        the equation is model dependent.
    !   06/24/2026      S.Y. Kwon           Kept active solution species amounts synchronized with
    !                                        model-projected endmember fractions after SUBL/SUBOM evaluation.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the chemical potentials of all solution phase
    !! constituents expected to be stable at equilibrium.
    !
    !
    ! Pertinent variables:
    ! ====================
    !
    !> \param[in]   lCompEverything     A logical scalar indicating whether everything should be computed, or
    !!                                   only what is necessary.  For the most part, it is only necessary to
    !!                                   compute chemical potentials of solution species that are expected to
    !!                                   be stable.
    !
    ! dChemicalPotential                A double real vector representing the chemical potential of a substance.
    ! dGibbsSolnPhase                   A double real vector representing the Gibbs energy of a solution phase.
    ! dPartialExcessGibbs               A double real vector representing the partial molar excess Gibbs energy
    !                                    of mixing of each solution phase constituent.
    ! nSpeciesPhase                     An integer vector representing the index of the last species in a
    !                                   particular solution phase.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
!
    implicit none
!
    integer :: i, j, k, m, n, iSlot
    real(8) :: dTemp
    logical :: lCompEverything
    logical :: lSyncProjectedMoles
!
!
    ! Initialize variables:
    dGibbsSolnPhase       = 0D0
    dChemicalPotential    = dStdGibbsEnergy
    dPartialEnthalpy      = dStdEnthalpy
    dPartialEntropy       = dStdEntropy
    dPartialHeatCapacity  = dStdHeatCapacity
    lSolnPhases           = .False.
!
    ! Compute the mole fractions of species in solution phases expected to be stable:
    do j = 1, nSolnPhases
!
        iSlot  = nElements - j + 1          ! Index of phase in iAssemblage
        dTemp  = dMolesPhase(iSlot)
        k      = -iAssemblage(iSlot)        ! Absolute index of solution phase
        m      = nSpeciesPhase(k-1) + 1
        n      = nSpeciesPhase(k)
        lSolnPhases(k) = .True.
!
        ! Compute the mole fractions of all solution phase constituents in the phase:
        do i = nSpeciesPhase(k-1) + 1, nSpeciesPhase(k)
            if (dTemp>1D-15) dMolFraction(i) = dMolesSpecies(i)/dTemp
            dPhasePotential(i) = 0d0
        end do
!  
        ! Compute excess terms:
        call CompExcessGibbsEnergy(k)
!
        ! SUBL/SUBOM Gibbs models project the caller's trial species fractions
        ! into a thermodynamically consistent endmember basis through site
        ! fractions.  Keep species moles in the same basis before the
        ! Lagrangian mass residual and next Newton matrix are assembled.  If a
        ! selected-component calculation has collapsed some endmembers onto
        ! omitted constituents, leave the original reduced-system moles in
        ! place; those endmembers need active-component filtering first.
        lSyncProjectedMoles = MINVAL(dSpeciesTotalAtoms(m:n)) > 1D-4
        if ((dTemp > 0D0).AND.lSyncProjectedMoles) then
            dMolesSpecies(m:n) = dTemp * dMolFraction(m:n)
            dMolesPhase(iSlot) = SUM(dMolesSpecies(m:n))
        end if
!
    end do
!
    ! Check if the chemical potentials for everything should be computed:
    if (lCompEverything) then
!
        ! Compute the chemical potentials for every solution phase in the system:
        LOOP_SolnPhasesSys: do j = 1, nSolnPhasesSys
!
            ! Skip this phase if it is already part of the system:
            if (lSolnPhases(j)) cycle LOOP_SolnPhasesSys
!
            ! Compute chemical potentials of unstable solution phase constituents:
            call CompMolFraction(j)
        end do LOOP_SolnPhasesSys
    end if
!
    return
!
end subroutine CompChemicalPotential
!
!
!
