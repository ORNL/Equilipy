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
    !                                       computed in the CompExcessGibbsEnergy.f90 subroutine because
    !                                       the equation is model dependent.
    !   07/20/2026      S.Y. Kwon           Kept active solution amounts synchronized and evaluated same-parent composition sets from canonical constitutions.
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
    logical :: lSyncProjectedMoles, lUseSlotFraction
!
!
    ! Initialize variables:
    dGibbsSolnPhase       = 0D0
    dChemicalPotential    = dStdGibbsEnergy
    dPartialEnthalpy      = dStdEnthalpy
    dPartialEntropy       = dStdEntropy
    dPartialHeatCapacity  = dStdHeatCapacity
    lSolnPhases           = .False.
    if (allocated(lActiveSlotPropValid)) lActiveSlotPropValid = .FALSE.
    if (allocated(dActiveSlotChemPot)) dActiveSlotChemPot = 0D0
    if (allocated(dActiveSlotPartialH)) dActiveSlotPartialH = 0D0
    if (allocated(dActiveSlotPartialS)) dActiveSlotPartialS = 0D0
    if (allocated(dActiveSlotPartialCp)) dActiveSlotPartialCp = 0D0
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
        ! A canonical order/disorder parent can own multiple composition-set
        ! slots.  Its shared species vector cannot represent both constitutions,
        ! so evaluate each repeated slot from the slot-local product fractions.
        lUseSlotFraction = CanonicalRepeatedParentSlot(iSlot, k)
        if (lUseSlotFraction) dMolFraction(m:n) = dActiveSlotMolFraction(iSlot,m:n)

        ! Compute the mole fractions of all solution phase constituents in the phase:
        do i = m, n
            if ((.NOT.lUseSlotFraction).AND.(dTemp>1D-15)) then
                dMolFraction(i) = dMolesSpecies(i)/dTemp
            end if
            dPhasePotential(i) = 0d0
        end do
!  
        ! Compute excess terms:
        call CompExcessGibbsEnergy(k)

        ! Preserve the properties immediately after this slot's evaluation.
        ! The phase-global arrays are retained for legacy callers and may be
        ! overwritten when another composition set shares the same parent.
        if (allocated(dActiveSlotMolFraction)) then
            dActiveSlotMolFraction(iSlot,m:n) = dMolFraction(m:n)
        end if
        if (allocated(dActiveSlotChemPot)) then
            dActiveSlotChemPot(iSlot,m:n) = dChemicalPotential(m:n)
            dActiveSlotPartialH(iSlot,m:n) = dPartialEnthalpy(m:n)
            dActiveSlotPartialS(iSlot,m:n) = dPartialEntropy(m:n)
            dActiveSlotPartialCp(iSlot,m:n) = dPartialHeatCapacity(m:n)
            lActiveSlotPropValid(iSlot) = .TRUE.
        end if
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
contains

    logical function CanonicalRepeatedParentSlot(iSlotIn, iPhaseIn)
        integer, intent(in) :: iSlotIn, iPhaseIn

        integer :: iSlotLocal, nMatchingSlots

        CanonicalRepeatedParentSlot = .FALSE.
        if (.NOT.lODPartitionUnifiedActive) return
        if (.NOT.allocated(iODCompanionPhase)) return
        if (.NOT.allocated(iODTopologyClass)) return
        if (.NOT.allocated(iActiveSlotThermoPhase)) return
        if (.NOT.allocated(iActiveSlotIdentityOrdinal)) return
        if (.NOT.allocated(dActiveSlotMolFraction)) return
        if ((iPhaseIn < 1).OR.(iPhaseIn > SIZE(iODCompanionPhase))) return
        if ((iODTopologyClass(iPhaseIn) < OD_TOPOLOGY_HELPER_STANDALONE).OR.&
            (iODTopologyClass(iPhaseIn) > OD_TOPOLOGY_HELPER_ONLY)) return
        if (iODCompanionPhase(iPhaseIn) <= 0) return
        if (iActiveSlotThermoPhase(iSlotIn) /= iPhaseIn) return
        if (iActiveSlotIdentityOrdinal(iSlotIn) <= 0) return

        nMatchingSlots = 0
        do iSlotLocal = 1, nElements
            if (iAssemblage(iSlotLocal) >= 0) cycle
            if (iActiveSlotThermoPhase(iSlotLocal) /= iPhaseIn) cycle
            nMatchingSlots = nMatchingSlots + 1
        end do
        CanonicalRepeatedParentSlot = nMatchingSlots > 1

        return
    end function CanonicalRepeatedParentSlot

end subroutine CompChemicalPotential
!
!
!
