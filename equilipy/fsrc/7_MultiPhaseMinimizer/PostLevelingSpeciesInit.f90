!> \brief Seed active solution endmembers after Leveling/PEA before Lagrangian GEM.
!!
!! \details Active solution phases may enter Lagrangian GEM from a single
!! Leveling pseudo-endmember.  This routine keeps nonzero selected endmembers,
!! applies a degeneracy-aware Eriksson exponential estimate only to zero
!! endmembers, replaces the degenerate exp(0) estimate with a trace seed, and
!! normalizes the final phase-local endmember vector.
!!
!! \param[out]  lInitialized  True if at least one active solution phase was initialized.



subroutine PostLevelingSpeciesInit(lInitialized)

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    PostLevelingSpeciesInit.f90
    !> \brief   Initialize active solution species amounts after Leveling/PEA.
    !> \author  S.Y. Kwon
    !> \date    Jun. 23, 2026
    !> \sa      Level2Lagrange.f90
    !> \sa      CompExcessGibbsEnergySUBL.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Initialized post-Leveling species amounts from mass-balanced phase constitutions and formula-atom scales.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to convert active Leveling/PEA
    !! solution candidates into Lagrangian GEM species amounts.  Nonzero
    !! selected endmembers keep their Leveling amount.  Zero endmembers receive
    !! a degeneracy-aware Eriksson-style exponential estimate from the
    !! Leveling element potentials and corrected endmember Gibbs energy.  When
    !! a phase already has selected species moles, those Eriksson estimates are
    !! used only as trace relative weights so the mass-balanced Leveling/PEA
    !! constitution is not overwritten before Lagrangian GEM.  If several zero
    !! endmembers have the same scaled composition and corrected Leveling
    !! intercept, they share one Eriksson class weight.  If the class estimate
    !! is degenerate, i.e., exp(0), the class is introduced only as a trace
    !! amount.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! iAssemblage                 Active assemblage after Leveling/PEA.
    ! dMolesPhase                 Active phase amounts in the Lagrangian slot layout.
    ! dMolesSpecies               Selected species amounts from Leveling/PEA candidate rows.
    ! dElementPotential           Leveling element potentials.
    ! dStoichSpecies              Elemental stoichiometry of each endmember.
    ! dLevelingChemicalPotential  Corrected Leveling Gibbs energy for each endmember.
    ! dLevelingSpeciesFormulaAtoms
    !                             Formula-atom scaling for the Leveling Gibbs energy.
    ! iParticlesPerMole           Formula-particle scaling for species potentials.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out]  lInitialized  True if at least one active solution phase was initialized.
    !
    ! dMolFraction                Active phase-local endmember fractions before CompChemicalPotential.
    ! dMolesSpecies               Active species amounts reconstructed from phase amount and mole fraction.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! EstimateZeroEndmemberFraction  Internal helper for Eriksson trace seeding.
    ! SplitDegenerateSeedClasses     Shares one Eriksson seed across equivalent zero endmembers.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! Level2Lagrange              Builds the Lagrangian starting point after Leveling/PEA.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - The Eriksson seed is applied only to zero endmembers when Leveling has
    !   selected a positive species amount for the phase.
    ! - A degenerate exp(0) estimate does not imply a mole fraction of one; it
    !   means the Leveling elemental potentials cannot distinguish that
    !   endmember, so a trace seed is used instead.
    ! - Site permutations that Leveling cannot distinguish share one Eriksson
    !   class weight, which is divided across the class members.
    ! - If selected species moles already exist, zero-endmember seeds are capped
    !   at trace scale and the selected nonzero fractions are rescaled only by
    !   that trace amount.  This preserves the Leveling/PEA mass-balanced
    !   constitution while avoiding exact-zero endmembers.
    ! - If no selected species moles exist, the whole phase-local endmember
    !   vector is normalized after Eriksson seeding so CompChemicalPotential
    !   receives a valid phase composition.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo

    implicit none

    logical, intent(out) :: lInitialized

    real(8), parameter :: dTraceZeroSeed      = 1D-12
    real(8), parameter :: dArgZeroTolerance   = 1D-12
    real(8), parameter :: dClassCompositionTol = 1D-12
    real(8), parameter :: dClassPotentialTol   = 1D-10
    real(8), parameter :: dLogOverflowLimit   = 709D0
    real(8), parameter :: dLogUnderflowLimit  = -745D0

    integer :: i, j, k, l, m, n, iSlot
    integer :: nSeeded
    real(8) :: dArg, dPhaseMoles, dParticleDenom, dSpeciesMolesSum, dSum
    real(8) :: dExistingFractionSum, dSeedSum, dSeedScale, dSeedTarget
    real(8), dimension(:), allocatable :: dSeedFraction
    logical, dimension(:), allocatable :: lNeedsSeed
    logical :: lHasExistingSpeciesMoles

    lInitialized = .FALSE.

    if (nSolnPhases <= 0) return

    do l = 1, nSolnPhases
        iSlot = nElements - l + 1
        k = -iAssemblage(iSlot)
        if ((k < 1).OR.(k > nSolnPhasesSys)) cycle

        m = nSpeciesPhase(k-1) + 1
        n = nSpeciesPhase(k)
        dPhaseMoles = DMAX1(dMolesPhase(iSlot), 0D0)
        dSpeciesMolesSum = SUM(dMolesSpecies(m:n))
        lHasExistingSpeciesMoles = (dPhaseMoles > 0D0).AND.(dSpeciesMolesSum > 0D0)

        if (allocated(dSeedFraction)) deallocate(dSeedFraction)
        if (allocated(lNeedsSeed)) deallocate(lNeedsSeed)
        allocate(dSeedFraction(m:n), lNeedsSeed(m:n))
        dSeedFraction = 0D0
        lNeedsSeed    = .FALSE.

        do i = m, n
            if (lHasExistingSpeciesMoles.AND.(dMolesSpecies(i) > 0D0)) then
                dMolFraction(i) = dMolesSpecies(i) / dPhaseMoles
            else
                call EstimateZeroEndmemberFraction(i, dSeedFraction(i))
                dMolFraction(i) = 0D0
                lNeedsSeed(i) = .TRUE.
            end if
        end do

        call SplitDegenerateSeedClasses(m, n, lNeedsSeed, dSeedFraction)

        if (lHasExistingSpeciesMoles) then
            dExistingFractionSum = SUM(dMolFraction(m:n), MASK=(.NOT.lNeedsSeed(m:n)))
            dSeedSum             = SUM(dMolFraction(m:n), MASK=lNeedsSeed(m:n))
            nSeeded              = COUNT(lNeedsSeed(m:n))
            dSeedTarget          = 0D0
            if ((nSeeded > 0).AND.(dSeedSum > 1D-300)) then
                dSeedTarget = DMIN1(dSeedSum, dTraceZeroSeed * DBLE(nSeeded))
                dSeedTarget = DMIN1(dSeedTarget, 5D-1)
                dSeedScale  = dSeedTarget / dSeedSum
                do i = m, n
                    if (lNeedsSeed(i)) dMolFraction(i) = dMolFraction(i) * dSeedScale
                end do
            else
                do i = m, n
                    if (lNeedsSeed(i)) dMolFraction(i) = 0D0
                end do
            end if

            if (dExistingFractionSum > 1D-300) then
                dSeedSum = SUM(dMolFraction(m:n), MASK=lNeedsSeed(m:n))
                do i = m, n
                    if (.NOT.lNeedsSeed(i)) then
                        dMolFraction(i) = dMolFraction(i) * (1D0 - dSeedSum) / dExistingFractionSum
                    end if
                end do
            end if
        else
            dSum = SUM(dMolFraction(m:n))
            if ((dSum > 1D-300).AND.(dSum < 1D300)) then
                dMolFraction(m:n) = dMolFraction(m:n) / dSum
            else
                dMolFraction(m:n) = 1D0 / DBLE(n - m + 1)
            end if
        end if

        do i = m, n
            dMolesSpecies(i) = dPhaseMoles * dMolFraction(i)
        end do

        lInitialized = .TRUE.
    end do

    if (allocated(dSeedFraction)) deallocate(dSeedFraction)
    if (allocated(lNeedsSeed)) deallocate(lNeedsSeed)

    do i = 1, nConPhases
        k = iAssemblage(i)
        if ((k >= 1).AND.(k <= nSpecies)) dMolesSpecies(k) = dMolesPhase(i)
    end do

    return

contains

    subroutine EstimateZeroEndmemberFraction(iSpecies, dFraction)
        implicit none

        integer, intent(in)  :: iSpecies
        real(8), intent(out) :: dFraction

        dArg = 0D0
        do j = 1, nElements
            dArg = dArg + dElementPotential(j) * dStoichSpecies(iSpecies,j)
        end do

        dParticleDenom = DBLE(iParticlesPerMole(iSpecies))
        if (dParticleDenom <= 0D0) dParticleDenom = 1D0

        dArg = dArg / dParticleDenom - &
            dLevelingChemicalPotential(iSpecies) * dLevelingSpeciesFormulaAtoms(iSpecies) / dParticleDenom

        if (ABS(dArg) <= dArgZeroTolerance) then
            dFraction = dTraceZeroSeed
        else
            dArg = DMAX1(DMIN1(dArg, dLogOverflowLimit), dLogUnderflowLimit)
            dFraction = DEXP(dArg)
        end if

        return
    end subroutine EstimateZeroEndmemberFraction

    subroutine SplitDegenerateSeedClasses(iFirst, iLast, lSeeded, dSeed)
        implicit none

        integer, intent(in) :: iFirst, iLast
        logical, intent(in), dimension(iFirst:iLast) :: lSeeded
        real(8), intent(in), dimension(iFirst:iLast) :: dSeed

        integer :: iClassCount, iMember, iRep
        real(8) :: dClassWeight
        logical, dimension(iFirst:iLast) :: lAssigned

        lAssigned = .FALSE.

        do iRep = iFirst, iLast
            if ((.NOT.lSeeded(iRep)).OR.lAssigned(iRep)) cycle

            iClassCount = 0
            dClassWeight = 0D0

            do iMember = iRep, iLast
                if ((.NOT.lSeeded(iMember)).OR.lAssigned(iMember)) cycle
                if (.NOT.AreLevelingSeedsEquivalent(iRep, iMember)) cycle

                iClassCount = iClassCount + 1
                dClassWeight = dClassWeight + dSeed(iMember)
            end do

            if (iClassCount <= 0) cycle

            dClassWeight = dClassWeight / DBLE(iClassCount)

            do iMember = iRep, iLast
                if ((.NOT.lSeeded(iMember)).OR.lAssigned(iMember)) cycle
                if (.NOT.AreLevelingSeedsEquivalent(iRep, iMember)) cycle

                dMolFraction(iMember) = dClassWeight / DBLE(iClassCount)
                lAssigned(iMember) = .TRUE.
            end do
        end do

        return
    end subroutine SplitDegenerateSeedClasses

    logical function AreLevelingSeedsEquivalent(iSpeciesA, iSpeciesB)
        implicit none

        integer, intent(in) :: iSpeciesA, iSpeciesB

        integer :: iElement
        real(8) :: dDenomA, dDenomB, dValueA, dValueB

        AreLevelingSeedsEquivalent = .FALSE.

        dDenomA = DBLE(iParticlesPerMole(iSpeciesA))
        dDenomB = DBLE(iParticlesPerMole(iSpeciesB))
        if (dDenomA <= 0D0) dDenomA = 1D0
        if (dDenomB <= 0D0) dDenomB = 1D0

        do iElement = 1, nElements
            dValueA = dStoichSpecies(iSpeciesA,iElement) / dDenomA
            dValueB = dStoichSpecies(iSpeciesB,iElement) / dDenomB
            if (ABS(dValueA - dValueB) > dClassCompositionTol) return
        end do

        dValueA = dLevelingSpeciesFormulaAtoms(iSpeciesA) / dDenomA
        dValueB = dLevelingSpeciesFormulaAtoms(iSpeciesB) / dDenomB
        if (ABS(dValueA - dValueB) > dClassCompositionTol) return

        dValueA = dLevelingChemicalPotential(iSpeciesA) * &
            dLevelingSpeciesFormulaAtoms(iSpeciesA) / dDenomA
        dValueB = dLevelingChemicalPotential(iSpeciesB) * &
            dLevelingSpeciesFormulaAtoms(iSpeciesB) / dDenomB
        if (ABS(dValueA - dValueB) > dClassPotentialTol) return

        AreLevelingSeedsEquivalent = .TRUE.

        return
    end function AreLevelingSeedsEquivalent

end subroutine PostLevelingSpeciesInit



subroutine SyncActiveSolutionFractionsFromSpecies

    USE ModuleThermo

    implicit none

    integer :: i, k, l, m, n, iSlot
    real(8) :: dPhaseMoles

    do l = 1, nSolnPhases
        iSlot = nElements - l + 1
        k = -iAssemblage(iSlot)
        if ((k < 1).OR.(k > nSolnPhasesSys)) cycle

        m = nSpeciesPhase(k-1) + 1
        n = nSpeciesPhase(k)
        dPhaseMoles = SUM(dMolesSpecies(m:n))
        dMolesPhase(iSlot) = dPhaseMoles

        if (dPhaseMoles > 0D0) then
            do i = m, n
                dMolFraction(i) = dMolesSpecies(i) / dPhaseMoles
            end do
        else
            dMolFraction(m:n) = 0D0
        end if
    end do

    do i = 1, nConPhases
        k = iAssemblage(i)
        if ((k >= 1).AND.(k <= nSpecies)) dMolesPhase(i) = dMolesSpecies(k)
    end do

    return

end subroutine SyncActiveSolutionFractionsFromSpecies
