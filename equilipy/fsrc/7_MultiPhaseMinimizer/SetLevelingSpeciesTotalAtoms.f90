!
!
!
!> \brief Define Leveling normalization atoms for one system species.
!! \details Separates the atoms used for composition constraints from the
!!          formula-atom scaling used for Gibbs-plane chemical potentials.
subroutine SetLevelingSpeciesTotalAtoms(iSystemSpecies, iCSPhase, iCSSpecies)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    SetLevelingSpeciesTotalAtoms.f90
    !> \brief   Define Leveling composition and formula-atom scaling for one species.
    !> \author  S.Y. Kwon
    !> \date    Jun. 23, 2026
    !> \sa      CheckSystem.f90
    !> \sa      LevelingSolver.f90
    !> \sa      ApplyOrderDisorderLevelingPotentials.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/23/2026      S.Y. Kwon           Original Leveling atom-normalization helper
    !   06/23/2026      S.Y. Kwon           Compute Leveling atoms from non-vacancy occupied sublattice sites
    !   06/23/2026      S.Y. Kwon           Harden sublattice bounds and store occupied-site count as FormulaAtoms
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details This subroutine stores two Leveling normalization factors for
    !! each species.  dLevelingSpeciesTotalAtoms is derived from the active
    !! system stoichiometry and is used by composition and mass-balance terms.
    !! dLevelingSpeciesFormulaAtoms is the formula-atom denominator used to
    !! convert endmember Gibbs energies into the Leveling Gibbs plane.  For
    !! sublattice phases, this denominator counts occupied non-vacancy
    !! sublattice sites, so explicit vacancy sublattices do not artificially
    !! rescale the Leveling Gibbs plane.
    !
    !
    ! Required input variables:
    ! =========================
    !
    !> \param[in] iSystemSpecies  System species index.
    !> \param[in] iCSPhase        ChemSage solution phase index for this species.
    !> \param[in] iCSSpecies      ChemSage species index for this species.
    !
    ! dStoichSpecies        Species stoichiometric coefficients in the active system.
    ! iPhaseSublatticeCS    Mapping from ChemSage phase index to sublattice-model index.
    ! iConstituentSublatticeCS  Constituent indexes for each species on each sublattice.
    ! nSublatticePhaseCS    Number of sublattices for each sublattice-model index.
    ! dStoichSublatticeCS   Sublattice stoichiometric coefficients.
    ! cConstituentNameSUBCS  Constituent names used to identify vacancy sites.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! dLevelingSpeciesTotalAtoms  Species atom count used by Leveling composition constraints.
    ! dLevelingSpeciesFormulaAtoms  Occupied non-vacancy formula-atom denominator for Leveling.
    !
    !
    ! Important called routines:
    ! ==========================
    !
    ! ComputeOccupiedSublatticeFormulaAtoms  Internal helper that counts non-vacancy occupied sites.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! CheckSystem  Initializes Leveling normalization when system species are constructed.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Composition normalization and formula-atom scaling are stored
    !   separately.  The composition count follows dStoichSpecies, while the
    !   formula-atom denominator can be recomputed from sublattice occupancy.
    ! - Vacancy constituents do not contribute to dLevelingSpeciesFormulaAtoms.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleParseCS

    implicit none

    integer, intent(in) :: iSystemSpecies, iCSPhase, iCSSpecies
    integer             :: iSublPhaseIndex
    real(8)             :: dLevelingAtoms, dFormulaAtoms, dOccupiedFormulaAtoms

    if ((iSystemSpecies < 1).OR.(iSystemSpecies > nSpecies)) return

    dLevelingAtoms = SUM(ABS(dStoichSpecies(iSystemSpecies,1:nElements)))
    dFormulaAtoms = dLevelingAtoms

    if ((iCSPhase < 1).OR.(iCSPhase > nSolnPhasesSysCS)) then
        dLevelingSpeciesTotalAtoms(iSystemSpecies) = dLevelingAtoms
        dLevelingSpeciesFormulaAtoms(iSystemSpecies) = dFormulaAtoms
        return
    end if

    if (.NOT.allocated(iPhaseSublatticeCS)) then
        dLevelingSpeciesTotalAtoms(iSystemSpecies) = dLevelingAtoms
        dLevelingSpeciesFormulaAtoms(iSystemSpecies) = dFormulaAtoms
        return
    end if
    if (.NOT.allocated(nSublatticePhaseCS)) then
        dLevelingSpeciesTotalAtoms(iSystemSpecies) = dLevelingAtoms
        dLevelingSpeciesFormulaAtoms(iSystemSpecies) = dFormulaAtoms
        return
    end if

    if (iCSPhase > SIZE(iPhaseSublatticeCS)) then
        dLevelingSpeciesTotalAtoms(iSystemSpecies) = dLevelingAtoms
        dLevelingSpeciesFormulaAtoms(iSystemSpecies) = dFormulaAtoms
        return
    end if
    iSublPhaseIndex = iPhaseSublatticeCS(iCSPhase)
    if ((iSublPhaseIndex <= 0).OR.(iSublPhaseIndex > SIZE(nSublatticePhaseCS))) then
        dLevelingSpeciesTotalAtoms(iSystemSpecies) = dLevelingAtoms
        dLevelingSpeciesFormulaAtoms(iSystemSpecies) = dFormulaAtoms
        return
    end if

    call ComputeOccupiedSublatticeFormulaAtoms(iCSPhase, iCSSpecies, iSublPhaseIndex, dOccupiedFormulaAtoms)
    if (dOccupiedFormulaAtoms > 0D0) then
        dFormulaAtoms = dOccupiedFormulaAtoms
    end if

    dLevelingSpeciesTotalAtoms(iSystemSpecies) = dLevelingAtoms
    dLevelingSpeciesFormulaAtoms(iSystemSpecies) = dFormulaAtoms

    return

contains

    subroutine ComputeOccupiedSublatticeFormulaAtoms(iPhaseIndex, iSpeciesIndex, iSublPhaseIndex, dOccupiedFormulaAtoms)
        integer, intent(in) :: iPhaseIndex
        integer, intent(in) :: iSpeciesIndex
        integer, intent(in) :: iSublPhaseIndex
        real(8), intent(out):: dOccupiedFormulaAtoms

        integer             :: iConstituentIndex, iRelativeSpecies, iSublattice, nSublattices
        character(8)        :: cConstituent

        dOccupiedFormulaAtoms = 0D0

        if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSysCS)) return
        if ((iSpeciesIndex <= nSpeciesPhaseCS(iPhaseIndex)).OR. &
            (iSpeciesIndex > nSpeciesPhaseCS(iPhaseIndex+1))) return
        if (.NOT.allocated(iConstituentSublatticeCS)) return
        if (.NOT.allocated(cConstituentNameSUBCS)) return
        if (.NOT.allocated(dStoichSublatticeCS)) return
        if (iSublPhaseIndex > SIZE(iConstituentSublatticeCS, DIM=1)) return
        if (iSublPhaseIndex > SIZE(cConstituentNameSUBCS, DIM=1)) return
        if (iSublPhaseIndex > SIZE(dStoichSublatticeCS, DIM=1)) return

        iRelativeSpecies = iSpeciesIndex - nSpeciesPhaseCS(iPhaseIndex)
        if (iRelativeSpecies > SIZE(iConstituentSublatticeCS, DIM=3)) return

        nSublattices = nSublatticePhaseCS(iSublPhaseIndex)
        if (nSublattices < 1) return
        if (nSublattices > SIZE(iConstituentSublatticeCS, DIM=2)) return
        if (nSublattices > SIZE(cConstituentNameSUBCS, DIM=2)) return
        if (nSublattices > SIZE(dStoichSublatticeCS, DIM=2)) return

        do iSublattice = 1, nSublattices
            iConstituentIndex = iConstituentSublatticeCS(iSublPhaseIndex, iSublattice, iRelativeSpecies)
            if (iConstituentIndex <= 0) cycle
            if (iConstituentIndex > SIZE(cConstituentNameSUBCS, DIM=3)) return

            cConstituent = ADJUSTL(cConstituentNameSUBCS(iSublPhaseIndex, iSublattice, iConstituentIndex))
            if (IsVacancyConstituent(cConstituent)) cycle

            dOccupiedFormulaAtoms = dOccupiedFormulaAtoms + dStoichSublatticeCS(iSublPhaseIndex, iSublattice)
        end do

        return
    end subroutine ComputeOccupiedSublatticeFormulaAtoms

    logical function IsVacancyConstituent(cConstituent)
        character(*), intent(in) :: cConstituent

        IsVacancyConstituent = (TRIM(cConstituent) == 'VA').OR. &
            (TRIM(cConstituent) == 'Va').OR.(TRIM(cConstituent) == 'vA').OR. &
            (TRIM(cConstituent) == 'va')

        return
    end function IsVacancyConstituent

end subroutine SetLevelingSpeciesTotalAtoms
