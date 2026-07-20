!> \brief Initialize Leveling candidate handoff storage.
!!
!! \details Allocates identity and constitution records for both immutable grid
!! rows and dynamic subminimized rows.  Dynamic rows start after the grid-row
!! offset so the two candidate classes cannot overwrite one another.
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Added typed candidate identity, immutable-grid ownership, and Leveling-row exclusion state.
    !
! Numerical assumptions:
! ======================
!
! - nSpeciesLevel already includes every row that may be registered.
! - The row offset is reset here and assigned after grid rows are restored.
subroutine InitLevelCandidatePool(nCapacity)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in) :: nCapacity
    integer :: nPoolRows, nLevelRows

    ! Step 1. Release the previous solve's candidate transport arrays.
    if (allocated(iLevelCandidatePhase)) deallocate(iLevelCandidatePhase)
    if (allocated(iLevelCandidateSource)) deallocate(iLevelCandidateSource)
    if (allocated(iLevelCandidateSubMinStatus)) deallocate(iLevelCandidateSubMinStatus)
    if (allocated(iLevelCandidateStaticRow)) deallocate(iLevelCandidateStaticRow)
    if (allocated(iLevelCandidateParentPhase)) deallocate(iLevelCandidateParentPhase)
    if (allocated(iLevelCandidateDisplayPhase)) deallocate(iLevelCandidateDisplayPhase)
    if (allocated(iLevelCandidateIdentityOrdinal)) deallocate(iLevelCandidateIdentityOrdinal)
    if (allocated(iLevelCandidateFromLevel)) deallocate(iLevelCandidateFromLevel)
    if (allocated(dLevelCandidateMolFraction)) deallocate(dLevelCandidateMolFraction)
    if (allocated(lLevelingRowExcluded)) deallocate(lLevelingRowExcluded)

    ! Step 2. Size candidate identity storage to the complete Leveling row table.
    nLevelRows = MAX(1, nSpeciesLevel)
    nLevelCandidate = 0
    nLevelCandidateCapacity = MAX(MAX(0, nCapacity), nLevelRows)
    nLevelCandidateRowOffset = 0
    nPoolRows = MAX(1, nLevelCandidateCapacity)

    ! Step 3. Allocate and initialize typed row identity and constitution state.
    allocate(&
        iLevelCandidatePhase(nPoolRows),&
        iLevelCandidateSource(nPoolRows),&
        iLevelCandidateSubMinStatus(nPoolRows),&
        iLevelCandidateStaticRow(nPoolRows),&
        iLevelCandidateParentPhase(nPoolRows),&
        iLevelCandidateDisplayPhase(nPoolRows),&
        iLevelCandidateIdentityOrdinal(nPoolRows),&
        iLevelCandidateFromLevel(nLevelRows),&
        dLevelCandidateMolFraction(nPoolRows,nSpecies),&
        lLevelingRowExcluded(nLevelRows))

    iLevelCandidatePhase = 0
    iLevelCandidateSource = 0
    iLevelCandidateSubMinStatus = SUBMIN_CANDIDATE_UNKNOWN
    iLevelCandidateStaticRow = 0
    iLevelCandidateParentPhase = 0
    iLevelCandidateDisplayPhase = 0
    iLevelCandidateIdentityOrdinal = 0
    iLevelCandidateFromLevel = 0
    dLevelCandidateMolFraction = 0D0
    lLevelingRowExcluded = .FALSE.

    return

end subroutine InitLevelCandidatePool


!> \brief Temporarily extend Leveling thermodynamic arrays for sampled candidates.
subroutine ExtendSampledLevelingThermoArrays

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    real(8), dimension(:), allocatable :: dChemicalPotentialBase, dLevelingChemicalPotentialBase
    real(8), dimension(:,:), allocatable :: dAtomFractionSpeciesBase, dLevelingCompositionSpeciesBase

    lSampledLevelingThermoExtended = .FALSE.
    if (nSpeciesLevel <= nSpecies) return

    if (allocated(dChemicalPotentialOld)) deallocate(dChemicalPotentialOld)
    if (allocated(dAtomFractionSpeciesOld)) deallocate(dAtomFractionSpeciesOld)
    if (allocated(dLevelingChemicalPotentialOld)) deallocate(dLevelingChemicalPotentialOld)
    if (allocated(dLevelingCompositionSpeciesOld)) deallocate(dLevelingCompositionSpeciesOld)
    allocate(dChemicalPotentialOld(nSpecies), dAtomFractionSpeciesOld(nSpecies,nElements))
    allocate(dLevelingChemicalPotentialOld(nSpecies), dLevelingCompositionSpeciesOld(nSpecies,nElements))
    dChemicalPotentialOld = dChemicalPotential(:nSpecies)
    dAtomFractionSpeciesOld = dAtomFractionSpecies(:nSpecies,:)
    dLevelingChemicalPotentialOld = dLevelingChemicalPotential(:nSpecies)
    dLevelingCompositionSpeciesOld = dLevelingCompositionSpecies(:nSpecies,:)

    allocate(dChemicalPotentialBase(nSpecies), dAtomFractionSpeciesBase(nSpecies,nElements))
    allocate(dLevelingChemicalPotentialBase(nSpecies), dLevelingCompositionSpeciesBase(nSpecies,nElements))
    dChemicalPotentialBase = dChemicalPotential(:nSpecies)
    dAtomFractionSpeciesBase = dAtomFractionSpecies(:nSpecies,:)
    dLevelingChemicalPotentialBase = dLevelingChemicalPotential(:nSpecies)
    dLevelingCompositionSpeciesBase = dLevelingCompositionSpecies(:nSpecies,:)

    if (allocated(dChemicalPotential)) deallocate(dChemicalPotential)
    if (allocated(dAtomFractionSpecies)) deallocate(dAtomFractionSpecies)
    if (allocated(dLevelingChemicalPotential)) deallocate(dLevelingChemicalPotential)
    if (allocated(dLevelingCompositionSpecies)) deallocate(dLevelingCompositionSpecies)
    allocate(dChemicalPotential(nSpeciesLevel), dAtomFractionSpecies(nSpeciesLevel,nElements))
    allocate(dLevelingChemicalPotential(nSpeciesLevel), dLevelingCompositionSpecies(nSpeciesLevel,nElements))

    dChemicalPotential = 5D9
    dAtomFractionSpecies = 0D0
    dLevelingChemicalPotential = 5D9
    dLevelingCompositionSpecies = 0D0
    dChemicalPotential(:nSpecies) = dChemicalPotentialBase
    dAtomFractionSpecies(:nSpecies,:) = dAtomFractionSpeciesBase
    dLevelingChemicalPotential(:nSpecies) = dLevelingChemicalPotentialBase
    dLevelingCompositionSpecies(:nSpecies,:) = dLevelingCompositionSpeciesBase

    deallocate(dChemicalPotentialBase, dAtomFractionSpeciesBase)
    deallocate(dLevelingChemicalPotentialBase, dLevelingCompositionSpeciesBase)
    lSampledLevelingThermoExtended = .TRUE.

    return

end subroutine ExtendSampledLevelingThermoArrays


!> \brief Restore physical-species Leveling arrays after candidate transport.
subroutine RestoreSampledLevelingThermoArrays

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    real(8), dimension(:), allocatable :: dPhasePotentialBase

    if (.NOT.lSampledLevelingThermoExtended) return

    allocate(dPhasePotentialBase(nSpecies))
    dPhasePotentialBase = dPhasePotential(:nSpecies)

    if (allocated(dChemicalPotential)) deallocate(dChemicalPotential)
    if (allocated(dAtomFractionSpecies)) deallocate(dAtomFractionSpecies)
    if (allocated(dLevelingChemicalPotential)) deallocate(dLevelingChemicalPotential)
    if (allocated(dLevelingCompositionSpecies)) deallocate(dLevelingCompositionSpecies)
    if (allocated(dPhasePotential)) deallocate(dPhasePotential)

    allocate(dChemicalPotential(nSpecies), dAtomFractionSpecies(nSpecies,nElements), dPhasePotential(nSpecies))
    allocate(dLevelingChemicalPotential(nSpecies), dLevelingCompositionSpecies(nSpecies,nElements))
    dChemicalPotential = dChemicalPotentialOld
    dAtomFractionSpecies = dAtomFractionSpeciesOld
    dLevelingChemicalPotential = dLevelingChemicalPotentialOld
    dLevelingCompositionSpecies = dLevelingCompositionSpeciesOld
    dPhasePotential = dPhasePotentialBase

    deallocate(dPhasePotentialBase)
    lSampledLevelingThermoExtended = .FALSE.

    return

end subroutine RestoreSampledLevelingThermoArrays


!> Copy a retained Leveling candidate composition into the GEM composition array.
subroutine CopyLevelCandidateMolFraction(iLevelRow, iTargetRow, iSolnPhaseIndex, lFound)

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in) :: iLevelRow, iTargetRow, iSolnPhaseIndex
    logical, intent(out) :: lFound
    integer :: iCand, m, n

    lFound = .FALSE.
    if (.NOT.allocated(iLevelCandidateFromLevel)) return
    if (iLevelRow < 1 .OR. iLevelRow > SIZE(iLevelCandidateFromLevel)) return

    iCand = iLevelCandidateFromLevel(iLevelRow)
    if (iCand < 1 .OR. iCand > nLevelCandidate) return
    if (iLevelCandidatePhase(iCand) /= iSolnPhaseIndex) return

    m = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    n = nSpeciesPhase(iSolnPhaseIndex)
    dMolFractionGEM(iTargetRow,m:n) = dLevelCandidateMolFraction(iCand,m:n)
    lFound = .TRUE.

    return

end subroutine CopyLevelCandidateMolFraction
