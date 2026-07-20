!> \brief Write one solution pseudo-compound row for PEA Leveling.
!!
!! \details Converts a subminimized solution constitution into the Leveling
!! composition and Gibbs-energy basis.  The caller must state whether the
!! candidate is valid; invalid candidates keep their composition diagnostics
!! but receive a large positive potential so they cannot become active.
!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    SetLevelingSolutionCandidateRow.f90
!> \brief   Write one solution pseudo-compound row for PEA Leveling.
!> \author  S.Y. Kwon
!> \date    Jun. 26, 2026
!> \sa      CompInitMinSolnPoint.f90
!> \sa      CompMinSolnPoint.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Registered normalized typed candidate rows with preserved grid identity, exclusion state, and atomic capacity failure.
    !
    !
! Purpose:
! ========
!
!> \details This routine updates the pseudo-compound row used by PEA Leveling
!! for one solution phase.  The row stores the phase composition on the
!! Leveling basis, the atom-normalized chemical potential used by legacy GEM
!! arrays, and the formula-normalized chemical potential used by Leveling.  If
!! Subminimization did not converge or produce an accepted CEF witness, the row
!! remains available for diagnostics but is assigned a large positive
!! potential so it cannot be selected as a stable pseudo-compound.
!
!
! Required input variables:
! =========================
!
!> \param[in] iLevelRow        Absolute Leveling row to update.
!> \param[in] iSolnPhaseIndex  Absolute solution phase index for the row.
!> \param[in] dCandidate       Phase-local endmember fractions.
!> \param[in] lCandidateValid  True only for converged or accepted witness candidates.
!
! nSpeciesPhase                 Maps the phase to its endmember range.
! dChemicalPotential            Endmember chemical potentials after subminimization.
! dStoichSpecies                Full stoichiometry of each endmember.
! dLevelingSpeciesTotalAtoms    Normalization used by Leveling composition rows.
! dLevelingSpeciesFormulaAtoms  Formula-atom normalization used by Leveling potentials.
!
!
! Output/updated variables:
! =========================
!
! dLevelingCompositionSpecies   Updated composition of the pseudo-compound row.
! dStoichSpeciesLevel           Updated stoichiometry of the pseudo-compound row.
! iPhaseLevel                   Updated solution phase id of the pseudo-compound row.
! dChemicalPotential            Updated atom-normalized potential for iLevelRow.
! dLevelingChemicalPotential    Updated formula-normalized Leveling potential for iLevelRow.
! iLevelCandidateParentPhase    Parent phase identity for the retained candidate row.
! iLevelCandidateDisplayPhase   Display/reporting phase identity for the retained candidate row.
! iLevelCandidateIdentityOrdinal Composition-set ordinal for the retained candidate row.
! lLevelingRowExcluded           True when Leveling must not admit this row as physical evidence.
!
!
! Called subroutines/functions:
! =============================
!
! None.
!
!
! Primary callers:
! ================
!
! CompInitMinSolnPoint          Registers initial PEA solution minima.
! CompMinSolnPoint              Refreshes PEA solution minima during CheckPhaseAssemblage.
!
!
! Numerical assumptions:
! ======================
!
! - dCandidate is normalized by the caller as a phase-local endmember-fraction vector.
! - Candidate Gibbs energies and compositions must be put on the same basis as
!   Leveling before they are compared with real stoichiometric compounds.
! - Unknown or unconverged submin candidates are not active-set evidence; they
!   are assigned a high potential rather than being silently trusted.
!
!-------------------------------------------------------------------------------------------------------------



subroutine SetLevelingSolutionCandidateRow(iLevelRow, iSolnPhaseIndex, dCandidate, lCandidateValid)

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in) :: iLevelRow, iSolnPhaseIndex
    real(8), intent(in) :: dCandidate(*)
    logical, intent(in) :: lCandidateValid
    integer             :: i, j, m, n, iCand
    real(8)             :: dLevelingDenom, dLevelingFormulaAtomDenom, dCandidateAtomDenom, dCandidateGibbs
    real(8)             :: dCandidateStoich(nElements)
    logical             :: lRowEligible

    if ((iLevelRow < 1).OR.(iLevelRow > nSpeciesLevel)) return
    if ((iSolnPhaseIndex < 1).OR.(iSolnPhaseIndex > nSolnPhasesSys)) return

    m = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    n = nSpeciesPhase(iSolnPhaseIndex)

    ! Refuse capacity exhaustion before changing either the row or pool.  The
    ! caller may then abandon the complete candidate generation atomically.
    if (allocated(iLevelCandidateFromLevel)) then
        iCand = iLevelCandidateFromLevel(iLevelRow)
        if (((iCand < 1).OR.(iCand > nLevelCandidate)).AND.&
            (nLevelCandidate >= nLevelCandidateCapacity)) then
            INFOThermo = 42
            return
        end if
    end if

    dLevelingCompositionSpecies(iLevelRow,:) = 0D0
    dLevelingDenom = 0D0
    dLevelingFormulaAtomDenom = 0D0
    dCandidateAtomDenom = 0D0
    dCandidateGibbs = 0D0
    dCandidateStoich = 0D0
    lRowEligible = .FALSE.

    do i = m, n
        dCandidateAtomDenom = dCandidateAtomDenom + &
            dCandidate(i-m+1) * dSpeciesTotalAtoms(i) / DBLE(iParticlesPerMole(i))
        dLevelingDenom = dLevelingDenom + &
            dCandidate(i-m+1) * dLevelingSpeciesTotalAtoms(i) / DBLE(iParticlesPerMole(i))
        dLevelingFormulaAtomDenom = dLevelingFormulaAtomDenom + &
            dCandidate(i-m+1) * dLevelingSpeciesFormulaAtoms(i) / DBLE(iParticlesPerMole(i))
        dCandidateGibbs = dCandidateGibbs + dCandidate(i-m+1) * dChemicalPotential(i)
        do j = 1, nElements
            dLevelingCompositionSpecies(iLevelRow,j) = dLevelingCompositionSpecies(iLevelRow,j) + &
                dCandidate(i-m+1) * dStoichSpecies(i,j) / DBLE(iParticlesPerMole(i))
            dCandidateStoich(j) = dCandidateStoich(j) + &
                dCandidate(i-m+1) * dStoichSpecies(i,j) / DBLE(iParticlesPerMole(i))
        end do
    end do

    iPhaseLevel(iLevelRow) = iSolnPhaseIndex
    dStoichSpeciesLevel(iLevelRow,:) = dCandidateStoich

    lRowEligible = lCandidateValid.AND.(dLevelingDenom > 1D-300).AND.&
        (dCandidateAtomDenom > 1D-300).AND.(dLevelingFormulaAtomDenom > 1D-300)
    if (dLevelingDenom > 1D-300) then
        dLevelingCompositionSpecies(iLevelRow,:) = dLevelingCompositionSpecies(iLevelRow,:) / dLevelingDenom
        if (dCandidateAtomDenom > 1D-300) then
            dChemicalPotential(iLevelRow) = dCandidateGibbs / dCandidateAtomDenom
        else
            dChemicalPotential(iLevelRow) = 5D9
        end if
        if (dLevelingFormulaAtomDenom > 1D-300) then
            dLevelingChemicalPotential(iLevelRow) = dCandidateGibbs / dLevelingFormulaAtomDenom
        else
            dLevelingChemicalPotential(iLevelRow) = 5D9
        end if
    else
        dLevelingCompositionSpecies(iLevelRow,:) = 0D0
        dChemicalPotential(iLevelRow) = 5D9
        dLevelingChemicalPotential(iLevelRow) = 5D9
    end if

    if (.NOT.lRowEligible) then
        dChemicalPotential(iLevelRow) = 5D9
        dLevelingChemicalPotential(iLevelRow) = 5D9
    end if
    if (allocated(lLevelingRowExcluded)) then
        if (iLevelRow <= SIZE(lLevelingRowExcluded)) &
            lLevelingRowExcluded(iLevelRow) = .NOT.lRowEligible
    end if

    if (allocated(iLevelCandidateFromLevel).AND.allocated(iLevelCandidatePhase).AND.&
        allocated(iLevelCandidateSource).AND.allocated(iLevelCandidateSubMinStatus).AND.&
        allocated(iLevelCandidateParentPhase).AND.&
        allocated(iLevelCandidateDisplayPhase).AND.allocated(iLevelCandidateIdentityOrdinal).AND.&
        allocated(dLevelCandidateMolFraction)) then
        if ((iLevelRow >= 1).AND.(iLevelRow <= SIZE(iLevelCandidateFromLevel))) then
            iCand = iLevelCandidateFromLevel(iLevelRow)
            if ((iCand < 1).OR.(iCand > nLevelCandidate)) then
                if (nLevelCandidate < nLevelCandidateCapacity) then
                    nLevelCandidate = nLevelCandidate + 1
                    iCand = nLevelCandidate
                    iLevelCandidateFromLevel(iLevelRow) = iCand
                else
                    INFOThermo = 42
                    iCand = 0
                end if
            end if

            if ((iCand >= 1).AND.(iCand <= SIZE(iLevelCandidatePhase))) then
                iLevelCandidatePhase(iCand) = iSolnPhaseIndex
                iLevelCandidateSource(iCand) = LEVEL_CANDIDATE_SOURCE_SUBMIN
                if (lRowEligible) then
                    iLevelCandidateSubMinStatus(iCand) = SUBMIN_CANDIDATE_CONVERGED
                else
                    iLevelCandidateSubMinStatus(iCand) = SUBMIN_CANDIDATE_REJECTED
                end if
                iLevelCandidateParentPhase(iCand) = iSolnPhaseIndex
                iLevelCandidateDisplayPhase(iCand) = iSolnPhaseIndex
                iLevelCandidateIdentityOrdinal(iCand) = 1
                dLevelCandidateMolFraction(iCand,:) = 0D0
                do i = m, n
                    dLevelCandidateMolFraction(iCand,i) = dCandidate(i-m+1)
                end do
            end if
        end if
    end if

    return

end subroutine SetLevelingSolutionCandidateRow
