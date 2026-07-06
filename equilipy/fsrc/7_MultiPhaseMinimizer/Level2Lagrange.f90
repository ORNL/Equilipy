!> \brief Translate selected Leveling/PEA candidates into Lagrangian GEM variables.
!!
!! \details Converts the selected Leveling/PEA rows into the fixed-active-set
!! Lagrangian layout used by RunLagrangianGEM.  The routine also records the
!! latest handoff identity, source, phase amount, and elemental-potential
!! snapshots so the embedded PEA-Lagrangian minimizer path can be audited.
!!
!-------------------------------------------------------------------------------------------------------------
!
!> \file    Level2Lagrange.f90
!> \brief   Translate selected Leveling/PEA candidates into Lagrangian GEM variables.
!> \date    Jun. 25, 2026
!> \sa      CheckPhaseAssemblage.f90
!> \sa      PostProcessPEA.f90
!> \sa      RunLagrangianGEM.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   06/25/2026      S.Y. Kwon           Added endmember-row compatibility checks before merging miscibility
!                                       copies into a base solution phase.
!   06/25/2026      S.Y. Kwon           Used the phase-mole numerical tolerance when pruning zero candidates
!                                       before Lagrangian GEM.
!   06/25/2026      S.Y. Kwon           Synchronized GEM active-set rows after Lagrangian translation so later
!                                       PEA repair starts from a consistent assemblage.
!   06/25/2026      S.Y. Kwon           Stored synchronized solution rows in PEA candidate-index form for
!                                       repair-stage CheckPhaseAssemblage.
!   06/26/2026      S.Y. Kwon           Recorded Leveling-to-Lagrangian handoff diagnostics without changing
!                                       candidate-selection logic.
!   06/26/2026      S.Y. Kwon           Preserved valid PEA phase-local solution minima when Leveling selects
!                                       an endmember row from that solution phase.
!   06/27/2026      S.Y. Kwon           Preserved PEA pseudo-solution candidates for disordered phases instead
!                                       of projecting true BCC_A2/FCC_A1 candidates into ordered partners.
!   06/27/2026      S.Y. Kwon           Projected only duplicate-composition DIS_PART pseudo-candidates while
!                                       preserving distinct disordered/ordered coexistence.
!   06/27/2026      S.Y. Kwon           Added explicit pseudo-compound merge policy and restricted immiscible
!                                       copies to database-provided miscibility phases.
!   06/28/2026      S.Y. Kwon           Preserved concrete Leveling endmember rows during solution-phase
!                                       merges so Level2Lagrange keeps the PEA mass-balanced row set.
!   07/01/2026      S.Y. Kwon           Rescaled solution-candidate amounts onto the actual Lagrangian
!                                       endmember-composition basis before storing species moles.
!   07/01/2026      S.Y. Kwon           Disabled disordered-to-ordered helper projection when the ordered
!                                       phase has no active ordering degree of freedom.
!   07/01/2026      S.Y. Kwon           Clarified that solution endmember vectors are handoff data, not
!                                       CEF Newton direction variables.
!   07/01/2026      S.Y. Kwon           Computed the initial handoff norm in the same CEF site-fraction
!                                       coordinate system used by RunLagrangianGEM.
!   07/02/2026      S.Y. Kwon           Filled passive active-slot model/display/composition-set diagnostics
!                                       during Lagrangian handoff synchronization.
!   07/02/2026      S.Y. Kwon           Stored passive per-active-slot constitutions and assigned identity
!                                       ordinals by display phase.
!   07/02/2026      S.Y. Kwon           Split passive active-slot thermodynamic parent identity from
!                                       display phase identity during Lagrangian handoff synchronization.
!   07/02/2026      S.Y. Kwon           Stored passive parent-model site fractions for active CEF
!                                       composition-set handoffs.
!   07/02/2026      S.Y. Kwon           Preserved distinct PEA solution rows as separate
!                                       Lagrangian composition-set slots during embedded polish.
!   07/02/2026      S.Y. Kwon           Limited active-slot parent remapping to actual DIS_PART helper
!                                       phases so real disordered phases remain independent.
!   07/03/2026      S.Y. Kwon           Added default-off c2 handoff logic that converts an unstable
!                                       DIS_PART helper candidate into a second composition-set slot
!                                       of the ordered SUBOM parent.
!   07/03/2026      S.Y. Kwon           Extended the switch-gated two-set handoff to the helper-first
!                                       ordering case where the DIS_PART candidate is selected before
!                                       the ordered SUBOM parent.
!   07/03/2026      S.Y. Kwon           Allowed unstable DIS_PART endpoint rows to seed a second
!                                       ordered-parent composition-set slot under the c2/c3 switch.
!   07/03/2026      S.Y. Kwon           Used the active order/disorder companion alias map for the
!                                       switch-gated two-set SUBOM handoff without changing switch-off
!                                       disordered-phase preservation.
!   07/03/2026      S.Y. Kwon           Prevented companion-to-ordered projection from assigning artificial
!                                       vacancy weight on non-fixed ordered substitutional sublattices.
!   07/03/2026      S.Y. Kwon           Routed switch-gated SUBOM two-set handoff through candidate-pool
!                                       identity rows instead of post-handoff slot injection.
!   07/03/2026      S.Y. Kwon           Counted direct same-parent SUBOM handoffs from synchronized active slots
!                                       instead of the retired post-handoff injection path.
!   07/03/2026      S.Y. Kwon           Removed retired post-handoff two-set helper routines after candidate-pool
!                                       handoff became the only active two-set path.
!
!
! Purpose:
! ========
!
!> \details Selected Leveling/PEA rows enter as stoichiometric species,
!! solution endmembers, or pseudo-compound solution candidates.  This routine
!! rebuilds the Lagrangian active set with stoichiometric phases in the front
!! of iAssemblage and solution phases as negative phase ids at the back.  For
!! CEF phases, the stored endmember vector is only the product-fraction
!! representation needed to hand off a candidate constitution.  The CEF
!! Lagrangian solver converts that vector to sublattice site fractions before
!! computing Newton directions.
!
!
! Required input variables:
! =========================
!
! iAssemblage(slot)       Species id, or nSpecies + solution phase id for a solution candidate.
! dMolesPhase(slot)       Candidate phase amount from Leveling/PEA.
! iPhaseGEM(slot)         Solution phase id associated with a selected candidate row.
! dMolFractionGEM(slot,:) Solution candidate composition in species/endmember coordinates.
!                         For CEF phases this is not a minimizer coordinate.
! dElementPotential       Elemental potentials from the incoming Leveling/PEA plane.
!
!
! Output/updated variables:
! =========================
!
! iAssemblage             Lagrangian assemblage; positive stoichiometric species in front and
!                         negative solution phase ids at the back.
! dMolesPhase             Phase amounts in the Lagrangian layout.
! dMolesSpecies           Species amounts for active stoichiometric and solution phases.
! dMolFraction            Phase-local solution mole fractions.
! dChemicalPotentialGEM   Synchronized GEM row potentials for later PEA repair.
! iLevel2Lagrange*        Latest handoff identity/source diagnostics.
! dLevel2Lagrange*        Latest handoff mole and elemental-potential diagnostics.
!
!
! Called subroutines/functions:
! =============================
!
! PostLevelingSpeciesInit       Applies zero-endmember protection before Lagrangian GEM.
! CompChemicalPotential         Recomputes chemical potentials for the translated active set.
! CompFunctionNorm              Computes the initial Lagrangian residual norm.
! SyncGEMRowsFromLagrangian     Converts the final Lagrangian layout back to GEM row arrays.
!
!
! Primary callers:
! ================
!
! PostProcessPEA                Converts PEA output into Lagrangian state.
! Python level2lagrange wrapper Used by focused handoff tests and diagnostics.
!
!
! Numerical assumptions:
! ======================
!
! - Zero or near-zero phase amounts are pruned before the Lagrangian layout is built.
! - Zero-endmember estimates are applied after the Leveling/PEA selection has been
!   translated into active solution phases.
! - Direct helper endmember rows are projected into their ordered partner only
!   when the ordered phase remains eligible in the screened active system.
!   PEA pseudo-solution helper rows are projected only when they duplicate an
!   already-active ordered candidate composition; distinct disordered/ordered
!   candidates remain separate phases.
! - Repeated pseudo-compound rows are merged unless the database provides
!   compatible miscibility-copy solution phases for that phase name.  During
!   embedded PEA polish, distinct solution rows are instead preserved
!   as separate composition-set slots because the PEA Leveling row set is the
!   mass-balanced active-set candidate.
! - Level2Lagrange may store CEF product endmember fractions in dMolFraction
!   and dMolFractionGEM for thermodynamic evaluation and later PEA handoff.
!   It must not construct CEF search directions; that is owned by
!   GEMNewtonCEF in site-fraction variables.
! - The diagnostic snapshots are passive; they must not control phase selection or
!   convergence decisions.
!
!-------------------------------------------------------------------------------------------------------------



subroutine Level2Lagrange
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
!
    integer                         :: i, j, k, m, n
    integer, dimension(nElements)    :: iSelectedAssemblage
    real(8), dimension(nElements)    :: dSelectedMolesPhase
    real(8), dimension(nSpecies)     :: dTempPhasePotential
    real(8)                         :: dPhasePruneTolerance, dTemp
    logical                         :: lCompEveryPhases, lPostLevelingInitialized
    logical                         :: lThermoArraysExtended
    logical                         :: lUseCEFSiteGEM
    logical                         :: OrderDisorderPhaseIsEligible
!
    lCompEveryPhases = .TRUE.
    lSolnPhases      = .FALSE.
    call ResetLevel2LagrangeCounters
!
    if (allocated(iAssemblageGEM)) deallocate(iAssemblageGEM)
    allocate(iAssemblageGEM(nElements))
    iAssemblageGEM = iAssemblage
!
!   Preserve the legacy small-negative cleanup for now.  A later active-set
!   repair pass should classify this instead of silently clipping it.
    if (MINVAL(dMolesPhase) < 0D0) then
        dMolesPhase(MINLOC(dMolesPhase, dim=1)) = 0D0
    end if
!
    dPhasePruneTolerance = DMAX1(1D-15, 10D0 * dTolerance(7))
    do i = 1, nElements
        if ((iAssemblage(i) /= 0).AND.(dMolesPhase(i) < dPhasePruneTolerance)) then
            nLevel2LagrangePruned = nLevel2LagrangePruned + 1
            iAssemblage(i) = 0
        end if
    end do
!
    iSelectedAssemblage = iAssemblage
    dSelectedMolesPhase = dMolesPhase
    call InitLevel2LagrangeDiagnostics(iSelectedAssemblage, dSelectedMolesPhase)
!
    if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)
    allocate(dMolesPhaseLast(nElements))
    dMolesPhaseLast = dSelectedMolesPhase
!
    dTempPhasePotential = 0D0
    nConPhases          = 0
    nSolnPhases         = 0
    iAssemblage         = 0
    dMolesPhase         = 0D0
    dMolesSpecies       = 0D0
    lSkipLagrange       = .FALSE.
!
    do i = 1, nElements
        call AddSelectedCandidate(i)
    end do
!
!   Restore arrays that sampled Leveling/PEA temporarily extended to
!   nSpeciesLevel rows.  Classical Leveling already has base-size arrays.
    lThermoArraysExtended = .FALSE.
    if (allocated(dChemicalPotential)) lThermoArraysExtended = SIZE(dChemicalPotential) > nSpecies
    if (lSampledLevelingThermoExtended.OR.lThermoArraysExtended) then
        dTempPhasePotential = dPhasePotential(:nSpecies)
        if (allocated(dAtomFractionSpecies)) deallocate(dAtomFractionSpecies)
        if (allocated(dLevelingCompositionSpecies)) deallocate(dLevelingCompositionSpecies)
        if (allocated(dPhasePotential))      deallocate(dPhasePotential)
        if (allocated(dChemicalPotential))   deallocate(dChemicalPotential)
        if (allocated(dLevelingChemicalPotential)) deallocate(dLevelingChemicalPotential)
        allocate(dAtomFractionSpecies(nSpecies,nElements), dLevelingCompositionSpecies(nSpecies,nElements),&
            dPhasePotential(nSpecies), dChemicalPotential(nSpecies), dLevelingChemicalPotential(nSpecies))
        dAtomFractionSpecies = dAtomFractionSpeciesOld
        dLevelingCompositionSpecies = dLevelingCompositionSpeciesOld
        dLevelingChemicalPotential = dLevelingChemicalPotentialOld
        dPhasePotential      = dTempPhasePotential
        lSampledLevelingThermoExtended = .FALSE.
    end if
!
!   Recalculate mole fractions from species moles.  This is the authoritative
!   Lagrangian starting composition and still only a handoff representation
!   for CEF phases; GEMNewtonCEF rebuilds site fractions before minimizing.
    do i = 1, nSolnPhases
        j = -iAssemblage(nElements - i + 1)
        m = nSpeciesPhase(j-1) + 1
        n = nSpeciesPhase(j)
        dTemp = SUM(dMolesSpecies(m:n))
        if (dTemp > 0D0) dMolFraction(m:n) = dMolesSpecies(m:n) / dTemp
    end do
!
!   Keep zero-endmember protection consistent with the direct leveling handoff.
    call PostLevelingSpeciesInit(lPostLevelingInitialized)
    do i = 1, nSolnPhases
        j = -iAssemblage(nElements - i + 1)
        m = nSpeciesPhase(j-1) + 1
        n = nSpeciesPhase(j)
        if (MINVAL(dMolesSpecies(m:n)) < 1D-100) lSkipLagrange = .TRUE.
    end do
!
    call CompChemicalPotential(lCompEveryPhases)
!
!   PEA repair and CEF residual evaluation use the GEM row arrays as the
!   active-slot handoff.  Synchronize before the first residual norm so
!   repeated CEF composition-set slots keep their own site fractions.
    call SyncGEMRowsFromLagrangian
    call UseCEFLagrangian(lUseCEFSiteGEM)
    lGEMCEFSiteLagrangianActive = lUseCEFSiteGEM
    call CompFunctionNorm
    dGEMFunctionNormLast = dGEMFunctionNorm
    call FinalizeLevel2LagrangeDiagnostics
!
    call ReleaseLevel2LagrangeWorkArrays
!
    return
!
contains
!
    subroutine ResetLevel2LagrangeCounters
!
        nLevel2LagrangeInput          = 0
        nLevel2LagrangePruned         = 0
        nLevel2LagrangeStoichSelected = 0
        nLevel2LagrangeSolnSelected   = 0
        nLevel2LagrangeSolnMerged     = 0
        nLevel2LagrangeSolnAdded      = 0
        nLevel2LagrangeOrderProjected = 0
        nLevel2LagrangeTwoSetCreated  = 0
        return
    end subroutine ResetLevel2LagrangeCounters
!
    subroutine InitLevel2LagrangeDiagnostics(iInputAssemblage, dInputMoles)
        integer, dimension(nElements), intent(in) :: iInputAssemblage
        real(8), dimension(nElements), intent(in) :: dInputMoles
!
        if (allocated(iLevel2LagrangeInputAssemblage))  deallocate(iLevel2LagrangeInputAssemblage)
        if (allocated(iLevel2LagrangeInputPhase))       deallocate(iLevel2LagrangeInputPhase)
        if (allocated(iLevel2LagrangeInputCandidate))   deallocate(iLevel2LagrangeInputCandidate)
        if (allocated(iLevel2LagrangeInputSource))      deallocate(iLevel2LagrangeInputSource)
        if (allocated(iLevel2LagrangeOutputAssemblage)) deallocate(iLevel2LagrangeOutputAssemblage)
        if (allocated(dLevel2LagrangeInputMoles))       deallocate(dLevel2LagrangeInputMoles)
        if (allocated(dLevel2LagrangeOutputMoles))      deallocate(dLevel2LagrangeOutputMoles)
        if (allocated(dLevel2LagrangeElementPotentialIn))  deallocate(dLevel2LagrangeElementPotentialIn)
        if (allocated(dLevel2LagrangeElementPotentialOut)) deallocate(dLevel2LagrangeElementPotentialOut)
!
        allocate(iLevel2LagrangeInputAssemblage(nElements),&
            iLevel2LagrangeInputPhase(nElements),&
            iLevel2LagrangeInputCandidate(nElements),&
            iLevel2LagrangeInputSource(nElements),&
            iLevel2LagrangeOutputAssemblage(nElements),&
            dLevel2LagrangeInputMoles(nElements),&
            dLevel2LagrangeOutputMoles(nElements),&
            dLevel2LagrangeElementPotentialIn(nElements),&
            dLevel2LagrangeElementPotentialOut(nElements))
!
        iLevel2LagrangeInputAssemblage  = iInputAssemblage
        iLevel2LagrangeInputPhase       = 0
        iLevel2LagrangeInputCandidate   = 0
        iLevel2LagrangeInputSource      = 0
        iLevel2LagrangeOutputAssemblage = 0
        dLevel2LagrangeInputMoles       = dInputMoles
        dLevel2LagrangeOutputMoles      = 0D0
        dLevel2LagrangeElementPotentialIn  = 0D0
        dLevel2LagrangeElementPotentialOut = 0D0
        nLevel2LagrangeInput = COUNT(iInputAssemblage /= 0)
!
        if (allocated(dElementPotential)) then
            if (SIZE(dElementPotential) >= nElements) then
                dLevel2LagrangeElementPotentialIn = dElementPotential(1:nElements)
            end if
        end if
!
        return
    end subroutine InitLevel2LagrangeDiagnostics
!
    subroutine RecordLevel2LagrangeInputCandidate(iSlot, iEntry, iSolnPhase)
        integer, intent(in) :: iSlot, iEntry, iSolnPhase
        integer             :: iCandidate
!
        if (.NOT.allocated(iLevel2LagrangeInputAssemblage)) return
        if ((iSlot < 1).OR.(iSlot > SIZE(iLevel2LagrangeInputAssemblage))) return
!
        iLevel2LagrangeInputPhase(iSlot) = iSolnPhase
!
        iCandidate = 0
        if (allocated(iLevelCandidateFromLevel)) then
            if ((iEntry >= 1).AND.(iEntry <= SIZE(iLevelCandidateFromLevel))) then
                iCandidate = iLevelCandidateFromLevel(iEntry)
            end if
        end if
!
        if ((iCandidate >= 1).AND.(iCandidate <= nLevelCandidate)) then
            iLevel2LagrangeInputCandidate(iSlot) = iCandidate
            iLevel2LagrangeInputSource(iSlot)    = iLevelCandidateSource(iCandidate)
        end if
!
        return
    end subroutine RecordLevel2LagrangeInputCandidate
!
    subroutine FinalizeLevel2LagrangeDiagnostics
!
        if (allocated(iLevel2LagrangeOutputAssemblage)) iLevel2LagrangeOutputAssemblage = iAssemblage
        if (allocated(dLevel2LagrangeOutputMoles))      dLevel2LagrangeOutputMoles      = dMolesPhase
!
        if (allocated(dLevel2LagrangeElementPotentialOut)) then
            if (allocated(dElementPotential)) then
                if (SIZE(dElementPotential) >= nElements) then
                    dLevel2LagrangeElementPotentialOut = dElementPotential(1:nElements)
                end if
            end if
        end if
!
        return
    end subroutine FinalizeLevel2LagrangeDiagnostics
!
    subroutine AddSelectedCandidate(iSlot)
        integer, intent(in) :: iSlot
        integer             :: iEntry, iSolnPhase
!
        iEntry = iSelectedAssemblage(iSlot)
        if (iEntry == 0) return
!
        iSolnPhase = CandidateSolutionPhase(iSlot, iEntry)
        call RecordLevel2LagrangeInputCandidate(iSlot, iEntry, iSolnPhase)
!
        if (iSolnPhase == 0) then
            call StoreStoichiometricCandidate(iEntry, dSelectedMolesPhase(iSlot))
        else if (iSolnPhase > 0) then
            call StoreSolutionCandidate(iSlot, iEntry, iSolnPhase, dSelectedMolesPhase(iSlot))
        end if
!
        return
    end subroutine AddSelectedCandidate

    integer function CandidateSolutionPhase(iSlot, iEntry)
        integer, intent(in) :: iSlot, iEntry
        integer             :: iCandidate
!
        CandidateSolutionPhase = iPhaseGEM(iSlot)
        iCandidate = LevelCandidateIndexForEntry(iEntry)
        if ((iCandidate >= 1).AND.(iCandidate <= nLevelCandidate)) then
            if (iLevelCandidateParentPhase(iCandidate) > 0) then
                CandidateSolutionPhase = iLevelCandidateParentPhase(iCandidate)
                return
            end if
        end if
        if ((CandidateSolutionPhase == 0) .AND. (iEntry > nSpecies)) then
            CandidateSolutionPhase = iEntry - nSpecies
        end if
        if ((CandidateSolutionPhase == 0) .AND. (iEntry > 0) .AND. (iEntry <= nSpecies)) then
            if (iPhase(iEntry) > 0) CandidateSolutionPhase = iPhase(iEntry)
        end if
!
        return
    end function CandidateSolutionPhase
!
    integer function LevelCandidateIndexForEntry(iEntry)
        integer, intent(in) :: iEntry
!
        LevelCandidateIndexForEntry = 0
        if (.NOT.allocated(iLevelCandidateFromLevel)) return
        if ((iEntry < 1).OR.(iEntry > SIZE(iLevelCandidateFromLevel))) return
        LevelCandidateIndexForEntry = iLevelCandidateFromLevel(iEntry)
        if ((LevelCandidateIndexForEntry < 1).OR.(LevelCandidateIndexForEntry > nLevelCandidate)) then
            LevelCandidateIndexForEntry = 0
        end if
!
        return
    end function LevelCandidateIndexForEntry
!
    logical function CandidateLevelRowHasExplicitCompositionSet(iEntry)
        integer, intent(in) :: iEntry
        integer             :: iCandidate
!
        CandidateLevelRowHasExplicitCompositionSet = .FALSE.
        iCandidate = LevelCandidateIndexForEntry(iEntry)
        if (iCandidate == 0) return
        if (.NOT.allocated(iLevelCandidateIdentityOrdinal)) return
        CandidateLevelRowHasExplicitCompositionSet = iLevelCandidateIdentityOrdinal(iCandidate) > 1
!
        return
    end function CandidateLevelRowHasExplicitCompositionSet
!
    subroutine StoreStoichiometricCandidate(iSpecies, dAmount)
        integer, intent(in) :: iSpecies
        real(8), intent(in) :: dAmount
!
        if (nConPhases + nSolnPhases >= nElements) return
!
        nConPhases                        = nConPhases + 1
        iAssemblage(nConPhases)           = iSpecies
        dMolesPhase(nConPhases)           = dAmount
        dMolesSpecies(iSpecies)           = dAmount
        nLevel2LagrangeStoichSelected     = nLevel2LagrangeStoichSelected + 1
!
        return
    end subroutine StoreStoichiometricCandidate
!
    subroutine StoreSolutionCandidate(iSlot, iEntry, iSolnPhase, dAmount)
        integer, intent(in) :: iSlot, iEntry, iSolnPhase
        real(8), intent(in) :: dAmount
        integer             :: iCandidateFirst, iCandidateLast
        integer             :: iExistingPhase, iExistingSlot, iExistingFirst, iExistingLast
        integer             :: iSoln, iTargetPhase, iBasePhase, iImmisciblePhase
        logical             :: lDuplicate, lMergePseudoCompounds, lPreservePEASolutionCandidate
!
        if (nConPhases + nSolnPhases >= nElements) return
        nLevel2LagrangeSolnSelected = nLevel2LagrangeSolnSelected + 1
!
        iCandidateFirst = nSpeciesPhase(iSolnPhase-1) + 1
        iCandidateLast  = nSpeciesPhase(iSolnPhase)
        iBasePhase      = MiscibleBasePhase(iSolnPhase)
        iTargetPhase    = 0
        lMergePseudoCompounds = .NOT.HasExplicitImmiscibilityCopies(iSolnPhase)
        lPreservePEASolutionCandidate = ShouldPreservePEASolutionCandidate().OR.&
            CandidateLevelRowHasExplicitCompositionSet(iEntry)
!
        iTargetPhase = OrderedPhaseForDisorderedHelper(iSolnPhase)
        if (iTargetPhase > 0) then
            if (OrderDisorderPhaseIsEligible(iTargetPhase)) then
                lDuplicate = .FALSE.
                if (iEntry > nSpecies) then
                    lDuplicate = DisorderedHelperDuplicatesOrderedCandidate(&
                        iSlot, iEntry, iSolnPhase, iTargetPhase, iCandidateFirst, iCandidateLast)
                end if
                if ((iEntry <= nSpecies).OR.lDuplicate) then
                    call StoreProjectedOrderDisorderCandidate(iSlot, iSolnPhase, iTargetPhase, dAmount)
                    return
                end if
            end if
        end if
!
        do iSoln = 1, nSolnPhases
            iExistingSlot  = nElements - iSoln + 1
            iExistingPhase = -iAssemblage(iExistingSlot)
!
            if (cSolnPhaseName(iExistingPhase) /= cSolnPhaseName(iSolnPhase)) cycle
!
            iExistingFirst = nSpeciesPhase(iExistingPhase-1) + 1
            iExistingLast  = nSpeciesPhase(iExistingPhase)
            if ((iExistingLast - iExistingFirst) /= (iCandidateLast - iCandidateFirst)) cycle
            if (.NOT.PhaseConstituentRowsCompatible(iSolnPhase, iExistingPhase)) cycle
            lDuplicate     = .FALSE.
!
            call AscertainMiscibilityGap(&
                iExistingPhase, iExistingFirst, iExistingLast, iExistingLast-iExistingFirst+1, &
                dMolFractionGEM(iSlot,iCandidateFirst:iCandidateLast), &
                dMolFraction(iExistingFirst:iExistingLast), &
                lDuplicate)
            if (lDuplicate) then
                if (SameParentCandidateConstitutionIsDistinct(&
                    iSlot, iEntry, iSolnPhase, iExistingPhase, &
                    iCandidateFirst, iCandidateLast, iExistingFirst, iExistingLast)) lDuplicate = .FALSE.
            end if
!
            if (lDuplicate) then
                call MergeIntoSolutionCandidate(&
                    iSlot, iExistingSlot, iEntry, dAmount, &
                    iCandidateFirst, iCandidateLast, iExistingFirst, iExistingLast)
            else if (lPreservePEASolutionCandidate) then
                call AddNewSolutionCandidate(&
                    iSlot, iEntry, iSolnPhase, dAmount, iCandidateFirst, iCandidateLast)
            else if (lMergePseudoCompounds) then
                call MergeIntoSolutionCandidate(&
                    iSlot, iExistingSlot, iEntry, dAmount, &
                    iCandidateFirst, iCandidateLast, iExistingFirst, iExistingLast)
            else
                iImmisciblePhase = NextImmisciblePhase(iBasePhase)
                if (iImmisciblePhase == 0) then
                    call MergeIntoSolutionCandidate(&
                        iSlot, iExistingSlot, iEntry, dAmount, &
                        iCandidateFirst, iCandidateLast, iExistingFirst, iExistingLast)
                    return
                end if
                iTargetPhase = iImmisciblePhase
                call AddNewSolutionCandidate(&
                    iSlot, iEntry, iTargetPhase, dAmount, iCandidateFirst, iCandidateLast)
            end if
            return
        end do
!
        if (lMiscibility(iSolnPhase).AND.&
            (cSolnPhaseName(iBasePhase) == cSolnPhaseName(iSolnPhase)).AND.&
            ((nSpeciesPhase(iBasePhase) - nSpeciesPhase(iBasePhase-1)) == &
             (iCandidateLast - iCandidateFirst + 1)).AND.&
            PhaseConstituentRowsCompatible(iSolnPhase, iBasePhase)) then
            iTargetPhase = iBasePhase
        else
            iTargetPhase = iSolnPhase
        end if
        call AddNewSolutionCandidate(iSlot, iEntry, iTargetPhase, dAmount, iCandidateFirst, iCandidateLast)
!
        return
    end subroutine StoreSolutionCandidate
!
    logical function ShouldPreservePEASolutionCandidate()
!
        ShouldPreservePEASolutionCandidate = .FALSE.
        if (.NOT.lPEALagrangianPolishActive) return
        if (nConPhases + nSolnPhases >= nElements) return
        ShouldPreservePEASolutionCandidate = .TRUE.
!
        return
    end function ShouldPreservePEASolutionCandidate
!
    subroutine MergeIntoSolutionCandidate(&
        iSlot, iTargetSlot, iEntry, dAmount, iCandidateFirst, iCandidateLast, iTargetFirst, iTargetLast)
        integer, intent(in) :: iSlot, iTargetSlot, iEntry, iCandidateFirst, iCandidateLast, iTargetFirst, iTargetLast
        real(8), intent(in) :: dAmount
        integer             :: iTargetSpecies, iCandidatePhase, nCandidateLocal
        logical             :: lUseCandidateFraction
        real(8)             :: dAdjustedAmount
        real(8), dimension(:), allocatable :: dCandidateFraction
!
        iCandidatePhase = 0
        if ((iEntry >= 1).AND.(iEntry <= nSpecies)) iCandidatePhase = iPhase(iEntry)
        if ((iCandidatePhase == 0).AND.(iCandidateFirst >= 1).AND.(iCandidateFirst <= nSpecies)) then
            iCandidatePhase = iPhase(iCandidateFirst)
        end if
        nCandidateLocal = iCandidateLast - iCandidateFirst + 1
        allocate(dCandidateFraction(nCandidateLocal))
        call GetSolutionCandidateFraction(&
            iSlot, iEntry, iCandidatePhase, iCandidateFirst, iCandidateLast, &
            dCandidateFraction, lUseCandidateFraction)
        dAdjustedAmount = RescaledSolutionCandidateAmount(&
            iSlot, iEntry, iCandidateFirst, iCandidateLast, &
            dCandidateFraction, lUseCandidateFraction, dAmount)
        dMolesPhase(iTargetSlot) = dMolesPhase(iTargetSlot) + dAdjustedAmount
!
        if (lUseCandidateFraction) then
            dMolesSpecies(iTargetFirst:iTargetLast) = dMolesSpecies(iTargetFirst:iTargetLast) + &
                dCandidateFraction * dAdjustedAmount
        else
            iTargetSpecies = MapCandidateSpeciesToTarget(iEntry, iCandidateFirst, iTargetFirst)
            dMolesSpecies(iTargetSpecies) = dMolesSpecies(iTargetSpecies) + dAdjustedAmount
        end if
!
        call UpdateSolutionMolFraction(iTargetSlot, iTargetFirst, iTargetLast)
        call StoreActiveSlotMolFractionFromGlobal(iTargetSlot, iTargetFirst, iTargetLast)
        deallocate(dCandidateFraction)
        nLevel2LagrangeSolnMerged = nLevel2LagrangeSolnMerged + 1
!
        return
    end subroutine MergeIntoSolutionCandidate
!
    subroutine AddNewSolutionCandidate(iSlot, iEntry, iTargetPhase, dAmount, iCandidateFirst, iCandidateLast)
        integer, intent(in) :: iSlot, iEntry, iTargetPhase, iCandidateFirst, iCandidateLast
        real(8), intent(in) :: dAmount
        integer             :: iResolvedTargetPhase, iTargetSlot, iTargetFirst, iTargetLast, iTargetSpecies
        integer             :: nCandidateLocal
        logical             :: lUseCandidateFraction
        real(8)             :: dAdjustedAmount
        real(8), dimension(:), allocatable :: dCandidateFraction
!
        if (nConPhases + nSolnPhases >= nElements) return
!
        iResolvedTargetPhase = iTargetPhase
        if (iEntry > nSpecies) then
            if ((nSpeciesPhase(iResolvedTargetPhase) - nSpeciesPhase(iResolvedTargetPhase-1)) /= &
                (iCandidateLast - iCandidateFirst + 1)) then
                if ((iCandidateFirst >= 1).AND.(iCandidateFirst <= nSpecies)) then
                    if (iPhase(iCandidateFirst) > 0) iResolvedTargetPhase = iPhase(iCandidateFirst)
                end if
            end if
        end if
!
        nSolnPhases = nSolnPhases + 1
        iTargetSlot = nElements - nSolnPhases + 1
        iTargetFirst = nSpeciesPhase(iResolvedTargetPhase-1) + 1
        iTargetLast  = nSpeciesPhase(iResolvedTargetPhase)
        nCandidateLocal = iCandidateLast - iCandidateFirst + 1
        allocate(dCandidateFraction(nCandidateLocal))
        call GetSolutionCandidateFraction(&
            iSlot, iEntry, iResolvedTargetPhase, iCandidateFirst, iCandidateLast, &
            dCandidateFraction, lUseCandidateFraction)
!
        iAssemblage(iTargetSlot) = -iResolvedTargetPhase
        nLevel2LagrangeSolnAdded = nLevel2LagrangeSolnAdded + 1
        dAdjustedAmount = RescaledSolutionCandidateAmount(&
            iSlot, iEntry, iCandidateFirst, iCandidateLast, &
            dCandidateFraction, lUseCandidateFraction, dAmount)
        dMolesPhase(iTargetSlot) = dAdjustedAmount
!
        if (lUseCandidateFraction) then
            dMolesSpecies(iTargetFirst:iTargetLast) = dCandidateFraction * dAdjustedAmount
        else
            iTargetSpecies = MapCandidateSpeciesToTarget(iEntry, iCandidateFirst, iTargetFirst)
            dMolesSpecies(iTargetSpecies) = dAdjustedAmount
        end if
!
        call UpdateSolutionMolFraction(iTargetSlot, iTargetFirst, iTargetLast)
        call StoreActiveSlotMolFractionFromCandidate(&
            iTargetSlot, iResolvedTargetPhase, iEntry, iCandidateFirst, iTargetFirst, iTargetLast, &
            dCandidateFraction, lUseCandidateFraction)
        if ((iResolvedTargetPhase >= 1) .AND. &
            (iResolvedTargetPhase <= SIZE(lSolnPhases))) lSolnPhases(iResolvedTargetPhase) = .TRUE.
        deallocate(dCandidateFraction)
!
        return
    end subroutine AddNewSolutionCandidate
!
    subroutine GetSolutionCandidateFraction(&
        iSlot, iEntry, iSolnPhase, iCandidateFirst, iCandidateLast, dFraction, lUseFraction)
        integer, intent(in) :: iSlot, iEntry, iSolnPhase, iCandidateFirst, iCandidateLast
        real(8), dimension(:), intent(out) :: dFraction
        logical, intent(out) :: lUseFraction
        integer             :: nCandidateLocal
        real(8)             :: dSum
!
        dFraction = 0D0
        lUseFraction = .FALSE.
        nCandidateLocal = iCandidateLast - iCandidateFirst + 1
        if (SIZE(dFraction) /= nCandidateLocal) return
!
        if (iEntry > nSpecies) then
            dFraction = MAX(dMolFractionGEM(iSlot,iCandidateFirst:iCandidateLast), 0D0)
            dSum = SUM(dFraction)
            if (dSum > 1D-300) then
                dFraction = dFraction / dSum
                lUseFraction = .TRUE.
            end if
            return
        end if
!
!       Concrete Leveling endmember/species rows must remain exact rows of the
!       mass-balanced Leveling solve.  Returning lUseFraction=.FALSE. makes the
!       caller map the selected row as one-hot instead of replacing it with the
!       latest phase-local subminimized composition.  If this row is CEF, the
!       one-hot vector is only the initial handoff from the Leveling plane.
!
        return
    end subroutine GetSolutionCandidateFraction
!
    logical function DisorderedHelperDuplicatesOrderedCandidate(&
        iSlot, iEntry, iHelperPhase, iOrderedPhase, iCandidateFirst, iCandidateLast)
        integer, intent(in) :: iSlot, iEntry, iHelperPhase, iOrderedPhase, iCandidateFirst, iCandidateLast
        integer             :: iSoln, iTargetSlot, iTargetFirst, iTargetLast, nCandidateLocal, nTargetLocal
        real(8)             :: dCompositionDiff
        real(8), dimension(nElements) :: dHelperComposition, dOrderedComposition
        real(8), dimension(:), allocatable :: dCandidateFraction, dOrderedFraction
        logical             :: lUseCandidateFraction
!
        DisorderedHelperDuplicatesOrderedCandidate = .FALSE.
        if (iEntry <= nSpecies) return
        if (iHelperPhase < 1 .OR. iOrderedPhase < 1) return
!
        iTargetSlot = 0
        do iSoln = 1, nSolnPhases
            if (-iAssemblage(nElements - iSoln + 1) == iOrderedPhase) then
                iTargetSlot = nElements - iSoln + 1
                exit
            end if
        end do
        if (iTargetSlot == 0) return
!
        iTargetFirst = nSpeciesPhase(iOrderedPhase-1) + 1
        iTargetLast  = nSpeciesPhase(iOrderedPhase)
        nCandidateLocal = iCandidateLast - iCandidateFirst + 1
        nTargetLocal    = iTargetLast - iTargetFirst + 1
        if ((nCandidateLocal <= 0).OR.(nTargetLocal <= 0)) return
!
        allocate(dCandidateFraction(nCandidateLocal), dOrderedFraction(nTargetLocal))
        call GetSolutionCandidateFraction(&
            iSlot, iEntry, iHelperPhase, iCandidateFirst, iCandidateLast, &
            dCandidateFraction, lUseCandidateFraction)
        if (.NOT.lUseCandidateFraction) then
            deallocate(dCandidateFraction, dOrderedFraction)
            return
        end if
!
        dOrderedFraction = MAX(dMolFraction(iTargetFirst:iTargetLast), 0D0)
        if (SUM(dOrderedFraction) <= 1D-300) then
            deallocate(dCandidateFraction, dOrderedFraction)
            return
        end if
        dOrderedFraction = dOrderedFraction / SUM(dOrderedFraction)
!
        call ComputeFractionComposition(iCandidateFirst, iCandidateLast, dCandidateFraction, dHelperComposition)
        call ComputeFractionComposition(iTargetFirst, iTargetLast, dOrderedFraction, dOrderedComposition)
        dCompositionDiff = SUM(DABS(dHelperComposition - dOrderedComposition)) / DFLOAT(MAX(1,nElements))
        if (dCompositionDiff <= 1D-8) DisorderedHelperDuplicatesOrderedCandidate = .TRUE.
!
        deallocate(dCandidateFraction, dOrderedFraction)
!
        return
    end function DisorderedHelperDuplicatesOrderedCandidate
!
    logical function SameParentCandidateConstitutionIsDistinct(&
        iSlot, iEntry, iSolnPhase, iExistingPhase, iCandidateFirst, iCandidateLast, &
        iExistingFirst, iExistingLast)
        integer, intent(in) :: iSlot, iEntry, iSolnPhase, iExistingPhase
        integer, intent(in) :: iCandidateFirst, iCandidateLast, iExistingFirst, iExistingLast
        integer             :: nCandidateLocal, nExistingLocal
        real(8)             :: dCandidateSum, dExistingSum, dConstitutionDiff
        real(8), dimension(:), allocatable :: dCandidateFraction, dExistingFraction
!
        SameParentCandidateConstitutionIsDistinct = .FALSE.
        if (.NOT.CandidateLevelRowHasExplicitCompositionSet(iEntry)) return
        if (iSolnPhase /= iExistingPhase) return
        if (iEntry <= nSpecies) return
!
        nCandidateLocal = iCandidateLast - iCandidateFirst + 1
        nExistingLocal  = iExistingLast - iExistingFirst + 1
        if ((nCandidateLocal <= 0).OR.(nCandidateLocal /= nExistingLocal)) return
!
        allocate(dCandidateFraction(nCandidateLocal), dExistingFraction(nExistingLocal))
        dCandidateFraction = DMAX1(dMolFractionGEM(iSlot,iCandidateFirst:iCandidateLast), 0D0)
        dExistingFraction = DMAX1(dMolFraction(iExistingFirst:iExistingLast), 0D0)
        dCandidateSum = SUM(dCandidateFraction)
        dExistingSum = SUM(dExistingFraction)
        if ((dCandidateSum <= 1D-300).OR.(dExistingSum <= 1D-300)) then
            deallocate(dCandidateFraction, dExistingFraction)
            return
        end if
        dCandidateFraction = dCandidateFraction / dCandidateSum
        dExistingFraction = dExistingFraction / dExistingSum
!
        dConstitutionDiff = SUM(DABS(dCandidateFraction - dExistingFraction)) / &
            DBLE(MAX(1,nCandidateLocal))
        if (dConstitutionDiff > 1D-8) SameParentCandidateConstitutionIsDistinct = .TRUE.
!
        deallocate(dCandidateFraction, dExistingFraction)
!
        return
    end function SameParentCandidateConstitutionIsDistinct
!
    subroutine ComputeFractionComposition(iFirst, iLast, dFraction, dComposition)
        integer, intent(in) :: iFirst, iLast
        real(8), dimension(:), intent(in) :: dFraction
        real(8), dimension(nElements), intent(out) :: dComposition
        integer :: iSpecies, iLocal, iElement
        real(8) :: dDenom
!
        dComposition = 0D0
        do iSpecies = iFirst, iLast
            iLocal = iSpecies - iFirst + 1
            if (iLocal > SIZE(dFraction)) exit
            do iElement = 1, nElements
                dComposition(iElement) = dComposition(iElement) + &
                    dFraction(iLocal) * dStoichSpecies(iSpecies,iElement) / &
                    DFLOAT(iParticlesPerMole(iSpecies))
            end do
        end do
!
        dDenom = SUM(DABS(dComposition))
        if (dDenom > 1D-300) dComposition = dComposition / dDenom
!
        return
    end subroutine ComputeFractionComposition
!
    real(8) function RescaledSolutionCandidateAmount(&
        iSlot, iEntry, iCandidateFirst, iCandidateLast, dFraction, lUseFraction, dAmount)
        integer, intent(in) :: iSlot, iEntry, iCandidateFirst, iCandidateLast
        real(8), dimension(:), intent(in) :: dFraction
        logical, intent(in) :: lUseFraction
        real(8), intent(in) :: dAmount
!
        integer :: iSpecies, iElement, iLocal
        real(8) :: dInputNorm, dStoredNorm
        real(8), dimension(nElements) :: dInputComposition, dStoredComposition
!
        RescaledSolutionCandidateAmount = dAmount
        if (dAmount <= 0D0) return
!
        call CandidateInputComposition(iSlot, iEntry, dInputComposition)
        dStoredComposition = 0D0
        if (lUseFraction) then
            do iSpecies = iCandidateFirst, iCandidateLast
                iLocal = iSpecies - iCandidateFirst + 1
                if (iLocal > SIZE(dFraction)) exit
                do iElement = 1, nElements
                    dStoredComposition(iElement) = dStoredComposition(iElement) + &
                        dFraction(iLocal) * dStoichSpecies(iSpecies,iElement) / &
                        DFLOAT(iParticlesPerMole(iSpecies))
                end do
            end do
        else if ((iEntry >= iCandidateFirst).AND.(iEntry <= iCandidateLast)) then
            do iElement = 1, nElements
                dStoredComposition(iElement) = dStoichSpecies(iEntry,iElement) / &
                    DFLOAT(iParticlesPerMole(iEntry))
            end do
        end if
!
        dInputNorm = SUM(DABS(dInputComposition))
        dStoredNorm = SUM(DABS(dStoredComposition))
        if ((dInputNorm > 1D-300).AND.(dStoredNorm > 1D-300)) then
            RescaledSolutionCandidateAmount = dAmount * dInputNorm / dStoredNorm
        end if
!
        return
    end function RescaledSolutionCandidateAmount
!
    real(8) function RescaledProjectedOrderDisorderAmount(iSlot, iTargetFirst, iTargetLast, &
        dFraction, dAmount)
        integer, intent(in) :: iSlot, iTargetFirst, iTargetLast
        real(8), dimension(:), intent(in) :: dFraction
        real(8), intent(in) :: dAmount
!
        integer :: iSpecies, iElement, iLocal
        real(8) :: dInputNorm, dStoredNorm
        real(8), dimension(nElements) :: dInputComposition, dStoredComposition
!
        RescaledProjectedOrderDisorderAmount = dAmount
        if (dAmount <= 0D0) return
!
        call CandidateInputComposition(iSlot, 0, dInputComposition)
        dStoredComposition = 0D0
        do iSpecies = iTargetFirst, iTargetLast
            iLocal = iSpecies - iTargetFirst + 1
            if (iLocal > SIZE(dFraction)) exit
            do iElement = 1, nElements
                dStoredComposition(iElement) = dStoredComposition(iElement) + &
                    dFraction(iLocal) * dStoichSpecies(iSpecies,iElement) / &
                    DFLOAT(iParticlesPerMole(iSpecies))
            end do
        end do
!
        dInputNorm = SUM(DABS(dInputComposition))
        dStoredNorm = SUM(DABS(dStoredComposition))
        if ((dInputNorm > 1D-300).AND.(dStoredNorm > 1D-300)) then
            RescaledProjectedOrderDisorderAmount = dAmount * dInputNorm / dStoredNorm
        end if
!
        return
    end function RescaledProjectedOrderDisorderAmount
!
    subroutine CandidateInputComposition(iSlot, iEntry, dComposition)
        integer, intent(in) :: iSlot, iEntry
        real(8), dimension(nElements), intent(out) :: dComposition
!
        dComposition = 0D0
        if ((iSlot >= 1).AND.(iSlot <= nElements)) then
            dComposition = dStoichSpeciesGEM(iSlot,1:nElements)
            if (SUM(DABS(dComposition)) > 1D-300) return
        end if
!
        if ((iEntry >= 1).AND.(iEntry <= SIZE(dLevelingCompositionSpecies,1))) then
            dComposition = dLevelingCompositionSpecies(iEntry,1:nElements)
        end if
!
        return
    end subroutine CandidateInputComposition
!
    subroutine StoreProjectedOrderDisorderCandidate(iSlot, iHelperPhase, iOrderedPhase, dAmount)
        integer, intent(in) :: iSlot, iHelperPhase, iOrderedPhase
        real(8), intent(in) :: dAmount
        integer             :: iSoln, iTargetSlot, iTargetFirst, iTargetLast, nTargetLocal
        logical             :: lCreated
        real(8)             :: dAdjustedAmount
        real(8), dimension(:), allocatable :: dProjectedFraction
!
        if (dAmount <= 0D0) return
!
        iTargetSlot = 0
        lCreated    = .FALSE.
        do iSoln = 1, nSolnPhases
            if (-iAssemblage(nElements - iSoln + 1) == iOrderedPhase) then
                iTargetSlot = nElements - iSoln + 1
                exit
            end if
        end do
!
        if (iTargetSlot == 0) then
            if (nConPhases + nSolnPhases >= nElements) return
            nSolnPhases = nSolnPhases + 1
            iTargetSlot = nElements - nSolnPhases + 1
            iAssemblage(iTargetSlot) = -iOrderedPhase
            dMolesPhase(iTargetSlot) = 0D0
            lCreated = .TRUE.
        end if
!
        iTargetFirst = nSpeciesPhase(iOrderedPhase-1) + 1
        iTargetLast  = nSpeciesPhase(iOrderedPhase)
        nTargetLocal = iTargetLast - iTargetFirst + 1
        if (nTargetLocal <= 0) return
!
        allocate(dProjectedFraction(nTargetLocal))
        call ProjectDisorderedHelperToOrderedFractions(iSlot, iHelperPhase, iOrderedPhase, dProjectedFraction)
        dAdjustedAmount = RescaledProjectedOrderDisorderAmount(&
            iSlot, iTargetFirst, iTargetLast, dProjectedFraction, dAmount)
!
        dMolesPhase(iTargetSlot) = dMolesPhase(iTargetSlot) + dAdjustedAmount
        dMolesSpecies(iTargetFirst:iTargetLast) = dMolesSpecies(iTargetFirst:iTargetLast) + &
            dAdjustedAmount * dProjectedFraction
        call UpdateSolutionMolFraction(iTargetSlot, iTargetFirst, iTargetLast)
        call StoreActiveSlotMolFractionFromFraction(&
            iTargetSlot, iOrderedPhase, iTargetFirst, iTargetLast, dProjectedFraction)
        if ((iOrderedPhase >= 1).AND.(iOrderedPhase <= SIZE(lSolnPhases))) lSolnPhases(iOrderedPhase) = .TRUE.
        nLevel2LagrangeOrderProjected = nLevel2LagrangeOrderProjected + 1
        if (lCreated) then
            nLevel2LagrangeSolnAdded = nLevel2LagrangeSolnAdded + 1
        else
            nLevel2LagrangeSolnMerged = nLevel2LagrangeSolnMerged + 1
        end if
!
        deallocate(dProjectedFraction)
!
        return
    end subroutine StoreProjectedOrderDisorderCandidate
!
    subroutine ProjectDisorderedHelperToOrderedFractions(iSlot, iHelperPhase, iOrderedPhase, dProjectedFraction)
        integer, intent(in) :: iSlot, iHelperPhase, iOrderedPhase
        real(8), dimension(:), intent(out) :: dProjectedFraction
!
        integer :: i, j, s, c, m, iHelperFirst, iHelperLast, iOrderedFirst, iOrderedLast
        integer :: iOrderedLocal, iSublPhase, iElement
        real(8) :: dHelperTotal, dNorm, dProduct, dX
        real(8), dimension(nElements) :: dHelperElementFraction
!
        dProjectedFraction = 0D0
        dHelperElementFraction = 0D0
!
        iHelperFirst = nSpeciesPhase(iHelperPhase-1) + 1
        iHelperLast  = nSpeciesPhase(iHelperPhase)
        do i = iHelperFirst, iHelperLast
            dX = dMolFractionGEM(iSlot,i)
            if (dX <= 0D0) dX = dMolFraction(i)
            do j = 1, nElements
                dHelperElementFraction(j) = dHelperElementFraction(j) + dX * DMAX1(dStoichSpecies(i,j), 0D0)
            end do
        end do
!
        dHelperTotal = SUM(dHelperElementFraction)
        if (dHelperTotal > 0D0) then
            dHelperElementFraction = dHelperElementFraction / dHelperTotal
        else
            dHelperElementFraction = 1D0 / DBLE(nElements)
        end if
!
        iSublPhase = iPhaseSublattice(iOrderedPhase)
        iOrderedFirst = nSpeciesPhase(iOrderedPhase-1) + 1
        iOrderedLast  = nSpeciesPhase(iOrderedPhase)
        do i = iOrderedFirst, iOrderedLast
            iOrderedLocal = i - iOrderedFirst + 1
            dProduct = 1D0
            do s = 1, nSublatticePhase(iSublPhase)
                m = iConstituentSublattice(iSublPhase,s,iOrderedLocal)
                if (m <= 0) cycle
                if (IsVacancyConstituent(cConstituentNameSUB(iSublPhase,s,m))) then
                    if (nConstituentSublattice(iSublPhase,s) > 1) dProduct = 0D0
                    cycle
                end if
!
                iElement = ConstituentElementIndex(cConstituentNameSUB(iSublPhase,s,m))
                if (iElement <= 0) then
                    dProduct = 0D0
                else
                    dProduct = dProduct * DMAX1(dHelperElementFraction(iElement), 0D0)
                end if
            end do
            dProjectedFraction(iOrderedLocal) = dProduct
        end do
!
        dNorm = SUM(dProjectedFraction)
        if (dNorm > 0D0) then
            dProjectedFraction = dProjectedFraction / dNorm
        else
            dProjectedFraction = 1D0 / DBLE(SIZE(dProjectedFraction))
        end if
!
        return
    end subroutine ProjectDisorderedHelperToOrderedFractions
!
    integer function OrderedPhaseForDisorderedHelper(iSolnPhase)
        integer, intent(in) :: iSolnPhase
        integer             :: iPhaseIndex
!
        OrderedPhaseForDisorderedHelper = 0
        if (.NOT.allocated(iDisorderedPhase)) return
!
        do iPhaseIndex = 1, nSolnPhasesSys
            if (iPhaseIndex > SIZE(iDisorderedPhase)) exit
            if (iDisorderedPhase(iPhaseIndex) == iSolnPhase) then
                if (TRIM(cSolnPhaseType(iPhaseIndex)) == 'SUBOM') then
                    OrderedPhaseForDisorderedHelper = iPhaseIndex
                    return
                end if
            end if
        end do
!
        return
    end function OrderedPhaseForDisorderedHelper
!
    integer function ConstituentElementIndex(cName)
        character(*), intent(in) :: cName
        integer                  :: iElement
        character(8)             :: cTarget, cElement
!
        ConstituentElementIndex = 0
        cTarget = UpperName(cName)
        do iElement = 1, nElements
            cElement = UpperName(cElementName(iElement))
            if (TRIM(cTarget) == TRIM(cElement)) then
                ConstituentElementIndex = iElement
                return
            end if
        end do
!
        return
    end function ConstituentElementIndex
!
    logical function IsVacancyConstituent(cName)
        character(*), intent(in) :: cName
        character(8)             :: cUpper
!
        cUpper = UpperName(cName)
        IsVacancyConstituent = TRIM(cUpper) == 'VA'
!
        return
    end function IsVacancyConstituent
!
    character(8) function UpperName(cName)
        character(*), intent(in) :: cName
        integer                  :: iChar, iCode, nChar
!
        UpperName = ' '
        nChar = MIN(LEN_TRIM(cName), LEN(UpperName))
        do iChar = 1, nChar
            iCode = IACHAR(cName(iChar:iChar))
            if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                UpperName(iChar:iChar) = ACHAR(iCode - 32)
            else
                UpperName(iChar:iChar) = cName(iChar:iChar)
            end if
        end do
!
        return
    end function UpperName
!
    integer function MapCandidateSpeciesToTarget(iSpecies, iCandidateFirst, iTargetFirst)
        integer, intent(in) :: iSpecies, iCandidateFirst, iTargetFirst
!
        MapCandidateSpeciesToTarget = iTargetFirst + iSpecies - iCandidateFirst
!
        return
    end function MapCandidateSpeciesToTarget
!
    subroutine UpdateSolutionMolFraction(iTargetSlot, iFirst, iLast)
        integer, intent(in) :: iTargetSlot, iFirst, iLast
        real(8)             :: dPhaseAmount
!
        dPhaseAmount = dMolesPhase(iTargetSlot)
        if (dPhaseAmount > 0D0) then
            dMolFraction(iFirst:iLast) = dMolesSpecies(iFirst:iLast) / dPhaseAmount
        end if
!
        return
    end subroutine UpdateSolutionMolFraction
!
    subroutine StoreActiveSlotMolFractionFromGlobal(iSlotIn, iFirstIn, iLastIn)
        integer, intent(in) :: iSlotIn, iFirstIn, iLastIn
        real(8)             :: dSumLocal
!
        if (.NOT.allocated(dActiveSlotMolFraction)) return
        if ((iSlotIn < 1).OR.(iSlotIn > SIZE(dActiveSlotMolFraction,1))) return
        if ((iFirstIn < 1).OR.(iLastIn > SIZE(dActiveSlotMolFraction,2))) return
!
        dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) = DMAX1(dMolFraction(iFirstIn:iLastIn), 0D0)
        dSumLocal = SUM(dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn))
        if (dSumLocal > 1D-300) dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) = &
            dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) / dSumLocal
!
        return
    end subroutine StoreActiveSlotMolFractionFromGlobal
!
    subroutine StoreActiveSlotMolFractionFromCandidate(&
        iSlotIn, iSolnPhaseIn, iEntryIn, iCandidateFirstIn, iFirstIn, iLastIn, dFractionIn, lUseFractionIn)
        integer, intent(in) :: iSlotIn, iSolnPhaseIn, iEntryIn, iCandidateFirstIn, iFirstIn, iLastIn
        real(8), dimension(:), intent(in) :: dFractionIn
        logical, intent(in) :: lUseFractionIn
        integer             :: iTargetSpecies
        real(8), dimension(:), allocatable :: dFractionLocal
!
        allocate(dFractionLocal(iLastIn-iFirstIn+1))
        dFractionLocal = 0D0
        if (lUseFractionIn) then
            dFractionLocal = dFractionIn
        else
            iTargetSpecies = MapCandidateSpeciesToTarget(iEntryIn, iCandidateFirstIn, iFirstIn)
            if ((iTargetSpecies >= iFirstIn).AND.(iTargetSpecies <= iLastIn)) then
                dFractionLocal(iTargetSpecies-iFirstIn+1) = 1D0
            end if
        end if
        call StoreActiveSlotMolFractionFromFraction(iSlotIn, iSolnPhaseIn, iFirstIn, iLastIn, dFractionLocal)
        deallocate(dFractionLocal)
!
        return
    end subroutine StoreActiveSlotMolFractionFromCandidate
!
    subroutine StoreActiveSlotMolFractionFromFraction(iSlotIn, iSolnPhaseIn, iFirstIn, iLastIn, dFractionIn)
        integer, intent(in) :: iSlotIn, iSolnPhaseIn, iFirstIn, iLastIn
        real(8), dimension(:), intent(in) :: dFractionIn
        real(8)             :: dSumLocal
!
        if (.NOT.allocated(dActiveSlotMolFraction)) return
        if ((iSlotIn < 1).OR.(iSlotIn > SIZE(dActiveSlotMolFraction,1))) return
        if ((iFirstIn < 1).OR.(iLastIn > SIZE(dActiveSlotMolFraction,2))) return
        if (SIZE(dFractionIn) /= (iLastIn - iFirstIn + 1)) return
        if ((iSolnPhaseIn < 1).OR.(iSolnPhaseIn > nSolnPhasesSys)) return
!
        dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) = DMAX1(dFractionIn, 0D0)
        dSumLocal = SUM(dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn))
        if (dSumLocal > 1D-300) then
            dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) = &
                dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) / dSumLocal
        else
            dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) = 0D0
        end if
!
        return
    end subroutine StoreActiveSlotMolFractionFromFraction
!
    subroutine SyncGEMRowsFromLagrangian
        integer :: iSlot, iEntry, iSolnPhase, iFirst, iLast, iSpecies, iElement
        integer :: iPrevSlot, iThermoPhase
        real(8) :: dStoichDenom, dFormulaDenom
        real(8) :: dSlotFractionSum
        logical :: lHadSlotMolFraction
        real(8), dimension(:,:), allocatable :: dSlotMolFractionSave
        real(8), dimension(:), allocatable :: dSlotFraction
!
        lHadSlotMolFraction = allocated(dActiveSlotMolFraction)
        if (lHadSlotMolFraction) dSlotMolFractionSave = dActiveSlotMolFraction
!
        if (allocated(iAssemblageGEM)) iAssemblageGEM = 0
!
        iPhaseGEM = 0
        dStoichSpeciesGEM = 0D0
        dAtomFractionSpeciesGEM = 0D0
        dChemicalPotentialGEM = 0D0
        dMolFractionGEM = 0D0
        if (allocated(dActiveSlotMolFraction)) dActiveSlotMolFraction = 0D0
        if (allocated(dActiveSlotSiteFraction)) dActiveSlotSiteFraction = 0D0
        if (allocated(iActiveSlotThermoPhase)) iActiveSlotThermoPhase = 0
        if (allocated(iActiveSlotDisplayPhase)) iActiveSlotDisplayPhase = 0
        if (allocated(iActiveSlotIdentityOrdinal)) iActiveSlotIdentityOrdinal = 0
!
        do iSlot = 1, nElements
            iEntry = iAssemblage(iSlot)
            if (iEntry == 0) cycle
!
            if (iEntry > 0) then
                if (allocated(iAssemblageGEM)) iAssemblageGEM(iSlot) = iEntry
                dStoichSpeciesGEM(iSlot,:) = dStoichSpecies(iEntry,:)
                dAtomFractionSpeciesGEM(iSlot,:) = dLevelingCompositionSpecies(iEntry,:)
                dChemicalPotentialGEM(iSlot) = dLevelingChemicalPotential(iEntry)
                cycle
            end if
!
            iSolnPhase = -iEntry
            iThermoPhase = ActiveSlotThermoPhaseForDisplay(iSolnPhase)
            if (allocated(iAssemblageGEM)) iAssemblageGEM(iSlot) = nSpecies + iSolnPhase
            iPhaseGEM(iSlot) = iSolnPhase
            if (allocated(iActiveSlotThermoPhase)) iActiveSlotThermoPhase(iSlot) = iThermoPhase
            if (allocated(iActiveSlotDisplayPhase)) iActiveSlotDisplayPhase(iSlot) = iSolnPhase
            if (allocated(iActiveSlotIdentityOrdinal)) then
                iActiveSlotIdentityOrdinal(iSlot) = 1
                do iPrevSlot = 1, iSlot - 1
                    if (iActiveSlotThermoPhase(iPrevSlot) == iThermoPhase) then
                        iActiveSlotIdentityOrdinal(iSlot) = iActiveSlotIdentityOrdinal(iSlot) + 1
                    end if
                end do
                if (lSUBOMTwoSetCandidateEnabled.AND.&
                    (TRIM(cSolnPhaseType(iThermoPhase)) == 'SUBOM').AND.&
                    (iActiveSlotIdentityOrdinal(iSlot) > 1)) then
                    nLevel2LagrangeTwoSetCreated = nLevel2LagrangeTwoSetCreated + 1
                end if
            end if
            iFirst = nSpeciesPhase(iSolnPhase-1) + 1
            iLast  = nSpeciesPhase(iSolnPhase)
            allocate(dSlotFraction(iLast-iFirst+1))
            dSlotFraction = DMAX1(dMolFraction(iFirst:iLast), 0D0)
            if (lHadSlotMolFraction) then
                if ((iSlot >= 1).AND.(iSlot <= SIZE(dSlotMolFractionSave,1)).AND.&
                    (iLast <= SIZE(dSlotMolFractionSave,2))) then
                    dSlotFractionSum = SUM(dSlotMolFractionSave(iSlot,iFirst:iLast))
                    if (dSlotFractionSum > 1D-300) then
                        dSlotFraction = DMAX1(dSlotMolFractionSave(iSlot,iFirst:iLast), 0D0)
                    end if
                end if
            end if
            dSlotFractionSum = SUM(dSlotFraction)
            if (dSlotFractionSum > 1D-300) dSlotFraction = dSlotFraction / dSlotFractionSum
            dMolFraction(iFirst:iLast) = dSlotFraction
            ! Store the product-fraction row for later PEA handoff.  CEF
            ! Lagrangian directions continue to be rebuilt from site fractions.
            dMolFractionGEM(iSlot,iFirst:iLast) = dSlotFraction
            if (allocated(dActiveSlotMolFraction)) then
                dActiveSlotMolFraction(iSlot,iFirst:iLast) = dSlotFraction
            end if
            call StoreActiveSlotSiteFraction(iSlot, iSolnPhase, iThermoPhase)
            if (lSUBOMTwoSetCandidateEnabled.AND.(TRIM(cSolnPhaseType(iThermoPhase)) == 'SUBOM')) then
                call RecordSUBOMTwoSetTrace(SUBOM_TRACE_HANDOFF, iThermoPhase, &
                    iActiveSlotIdentityOrdinal(iSlot), iSlot, SIZE(dSlotFraction), dSlotFraction)
            end if
!
            do iSpecies = iFirst, iLast
                do iElement = 1, nElements
                    dStoichSpeciesGEM(iSlot,iElement) = dStoichSpeciesGEM(iSlot,iElement) + &
                        dSlotFraction(iSpecies-iFirst+1) * dStoichSpecies(iSpecies,iElement) / &
                        DFLOAT(iParticlesPerMole(iSpecies))
                end do
            end do
!
            dStoichDenom = SUM(dStoichSpeciesGEM(iSlot,:))
            if (dStoichDenom > 1D-300) then
                dAtomFractionSpeciesGEM(iSlot,:) = dStoichSpeciesGEM(iSlot,:) / dStoichDenom
            end if
!
            dFormulaDenom = SUM(dSlotFraction * &
                dLevelingSpeciesFormulaAtoms(iFirst:iLast) / DFLOAT(iParticlesPerMole(iFirst:iLast)))
            if (dFormulaDenom > 1D-300) then
                dChemicalPotentialGEM(iSlot) = &
                    DOT_PRODUCT(dSlotFraction, dChemicalPotential(iFirst:iLast)) / dFormulaDenom
            else
                dChemicalPotentialGEM(iSlot) = 5D9
            end if
            deallocate(dSlotFraction)
        end do
!
        if (allocated(dSlotMolFractionSave)) deallocate(dSlotMolFractionSave)
!
        return
    end subroutine SyncGEMRowsFromLagrangian

    integer function ActiveSlotThermoPhaseForDisplay(iDisplayPhase)
        integer, intent(in) :: iDisplayPhase
        integer             :: iOrderedPhase
!
        ActiveSlotThermoPhaseForDisplay = iDisplayPhase
        if (.NOT.allocated(iDisorderedPhase)) return
!
        do iOrderedPhase = 1, MIN(nSolnPhasesSys, SIZE(iDisorderedPhase))
            if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') cycle
            if (iDisorderedPhase(iOrderedPhase) == iDisplayPhase) then
                ActiveSlotThermoPhaseForDisplay = iOrderedPhase
                return
            end if
        end do
!
        return
    end function ActiveSlotThermoPhaseForDisplay

    subroutine StoreActiveSlotSiteFraction(iSlotIn, iDisplayPhaseIn, iThermoPhaseIn)
        integer, intent(in) :: iSlotIn, iDisplayPhaseIn, iThermoPhaseIn
        integer             :: iDisplayPhaseID, iThermoPhaseID
        real(8)             :: dDisplaySite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8)             :: dParentSite(nMaxSublatticeSys,nMaxConstituentSys)
!
        if (.NOT.allocated(dActiveSlotSiteFraction)) return
        if (.NOT.allocated(iPhaseSublattice)) return
        if ((iDisplayPhaseIn <= 0).OR.(iDisplayPhaseIn > nSolnPhasesSys)) return
        if ((iThermoPhaseIn <= 0).OR.(iThermoPhaseIn > nSolnPhasesSys)) return
!
        iDisplayPhaseID = iPhaseSublattice(iDisplayPhaseIn)
        iThermoPhaseID = iPhaseSublattice(iThermoPhaseIn)
        if ((iDisplayPhaseID <= 0).OR.(iThermoPhaseID <= 0)) return
!
        dDisplaySite = 0D0
        dParentSite = 0D0
        call BuildSiteFractionForPhase(iDisplayPhaseIn, iDisplayPhaseID, dDisplaySite)
        if (iDisplayPhaseIn == iThermoPhaseIn) then
            dParentSite = dDisplaySite
        else
            call MapDisplaySiteToParent(iDisplayPhaseID, iThermoPhaseID, dDisplaySite, dParentSite)
        end if
        call NormalizeSiteFractionForPhase(iThermoPhaseID, dParentSite)
        dActiveSlotSiteFraction(iSlotIn,:,:) = dParentSite
!
        return
    end subroutine StoreActiveSlotSiteFraction

    subroutine BuildSiteFractionForPhase(iSolnPhaseIn, iPhaseIDIn, dSiteOut)
        integer, intent(in)  :: iSolnPhaseIn, iPhaseIDIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)
        integer              :: iFirstLocal, iLastLocal, iSpeciesLocal, iLocal
        integer              :: iSub, iCon
!
        dSiteOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIn)
        do iSpeciesLocal = iFirstLocal, iLastLocal
            iLocal = iSpeciesLocal - iFirstLocal + 1
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                    DMAX1(dMolFraction(iSpeciesLocal), 0D0)
            end do
        end do
        call NormalizeSiteFractionForPhase(iPhaseIDIn, dSiteOut)
!
        return
    end subroutine BuildSiteFractionForPhase

    subroutine MapDisplaySiteToParent(iDisplayPhaseIDIn, iParentPhaseIDIn, dDisplaySiteIn, dParentSiteOut)
        integer, intent(in)  :: iDisplayPhaseIDIn, iParentPhaseIDIn
        real(8), intent(in)  :: dDisplaySiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dParentSiteOut(nMaxSublatticeSys,nMaxConstituentSys)
        integer              :: iParentSub, iParentCon, iDisplaySub, iDisplayCon
        integer              :: iBestDisplaySub, iBestCount
        logical              :: lSubset, lFound
!
        dParentSiteOut = 0D0
        do iParentSub = 1, nSublatticePhase(iParentPhaseIDIn)
            iBestDisplaySub = 0
            iBestCount = nMaxConstituentSys + 1
            do iDisplaySub = 1, nSublatticePhase(iDisplayPhaseIDIn)
                lSubset = .TRUE.
                do iParentCon = 1, nConstituentSublattice(iParentPhaseIDIn,iParentSub)
                    lFound = .FALSE.
                    do iDisplayCon = 1, nConstituentSublattice(iDisplayPhaseIDIn,iDisplaySub)
                        if (TRIM(cConstituentNameSUB(iParentPhaseIDIn,iParentSub,iParentCon)) == &
                            TRIM(cConstituentNameSUB(iDisplayPhaseIDIn,iDisplaySub,iDisplayCon))) then
                            lFound = .TRUE.
                            exit
                        end if
                    end do
                    if (.NOT.lFound) then
                        lSubset = .FALSE.
                        exit
                    end if
                end do
                if (lSubset.AND.(nConstituentSublattice(iDisplayPhaseIDIn,iDisplaySub) < iBestCount)) then
                    iBestDisplaySub = iDisplaySub
                    iBestCount = nConstituentSublattice(iDisplayPhaseIDIn,iDisplaySub)
                end if
            end do
            if (iBestDisplaySub <= 0) return
!
            do iParentCon = 1, nConstituentSublattice(iParentPhaseIDIn,iParentSub)
                do iDisplayCon = 1, nConstituentSublattice(iDisplayPhaseIDIn,iBestDisplaySub)
                    if (TRIM(cConstituentNameSUB(iParentPhaseIDIn,iParentSub,iParentCon)) == &
                        TRIM(cConstituentNameSUB(iDisplayPhaseIDIn,iBestDisplaySub,iDisplayCon))) then
                        dParentSiteOut(iParentSub,iParentCon) = dDisplaySiteIn(iBestDisplaySub,iDisplayCon)
                        exit
                    end if
                end do
            end do
        end do
!
        return
    end subroutine MapDisplaySiteToParent

    subroutine NormalizeSiteFractionForPhase(iPhaseIDIn, dSiteInOut)
        integer, intent(in)    :: iPhaseIDIn
        real(8), intent(inout) :: dSiteInOut(nMaxSublatticeSys,nMaxConstituentSys)
        integer                :: iSub, iCon
        real(8)                :: dSumLocal
!
        do iSub = 1, nSublatticePhase(iPhaseIDIn)
            dSumLocal = 0D0
            do iCon = 1, nConstituentSublattice(iPhaseIDIn,iSub)
                dSumLocal = dSumLocal + dSiteInOut(iSub,iCon)
            end do
            if (dSumLocal <= 0D0) cycle
            do iCon = 1, nConstituentSublattice(iPhaseIDIn,iSub)
                dSiteInOut(iSub,iCon) = dSiteInOut(iSub,iCon) / dSumLocal
            end do
        end do
!
        return
    end subroutine NormalizeSiteFractionForPhase

    integer function MiscibleBasePhase(iSolnPhase)
        integer, intent(in) :: iSolnPhase
        integer             :: iPhaseIndex
!
        MiscibleBasePhase = iSolnPhase
        do iPhaseIndex = iSolnPhase, 1, -1
            if (.NOT. lMiscibility(iPhaseIndex)) then
                MiscibleBasePhase = iPhaseIndex
                return
            end if
        end do
!
        return
    end function MiscibleBasePhase
!
    logical function HasExplicitImmiscibilityCopies(iSolnPhase)
        integer, intent(in) :: iSolnPhase
        integer             :: iBasePhase, iPhaseIndex
!
        HasExplicitImmiscibilityCopies = .FALSE.
        if ((iSolnPhase < 1).OR.(iSolnPhase > nSolnPhasesSys)) return
!
        iBasePhase = MiscibleBasePhase(iSolnPhase)
        do iPhaseIndex = 1, nSolnPhasesSys
            if (iPhaseIndex == iBasePhase) cycle
            if (.NOT.lMiscibility(iPhaseIndex)) cycle
            if (cSolnPhaseName(iPhaseIndex) /= cSolnPhaseName(iBasePhase)) cycle
            if (.NOT.PhaseConstituentRowsCompatible(iPhaseIndex, iBasePhase)) cycle
            HasExplicitImmiscibilityCopies = .TRUE.
            return
        end do
!
        return
    end function HasExplicitImmiscibilityCopies
!
    integer function NextImmisciblePhase(iBasePhase)
        integer, intent(in) :: iBasePhase
        integer             :: iPhaseIndex
!
        NextImmisciblePhase = 0
        if ((iBasePhase < 1).OR.(iBasePhase > nSolnPhasesSys)) return
!
        do iPhaseIndex = 1, nSolnPhasesSys
            if (iPhaseIndex == iBasePhase) cycle
            if (.NOT.lMiscibility(iPhaseIndex)) cycle
            if (cSolnPhaseName(iPhaseIndex) /= cSolnPhaseName(iBasePhase)) cycle
            if (.NOT.PhaseConstituentRowsCompatible(iPhaseIndex, iBasePhase)) cycle
            if (SolutionPhaseActive(iPhaseIndex)) cycle
            NextImmisciblePhase = iPhaseIndex
            return
        end do
!
        return
    end function NextImmisciblePhase
!
    logical function SolutionPhaseActive(iSolnPhase)
        integer, intent(in) :: iSolnPhase
        integer             :: iSoln, iSlot
!
        SolutionPhaseActive = .FALSE.
        do iSoln = 1, nSolnPhases
            iSlot = nElements - iSoln + 1
            if (iAssemblage(iSlot) == -iSolnPhase) then
                SolutionPhaseActive = .TRUE.
                return
            end if
        end do
!
        return
    end function SolutionPhaseActive
!
    logical function PhaseConstituentRowsCompatible(iCandidatePhase, iTargetPhase)
        integer, intent(in) :: iCandidatePhase, iTargetPhase
        integer             :: iCandidateFirst, iCandidateLast
        integer             :: iTargetFirst, iTargetLast
        integer             :: iLocal, nLocal
!
        PhaseConstituentRowsCompatible = .FALSE.
        if ((iCandidatePhase < 1).OR.(iCandidatePhase > nSolnPhasesSys)) return
        if ((iTargetPhase < 1).OR.(iTargetPhase > nSolnPhasesSys)) return
!
        iCandidateFirst = nSpeciesPhase(iCandidatePhase-1) + 1
        iCandidateLast  = nSpeciesPhase(iCandidatePhase)
        iTargetFirst    = nSpeciesPhase(iTargetPhase-1) + 1
        iTargetLast     = nSpeciesPhase(iTargetPhase)
        nLocal          = iCandidateLast - iCandidateFirst
!
        if ((iTargetLast - iTargetFirst) /= nLocal) return
!
        do iLocal = 0, nLocal
            if (TRIM(cSpeciesName(iCandidateFirst+iLocal)) /= &
                TRIM(cSpeciesName(iTargetFirst+iLocal))) return
        end do
!
        PhaseConstituentRowsCompatible = .TRUE.
!
        return
    end function PhaseConstituentRowsCompatible
!
    subroutine ReleaseLevel2LagrangeWorkArrays
!
        if (allocated(dMolesPhaseHistory))  deallocate(dMolesPhaseHistory)
        if (allocated(dMolesPhaseLast))     deallocate(dMolesPhaseLast)
        if (allocated(dMolFractionOld))     deallocate(dMolFractionOld)
        if (allocated(iPhaseLevel))         deallocate(iPhaseLevel)
        if (allocated(dStoichSpeciesLevel)) deallocate(dStoichSpeciesLevel)
!
        return
    end subroutine ReleaseLevel2LagrangeWorkArrays
!
end subroutine Level2Lagrange
