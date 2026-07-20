!> \brief Legacy one-shot static-grid workflow retained for API/source compatibility.
!!
!! \details This module is unreachable from the production minimizer as of
!! 07/19/2026.  Grid discovery now contributes starts to standard PEA; the
!! coupled one-shot, certificate, recovery, and fallback path below is legacy.
!
module GridFrontEndWorkflow
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GridFrontEndWorkflow.f90
    !> \brief   Orchestrate default-off static-grid minimization.
    !> \author  S.Y. Kwon
    !> \date    Jul. 17, 2026
    !> \sa      MultiPhaseMinimizer.f90
    !> \sa      CompDrivingForceAll.f90
    !> \sa      CheckPhaseAssemblage.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Implemented the typed static-grid Newton and certification workflow and retained it as unreachable legacy after standard PEA integration.
    !
    ! Purpose:
    ! ========
    !
    !> \details A static grid is useful only if it can replace repeated
    !! plane-dependent discovery.  After the finite grid exchange, this routine
    !! uses the existing subminimizer once at the settled grid plane and performs
    !! one finite exchange over those refined rows.  It then hands that active
    !! set to the coupled Lagrangian solve and subminimizes every eligible phase
    !! once against the converged plane for certification.  A favorable inactive
    !! witness or an unconverged Newton solve triggers at most one bounded
    !! recovery.  The failing phases are subminimized once more with the existing
    !! endpoint-ranked and symmetry-partner starts, the original sampled support
    !! rows are restored, and finite Leveling supplies lever amounts before one
    !! more coupled solve.  A second failure records the typed reason and fails
    !! closed.  The workflow never recursively enters the production dynamic path.
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - The certificate uses the existing PEA strict Leveling tolerance of
    !   -1D-8.  Static-grid discovery introduces no new thermodynamic tolerance.
    ! - Certificate subminimization is transactional.  All final active-state
    !   thermodynamic arrays are restored before accepting or falling back.
    ! - CompDrivingForceAll is used only as the fresh all-phase sweep.  It does
    !   not change the elemental-potential plane during certification.
    ! - Repeated refined composition-set rows and deliberately empty active-set
    !   slots can make the temporary Leveling plane singular.  Those events are
    !   typed and handed to Level2Lagrange, whose existing duplicate/under-filled
    !   handling constructs the Newton set.
    ! - INFOThermo=42 is the generic grid fail-closed return.  Every such return
    !   must set nGridFallback and a non-NONE iGridFallbackReason; the reason is
    !   the authoritative physical or workflow classification.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE GridDiscovery, ONLY: CapturePreviousGridSeeds, ApplySelectedGridCandidateFractions, &
        CopyGridPointMolFraction, GridTangentRowIndex, GridPartitionRowIndex

    implicit none

    private
    public :: RunStaticGridFrontEnd

contains

!> \brief Refine static-grid winners, solve their active set, and certify the plane.
!!
!! \details Preserves the sampled LP support, runs one fixed-plane refinement,
!! transports typed composition-set evidence through Level2, and permits only
!! one bounded recovery before an uncertified solve fails closed.
subroutine RunStaticGridFrontEnd
    implicit none

    real(8) :: dTimerStart, dTimerStop
    logical :: lRecoverySettled
    integer :: iGridSupportBasis(nElements)
    real(8) :: dGridSupportMoles(nElements)

    ! Step 1. Retain the sampled LP support rows before Level2 decides how many
    ! physical composition sets they represent.  A bounded recovery needs
    ! those rows for a finite lever seed when a newly admitted phase lies on
    ! one side of the converged active composition.
    iGridSupportBasis = iAssemblage
    dGridSupportMoles = dMolesPhase

    ! Step 2. Require a valid, nonempty finite row problem.
    if (INFOThermo /= 0) then
        nGridFallback = 1
        iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
        call RecordSwallowedGridInfo(GRID_INFO_STAGE_ENTRY, INFOThermo)
        INFOThermo = 0
        call StopUncertifiedGridSolve
        call CapturePreviousGridSeeds
        return
    end if

    if (nGridPoint <= 0) then
        nGridFallback = 1
        iGridFallbackReason = GRID_FALLBACK_EMPTY_GRID
        call StopUncertifiedGridSolve
        call CapturePreviousGridSeeds
        return
    end if

    ! Step 3. Subminimize the selected row constitutions once at the static
    ! Leveling plane, then solve one more finite exchange.
    call RefineStaticGridPlaneOnce
    if (INFOThermo /= 0) then
        nGridFallback = 1
        iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
        call RecordSwallowedGridInfo(GRID_INFO_STAGE_REFINEMENT, INFOThermo)
        INFOThermo = 0
        call StopUncertifiedGridSolve
        call CapturePreviousGridSeeds
        return
    end if

    ! Step 4. Transport the physically distinct refined sets into Level2.
    call cpu_time(dTimerStart)
    call ApplySelectedGridCandidateFractions
    call Level2Lagrange
    call cpu_time(dTimerStop)
    dGEMTimingHandoff = dGEMTimingHandoff + MAX(0D0,dTimerStop-dTimerStart)
    if (INFOThermo /= 0) then
        nGridFallback = 1
        iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
        call RecordSwallowedGridInfo(GRID_INFO_STAGE_HANDOFF, INFOThermo)
        INFOThermo = 0
        call StopUncertifiedGridSolve
        call CapturePreviousGridSeeds
        return
    end if

    ! Step 5. Run one coupled fixed-active-set Newton solve.
    call cpu_time(dTimerStart)
    call RunLagrangianGEM
    call cpu_time(dTimerStop)
    dGEMTimingLagrangian = dGEMTimingLagrangian + MAX(0D0,dTimerStop-dTimerStart)

    ! Step 6. Certify the converged plane.  One failed certificate may use the
    ! bounded static-grid recovery; no dynamic discovery recursion is allowed.
    if (lConverged) then
        call CertifyStaticGridPlane
        if (INFOThermo /= 0) then
            call StopUncertifiedGridSolve
            call CapturePreviousGridSeeds
            return
        end if
        if (iPEAExitStatus == PEA_EXIT_STATUS_CERTIFIED_SETTLED) then
            call CapturePreviousGridSeeds
            return
        end if
        call RunBoundedGridRecovery(iGridSupportBasis, dGridSupportMoles, lRecoverySettled)
        if (lRecoverySettled) then
            call CapturePreviousGridSeeds
            return
        end if
        nGridFallback = 1
        if (iGridFallbackReason == GRID_FALLBACK_NONE) then
            if (lConverged) then
                iGridFallbackReason = GRID_FALLBACK_CERTIFICATE
            else
                iGridFallbackReason = GRID_FALLBACK_LAGRANGIAN
            end if
        end if
    else
        nGridFallback = 1
        iGridFallbackReason = GRID_FALLBACK_LAGRANGIAN
    end if

    call StopUncertifiedGridSolve
    call CapturePreviousGridSeeds

    return
end subroutine RunStaticGridFrontEnd


!> \brief Retry one uncertified grid active set without dynamic rediscovery.
!!
!! \details Reuses the original sampled support, production multistart minima,
!! and finite Leveling lever amounts before one final coupled Newton solve.
subroutine RunBoundedGridRecovery(iGridSupportBasis, dGridSupportMoles, lSettled)
    implicit none

    integer, intent(in) :: iGridSupportBasis(nElements)
    real(8), intent(in) :: dGridSupportMoles(nElements)
    logical, intent(out) :: lSettled
    real(8) :: dTimerStart, dTimerStop

    ! Step 1. Arm the recovery exactly once.
    lSettled = .FALSE.
    if (nGridRecoveryAttempt /= 0) return
    nGridRecoveryAttempt = 1
    lGridRecoveryPassActive = .TRUE.

    ! Step 2. Re-subminimize certificate-failing phases and solve their lever
    ! amounts against the original sampled support rows.
    call RefineStaticGridPlaneOnce(iGridSupportBasis, dGridSupportMoles)
    if (INFOThermo /= 0) then
        lGridRecoveryPassActive = .FALSE.
        return
    end if

    ! Step 3. Transport the recovered physical composition sets through Level2.
    call cpu_time(dTimerStart)
    call ApplySelectedGridCandidateFractions
    call Level2Lagrange
    call cpu_time(dTimerStop)
    dGEMTimingHandoff = dGEMTimingHandoff + MAX(0D0,dTimerStop-dTimerStart)
    if (INFOThermo /= 0) then
        lGridRecoveryPassActive = .FALSE.
        return
    end if

    ! Step 4. Run and certify one final coupled Newton solve.
    call cpu_time(dTimerStart)
    call RunLagrangianGEM
    call cpu_time(dTimerStop)
    dGEMTimingLagrangian = dGEMTimingLagrangian + MAX(0D0,dTimerStop-dTimerStart)
    if (lConverged) then
        call CertifyStaticGridPlane
        lSettled = iPEAExitStatus == PEA_EXIT_STATUS_CERTIFIED_SETTLED
    end if
    lGridRecoveryPassActive = .FALSE.

    return
end subroutine RunBoundedGridRecovery


!> \brief Refine selected constitutions once at a fixed elemental-potential plane.
subroutine RefineStaticGridPlaneOnce(iGridSupportBasis, dGridSupportMoles)
    implicit none

    integer, intent(in), optional :: iGridSupportBasis(nElements)
    real(8), intent(in), optional :: dGridSupportMoles(nElements)
    real(8) :: dTimerStart, dTimerStop

    ! Step 1. Run the production endpoint-ranked and symmetry-orbit multistart.
    call cpu_time(dTimerStart)
    call InitCheckPhaseAssemblage
    ! Step 2. During recovery, restore the original sampled LP support and keep
    ! only tangent rows belonging to certificate-failing phases.
    ! Step 3. Reprice the finite rows and perform one bounded Leveling exchange.
    if (INFOThermo == 0) then
        if (lGridRecoveryPassActive.AND.PRESENT(iGridSupportBasis).AND.&
            PRESENT(dGridSupportMoles)) then
            call RestoreGridSupportBasis(iGridSupportBasis, dGridSupportMoles)
            call SuppressNonfailingRecoveryCandidates
        end if
    end if
    if (INFOThermo == 0) then
        iGridLevelingRepeatedAssemblage = 0
        lGridRefinementSweepActive = .TRUE.
        call LevelingSolver
        lGridRefinementSweepActive = .FALSE.
        if (INFOThermo == 0) then
            iGridRefinementStatus = GRID_REFINEMENT_STATUS_SETTLED
        else if ((INFOThermo == 10).AND.&
            (iGridLevelingRepeatedAssemblage == 1)) then
            iGridRefinementStatus = &
                GRID_REFINEMENT_STATUS_REPEATED_SINGULAR_HANDOFF
            call RecordSwallowedGridInfo(GRID_INFO_STAGE_REFINEMENT, INFOThermo)
            INFOThermo = 0
        else if ((INFOThermo == 10).AND.ANY(iAssemblage == 0)) then
            iGridRefinementStatus = &
                GRID_REFINEMENT_STATUS_UNDERFILLED_SINGULAR_HANDOFF
            call RecordSwallowedGridInfo(GRID_INFO_STAGE_REFINEMENT, INFOThermo)
            INFOThermo = 0
        else
            iGridRefinementStatus = GRID_REFINEMENT_STATUS_ERROR
        end if
    end if
    ! Step 4. Record the typed refinement outcome and elapsed CPU time.
    call cpu_time(dTimerStop)
    dGEMTimingGridRefinement = dGEMTimingGridRefinement + &
        MAX(0D0,dTimerStop-dTimerStart)
    lGridRefinementSweepActive = .FALSE.

    return
end subroutine RefineStaticGridPlaneOnce


!> \brief Restore the original sampled LP support as recovery lever endpoints.
subroutine RestoreGridSupportBasis(iGridSupportBasis, dGridSupportMoles)
    implicit none

    integer, intent(in) :: iGridSupportBasis(nElements)
    real(8), intent(in) :: dGridSupportMoles(nElements)

    integer :: iSlot, iLevelRow, iPhaseIndex, iFirst, iLast
    logical :: lCandidateFound

    ! Step 1. Restore the sampled support row identities and phase amounts.
    iAssemblage = iGridSupportBasis
    dMolesPhase = dGridSupportMoles
    iPhaseGEM = 0
    dChemicalPotentialGEM = 0D0
    dStoichSpeciesGEM = 0D0
    dAtomFractionSpeciesGEM = 0D0
    dMolFractionGEM = 0D0

    ! Step 2. Reconstruct each support row's phase basis and constitution.
    do iSlot = 1, nElements
        iLevelRow = iAssemblage(iSlot)
        if ((iLevelRow < 1).OR.(iLevelRow > nSpeciesLevel)) then
            nGridFallback = 1
            iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
            INFOThermo = 42
            return
        end if

        iPhaseIndex = iPhaseLevel(iLevelRow)
        iPhaseGEM(iSlot) = iPhaseIndex
        dChemicalPotentialGEM(iSlot) = dLevelingChemicalPotential(iLevelRow)
        dStoichSpeciesGEM(iSlot,:) = dStoichSpeciesLevel(iLevelRow,:)
        dAtomFractionSpeciesGEM(iSlot,:) = dLevelingCompositionSpecies(iLevelRow,:)
        if (iPhaseIndex <= 0) cycle

        iFirst = nSpeciesPhase(iPhaseIndex-1) + 1
        iLast = nSpeciesPhase(iPhaseIndex)
        lCandidateFound = .FALSE.
        call CopyGridPointMolFraction(iLevelRow, iSlot, iPhaseIndex, lCandidateFound)
        if (.NOT.lCandidateFound) then
            call CopyLevelCandidateMolFraction(iLevelRow, iSlot, iPhaseIndex, lCandidateFound)
        end if
        if (.NOT.lCandidateFound) then
            nGridFallback = 1
            iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
            INFOThermo = 42
            return
        end if
        if (SUM(dMolFractionGEM(iSlot,iFirst:iLast)) <= 1D-300) then
            nGridFallback = 1
            iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
            INFOThermo = 42
            return
        end if
    end do

    return
end subroutine RestoreGridSupportBasis


!> \brief Keep recovery tangent rows only for certificate-failing phases.
subroutine SuppressNonfailingRecoveryCandidates
    implicit none

    integer :: iPhaseIndex, iLevelRow
    logical :: lFailingPhase

    if (.NOT.allocated(lLevelingRowExcluded)) then
        nGridFallback = 1
        iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
        INFOThermo = 42
        return
    end if
    if (SIZE(lLevelingRowExcluded) /= nSpeciesLevel) then
        nGridFallback = 1
        iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
        INFOThermo = 42
        return
    end if

    ! Step 1. Walk every solution phase and retain its primary tangent row only
    ! when the fresh certificate marked that phase for recovery.
    do iPhaseIndex = 1, nSolnPhasesSys
        lFailingPhase = .FALSE.
        if (allocated(lGridRecoveryPhase)) then
            if (iPhaseIndex <= SIZE(lGridRecoveryPhase)) then
                lFailingPhase = lGridRecoveryPhase(iPhaseIndex)
            end if
        end if
        ! The secondary partition row represents the disordered manifold at
        ! the current bulk composition.  It is valid discovery evidence, but
        ! retaining it beside a newly admitted ordered minimum would again
        ! force the ordered lever amount to zero.
        iLevelRow = GridPartitionRowIndex(iPhaseIndex)
        if ((iLevelRow >= 1).AND.(iLevelRow <= nSpeciesLevel)) then
            lLevelingRowExcluded(iLevelRow) = .TRUE.
        end if
        if (lFailingPhase) cycle

        iLevelRow = GridTangentRowIndex(iPhaseIndex)
        if ((iLevelRow >= 1).AND.(iLevelRow <= nSpeciesLevel)) then
            lLevelingRowExcluded(iLevelRow) = .TRUE.
        end if
    end do

    return
end subroutine SuppressNonfailingRecoveryCandidates


!> \brief Record one nonzero solver code translated into a typed grid outcome.
subroutine RecordSwallowedGridInfo(iStage, iInfo)
    implicit none

    integer, intent(in) :: iStage, iInfo

    if (iInfo == 0) return
    if (nGridSwallowedInfo >= GRID_SWALLOWED_INFO_CAP) then
        nGridSwallowedInfoDropped = nGridSwallowedInfoDropped + 1
        return
    end if

    nGridSwallowedInfo = nGridSwallowedInfo + 1
    iGridSwallowedInfo(nGridSwallowedInfo) = iInfo
    iGridSwallowedInfoStage(nGridSwallowedInfo) = iStage

    return
end subroutine RecordSwallowedGridInfo


!> \brief Stop an uncertified grid solve without entering dynamic discovery.
subroutine StopUncertifiedGridSolve
    implicit none

    nGridFallback = 1
    if (iGridFallbackReason == GRID_FALLBACK_NONE) &
        iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
    lConverged = .FALSE.
    INFOThermo = 42

    return
end subroutine StopUncertifiedGridSolve


!> \brief Restore the pre-certificate site-fraction allocation and shape.
!!
!! \details Shape drift is reported after the original state is reconstructed,
!! so certification remains transactional even on failure.
subroutine RestoreCertifiedSiteFractionState(&
    dSiteFractionSave, lWasAllocated, lShapeMatches)
    implicit none

    real(8), dimension(:,:,:), allocatable, intent(inout) :: dSiteFractionSave
    logical, intent(in) :: lWasAllocated
    logical, intent(out) :: lShapeMatches

    lShapeMatches = .TRUE.
    if (lWasAllocated) then
        if (.NOT.allocated(dSiteFractionSave)) then
            lShapeMatches = .FALSE.
            return
        end if
        if (.NOT.allocated(dSiteFraction)) then
            lShapeMatches = .FALSE.
        else
            lShapeMatches = (SIZE(dSiteFraction,1) == SIZE(dSiteFractionSave,1)).AND.&
                (SIZE(dSiteFraction,2) == SIZE(dSiteFractionSave,2)).AND.&
                (SIZE(dSiteFraction,3) == SIZE(dSiteFractionSave,3))
        end if
        if (.NOT.lShapeMatches) then
            if (allocated(dSiteFraction)) deallocate(dSiteFraction)
            allocate(dSiteFraction(&
                SIZE(dSiteFractionSave,1), SIZE(dSiteFractionSave,2), &
                SIZE(dSiteFractionSave,3)))
        end if
        dSiteFraction = dSiteFractionSave
        deallocate(dSiteFractionSave)
    else if (allocated(dSiteFraction)) then
        lShapeMatches = .FALSE.
        deallocate(dSiteFraction)
    end if

    return
end subroutine RestoreCertifiedSiteFractionState


!> \brief Certify a converged grid plane with one transactional all-phase sweep.
!!
!! \details Runs the existing fresh driving-force sweep, records phases with
!! favorable or unknown candidates for a possible bounded recovery, and restores
!! every active thermodynamic array before returning the typed certificate.
subroutine CertifyStaticGridPlane
    implicit none

    integer :: iFirstCompound, iLastCompound, iInfoSave, iCertificateInfo
    integer :: iPhaseIndex
    real(8) :: dMinimumPotential, dTimerStart, dTimerStop
    real(8), dimension(:), allocatable :: dMolFractionSave, dMolesSpeciesSave
    real(8), dimension(:), allocatable :: dChemicalPotentialSave, dActivitySave
    real(8), dimension(:), allocatable :: dPartialHSave, dPartialSSave, dPartialCpSave
    real(8), dimension(:), allocatable :: dPartialExcessSave, dPartialHXSSave
    real(8), dimension(:), allocatable :: dPartialSXSSave, dPartialCpXSSave
    real(8), dimension(:), allocatable :: dMagGSave, dMagHSave, dMagSSave, dMagCpSave
    real(8), dimension(:), allocatable :: dGibbsPhaseSave, dPhasePotentialSave
    real(8), dimension(:), allocatable :: dDrivingForceSave, dSumFractionSave
    real(8), dimension(:,:), allocatable :: dEffStoichSave
    real(8), dimension(:,:,:), allocatable :: dSiteFractionSave
    logical :: lSiteFractionWasAllocated, lSiteFractionShapeMatches

    ! Step 1. Snapshot every active thermodynamic array touched by certification.
    allocate(dMolFractionSave(nSpecies), dMolesSpeciesSave(nSpecies))
    allocate(dChemicalPotentialSave(nSpecies), dActivitySave(nSpecies))
    allocate(dPartialHSave(nSpecies), dPartialSSave(nSpecies), dPartialCpSave(nSpecies))
    allocate(dPartialExcessSave(nSpecies), dPartialHXSSave(nSpecies))
    allocate(dPartialSXSSave(nSpecies), dPartialCpXSSave(nSpecies))
    allocate(dMagGSave(nSpecies), dMagHSave(nSpecies), dMagSSave(nSpecies), dMagCpSave(nSpecies))
    allocate(dGibbsPhaseSave(nSolnPhasesSys), dPhasePotentialSave(nSpecies))
    allocate(dDrivingForceSave(nSolnPhasesSys), dSumFractionSave(nSolnPhasesSys))
    allocate(dEffStoichSave(nSolnPhasesSys,nElements))

    dMolFractionSave = dMolFraction
    dMolesSpeciesSave = dMolesSpecies
    dChemicalPotentialSave = dChemicalPotential
    dActivitySave = dActivity
    dPartialHSave = dPartialEnthalpy
    dPartialSSave = dPartialEntropy
    dPartialCpSave = dPartialHeatCapacity
    dPartialExcessSave = dPartialExcessGibbs
    dPartialHXSSave = dPartialEnthalpyXS
    dPartialSXSSave = dPartialEntropyXS
    dPartialCpXSSave = dPartialHeatCapacityXS
    dMagGSave = dMagGibbsEnergy
    dMagHSave = dMagEnthalpy
    dMagSSave = dMagEntropy
    dMagCpSave = dMagHeatCapacity
    dGibbsPhaseSave = dGibbsSolnPhase
    dPhasePotentialSave = dPhasePotential
    dDrivingForceSave = dDrivingForceSoln
    dSumFractionSave = dSumMolFractionSoln
    dEffStoichSave = dEffStoichSolnPhase
    lSiteFractionWasAllocated = allocated(dSiteFraction)
    if (lSiteFractionWasAllocated) then
        allocate(dSiteFractionSave(SIZE(dSiteFraction,1),SIZE(dSiteFraction,2),SIZE(dSiteFraction,3)))
        dSiteFractionSave = dSiteFraction
    end if

    ! Step 2. Subminimize every eligible phase once against the converged plane.
    iInfoSave = INFOThermo
    call cpu_time(dTimerStart)
    call CompDrivingForceAll
    iCertificateInfo = INFOThermo
    call cpu_time(dTimerStop)
    dGEMTimingCertification = dGEMTimingCertification + MAX(0D0,dTimerStop-dTimerStart)
    nGEMCertificationSweep = nGEMCertificationSweep + 1
    nGridCertificateSweep = nGridCertificateSweep + 1

    ! Step 3. Find the most favorable fresh solution or compound candidate.
    dMinimumPotential = HUGE(1D0)
    if (nSolnPhasesSys > 0) dMinimumPotential = MIN(dMinimumPotential,MINVAL(dDrivingForceSoln))
    iFirstCompound = nSpeciesPhase(nSolnPhasesSys) + 1
    iLastCompound = nSpecies - nDummySpecies
    if (iFirstCompound <= iLastCompound) then
        dMinimumPotential = MIN(dMinimumPotential, &
            MINVAL(dPhasePotential(iFirstCompound:iLastCompound)))
    end if

    iPEAExitFreshMinPointSweep = 1
    dPEAExitMinPhasePotential = dMinimumPotential
    dPEAExitTolerance = -1D-8
    ! Step 4. Mark phases carrying favorable or unknown candidate evidence.
    if (allocated(lGridRecoveryPhase)) deallocate(lGridRecoveryPhase)
    allocate(lGridRecoveryPhase(MAX(1,nSolnPhasesSys)))
    lGridRecoveryPhase = .FALSE.
    do iPhaseIndex = 1, nSolnPhasesSys
        if (dDrivingForceSoln(iPhaseIndex) < dPEAExitTolerance) then
            lGridRecoveryPhase(iPhaseIndex) = .TRUE.
        else if (allocated(iSubMinCandidateStatusSoln)) then
            if ((iSubMinCandidateStatusSoln(iPhaseIndex) /= SUBMIN_CANDIDATE_CONVERGED).AND.&
                (iSubMinCandidateStatusSoln(iPhaseIndex) /= SUBMIN_CANDIDATE_NEGATIVE_WITNESS).AND.&
                (iSubMinCandidateStatusSoln(iPhaseIndex) /= SUBMIN_CANDIDATE_DUPLICATE)) then
                lGridRecoveryPhase(iPhaseIndex) = .TRUE.
            end if
        end if
    end do
    ! Step 5. Issue the same settled/unsettled certificate used by production PEA.
    if ((INFOThermo == 0).AND.(dMinimumPotential >= dPEAExitTolerance)) then
        iPEAExitStatus = PEA_EXIT_STATUS_CERTIFIED_SETTLED
        iPEAExitReason = PEA_EXIT_REASON_SETTLED
    else
        iPEAExitStatus = PEA_EXIT_STATUS_UNSETTLED
        iPEAExitReason = PEA_EXIT_REASON_NONE
    end if

    ! Step 6. Restore the converged active state before the caller accepts or
    ! retries the certificate.
    dMolFraction = dMolFractionSave
    dMolesSpecies = dMolesSpeciesSave
    dChemicalPotential = dChemicalPotentialSave
    dActivity = dActivitySave
    dPartialEnthalpy = dPartialHSave
    dPartialEntropy = dPartialSSave
    dPartialHeatCapacity = dPartialCpSave
    dPartialExcessGibbs = dPartialExcessSave
    dPartialEnthalpyXS = dPartialHXSSave
    dPartialEntropyXS = dPartialSXSSave
    dPartialHeatCapacityXS = dPartialCpXSSave
    dMagGibbsEnergy = dMagGSave
    dMagEnthalpy = dMagHSave
    dMagEntropy = dMagSSave
    dMagHeatCapacity = dMagCpSave
    dGibbsSolnPhase = dGibbsPhaseSave
    dPhasePotential = dPhasePotentialSave
    dDrivingForceSoln = dDrivingForceSave
    dSumMolFractionSoln = dSumFractionSave
    dEffStoichSolnPhase = dEffStoichSave
    call RestoreCertifiedSiteFractionState(&
        dSiteFractionSave, lSiteFractionWasAllocated, lSiteFractionShapeMatches)
    if ((iInfoSave == 0).AND.(iCertificateInfo /= 0)) &
        call RecordSwallowedGridInfo(GRID_INFO_STAGE_CERTIFICATE, iCertificateInfo)
    if (.NOT.lSiteFractionShapeMatches) then
        nGridFallback = 1
        iGridFallbackReason = GRID_FALLBACK_STATE_SHAPE
        INFOThermo = 42
    else if (iPEAExitStatus == PEA_EXIT_STATUS_CERTIFIED_SETTLED) then
        INFOThermo = iInfoSave
    else if (iInfoSave == 0) then
        ! The bounded census above preserves a certificate-only evaluator code;
        ! the caller sees the unsettled certificate and may use its one recovery.
        INFOThermo = 0
    end if

    deallocate(dMolFractionSave, dMolesSpeciesSave, dChemicalPotentialSave, dActivitySave)
    deallocate(dPartialHSave, dPartialSSave, dPartialCpSave, dPartialExcessSave)
    deallocate(dPartialHXSSave, dPartialSXSSave, dPartialCpXSSave)
    deallocate(dMagGSave, dMagHSave, dMagSSave, dMagCpSave)
    deallocate(dGibbsPhaseSave, dPhasePotentialSave, dDrivingForceSave, dSumFractionSave)
    deallocate(dEffStoichSave)

    return
end subroutine CertifyStaticGridPlane

end module GridFrontEndWorkflow
