!> \brief Build and price the default-off static constitution grid.
!!
!! \details Generates deterministic phase constitutions once at fixed T and P,
!! evaluates their Gibbs energies without subminimization, and appends them to
!! the finite Leveling candidate set.
!
module GridDiscovery
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GridDiscovery.f90
    !> \brief   Deterministic static-grid discovery front end for Leveling.
    !> \author  S.Y. Kwon
    !> \date    Jul. 17, 2026
    !> \sa      RunLeveling.f90
    !> \sa      LevelCandidatePool.f90
    !> \sa      ProjectOrderDisorderCompanionFraction.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Added deterministic immutable static-grid generation, direct batch pricing, typed refinement and failure accounting, and collision-safe selected-row handoff.
    !
    ! Purpose:
    ! ========
    !
    !> \details The dynamic PEA rows move whenever the elemental-potential plane
    !! moves.  This module supplies a plane-independent discovery set: simplex
    !! ladders for ordinary solutions, OpenCalphad-style endmember mixtures for
    !! multi-sublattice models, and disordered-manifold projections for declared
    !! order/disorder partitions.  Grid rows remain immutable during Leveling.
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Sobol coordinates are mapped to a simplex by sorted spacings.  Direct
    !   normalization of hypercube coordinates is intentionally not used.
    ! - Generation order is phase, face, and sequence-index order; no random
    !   state or compiler-dependent ordering participates.
    ! - Grid evaluation calls only the fixed-constitution Gibbs evaluator.  It
    !   does not call Subminimization, a Hessian, or a composition projection.
    ! - The 4096-point phase budget is enforced by typed layer coarsening.  A
    !   vertex layer that cannot fit is rejected instead of silently truncated.
    ! - Dynamic tangent rows follow all immutable rows; partition rows follow
    !   all primary tangent rows.  Named row-index functions own this layout.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    USE GridBatchKernels, ONLY: EvaluateGridPhaseGibbsBatch, EvaluateGridPhaseGibbsScalar

    implicit none

    private
    public :: GenerateStaticGrid
    public :: CopyGridPointMolFraction
    public :: RegisterSelectedGridSeeds
    public :: RefineSelectedGridSeeds
    public :: ApplySelectedGridCandidateFractions
    public :: RestoreStaticGridLevelingRows
    public :: CapturePreviousGridSeeds
    public :: GridTangentRowIndex
    public :: GridPartitionRowIndex

    integer, parameter :: nGridLayer = 5
    integer, parameter :: nSobolBit = 32

contains

!> \brief Generate, store, and price one deterministic constitution grid.
!!
!! \details Plans every eligible phase under the per-phase point cap, emits
!! each fixed constitution in canonical order, and evaluates its Gibbs energy
!! once before Leveling changes the elemental-potential plane.
subroutine GenerateStaticGrid
    implicit none

    integer :: iPhaseIndex, iTargetPhase, iDefinitionPhase
    integer :: iClass, iSubstitutionalSublattice, iInterstitialSublattice
    integer :: nSimplexConstituent, nPhasePointLocal, nFractionLocal
    integer :: nPreviousPointLocal, nPreviousFractionLocal
    integer :: iPointCursor, iFractionCursor
    integer(8) :: nFractionTotal
    logical :: lSupported, lHelper

    ! Step 1. Reset per-solve grid state and leave the production OFF path empty.
    call ResetStaticGridSolveState
    if (.NOT.lGridFrontEndActive) return
    if (nSolnPhasesSys <= 0) return

    ! Step 2. Classify every phase and choose the finest complete dimensional
    ! layers that fit its point budget.
    allocate(nGridPhasePoint(nSolnPhasesSys), nGridPhaseProjectionFailure(nSolnPhasesSys))
    allocate(iGridPhaseClass(nSolnPhasesSys), iGridPhaseEvaluationMode(nSolnPhasesSys))
    allocate(iGridPhaseStatus(nSolnPhasesSys), iGridPhaseLayerRung(nSolnPhasesSys,nGridLayer))
    nGridPhasePoint = 0
    nGridPhaseProjectionFailure = 0
    iGridPhaseClass = GRID_PHASE_CLASS_NONE
    iGridPhaseEvaluationMode = GRID_EVALUATION_NOT_RUN
    iGridPhaseStatus = GRID_PHASE_STATUS_NOT_EVALUATED
    iGridPhaseLayerRung = 0

    nFractionTotal = 0_8
    do iPhaseIndex = 1, nSolnPhasesSys
        call ClassifyGridPhase(iPhaseIndex, iTargetPhase, iDefinitionPhase, iClass, &
            iSubstitutionalSublattice, iInterstitialSublattice, nSimplexConstituent, &
            lSupported, lHelper)
        if (lHelper) then
            iGridPhaseStatus(iPhaseIndex) = GRID_PHASE_STATUS_HELPER_SKIPPED
            cycle
        end if
        if (.NOT.lSupported) then
            iGridPhaseStatus(iPhaseIndex) = GRID_PHASE_STATUS_UNSUPPORTED
            cycle
        end if

        call PlanGridPhase(iTargetPhase, iDefinitionPhase, iClass, &
            iSubstitutionalSublattice, iInterstitialSublattice, nSimplexConstituent, &
            nPhasePointLocal)
        nGridPhasePoint(iTargetPhase) = nPhasePointLocal
        iGridPhaseClass(iTargetPhase) = iClass
        if (nPhasePointLocal <= 0) cycle

        nFractionLocal = nSpeciesPhase(iTargetPhase) - nSpeciesPhase(iTargetPhase-1)
        nFractionTotal = nFractionTotal + INT(nPhasePointLocal,8) * INT(nFractionLocal,8)
    end do

    ! Step 3. Admit compatible previous-solution constitutions without letting
    ! them displace any required grid layer.
    call CountCompatiblePreviousSeeds(nPreviousPointLocal,nPreviousFractionLocal)
    nGridPoint = SUM(nGridPhasePoint) + nPreviousPointLocal
    if (nGridPoint <= 0) return
    nFractionTotal = nFractionTotal + INT(nPreviousFractionLocal,8)
    if (nFractionTotal > INT(HUGE(nFractionLocal),8)) then
        call SetStaticGridWorkflowError
        nGridPoint = 0
        return
    end if
    nGridFraction = INT(nFractionTotal)

    ! Step 4. Allocate the exact ragged constitution store determined above.
    allocate(iGridPointPhase(nGridPoint), iGridPointDisplayPhase(nGridPoint))
    allocate(iGridPointStatus(nGridPoint), iGridPointIdentityOrdinal(nGridPoint))
    allocate(iGridPointLevelRow(nGridPoint), iGridPointFractionOffset(nGridPoint))
    allocate(iGridPointFractionSize(nGridPoint), dGridPointFraction(nGridFraction))
    iGridPointPhase = 0
    iGridPointDisplayPhase = 0
    iGridPointStatus = GRID_PHASE_STATUS_NOT_EVALUATED
    iGridPointIdentityOrdinal = 0
    iGridPointLevelRow = 0
    iGridPointFractionOffset = 0
    iGridPointFractionSize = 0
    dGridPointFraction = 0D0

    ! Step 5. Emit phase points in phase, face, and sequence-index order.
    iPointCursor = 0
    iFractionCursor = 1
    do iPhaseIndex = 1, nSolnPhasesSys
        if (nGridPhasePoint(iPhaseIndex) <= 0) cycle
        call ClassifyGridPhase(iPhaseIndex, iTargetPhase, iDefinitionPhase, iClass, &
            iSubstitutionalSublattice, iInterstitialSublattice, nSimplexConstituent, &
            lSupported, lHelper)
        if ((.NOT.lSupported).OR.lHelper) cycle
        if (iClass == GRID_PHASE_CLASS_ENDMEMBER_MIX) then
            call FillEndmemberMixGrid(iTargetPhase, iPointCursor, iFractionCursor)
        else
            call FillSimplexGrid(iTargetPhase, iDefinitionPhase, &
                iSubstitutionalSublattice, iInterstitialSublattice, nSimplexConstituent, &
                iPointCursor, iFractionCursor)
        end if
        if (INFOThermo /= 0) return
    end do
    call AppendPreviousGridSeeds(iPointCursor,iFractionCursor)

    if ((iPointCursor /= nGridPoint).OR.(iFractionCursor /= nGridFraction + 1)) then
        call SetStaticGridWorkflowError
        return
    end if

    ! Step 6. Append the immutable rows to Leveling and price each row at its
    ! fixed constitution.
    call ExtendGridLevelingArrays
    if (INFOThermo /= 0) return
    call EvaluateStaticGridRows
    if (INFOThermo == 0) then
        iGridPoolGeneration = iGridPoolGeneration + 1
        lGridPoolValid = .TRUE.
    end if

    return
end subroutine GenerateStaticGrid


!> \brief Reset all grid state owned by the current minimization.
!!
!! \details Clears generated rows and per-solve census values while leaving
!! optional previous-solution seed storage available for an enabled continuation.
subroutine ResetStaticGridSolveState
    implicit none

    nGridPoint = 0
    nGridFraction = 0
    nGridLevelingPass = 0
    nGridRefinementLevelingPass = 0
    iGridRefinementStatus = GRID_REFINEMENT_STATUS_NONE
    iGridLevelingRepeatedAssemblage = 0
    nGridBudgetCoarsening = 0
    nGridBudgetDroppedPoint = 0
    nGridPreviousSeedAccepted = 0
    nGridPreviousSeedDropped = 0
    nGridFallback = 0
    iGridFallbackReason = GRID_FALLBACK_NONE
    nGridCertificateSweep = 0
    nGridRecoveryAttempt = 0
    nGridSwallowedInfo = 0
    nGridSwallowedInfoDropped = 0
    nGridProjectionFailure = 0
    nGridBatchPhase = 0
    nGridScalarFallbackPhase = 0
    nGridProjectionFailureStatus = 0
    iGridSwallowedInfo = 0
    iGridSwallowedInfoStage = 0
    dGEMTimingGridGeneration = 0D0
    dGEMTimingGridRefinement = 0D0
    lGridRefinementSweepActive = .FALSE.
    lGridRecoveryPassActive = .FALSE.
    lGridPoolValid = .FALSE.

    if (allocated(iGridPointPhase)) deallocate(iGridPointPhase)
    if (allocated(iGridPointDisplayPhase)) deallocate(iGridPointDisplayPhase)
    if (allocated(iGridPointStatus)) deallocate(iGridPointStatus)
    if (allocated(iGridPointIdentityOrdinal)) deallocate(iGridPointIdentityOrdinal)
    if (allocated(iGridPointLevelRow)) deallocate(iGridPointLevelRow)
    if (allocated(iGridPointFractionOffset)) deallocate(iGridPointFractionOffset)
    if (allocated(iGridPointFractionSize)) deallocate(iGridPointFractionSize)
    if (allocated(iGridPointFromLevel)) deallocate(iGridPointFromLevel)
    if (allocated(nGridPhasePoint)) deallocate(nGridPhasePoint)
    if (allocated(nGridPhaseProjectionFailure)) deallocate(nGridPhaseProjectionFailure)
    if (allocated(iGridPhaseClass)) deallocate(iGridPhaseClass)
    if (allocated(iGridPhaseStatus)) deallocate(iGridPhaseStatus)
    if (allocated(iGridPhaseEvaluationMode)) deallocate(iGridPhaseEvaluationMode)
    if (allocated(iGridPhaseLayerRung)) deallocate(iGridPhaseLayerRung)
    if (allocated(dGridPointFraction)) deallocate(dGridPointFraction)
    if (allocated(dGridPointGibbs)) deallocate(dGridPointGibbs)
    if (allocated(dGridPointComposition)) deallocate(dGridPointComposition)
    if (allocated(dGridPointStoich)) deallocate(dGridPointStoich)
    if (allocated(lGridRecoveryPhase)) deallocate(lGridRecoveryPhase)

    return
end subroutine ResetStaticGridSolveState


!> \brief Mark a grid construction or row-layout invariant failure.
subroutine SetStaticGridWorkflowError
    implicit none

    nGridFallback = 1
    iGridFallbackReason = GRID_FALLBACK_WORKFLOW_ERROR
    INFOThermo = 42

    return
end subroutine SetStaticGridWorkflowError


!> \brief Count previous active constitutions that remain valid grid seeds.
subroutine CountCompatiblePreviousSeeds(nAccepted, nFractionAccepted)
    implicit none

    integer, intent(out) :: nAccepted, nFractionAccepted
    integer :: iSeed, iPhaseIndex, nLocal
    integer, dimension(:), allocatable :: nAcceptedPhase

    nAccepted = 0
    nFractionAccepted = 0
    nGridPreviousSeedAccepted = 0
    nGridPreviousSeedDropped = 0
    if (.NOT.lGridPreviousSolutionSeedActive) return
    if (nGridPreviousSeed <= 0) return
    if (.NOT.allocated(cGridPreviousSeedPhaseName)) return
    if (.NOT.allocated(iGridPreviousSeedFractionOffset)) return
    if (.NOT.allocated(iGridPreviousSeedFractionSize)) return
    if (.NOT.allocated(dGridPreviousSeedFraction)) return

    ! Step 1. Validate each stored phase name and constitution against the
    ! current system and per-phase point budget.
    allocate(nAcceptedPhase(nSolnPhasesSys))
    nAcceptedPhase = 0
    do iSeed = 1, nGridPreviousSeed
        iPhaseIndex = FindGridPhaseByName(cGridPreviousSeedPhaseName(iSeed))
        if (.NOT.PreviousSeedIsCompatible(iSeed,iPhaseIndex)) then
            nGridPreviousSeedDropped = nGridPreviousSeedDropped + 1
            cycle
        end if
        if (nGridPhasePoint(iPhaseIndex)+nAcceptedPhase(iPhaseIndex) >= &
            GRID_PHASE_POINT_CAP) then
            nGridPreviousSeedDropped = nGridPreviousSeedDropped + 1
            cycle
        end if
        nLocal = iGridPreviousSeedFractionSize(iSeed)
        nAcceptedPhase(iPhaseIndex) = nAcceptedPhase(iPhaseIndex) + 1
        nAccepted = nAccepted + 1
        nFractionAccepted = nFractionAccepted + nLocal
    end do
    deallocate(nAcceptedPhase)
    nGridPreviousSeedAccepted = nAccepted

    return
end subroutine CountCompatiblePreviousSeeds


!> \brief Append accepted previous active constitutions after generated rows.
subroutine AppendPreviousGridSeeds(iPointCursor, iFractionCursor)
    implicit none

    integer, intent(inout) :: iPointCursor, iFractionCursor
    integer :: iSeed, iPhaseIndex, iOffset, nLocal
    integer, dimension(:), allocatable :: nAcceptedPhase

    if (nGridPreviousSeedAccepted <= 0) return
    allocate(nAcceptedPhase(nSolnPhasesSys))
    nAcceptedPhase = 0
    do iSeed = 1, nGridPreviousSeed
        iPhaseIndex = FindGridPhaseByName(cGridPreviousSeedPhaseName(iSeed))
        if (.NOT.PreviousSeedIsCompatible(iSeed,iPhaseIndex)) cycle
        if (nGridPhasePoint(iPhaseIndex)+nAcceptedPhase(iPhaseIndex) >= &
            GRID_PHASE_POINT_CAP) cycle
        iOffset = iGridPreviousSeedFractionOffset(iSeed)
        nLocal = iGridPreviousSeedFractionSize(iSeed)
        call StoreGridPoint(iPhaseIndex, GridDisplayPhase(iPhaseIndex), &
            dGridPreviousSeedFraction(iOffset:iOffset+nLocal-1), &
            iPointCursor, iFractionCursor)
        nAcceptedPhase(iPhaseIndex) = nAcceptedPhase(iPhaseIndex) + 1
    end do
    deallocate(nAcceptedPhase)

    return
end subroutine AppendPreviousGridSeeds


!> \brief Test whether one stored constitution belongs to the current phase model.
logical function PreviousSeedIsCompatible(iSeed, iPhaseIndex)
    implicit none

    integer, intent(in) :: iSeed, iPhaseIndex
    integer :: nLocal, iOffset
    integer :: OrderDisorderProjectionTopologyClass

    PreviousSeedIsCompatible = .FALSE.
    if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) return
    if ((iSeed < 1).OR.(iSeed > nGridPreviousSeed)) return
    if (allocated(iODCompanionPhase)) then
        if ((iPhaseIndex <= SIZE(iODCompanionPhase)).AND.&
            (iODCompanionPhase(iPhaseIndex) > 0).AND.&
            (OrderDisorderProjectionTopologyClass(iPhaseIndex) /= &
                OD_PROJECTION_TOPOLOGY_UNSUPPORTED)) return
    end if
    nLocal = iGridPreviousSeedFractionSize(iSeed)
    if (nLocal /= nSpeciesPhase(iPhaseIndex)-nSpeciesPhase(iPhaseIndex-1)) return
    iOffset = iGridPreviousSeedFractionOffset(iSeed)
    if ((iOffset < 1).OR.(iOffset+nLocal-1 > SIZE(dGridPreviousSeedFraction))) return
    if (SUM(dGridPreviousSeedFraction(iOffset:iOffset+nLocal-1)) <= 1D-300) return
    PreviousSeedIsCompatible = .TRUE.

    return
end function PreviousSeedIsCompatible


!> \brief Resolve a stored phase name to the current system phase index.
integer function FindGridPhaseByName(cPhaseName)
    implicit none

    character(*), intent(in) :: cPhaseName
    integer :: iPhaseIndex

    FindGridPhaseByName = 0
    do iPhaseIndex = 1, nSolnPhasesSys
        if (TRIM(cSolnPhaseName(iPhaseIndex)) /= TRIM(cPhaseName)) cycle
        FindGridPhaseByName = iPhaseIndex
        return
    end do

    return
end function FindGridPhaseByName


!> \brief Capture converged active solution constitutions for a later solve.
!!
!! \details Stores phase names rather than relying on transient system indices,
!! and preserves one normalized constitution for every positive active set.
subroutine CapturePreviousGridSeeds
    implicit none

    integer :: iSlot, iSeed, iPhaseIndex, iFirst, iLast, iOffset, nFractionTotal
    real(8) :: dNorm

    if (.NOT.lGridPreviousSolutionSeedActive) return
    if (INFOThermo /= 0) return
    if (.NOT.lConverged) return

    if (allocated(iGridPreviousSeedPhase)) deallocate(iGridPreviousSeedPhase)
    if (allocated(iGridPreviousSeedFractionOffset)) deallocate(iGridPreviousSeedFractionOffset)
    if (allocated(iGridPreviousSeedFractionSize)) deallocate(iGridPreviousSeedFractionSize)
    if (allocated(dGridPreviousSeedFraction)) deallocate(dGridPreviousSeedFraction)
    if (allocated(cGridPreviousSeedPhaseName)) deallocate(cGridPreviousSeedPhaseName)
    nGridPreviousSeed = 0
    nFractionTotal = 0

    ! Step 1. Count active solution sets and their ragged constitution sizes.
    do iSlot = 1, nElements
        if (iAssemblage(iSlot) >= 0) cycle
        if (dMolesPhase(iSlot) <= 0D0) cycle
        iPhaseIndex = -iAssemblage(iSlot)
        if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) cycle
        iFirst = nSpeciesPhase(iPhaseIndex-1) + 1
        iLast = nSpeciesPhase(iPhaseIndex)
        dNorm = 0D0
        if (allocated(dActiveSlotMolFraction)) &
            dNorm = SUM(dActiveSlotMolFraction(iSlot,iFirst:iLast))
        if (dNorm <= 1D-300) dNorm = SUM(dMolFraction(iFirst:iLast))
        if (dNorm <= 1D-300) cycle
        nGridPreviousSeed = nGridPreviousSeed + 1
        nFractionTotal = nFractionTotal + iLast-iFirst+1
    end do
    if (nGridPreviousSeed <= 0) return

    ! Step 2. Allocate an exact-size snapshot.
    allocate(iGridPreviousSeedPhase(nGridPreviousSeed))
    allocate(iGridPreviousSeedFractionOffset(nGridPreviousSeed))
    allocate(iGridPreviousSeedFractionSize(nGridPreviousSeed))
    allocate(cGridPreviousSeedPhaseName(nGridPreviousSeed))
    allocate(dGridPreviousSeedFraction(nFractionTotal))
    iGridPreviousSeedPhase = 0
    iGridPreviousSeedFractionOffset = 0
    iGridPreviousSeedFractionSize = 0
    cGridPreviousSeedPhaseName = ' '
    dGridPreviousSeedFraction = 0D0

    ! Step 3. Store normalized slot-local constitutions in active-slot order.
    iSeed = 0
    iOffset = 1
    do iSlot = 1, nElements
        if (iAssemblage(iSlot) >= 0) cycle
        if (dMolesPhase(iSlot) <= 0D0) cycle
        iPhaseIndex = -iAssemblage(iSlot)
        if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) cycle
        iFirst = nSpeciesPhase(iPhaseIndex-1) + 1
        iLast = nSpeciesPhase(iPhaseIndex)
        dNorm = 0D0
        if (allocated(dActiveSlotMolFraction)) &
            dNorm = SUM(dActiveSlotMolFraction(iSlot,iFirst:iLast))
        if (dNorm <= 1D-300) dNorm = SUM(dMolFraction(iFirst:iLast))
        if (dNorm <= 1D-300) cycle

        iSeed = iSeed + 1
        iGridPreviousSeedPhase(iSeed) = iPhaseIndex
        iGridPreviousSeedFractionOffset(iSeed) = iOffset
        iGridPreviousSeedFractionSize(iSeed) = iLast-iFirst+1
        cGridPreviousSeedPhaseName(iSeed) = cSolnPhaseName(iPhaseIndex)
        if (allocated(dActiveSlotMolFraction).AND.&
            (SUM(dActiveSlotMolFraction(iSlot,iFirst:iLast)) > 1D-300)) then
            dGridPreviousSeedFraction(iOffset:iOffset+iLast-iFirst) = &
                dActiveSlotMolFraction(iSlot,iFirst:iLast) / &
                SUM(dActiveSlotMolFraction(iSlot,iFirst:iLast))
        else
            dGridPreviousSeedFraction(iOffset:iOffset+iLast-iFirst) = &
                dMolFraction(iFirst:iLast) / SUM(dMolFraction(iFirst:iLast))
        end if
        iOffset = iOffset + iLast-iFirst+1
    end do

    return
end subroutine CapturePreviousGridSeeds


!> \brief Classify one solution phase for constitution-space sampling.
!!
!! \details Distinguishes ordinary simplexes, multi-sublattice endmember
!! mixtures, and supported order/disorder manifolds without phase-name rules.
subroutine ClassifyGridPhase(iPhaseIndex, iTargetPhase, iDefinitionPhase, iClass, &
    iSubstitutionalSublattice, iInterstitialSublattice, nSimplexConstituent, &
    lSupported, lHelper)
    implicit none

    integer, intent(in) :: iPhaseIndex
    integer, intent(out) :: iTargetPhase, iDefinitionPhase, iClass
    integer, intent(out) :: iSubstitutionalSublattice, iInterstitialSublattice
    integer, intent(out) :: nSimplexConstituent
    logical, intent(out) :: lSupported, lHelper

    integer :: iPhaseLocal, iSublatticePhaseLocal, iSublattice
    integer :: nMixed, iMixedFirst, iMixedSecond
    logical :: lFirstInterstitial, lSecondInterstitial
    integer :: OrderDisorderProjectionTopologyClass

    iTargetPhase = iPhaseIndex
    iDefinitionPhase = iPhaseIndex
    iClass = GRID_PHASE_CLASS_NONE
    iSubstitutionalSublattice = 0
    iInterstitialSublattice = 0
    nSimplexConstituent = 0
    lSupported = .FALSE.
    lHelper = .FALSE.

    if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) return

    ! Step 1. Skip structural helper phases that are sampled through their parent.
    if (allocated(iODCompanionPhase)) then
        do iPhaseLocal = 1, MIN(nSolnPhasesSys, SIZE(iODCompanionPhase))
            if (iODCompanionPhase(iPhaseLocal) /= iPhaseIndex) cycle
            if (iPhaseLocal == iPhaseIndex) cycle
            if (allocated(iODStandalonePhase)) then
                if (iPhaseLocal <= SIZE(iODStandalonePhase)) then
                    if (iODStandalonePhase(iPhaseLocal) == iPhaseIndex) cycle
                end if
            end if
            lHelper = .TRUE.
            return
        end do
    end if

    ! Step 2. Map a supported ordered parent onto its disordered definition.
    if (allocated(iODCompanionPhase)) then
        if (iPhaseIndex <= SIZE(iODCompanionPhase)) then
            if ((iODCompanionPhase(iPhaseIndex) > 0).AND.&
                (OrderDisorderProjectionTopologyClass(iPhaseIndex) /= &
                    OD_PROJECTION_TOPOLOGY_UNSUPPORTED)) then
                iDefinitionPhase = iODCompanionPhase(iPhaseIndex)
                iClass = GRID_PHASE_CLASS_DISORDERED_MANIFOLD
            end if
        end if
    end if

    ! Step 3. Classify non-sublattice models directly from their endmembers.
    if ((TRIM(cSolnPhaseType(iDefinitionPhase)) /= 'SUBL').AND.&
        (TRIM(cSolnPhaseType(iDefinitionPhase)) /= 'SUBLM').AND.&
        (TRIM(cSolnPhaseType(iDefinitionPhase)) /= 'SUBOM')) then
        if (iClass == GRID_PHASE_CLASS_DISORDERED_MANIFOLD) then
            nSimplexConstituent = nSpeciesPhase(iDefinitionPhase) - &
                nSpeciesPhase(iDefinitionPhase-1)
            lSupported = (nSimplexConstituent > 0)
        else if ((TRIM(cSolnPhaseType(iDefinitionPhase)) == 'IDMX').OR.&
            (TRIM(cSolnPhaseType(iDefinitionPhase)) == 'QKTO').OR.&
            (TRIM(cSolnPhaseType(iDefinitionPhase)) == 'QKTOM').OR.&
            (TRIM(cSolnPhaseType(iDefinitionPhase)) == 'RKMP').OR.&
            (TRIM(cSolnPhaseType(iDefinitionPhase)) == 'RKMPM')) then
            iClass = GRID_PHASE_CLASS_SIMPLEX
            nSimplexConstituent = nSpeciesPhase(iDefinitionPhase) - &
                nSpeciesPhase(iDefinitionPhase-1)
            lSupported = (nSimplexConstituent > 0)
        else
            iClass = GRID_PHASE_CLASS_ENDMEMBER_MIX
            lSupported = (nSpeciesPhase(iDefinitionPhase) > nSpeciesPhase(iDefinitionPhase-1))
        end if
        return
    end if

    ! Step 4. Identify the substitutional simplex and optional interstitial axis.
    iSublatticePhaseLocal = iPhaseSublattice(iDefinitionPhase)
    if ((iSublatticePhaseLocal < 1).OR.(iSublatticePhaseLocal > nCountSublattice)) return

    nMixed = 0
    iMixedFirst = 0
    iMixedSecond = 0
    do iSublattice = 1, nSublatticePhase(iSublatticePhaseLocal)
        if (nConstituentSublattice(iSublatticePhaseLocal,iSublattice) <= 1) cycle
        nMixed = nMixed + 1
        if (nMixed == 1) iMixedFirst = iSublattice
        if (nMixed == 2) iMixedSecond = iSublattice
    end do

    if (nMixed == 1) then
        iSubstitutionalSublattice = iMixedFirst
    else if (nMixed == 2) then
        lFirstInterstitial = IsInterstitialSublattice(iSublatticePhaseLocal, iMixedFirst)
        lSecondInterstitial = IsInterstitialSublattice(iSublatticePhaseLocal, iMixedSecond)
        if (lFirstInterstitial.EQV.lSecondInterstitial) then
            if (iClass /= GRID_PHASE_CLASS_DISORDERED_MANIFOLD) then
                iClass = GRID_PHASE_CLASS_ENDMEMBER_MIX
                lSupported = .TRUE.
            end if
            return
        end if
        if (lFirstInterstitial) then
            iInterstitialSublattice = iMixedFirst
            iSubstitutionalSublattice = iMixedSecond
        else
            iInterstitialSublattice = iMixedSecond
            iSubstitutionalSublattice = iMixedFirst
        end if
    else
        if (iClass /= GRID_PHASE_CLASS_DISORDERED_MANIFOLD) then
            iClass = GRID_PHASE_CLASS_ENDMEMBER_MIX
            lSupported = .TRUE.
        end if
        return
    end if

    nSimplexConstituent = GridSimplexConstituentCount(iTargetPhase, iDefinitionPhase, &
        iSublatticePhaseLocal, iSubstitutionalSublattice)
    if ((nSimplexConstituent <= 0).AND.&
        (iClass /= GRID_PHASE_CLASS_DISORDERED_MANIFOLD)) &
        nSimplexConstituent = nConstituentSublattice(&
            iSublatticePhaseLocal,iSubstitutionalSublattice)
    if (iClass /= GRID_PHASE_CLASS_DISORDERED_MANIFOLD) iClass = GRID_PHASE_CLASS_SIMPLEX
    lSupported = (nSimplexConstituent > 0)

    return
end subroutine ClassifyGridPhase


!> \brief Identify a binary vacancy/non-vacancy interstitial sublattice.
logical function IsInterstitialSublattice(iSublatticePhase, iSublattice)
    implicit none

    integer, intent(in) :: iSublatticePhase, iSublattice
    integer :: iConstituent
    character(8) :: cName

    IsInterstitialSublattice = .FALSE.
    if (nConstituentSublattice(iSublatticePhase,iSublattice) /= 2) return
    do iConstituent = 1, 2
        cName = UpperGridName(cConstituentNameSUB(&
            iSublatticePhase,iSublattice,iConstituent))
        if ((TRIM(cName) == 'VA').OR.(TRIM(cName) == 'VACANCY')) then
            IsInterstitialSublattice = .TRUE.
            return
        end if
    end do

    return
end function IsInterstitialSublattice


!> \brief Select complete grid layers under the per-phase point cap.
!!
!! \details Coarsens an entire dimensional layer to the next prescribed rung;
!! it never silently truncates individual faces or sequence points.
subroutine PlanGridPhase(iTargetPhase, iDefinitionPhase, iClass, &
    iSubstitutionalSublattice, iInterstitialSublattice, nSimplexConstituent, &
    nPhasePointOut)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase, iClass
    integer, intent(in) :: iSubstitutionalSublattice, iInterstitialSublattice
    integer, intent(in) :: nSimplexConstituent
    integer, intent(out) :: nPhasePointOut

    integer :: nEndmember, nDepth, nCross, nRemaining, nLayerNominal, nLayerSelected
    integer :: nFace, nInteractionFace, nFreeFace, iRung
    integer, dimension(5) :: nEdgeRung
    integer, dimension(4) :: nTernaryMeshRung, nTernarySobolRung
    integer, dimension(7) :: nSobolRung

    nPhasePointOut = 0
    iGridPhaseLayerRung(iTargetPhase,:) = 0
    ! Step 1. Apply the fixed endmember-mixing depth schedule to general
    ! multi-sublattice phases.
    if (iClass == GRID_PHASE_CLASS_ENDMEMBER_MIX) then
        nEndmember = nSpeciesPhase(iTargetPhase) - nSpeciesPhase(iTargetPhase-1)
        if (nEndmember <= 4) then
            nDepth = 5
        else if (nEndmember <= 7) then
            nDepth = 4
        else if (nEndmember <= 13) then
            nDepth = 3
        else if (nEndmember <= 50) then
            nDepth = 2
        else
            nDepth = 1
        end if
        nPhasePointOut = nEndmember**nDepth
        iGridPhaseLayerRung(iTargetPhase,GRID_LAYER_VERTICES) = nDepth
        iGridPhaseStatus(iTargetPhase) = GRID_PHASE_STATUS_GENERATED
        return
    end if

    ! Step 2. Cross every substitutional point with the five interstitial levels.
    nCross = 1
    if (iInterstitialSublattice > 0) nCross = 5
    nRemaining = GRID_PHASE_POINT_CAP
    ! Step 3. Preserve all vertices or reject the phase as a typed overflow.
    nLayerSelected = nSimplexConstituent * nCross
    if (nLayerSelected > nRemaining) then
        iGridPhaseStatus(iTargetPhase) = GRID_PHASE_STATUS_VERTEX_OVERFLOW
        nGridBudgetDroppedPoint = nGridBudgetDroppedPoint + nLayerSelected
        return
    end if
    nPhasePointOut = nLayerSelected
    nRemaining = nRemaining - nLayerSelected
    iGridPhaseLayerRung(iTargetPhase,GRID_LAYER_VERTICES) = 1

    ! Step 4. Add binary edges, including the dilute near-vertex points.
    nEdgeRung = (/17, 9, 4, 1, 0/)
    nFace = BinomialCount(nSimplexConstituent,2)
    nLayerNominal = nFace * nEdgeRung(1) * nCross
    call SelectUniformLayerRung(nFace, nCross, nEdgeRung, nRemaining, iRung, nLayerSelected)
    call RecordLayerPlan(iTargetPhase, GRID_LAYER_EDGES, iRung, nLayerNominal, &
        nLayerSelected, nRemaining, nPhasePointOut)

    ! Step 5. Use meshes on interaction-bearing ternary faces and low-discrepancy
    ! points on parameter-free faces.
    nInteractionFace = CountTernaryInteractionFaces(iTargetPhase, iDefinitionPhase, &
        iSubstitutionalSublattice, nSimplexConstituent)
    nFace = BinomialCount(nSimplexConstituent,3)
    nFreeFace = MAX(0, nFace - nInteractionFace)
    nTernaryMeshRung = (/36, 6, 1, 0/)
    nTernarySobolRung = (/32, 16, 8, 0/)
    nLayerNominal = (nInteractionFace*nTernaryMeshRung(1) + &
        nFreeFace*nTernarySobolRung(1))*nCross
    iRung = 4
    nLayerSelected = 0
    do nDepth = 1, 4
        nLayerSelected = (nInteractionFace*nTernaryMeshRung(nDepth) + &
            nFreeFace*nTernarySobolRung(nDepth))*nCross
        if (nLayerSelected <= nRemaining) then
            iRung = nDepth
            exit
        end if
    end do
    call RecordLayerPlan(iTargetPhase, GRID_LAYER_TERNARY, iRung, nLayerNominal, &
        nLayerSelected, nRemaining, nPhasePointOut)

    ! Step 6. Add low-discrepancy quaternary faces.
    nSobolRung = (/32, 16, 8, 4, 2, 1, 0/)
    nFace = BinomialCount(nSimplexConstituent,4)
    nLayerNominal = nFace*nSobolRung(1)*nCross
    call SelectUniformLayerRung(nFace, nCross, nSobolRung, nRemaining, iRung, nLayerSelected)
    call RecordLayerPlan(iTargetPhase, GRID_LAYER_QUATERNARY, iRung, nLayerNominal, &
        nLayerSelected, nRemaining, nPhasePointOut)

    ! Step 7. Sprinkle the full simplex only when five or more constituents exist.
    if (nSimplexConstituent >= 5) then
        nSobolRung = (/256, 128, 64, 32, 16, 8, 0/)
        nLayerNominal = nSobolRung(1)*nCross
        call SelectUniformLayerRung(1, nCross, nSobolRung, nRemaining, iRung, nLayerSelected)
        call RecordLayerPlan(iTargetPhase, GRID_LAYER_INTERIOR, iRung, nLayerNominal, &
            nLayerSelected, nRemaining, nPhasePointOut)
    end if

    iGridPhaseStatus(iTargetPhase) = GRID_PHASE_STATUS_GENERATED

    return
end subroutine PlanGridPhase


!> \brief Choose the finest uniform point-count rung that fits one layer.
subroutine SelectUniformLayerRung(nFace, nCross, nRungPoint, nRemaining, iRung, nSelected)
    implicit none

    integer, intent(in) :: nFace, nCross, nRemaining
    integer, dimension(:), intent(in) :: nRungPoint
    integer, intent(out) :: iRung, nSelected
    integer :: i

    iRung = SIZE(nRungPoint)
    nSelected = 0
    do i = 1, SIZE(nRungPoint)
        nSelected = nFace*nRungPoint(i)*nCross
        if (nSelected <= nRemaining) then
            iRung = i
            return
        end if
    end do

    return
end subroutine SelectUniformLayerRung


!> \brief Record one selected dimensional layer and its typed coarsening census.
subroutine RecordLayerPlan(iPhaseIndex, iLayer, iRung, nNominal, nSelected, &
    nRemaining, nTotal)
    implicit none

    integer, intent(in) :: iPhaseIndex, iLayer, iRung, nNominal, nSelected
    integer, intent(inout) :: nRemaining, nTotal

    iGridPhaseLayerRung(iPhaseIndex,iLayer) = iRung
    nTotal = nTotal + nSelected
    nRemaining = nRemaining - nSelected
    if (nSelected < nNominal) then
        nGridBudgetCoarsening = nGridBudgetCoarsening + 1
        nGridBudgetDroppedPoint = nGridBudgetDroppedPoint + nNominal - nSelected
    end if

    return
end subroutine RecordLayerPlan


!> \brief Emit the dimensional simplex ladder selected for one phase.
!!
!! \details Keeps the constituent and face loops visible so vertices, edges,
!! ternary faces, quaternary faces, and full-simplex points can be audited in
!! their physical generation order.
subroutine FillSimplexGrid(iTargetPhase, iDefinitionPhase, iSubstitutionalSublattice, &
    iInterstitialSublattice, nConstituent, iPointCursor, iFractionCursor)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase
    integer, intent(in) :: iSubstitutionalSublattice, iInterstitialSublattice, nConstituent
    integer, intent(inout) :: iPointCursor, iFractionCursor

    integer :: i, j, k, l, iStep, iSobol, iRung, nPointLayer, nMeshDenom
    real(8), dimension(nConstituent) :: dSimplex
    real(8), dimension(4) :: dFaceFraction
    real(8), dimension(17), parameter :: dDenseEdgeFraction = (/&
        0.01D0, 0.02D0, 0.05D0, 0.08D0, 0.10D0, 0.20D0, 0.30D0, 0.40D0, 0.50D0, &
        0.60D0, 0.70D0, 0.80D0, 0.90D0, 0.92D0, 0.95D0, 0.98D0, 0.99D0/)
    logical :: lInteraction

    ! Step 1. Emit every pure substitutional vertex.
    do i = 1, nConstituent
        dSimplex = 0D0
        dSimplex(i) = 1D0
        call EmitSimplexPoint(iTargetPhase, iDefinitionPhase, iSubstitutionalSublattice, &
            iInterstitialSublattice, dSimplex, iPointCursor, iFractionCursor)
    end do

    ! Step 2. Walk every binary edge with the selected fixed weight ladder.
    iRung = iGridPhaseLayerRung(iTargetPhase,GRID_LAYER_EDGES)
    do i = 1, nConstituent - 1
        do j = i + 1, nConstituent
            select case (iRung)
                case (1)
                    nPointLayer = 17
                case (2)
                    nPointLayer = 9
                case (3)
                    nPointLayer = 4
                case (4)
                    nPointLayer = 1
                case default
                    nPointLayer = 0
            end select
            do iStep = 1, nPointLayer
                dSimplex = 0D0
                if (iRung == 1) then
                    dSimplex(j) = dDenseEdgeFraction(iStep)
                else if (iRung == 2) then
                    dSimplex(j) = DBLE(iStep)/10D0
                else if (iRung == 3) then
                    dSimplex(j) = DBLE(2*iStep)/10D0
                else
                    dSimplex(j) = 0.5D0
                end if
                dSimplex(i) = 1D0 - dSimplex(j)
                call EmitSimplexPoint(iTargetPhase, iDefinitionPhase, &
                    iSubstitutionalSublattice, iInterstitialSublattice, dSimplex, &
                    iPointCursor, iFractionCursor)
            end do
        end do
    end do

    ! Step 3. Walk every ternary face, choosing its mesh from database evidence.
    iRung = iGridPhaseLayerRung(iTargetPhase,GRID_LAYER_TERNARY)
    do i = 1, nConstituent - 2
        do j = i + 1, nConstituent - 1
            do k = j + 1, nConstituent
                lInteraction = HasTernaryInteraction(iTargetPhase, iDefinitionPhase, &
                    iSubstitutionalSublattice, i, j, k)
                if (lInteraction) then
                    select case (iRung)
                        case (1)
                            nMeshDenom = 10
                        case (2)
                            nMeshDenom = 5
                        case (3)
                            nMeshDenom = 3
                        case default
                            nMeshDenom = 0
                    end select
                    do iStep = 1, MAX(0,nMeshDenom-2)
                        do l = 1, nMeshDenom-iStep-1
                            dSimplex = 0D0
                            dSimplex(i) = DBLE(iStep)/DBLE(nMeshDenom)
                            dSimplex(j) = DBLE(l)/DBLE(nMeshDenom)
                            dSimplex(k) = 1D0 - dSimplex(i) - dSimplex(j)
                            call EmitSimplexPoint(iTargetPhase, iDefinitionPhase, &
                                iSubstitutionalSublattice, iInterstitialSublattice, dSimplex, &
                                iPointCursor, iFractionCursor)
                        end do
                    end do
                else
                    select case (iRung)
                        case (1)
                            nPointLayer = 32
                        case (2)
                            nPointLayer = 16
                        case (3)
                            nPointLayer = 8
                        case default
                            nPointLayer = 0
                    end select
                    do iSobol = 1, nPointLayer
                        call SobolSimplexPoint(iSobol, 3, dFaceFraction(1:3))
                        dSimplex = 0D0
                        dSimplex(i) = dFaceFraction(1)
                        dSimplex(j) = dFaceFraction(2)
                        dSimplex(k) = dFaceFraction(3)
                        call EmitSimplexPoint(iTargetPhase, iDefinitionPhase, &
                            iSubstitutionalSublattice, iInterstitialSublattice, dSimplex, &
                            iPointCursor, iFractionCursor)
                    end do
                end if
            end do
        end do
    end do

    ! Step 4. Walk every quaternary face with sorted-spacings Sobol points.
    iRung = iGridPhaseLayerRung(iTargetPhase,GRID_LAYER_QUATERNARY)
    nPointLayer = SobolRungPointCount(iRung, 32)
    do i = 1, nConstituent - 3
        do j = i + 1, nConstituent - 2
            do k = j + 1, nConstituent - 1
                do l = k + 1, nConstituent
                    do iSobol = 1, nPointLayer
                        call SobolSimplexPoint(iSobol, 4, dFaceFraction)
                        dSimplex = 0D0
                        dSimplex(i) = dFaceFraction(1)
                        dSimplex(j) = dFaceFraction(2)
                        dSimplex(k) = dFaceFraction(3)
                        dSimplex(l) = dFaceFraction(4)
                        call EmitSimplexPoint(iTargetPhase, iDefinitionPhase, &
                            iSubstitutionalSublattice, iInterstitialSublattice, dSimplex, &
                            iPointCursor, iFractionCursor)
                    end do
                end do
            end do
        end do
    end do

    ! Step 5. Cover the full five-plus simplex with the selected Sobol sprinkle.
    iRung = iGridPhaseLayerRung(iTargetPhase,GRID_LAYER_INTERIOR)
    nPointLayer = SobolRungPointCount(iRung, 256)
    do iSobol = 1, nPointLayer
        call SobolSimplexPoint(iSobol, nConstituent, dSimplex)
        call EmitSimplexPoint(iTargetPhase, iDefinitionPhase, &
            iSubstitutionalSublattice, iInterstitialSublattice, dSimplex, &
            iPointCursor, iFractionCursor)
    end do

    return
end subroutine FillSimplexGrid


!> \brief Convert a stored Sobol budget rung to its fixed point count.
integer function SobolRungPointCount(iRung, nFull)
    implicit none

    integer, intent(in) :: iRung, nFull

    SobolRungPointCount = 0
    if (iRung <= 0) return
    if (nFull == 32) then
        select case (iRung)
            case (1); SobolRungPointCount = 32
            case (2); SobolRungPointCount = 16
            case (3); SobolRungPointCount = 8
            case (4); SobolRungPointCount = 4
            case (5); SobolRungPointCount = 2
            case (6); SobolRungPointCount = 1
        end select
    else
        select case (iRung)
            case (1); SobolRungPointCount = 256
            case (2); SobolRungPointCount = 128
            case (3); SobolRungPointCount = 64
            case (4); SobolRungPointCount = 32
            case (5); SobolRungPointCount = 16
            case (6); SobolRungPointCount = 8
        end select
    end if

    return
end function SobolRungPointCount


!> \brief Emit the fixed-depth endmember mixtures for a general sublattice phase.
subroutine FillEndmemberMixGrid(iPhaseIndex, iPointCursor, iFractionCursor)
    implicit none

    integer, intent(in) :: iPhaseIndex
    integer, intent(inout) :: iPointCursor, iFractionCursor

    integer :: nEndmember, nDepth, nCombination, iCombination, iDigit, iIndex
    integer, dimension(5) :: iEndmember
    real(8), dimension(5) :: dWeight
    real(8), dimension(:), allocatable :: dFraction

    ! Step 1. Select the fixed mixture depth from the endmember count.
    nEndmember = nSpeciesPhase(iPhaseIndex) - nSpeciesPhase(iPhaseIndex-1)
    if (nEndmember <= 4) then
        nDepth = 5
        dWeight = (/0.07D0, 0.28D0, 0.16D0, 0.45D0, 0.04D0/)
    else if (nEndmember <= 7) then
        nDepth = 4
        dWeight = (/0.07D0, 0.28D0, 0.16D0, 0.49D0, 0D0/)
    else if (nEndmember <= 13) then
        nDepth = 3
        dWeight = (/0.07D0, 0.28D0, 0.65D0, 0D0, 0D0/)
    else if (nEndmember <= 50) then
        nDepth = 2
        dWeight = (/0.07D0, 0.93D0, 0D0, 0D0, 0D0/)
    else
        nDepth = 1
        dWeight = (/1D0, 0D0, 0D0, 0D0, 0D0/)
    end if

    ! Step 2. Enumerate every ordered endmember tuple at that depth and combine
    ! repeated tuple members into one normalized constitution.
    nCombination = nEndmember**nDepth
    allocate(dFraction(nEndmember))
    do iCombination = 0, nCombination - 1
        iIndex = iCombination
        iEndmember = 1
        do iDigit = nDepth, 1, -1
            iEndmember(iDigit) = MOD(iIndex,nEndmember) + 1
            iIndex = iIndex / nEndmember
        end do
        dFraction = 0D0
        do iDigit = 1, nDepth
            dFraction(iEndmember(iDigit)) = dFraction(iEndmember(iDigit)) + dWeight(iDigit)
        end do
        call StoreGridPoint(iPhaseIndex, iPhaseIndex, dFraction, iPointCursor, iFractionCursor)
    end do
    deallocate(dFraction)

    return
end subroutine FillEndmemberMixGrid


!> \brief Cross one substitutional simplex point with interstitial occupancy levels.
subroutine EmitSimplexPoint(iTargetPhase, iDefinitionPhase, iSubstitutionalSublattice, &
    iInterstitialSublattice, dSimplex, iPointCursor, iFractionCursor)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase
    integer, intent(in) :: iSubstitutionalSublattice, iInterstitialSublattice
    real(8), dimension(:), intent(in) :: dSimplex
    integer, intent(inout) :: iPointCursor, iFractionCursor

    integer :: iLevel, nLevel
    real(8), dimension(5) :: dInterstitialLevel
    real(8), dimension(:), allocatable :: dTargetFraction
    logical :: lMapped
    integer :: iProjectionStatus

    if (INFOThermo /= 0) return
    dInterstitialLevel = (/1D-12, 0.25D0, 0.5D0, 0.75D0, 1D0-1D-12/)
    nLevel = 1
    if (iInterstitialSublattice > 0) nLevel = 5
    allocate(dTargetFraction(nSpeciesPhase(iTargetPhase)-nSpeciesPhase(iTargetPhase-1)))
    do iLevel = 1, nLevel
        call MapSimplexToTargetFraction(iTargetPhase, iDefinitionPhase, &
            iSubstitutionalSublattice, iInterstitialSublattice, dSimplex, &
            dInterstitialLevel(iLevel), dTargetFraction, lMapped, iProjectionStatus)
        if (.NOT.lMapped) then
            call RecordGridProjectionFailure(iTargetPhase, iProjectionStatus)
            call SetStaticGridWorkflowError
            deallocate(dTargetFraction)
            return
        end if
        call StoreGridPoint(iTargetPhase, GridDisplayPhase(iTargetPhase), dTargetFraction, &
            iPointCursor, iFractionCursor)
    end do
    deallocate(dTargetFraction)

    return
end subroutine EmitSimplexPoint


!> \brief Map a simplex point into the target phase endmember basis.
!!
!! \details Builds product endmember fractions for sublattice models and uses
!! the structural order/disorder projection when the sampled definition is a
!! disordered companion of the thermodynamic parent.
subroutine MapSimplexToTargetFraction(iTargetPhase, iDefinitionPhase, &
    iSubstitutionalSublattice, iInterstitialSublattice, dSimplex, &
    dInterstitialOccupancy, dTargetFraction, lMapped, iProjectionStatus)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase
    integer, intent(in) :: iSubstitutionalSublattice, iInterstitialSublattice
    real(8), dimension(:), intent(in) :: dSimplex
    real(8), intent(in) :: dInterstitialOccupancy
    real(8), dimension(:), intent(out) :: dTargetFraction
    logical, intent(out) :: lMapped
    integer, intent(out) :: iProjectionStatus

    integer :: iFirst, iLast, iLocal, iSublatticePhaseLocal, iSublattice
    integer :: iConstituent, iSimplexConstituent
    integer :: iInterstitialVacancy, iInterstitialNonVacancy
    real(8) :: dProduct, dNorm
    real(8), dimension(:,:), allocatable :: dSite
    logical :: lProjected

    dTargetFraction = 0D0
    lMapped = .FALSE.
    iProjectionStatus = OD_PROJECTION_INVALID_INPUT
    iFirst = nSpeciesPhase(iDefinitionPhase-1) + 1
    iLast = nSpeciesPhase(iDefinitionPhase)

    ! Step 1. Express the sampled constitution in the definition phase basis.
    if ((TRIM(cSolnPhaseType(iDefinitionPhase)) /= 'SUBL').AND.&
        (TRIM(cSolnPhaseType(iDefinitionPhase)) /= 'SUBLM').AND.&
        (TRIM(cSolnPhaseType(iDefinitionPhase)) /= 'SUBOM')) then
        if (SIZE(dSimplex) /= iLast-iFirst+1) return
        dMolFraction(iFirst:iLast) = dSimplex
    else
        iSublatticePhaseLocal = iPhaseSublattice(iDefinitionPhase)
        allocate(dSite(nSublatticePhase(iSublatticePhaseLocal),nMaxConstituentSys))
        dSite = 0D0
        do iSublattice = 1, nSublatticePhase(iSublatticePhaseLocal)
            dSite(iSublattice,1) = 1D0
        end do
        dSite(iSubstitutionalSublattice,:) = 0D0
        do iConstituent = 1, SIZE(dSimplex)
            iSimplexConstituent = GridSimplexConstituentID(iTargetPhase, &
                iDefinitionPhase, iSublatticePhaseLocal, iSubstitutionalSublattice, &
                iConstituent)
            if (iSimplexConstituent <= 0) then
                deallocate(dSite)
                return
            end if
            dSite(iSubstitutionalSublattice,iSimplexConstituent) = dSimplex(iConstituent)
        end do
        if (iInterstitialSublattice > 0) then
            iInterstitialVacancy = 0
            iInterstitialNonVacancy = 0
            do iConstituent = 1, nConstituentSublattice(&
                iSublatticePhaseLocal,iInterstitialSublattice)
                if (IsVacancyName(cConstituentNameSUB(iSublatticePhaseLocal, &
                    iInterstitialSublattice,iConstituent))) then
                    iInterstitialVacancy = iConstituent
                else
                    iInterstitialNonVacancy = iConstituent
                end if
            end do
            if ((iInterstitialVacancy <= 0).OR.(iInterstitialNonVacancy <= 0)) then
                deallocate(dSite)
                return
            end if
            dSite(iInterstitialSublattice,:) = 0D0
            dSite(iInterstitialSublattice,iInterstitialVacancy) = 1D0-dInterstitialOccupancy
            dSite(iInterstitialSublattice,iInterstitialNonVacancy) = dInterstitialOccupancy
        end if

        do iLocal = 1, iLast-iFirst+1
            dProduct = 1D0
            do iSublattice = 1, nSublatticePhase(iSublatticePhaseLocal)
                iConstituent = iConstituentSublattice(&
                    iSublatticePhaseLocal,iSublattice,iLocal)
                dProduct = dProduct*dSite(iSublattice,iConstituent)
            end do
            dMolFraction(iFirst+iLocal-1) = dProduct
        end do
        dNorm = SUM(dMolFraction(iFirst:iLast))
        if (dNorm <= 1D-300) then
            deallocate(dSite)
            return
        end if
        dMolFraction(iFirst:iLast) = dMolFraction(iFirst:iLast)/dNorm
        deallocate(dSite)
    end if

    ! Step 2. Copy an ordinary phase directly or project the disordered
    ! constitution into its ordered parent representation.
    if (iTargetPhase == iDefinitionPhase) then
        if (SIZE(dTargetFraction) /= iLast-iFirst+1) return
        dTargetFraction = dMolFraction(iFirst:iLast)
        lMapped = .TRUE.
        iProjectionStatus = OD_PROJECTION_SUCCESS
    else
        call ProjectOrderDisorderCompanionFraction(iDefinitionPhase, iTargetPhase, &
            SIZE(dTargetFraction), dTargetFraction, lProjected, iProjectionStatus)
        lMapped = lProjected.AND.(iProjectionStatus == OD_PROJECTION_SUCCESS)
    end if

    return
end subroutine MapSimplexToTargetFraction


!> \brief Count and type one unexpected companion-projection failure.
subroutine RecordGridProjectionFailure(iTargetPhase, iProjectionStatus)
    implicit none

    integer, intent(in) :: iTargetPhase, iProjectionStatus

    nGridProjectionFailure = nGridProjectionFailure + 1
    if (allocated(nGridPhaseProjectionFailure)) then
        if ((iTargetPhase >= 1).AND.(iTargetPhase <= SIZE(nGridPhaseProjectionFailure))) then
            nGridPhaseProjectionFailure(iTargetPhase) = &
                nGridPhaseProjectionFailure(iTargetPhase) + 1
        end if
    end if
    if ((iProjectionStatus >= OD_PROJECTION_INVALID_INPUT).AND.&
        (iProjectionStatus <= OD_PROJECTION_AMBIGUOUS_TOPOLOGY)) then
        nGridProjectionFailureStatus(iProjectionStatus) = &
            nGridProjectionFailureStatus(iProjectionStatus) + 1
    end if
    if (allocated(iGridPhaseStatus)) then
        if ((iTargetPhase >= 1).AND.(iTargetPhase <= SIZE(iGridPhaseStatus))) &
            iGridPhaseStatus(iTargetPhase) = GRID_PHASE_STATUS_PROJECTION_FAILED
    end if

    return
end subroutine RecordGridProjectionFailure


!> \brief Store one phase constitution in the deterministic ragged grid arrays.
subroutine StoreGridPoint(iPhaseIndex, iDisplayPhase, dFraction, iPointCursor, iFractionCursor)
    implicit none

    integer, intent(in) :: iPhaseIndex, iDisplayPhase
    real(8), dimension(:), intent(in) :: dFraction
    integer, intent(inout) :: iPointCursor, iFractionCursor

    integer :: nLocal

    nLocal = SIZE(dFraction)
    iPointCursor = iPointCursor + 1
    if ((iPointCursor > nGridPoint).OR.(iFractionCursor+nLocal-1 > nGridFraction)) then
        call SetStaticGridWorkflowError
        return
    end if
    iGridPointPhase(iPointCursor) = iPhaseIndex
    iGridPointDisplayPhase(iPointCursor) = iDisplayPhase
    iGridPointFractionOffset(iPointCursor) = iFractionCursor
    iGridPointFractionSize(iPointCursor) = nLocal
    dGridPointFraction(iFractionCursor:iFractionCursor+nLocal-1) = dFraction
    iFractionCursor = iFractionCursor + nLocal

    return
end subroutine StoreGridPoint


!> \brief Extend Leveling storage by exactly the number of static grid rows.
subroutine ExtendGridLevelingArrays
    implicit none

    integer :: nBaseLevel
    integer, dimension(:), allocatable :: iPhaseBase
    real(8), dimension(:), allocatable :: dPhasePotentialBase
    real(8), dimension(:,:), allocatable :: dStoichBase

    nBaseLevel = nSpeciesLevel
    if (nBaseLevel /= nSpecies) then
        call SetStaticGridWorkflowError
        return
    end if

    allocate(iPhaseBase(nSpecies), dPhasePotentialBase(nSpecies))
    allocate(dStoichBase(nSpecies,nElements))
    iPhaseBase = iPhaseLevel(:nSpecies)
    dPhasePotentialBase = dPhasePotential(:nSpecies)
    dStoichBase = dStoichSpeciesLevel(:nSpecies,:)

    nSpeciesLevel = nSpecies + nGridPoint
    call ExtendSampledLevelingThermoArrays

    if (allocated(iPhaseLevel)) deallocate(iPhaseLevel)
    if (allocated(dPhasePotential)) deallocate(dPhasePotential)
    if (allocated(dStoichSpeciesLevel)) deallocate(dStoichSpeciesLevel)
    if (allocated(lLevelingRowExcluded)) deallocate(lLevelingRowExcluded)
    allocate(iPhaseLevel(nSpeciesLevel), dPhasePotential(nSpeciesLevel))
    allocate(dStoichSpeciesLevel(nSpeciesLevel,nElements))
    allocate(lLevelingRowExcluded(nSpeciesLevel))
    iPhaseLevel = 0
    dPhasePotential = 5D9
    dStoichSpeciesLevel = 0D0
    lLevelingRowExcluded = .FALSE.
    iPhaseLevel(:nSpecies) = iPhaseBase
    dPhasePotential(:nSpecies) = dPhasePotentialBase
    dStoichSpeciesLevel(:nSpecies,:) = dStoichBase

    allocate(iGridPointFromLevel(nSpeciesLevel))
    iGridPointFromLevel = 0

    deallocate(iPhaseBase, dPhasePotentialBase, dStoichBase)

    return
end subroutine ExtendGridLevelingArrays


!> \brief Rebuild immutable row mappings after temporary dynamic rows are removed.
subroutine RestoreStaticGridLevelingRows
    implicit none

    integer :: iPoint, iLevelRow

    if (.NOT.lGridFrontEndActive) return
    if (nGridPoint <= 0) return
    if (nSpeciesLevel < nSpecies+nGridPoint) then
        call SetStaticGridWorkflowError
        return
    end if

    if (allocated(iGridPointFromLevel)) deallocate(iGridPointFromLevel)
    allocate(iGridPointFromLevel(nSpeciesLevel))
    iGridPointFromLevel = 0
    do iPoint = 1, nGridPoint
        iLevelRow = nSpecies + iPoint
        iGridPointLevelRow(iPoint) = iLevelRow
        iGridPointFromLevel(iLevelRow) = iPoint
    end do
    if (.NOT.lGridPoolValid) then
        call SetStaticGridWorkflowError
        return
    end if
    do iPoint = 1, nGridPoint
        iLevelRow = nSpecies + iPoint
        iGridPointLevelRow(iPoint) = iLevelRow
        iGridPointFromLevel(iLevelRow) = iPoint
        iPhaseLevel(iLevelRow) = iGridPointPhase(iPoint)
        dChemicalPotential(iLevelRow) = dGridPointGibbs(iPoint)
        dLevelingChemicalPotential(iLevelRow) = dGridPointGibbs(iPoint)
        dLevelingCompositionSpecies(iLevelRow,:) = dGridPointComposition(iPoint,:)
        dStoichSpeciesLevel(iLevelRow,:) = dGridPointStoich(iPoint,:)
    end do

    return
end subroutine RestoreStaticGridLevelingRows


!> \brief Return the primary dynamic tangent-row index for one solution phase.
integer function GridTangentRowIndex(iPhaseIndex)
    implicit none

    integer, intent(in) :: iPhaseIndex

    GridTangentRowIndex = 0
    if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) return
    GridTangentRowIndex = nSpecies + nLevelCandidateRowOffset + iPhaseIndex

    return
end function GridTangentRowIndex


!> \brief Return the secondary order/disorder partition-row index for one phase.
integer function GridPartitionRowIndex(iPhaseIndex)
    implicit none

    integer, intent(in) :: iPhaseIndex

    GridPartitionRowIndex = 0
    if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) return
    GridPartitionRowIndex = GridTangentRowIndex(iPhaseIndex) + nSolnPhasesSys

    return
end function GridPartitionRowIndex


!> \brief Evaluate Gibbs energy for every fixed grid constitution transactionally.
!!
!! \details The loop calls only the shared fixed-constitution phase evaluator.
!! All mutable thermodynamic work arrays are restored after the final row.
subroutine EvaluateStaticGridRows
    implicit none

    integer :: iPoint, iPhaseIndex, iFirst, iLast, iOffset, iLevelRow
    integer :: nPhasePoint, nPhaseSpecies, iRow, iMode, iKernelInfo
    real(8), dimension(:,:), allocatable :: dPhaseFraction
    real(8), dimension(:), allocatable :: dPhaseGibbs
    integer, dimension(:), allocatable :: iPhasePoint
    real(8), dimension(:), allocatable :: dMolFractionSave
    real(8), dimension(:), allocatable :: dChemicalPotentialSave
    real(8), dimension(:), allocatable :: dPartialTotalHSave, dPartialTotalSSave
    real(8), dimension(:), allocatable :: dPartialTotalCpSave
    real(8), dimension(:), allocatable :: dPartialExcessSave, dPartialHSave
    real(8), dimension(:), allocatable :: dPartialSSave, dPartialCpSave
    real(8), dimension(:), allocatable :: dMagGSave, dMagHSave, dMagSSave, dMagCpSave
    real(8), dimension(:), allocatable :: dGibbsPhaseSave
    real(8), dimension(:,:,:), allocatable :: dSiteFractionSave

    if (allocated(dGridPointGibbs)) deallocate(dGridPointGibbs)
    if (allocated(dGridPointComposition)) deallocate(dGridPointComposition)
    if (allocated(dGridPointStoich)) deallocate(dGridPointStoich)
    allocate(dGridPointGibbs(nGridPoint))
    allocate(dGridPointComposition(nGridPoint,nElements))
    allocate(dGridPointStoich(nGridPoint,nElements))
    dGridPointGibbs = 5D9
    dGridPointComposition = 0D0
    dGridPointStoich = 0D0

    ! Step 1. Snapshot every shared thermodynamic array modified by phase pricing.
    allocate(dMolFractionSave(nSpecies), dChemicalPotentialSave(nSpecies))
    allocate(dPartialTotalHSave(nSpecies), dPartialTotalSSave(nSpecies))
    allocate(dPartialTotalCpSave(nSpecies))
    allocate(dPartialExcessSave(nSpecies), dPartialHSave(nSpecies))
    allocate(dPartialSSave(nSpecies), dPartialCpSave(nSpecies))
    allocate(dMagGSave(nSpecies), dMagHSave(nSpecies), dMagSSave(nSpecies), dMagCpSave(nSpecies))
    allocate(dGibbsPhaseSave(nSolnPhasesSys))
    dMolFractionSave = dMolFraction
    dChemicalPotentialSave = dChemicalPotential(:nSpecies)
    dPartialTotalHSave = dPartialEnthalpy
    dPartialTotalSSave = dPartialEntropy
    dPartialTotalCpSave = dPartialHeatCapacity
    dPartialExcessSave = dPartialExcessGibbs
    dPartialHSave = dPartialEnthalpyXS
    dPartialSSave = dPartialEntropyXS
    dPartialCpSave = dPartialHeatCapacityXS
    dMagGSave = dMagGibbsEnergy
    dMagHSave = dMagEnthalpy
    dMagSSave = dMagEntropy
    dMagCpSave = dMagHeatCapacity
    dGibbsPhaseSave = dGibbsSolnPhase
    if (allocated(dSiteFraction)) then
        allocate(dSiteFractionSave(SIZE(dSiteFraction,1),SIZE(dSiteFraction,2),SIZE(dSiteFraction,3)))
        dSiteFractionSave = dSiteFraction
    end if

    ! Step 2. Evaluate each phase block with a direct-G kernel. A phase whose
    ! family or parameter topology is uncovered takes the typed scalar path.
    do iPhaseIndex = 1, nSolnPhasesSys
        nPhasePoint = COUNT(iGridPointPhase == iPhaseIndex)
        if (nPhasePoint <= 0) cycle
        iFirst = nSpeciesPhase(iPhaseIndex-1)+1
        iLast = nSpeciesPhase(iPhaseIndex)
        nPhaseSpecies = iLast-iFirst+1
        allocate(iPhasePoint(nPhasePoint), dPhaseFraction(nPhasePoint,nPhaseSpecies))
        allocate(dPhaseGibbs(nPhasePoint))
        iRow = 0
        do iPoint = 1, nGridPoint
            if (iGridPointPhase(iPoint) /= iPhaseIndex) cycle
            iRow = iRow+1
            iPhasePoint(iRow) = iPoint
            iOffset = iGridPointFractionOffset(iPoint)
            dPhaseFraction(iRow,:) = dGridPointFraction(iOffset:iOffset+nPhaseSpecies-1)
        end do

        if (lGridBatchKernelActive) then
            call EvaluateGridPhaseGibbsBatch(iPhaseIndex, dPhaseFraction, dPhaseGibbs, iMode, iKernelInfo)
        else
            iMode = GRID_EVALUATION_SCALAR_REFERENCE
            iKernelInfo = 0
            do iRow = 1, nPhasePoint
                call EvaluateGridPhaseGibbsScalar(iPhaseIndex, dPhaseFraction(iRow,:), &
                    dPhaseGibbs(iRow), iKernelInfo)
                if (iKernelInfo /= 0) exit
            end do
        end if

        if (iKernelInfo == 0) then
            iGridPhaseEvaluationMode(iPhaseIndex) = iMode
            if (lGridBatchKernelActive) then
                nGridBatchPhase = nGridBatchPhase+1
            else
                nGridScalarFallbackPhase = nGridScalarFallbackPhase+1
            end if
            do iRow = 1, nPhasePoint
                iPoint = iPhasePoint(iRow)
                iLevelRow = nSpecies+iPoint
                call RegisterGridLevelRow(iPoint,iLevelRow)
                call SetStaticGridLevelingRow(iPoint,iLevelRow,dPhaseGibbs(iRow))
                dGridPointGibbs(iPoint) = dPhaseGibbs(iRow)
                dGridPointComposition(iPoint,:) = dLevelingCompositionSpecies(iLevelRow,:)
                dGridPointStoich(iPoint,:) = dStoichSpeciesLevel(iLevelRow,:)
                iGridPointStatus(iPoint) = GRID_PHASE_STATUS_GENERATED
                iGridPointIdentityOrdinal(iPoint) = iPoint
            end do
        else
            if (iKernelInfo == 1) then
                iGridPhaseEvaluationMode(iPhaseIndex) = GRID_EVALUATION_SCALAR_UNSUPPORTED
            else
                iGridPhaseEvaluationMode(iPhaseIndex) = GRID_EVALUATION_SCALAR_KERNEL_ERROR
            end if
            nGridScalarFallbackPhase = nGridScalarFallbackPhase+1
            do iRow = 1, nPhasePoint
                call EvaluateScalarGridPoint(iPhasePoint(iRow))
                if (INFOThermo /= 0) exit
            end do
        end if
        deallocate(iPhasePoint,dPhaseFraction,dPhaseGibbs)
        if (INFOThermo /= 0) exit
    end do

    ! Step 3. Restore the incoming thermodynamic state exactly.
    dMolFraction = dMolFractionSave
    dChemicalPotential(:nSpecies) = dChemicalPotentialSave
    dPartialEnthalpy = dPartialTotalHSave
    dPartialEntropy = dPartialTotalSSave
    dPartialHeatCapacity = dPartialTotalCpSave
    dPartialExcessGibbs = dPartialExcessSave
    dPartialEnthalpyXS = dPartialHSave
    dPartialEntropyXS = dPartialSSave
    dPartialHeatCapacityXS = dPartialCpSave
    dMagGibbsEnergy = dMagGSave
    dMagEnthalpy = dMagHSave
    dMagEntropy = dMagSSave
    dMagHeatCapacity = dMagCpSave
    dGibbsSolnPhase = dGibbsPhaseSave
    if (allocated(dSiteFractionSave)) then
        dSiteFraction = dSiteFractionSave
        deallocate(dSiteFractionSave)
    end if
    deallocate(dMolFractionSave, dChemicalPotentialSave)
    deallocate(dPartialTotalHSave, dPartialTotalSSave, dPartialTotalCpSave)
    deallocate(dPartialExcessSave, dPartialHSave, dPartialSSave, dPartialCpSave)
    deallocate(dMagGSave, dMagHSave, dMagSSave, dMagCpSave)
    deallocate(dGibbsPhaseSave)

    return
end subroutine EvaluateStaticGridRows


!> \brief Evaluate one row through the unchanged shared scalar thermodynamic path.
subroutine EvaluateScalarGridPoint(iPoint)
    implicit none

    integer, intent(in) :: iPoint
    integer :: iPhaseIndex, iFirst, iLast, iOffset, iLevelRow, i
    real(8) :: dScalarFallbackG

    iPhaseIndex = iGridPointPhase(iPoint)
    iFirst = nSpeciesPhase(iPhaseIndex-1)+1
    iLast = nSpeciesPhase(iPhaseIndex)
    iOffset = iGridPointFractionOffset(iPoint)
    dMolFraction(iFirst:iLast) = dGridPointFraction(iOffset:iOffset+iLast-iFirst)
    dGibbsSolnPhase(iPhaseIndex) = 0D0
    call CompExcessGibbsEnergy(iPhaseIndex)
    if (INFOThermo /= 0) return
    dScalarFallbackG = 0D0
    do i = iFirst, iLast
        dScalarFallbackG = dScalarFallbackG+ &
            dGridPointFraction(iOffset+i-iFirst)*dChemicalPotential(i)
    end do
    iLevelRow = nSpecies+iPoint
    call RegisterGridLevelRow(iPoint,iLevelRow)
    call SetStaticGridLevelingRow(iPoint,iLevelRow,dScalarFallbackG)

    return
end subroutine EvaluateScalarGridPoint


!> \brief Preserve the immutable point-to-Leveling-row identity mapping.
subroutine RegisterGridLevelRow(iPoint,iLevelRow)
    implicit none

    integer, intent(in) :: iPoint,iLevelRow

    iGridPointLevelRow(iPoint) = iLevelRow
    iGridPointFromLevel(iLevelRow) = iPoint

    return
end subroutine RegisterGridLevelRow


!> \brief Convert one priced constitution to the Leveling composition basis.
subroutine SetStaticGridLevelingRow(iPoint, iLevelRow, dGibbsInput)
    implicit none

    integer, intent(in) :: iPoint, iLevelRow
    real(8), intent(in) :: dGibbsInput
    integer :: iPhaseIndex, iFirst, iLast, i, j, iOffset
    real(8) :: dLevelingDenom, dFormulaDenom, dAtomDenom, dGibbs
    real(8), dimension(nElements) :: dStoich
    real(8) :: dX

    iPhaseIndex = iGridPointPhase(iPoint)
    iFirst = nSpeciesPhase(iPhaseIndex-1) + 1
    iLast = nSpeciesPhase(iPhaseIndex)
    iOffset = iGridPointFractionOffset(iPoint)
    dLevelingDenom = 0D0
    dFormulaDenom = 0D0
    dAtomDenom = 0D0
    dGibbs = 0D0
    dStoich = 0D0
    dLevelingCompositionSpecies(iLevelRow,:) = 0D0

    do i = iFirst, iLast
        dX = dGridPointFraction(iOffset+i-iFirst)
        dAtomDenom = dAtomDenom + dX*dSpeciesTotalAtoms(i)/DBLE(iParticlesPerMole(i))
        dLevelingDenom = dLevelingDenom + &
            dX*dLevelingSpeciesTotalAtoms(i)/DBLE(iParticlesPerMole(i))
        dFormulaDenom = dFormulaDenom + &
            dX*dLevelingSpeciesFormulaAtoms(i)/DBLE(iParticlesPerMole(i))
        do j = 1, nElements
            dLevelingCompositionSpecies(iLevelRow,j) = &
                dLevelingCompositionSpecies(iLevelRow,j) + &
                dX*dStoichSpecies(i,j)/DBLE(iParticlesPerMole(i))
            dStoich(j) = dStoich(j) + dX*dStoichSpecies(i,j)/DBLE(iParticlesPerMole(i))
        end do
    end do

    dGibbs = dGibbsInput

    iPhaseLevel(iLevelRow) = iPhaseIndex
    dStoichSpeciesLevel(iLevelRow,:) = dStoich
    if ((dLevelingDenom <= 1D-300).OR.(dFormulaDenom <= 1D-300).OR.&
        (dAtomDenom <= 1D-300)) then
        dChemicalPotential(iLevelRow) = 5D9
        dLevelingChemicalPotential(iLevelRow) = 5D9
        return
    end if
    dLevelingCompositionSpecies(iLevelRow,:) = &
        dLevelingCompositionSpecies(iLevelRow,:)/dLevelingDenom
    dChemicalPotential(iLevelRow) = dGibbs/dAtomDenom
    dLevelingChemicalPotential(iLevelRow) = dGibbs/dFormulaDenom

    return
end subroutine SetStaticGridLevelingRow


!> \brief Copy a selected grid constitution into one Lagrangian candidate row.
subroutine CopyGridPointMolFraction(iLevelRow, iTargetRow, iSolnPhaseIndex, lFound)
    implicit none

    integer, intent(in) :: iLevelRow, iTargetRow, iSolnPhaseIndex
    logical, intent(out) :: lFound
    integer :: iPoint, iFirst, iLast, iOffset

    lFound = .FALSE.
    if (.NOT.allocated(iGridPointFromLevel)) return
    if ((iLevelRow < 1).OR.(iLevelRow > SIZE(iGridPointFromLevel))) return
    iPoint = iGridPointFromLevel(iLevelRow)
    if ((iPoint < 1).OR.(iPoint > nGridPoint)) return
    if (iGridPointPhase(iPoint) /= iSolnPhaseIndex) return

    iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLast = nSpeciesPhase(iSolnPhaseIndex)
    iOffset = iGridPointFractionOffset(iPoint)
    dMolFractionGEM(iTargetRow,iFirst:iLast) = &
        dGridPointFraction(iOffset:iOffset+iLast-iFirst)
    lFound = .TRUE.

    return
end subroutine CopyGridPointMolFraction


!> \brief Register selected grid rows as PEA starts and retire unselected rows.
!!
!! \details Copies every selected immutable grid constitution into the PEA
!! transport pool without subminimization.  These rows are starts for PEALG,
!! never fabricated evidence that a continuous phase minimum converged.  Rows
!! not selected by the initial grid-plus-stoichiometric Leveling solve are
!! excluded from later PEA Leveling passes; the authoritative DF sweep supplies
!! fresh continuous-minimum candidates for those passes.
subroutine RegisterSelectedGridSeeds
    implicit none

    integer :: iSlot, iLevelRow, iPoint, iPhaseIndex, iFirst, iLast
    integer :: iCandidate, nRequired

    if (.NOT.lGridFrontEndActive) return
    if (.NOT.lGridPoolValid) return
    if (.NOT.allocated(iGridPointFromLevel)) return
    if (.NOT.allocated(iLevelCandidateFromLevel)) return
    if (.NOT.allocated(iLevelCandidateSubMinStatus)) return

    nRequired = 0
    do iSlot = 1, nElements
        iLevelRow = iAssemblage(iSlot)
        if ((iLevelRow < 1).OR.(iLevelRow > SIZE(iGridPointFromLevel))) cycle
        if (iGridPointFromLevel(iLevelRow) > 0) nRequired = nRequired+1
    end do
    if (nLevelCandidate+nRequired > nLevelCandidateCapacity) then
        INFOThermo = 42
        return
    end if

    ! The static pool has completed its only Leveling role.  Only the selected
    ! support rows remain available as initial PEALG starts; allowing every
    ! sampled row into later PEA Leveling passes can replace the selected basis
    ! with unregistered fixed points and create an under-filled handoff.
    if (allocated(lLevelingRowExcluded)) then
        lLevelingRowExcluded(nSpecies+1:nSpecies+nGridPoint) = .TRUE.
    end if

    ! Register every selected static row by exact row identity and bytes.
    do iSlot = 1, nElements
        iLevelRow = iAssemblage(iSlot)
        if ((iLevelRow < 1).OR.(iLevelRow > SIZE(iLevelCandidateFromLevel))) cycle
        iPhaseIndex = iPhaseGEM(iSlot)
        iPoint = 0
        if ((iLevelRow >= 1).AND.(iLevelRow <= SIZE(iGridPointFromLevel))) then
            iPoint = iGridPointFromLevel(iLevelRow)
        end if
        if ((iPhaseIndex <= 0).AND.(iPoint >= 1).AND.(iPoint <= nGridPoint)) then
            iPhaseIndex = iGridPointPhase(iPoint)
        end if
        if ((iPhaseIndex <= 0).AND.(iLevelRow <= nSpecies)) iPhaseIndex = iPhase(iLevelRow)
        if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) cycle

        iFirst = nSpeciesPhase(iPhaseIndex-1) + 1
        iLast = nSpeciesPhase(iPhaseIndex)
        if ((iPoint < 1).OR.(iPoint > nGridPoint)) cycle
        if (allocated(lLevelingRowExcluded)) lLevelingRowExcluded(iLevelRow) = .FALSE.

        ! Step 3. Register the refined constitution as typed transport evidence.
        iCandidate = iLevelCandidateFromLevel(iLevelRow)
        if ((iCandidate < 1).OR.(iCandidate > nLevelCandidate)) then
            if (nLevelCandidate >= nLevelCandidateCapacity) then
                INFOThermo = 42
                return
            end if
            nLevelCandidate = nLevelCandidate + 1
            iCandidate = nLevelCandidate
        end if
        iLevelCandidateFromLevel(iLevelRow) = iCandidate
        iLevelCandidatePhase(iCandidate) = iPhaseIndex
        iLevelCandidateSource(iCandidate) = LEVEL_CANDIDATE_SOURCE_GRID
        iLevelCandidateSubMinStatus(iCandidate) = PEA_CANDIDATE_PROVISIONAL_GRID_SEED
        iLevelCandidateStaticRow(iCandidate) = iPoint
        iLevelCandidateParentPhase(iCandidate) = iPhaseIndex
        iLevelCandidateDisplayPhase(iCandidate) = iPhaseIndex
        if ((iPoint >= 1).AND.(iPoint <= nGridPoint)) then
            iLevelCandidateDisplayPhase(iCandidate) = iGridPointDisplayPhase(iPoint)
        end if
        ! A selected-row id is transport evidence only.  Level2 assigns a
        ! composition-set ordinal after the physical duplicate test.
        iLevelCandidateIdentityOrdinal(iCandidate) = iGridPointIdentityOrdinal(iPoint)
        dLevelCandidateMolFraction(iCandidate,:) = 0D0
        dLevelCandidateMolFraction(iCandidate,iFirst:iLast) = &
            dGridPointFraction(iGridPointFractionOffset(iPoint):&
            iGridPointFractionOffset(iPoint)+iLast-iFirst)
    end do

    return
end subroutine RegisterSelectedGridSeeds


!> \brief Subminimize every selected grid support row exactly once for init-PEA.
!!
!! \details Consumes provisional grid candidates as phase-local starts.  Only
!! converged minima or negative-driving-force witnesses are written into the
!! reserved nElements tail rows.  The immutable sampled row is then retired
!! from PEA Leveling, so status 6 remains transport metadata only.
subroutine RefineSelectedGridSeeds
    implicit none

    integer :: iSlot, iLevelRow, iCandidate, iPhaseIndex, iFirst, iLast
    integer :: iTailRow, iTailCandidate, iPoint, nRequired
    integer, dimension(nElements) :: iSelectedRow
    real(8), dimension(nSpecies) :: dChemicalPotentialSave
    logical :: lAddPhase, lAccepted, lRowAlreadyInstalled

    if (.NOT.lGridFrontEndActive) return
    if (.NOT.allocated(iLevelCandidateFromLevel)) return
    if (.NOT.allocated(iLevelCandidateSubMinStatus)) return

    iSelectedRow = iAssemblage
    nRequired = 0
    do iSlot = 1, nElements
        iLevelRow = iSelectedRow(iSlot)
        if ((iLevelRow < 1).OR.(iLevelRow > SIZE(iLevelCandidateFromLevel))) cycle
        iCandidate = iLevelCandidateFromLevel(iLevelRow)
        if ((iCandidate < 1).OR.(iCandidate > nLevelCandidate)) cycle
        if (iLevelCandidateSubMinStatus(iCandidate) == PEA_CANDIDATE_PROVISIONAL_GRID_SEED) &
            nRequired = nRequired + 1
    end do
    if (nLevelCandidate+nRequired > nLevelCandidateCapacity) then
        INFOThermo = 42
        return
    end if

    dChemicalPotentialSave = dChemicalPotential(:nSpecies)
    do iSlot = 1, nElements
        iLevelRow = iSelectedRow(iSlot)
        if ((iLevelRow < 1).OR.(iLevelRow > SIZE(iLevelCandidateFromLevel))) cycle
        iCandidate = iLevelCandidateFromLevel(iLevelRow)
        if ((iCandidate < 1).OR.(iCandidate > nLevelCandidate)) cycle
        if (iLevelCandidateSubMinStatus(iCandidate) /= PEA_CANDIDATE_PROVISIONAL_GRID_SEED) cycle

        iPhaseIndex = iLevelCandidateParentPhase(iCandidate)
        if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) cycle
        iFirst = nSpeciesPhase(iPhaseIndex-1)+1
        iLast = nSpeciesPhase(iPhaseIndex)
        dMolFraction(iFirst:iLast) = dLevelCandidateMolFraction(iCandidate,iFirst:iLast)
        lAddPhase = .FALSE.
        call Subminimization(iPhaseIndex,lAddPhase)
        lAccepted = (iSubMinCandidateStatusSoln(iPhaseIndex) == SUBMIN_CANDIDATE_CONVERGED).OR.&
            (iSubMinCandidateStatusSoln(iPhaseIndex) == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)

        if (allocated(lLevelingRowExcluded)) lLevelingRowExcluded(iLevelRow) = .TRUE.
        if (.NOT.lAccepted) then
            iLevelCandidateSubMinStatus(iCandidate) = SUBMIN_CANDIDATE_REJECTED
            cycle
        end if

        iTailRow = nSpecies+nGridPoint+2*nSolnPhasesSys+iSlot
        call SetLevelingSolutionCandidateRow(&
            iTailRow,iPhaseIndex,dMolFraction(iFirst:iLast),.TRUE.)
        if (INFOThermo /= 0) exit
        iTailCandidate = iLevelCandidateFromLevel(iTailRow)
        if ((iTailCandidate < 1).OR.(iTailCandidate > nLevelCandidate)) then
            INFOThermo = 42
            exit
        end if
        iPoint = iLevelCandidateStaticRow(iCandidate)
        iLevelCandidateSource(iTailCandidate) = LEVEL_CANDIDATE_SOURCE_GRID
        iLevelCandidateSubMinStatus(iTailCandidate) = iSubMinCandidateStatusSoln(iPhaseIndex)
        iLevelCandidateStaticRow(iTailCandidate) = iPoint
        iLevelCandidateParentPhase(iTailCandidate) = iPhaseIndex
        iLevelCandidateDisplayPhase(iTailCandidate) = iLevelCandidateDisplayPhase(iCandidate)
        iLevelCandidateIdentityOrdinal(iTailCandidate) = iLevelCandidateIdentityOrdinal(iCandidate)
        lRowAlreadyInstalled = .FALSE.
        do iPoint = 1, iSlot-1
            if (iAssemblage(iPoint) == iTailRow) then
                lRowAlreadyInstalled = .TRUE.
                exit
            end if
        end do
        if (lRowAlreadyInstalled) cycle
        iAssemblage(iSlot) = iTailRow
        iPhaseGEM(iSlot) = iPhaseIndex
        dMolFractionGEM(iSlot,:) = 0D0
        dMolFractionGEM(iSlot,iFirst:iLast) = dMolFraction(iFirst:iLast)
        dStoichSpeciesGEM(iSlot,:) = dStoichSpeciesLevel(iTailRow,:)
        dAtomFractionSpeciesGEM(iSlot,:) = dLevelingCompositionSpecies(iTailRow,:)
        dChemicalPotentialGEM(iSlot) = dLevelingChemicalPotential(iTailRow)
    end do
    dChemicalPotential(:nSpecies) = dChemicalPotentialSave

    return
end subroutine RefineSelectedGridSeeds


!> \brief Apply refined selected-row constitutions before Level2 transport.
subroutine ApplySelectedGridCandidateFractions
    implicit none

    integer :: iSlot, iLevelRow, iCandidate, iPhaseIndex, iFirst, iLast

    if (.NOT.lGridFrontEndActive) return
    if (.NOT.allocated(iLevelCandidateFromLevel)) return
    if (.NOT.allocated(iLevelCandidateSubMinStatus)) return

    do iSlot = 1, nElements
        iLevelRow = iAssemblage(iSlot)
        if ((iLevelRow < 1).OR.(iLevelRow > SIZE(iLevelCandidateFromLevel))) cycle
        iCandidate = iLevelCandidateFromLevel(iLevelRow)
        if ((iCandidate < 1).OR.(iCandidate > nLevelCandidate)) cycle
        if ((iLevelCandidateSubMinStatus(iCandidate) /= SUBMIN_CANDIDATE_CONVERGED).AND.&
            (iLevelCandidateSubMinStatus(iCandidate) /= SUBMIN_CANDIDATE_NEGATIVE_WITNESS).AND.&
            (iLevelCandidateSubMinStatus(iCandidate) /= PEA_CANDIDATE_PROVISIONAL_GRID_SEED)) cycle
        iPhaseIndex = iLevelCandidateParentPhase(iCandidate)
        if ((iPhaseIndex < 1).OR.(iPhaseIndex > nSolnPhasesSys)) cycle
        iFirst = nSpeciesPhase(iPhaseIndex-1) + 1
        iLast = nSpeciesPhase(iPhaseIndex)
        dMolFractionGEM(iSlot,iFirst:iLast) = &
            dLevelCandidateMolFraction(iCandidate,iFirst:iLast)
        iPhaseGEM(iSlot) = iPhaseIndex
    end do

    return
end subroutine ApplySelectedGridCandidateFractions


!> \brief Count ternary simplex faces carrying database interaction parameters.
integer function CountTernaryInteractionFaces(iTargetPhase, iDefinitionPhase, &
    iSubstitutionalSublattice, nConstituent)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase
    integer, intent(in) :: iSubstitutionalSublattice, nConstituent
    integer :: i, j, k

    CountTernaryInteractionFaces = 0
    do i = 1, nConstituent - 2
        do j = i + 1, nConstituent - 1
            do k = j + 1, nConstituent
                if (HasTernaryInteraction(iTargetPhase, iDefinitionPhase, &
                    iSubstitutionalSublattice, i, j, k)) &
                    CountTernaryInteractionFaces = CountTernaryInteractionFaces + 1
            end do
        end do
    end do

    return
end function CountTernaryInteractionFaces


!> \brief Test one ternary constituent face for a phase interaction parameter.
logical function HasTernaryInteraction(iTargetPhase, iDefinitionPhase, &
    iSubstitutionalSublattice, iA, iB, iC)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase
    integer, intent(in) :: iSubstitutionalSublattice, iA, iB, iC
    integer :: iParam, iEntry, iSublattice, iConstituent
    integer, dimension(3) :: iFace, iParamFace

    HasTernaryInteraction = .FALSE.
    iFace = (/iA, iB, iC/)
    if ((TRIM(cSolnPhaseType(iDefinitionPhase)) == 'SUBL').OR.&
        (TRIM(cSolnPhaseType(iDefinitionPhase)) == 'SUBLM').OR.&
        (TRIM(cSolnPhaseType(iDefinitionPhase)) == 'SUBOM')) then
        iFace(1) = GridSimplexConstituentID(iTargetPhase, iDefinitionPhase, &
            iPhaseSublattice(iDefinitionPhase), iSubstitutionalSublattice, iA)
        iFace(2) = GridSimplexConstituentID(iTargetPhase, iDefinitionPhase, &
            iPhaseSublattice(iDefinitionPhase), iSubstitutionalSublattice, iB)
        iFace(3) = GridSimplexConstituentID(iTargetPhase, iDefinitionPhase, &
            iPhaseSublattice(iDefinitionPhase), iSubstitutionalSublattice, iC)
        call SortInteger3(iFace)
        do iParam = nParamPhase(iDefinitionPhase-1)+1, nParamPhase(iDefinitionPhase)
            if (iRegularParam(iParam,1) /= 3) cycle
            if ((iSUBLParamData(iParam,1) /= 1).OR.(iSUBLParamData(iParam,3) /= 3)) cycle
            do iEntry = 1, 3
                iConstituent = MOD(iRegularParam(iParam,iEntry+1),10000)
                iSublattice = (iRegularParam(iParam,iEntry+1)-iConstituent)/10000
                if (iSublattice /= iSubstitutionalSublattice) exit
                iParamFace(iEntry) = iConstituent
            end do
            if (iEntry <= 3) cycle
            call SortInteger3(iParamFace)
            if (ALL(iParamFace == iFace)) then
                HasTernaryInteraction = .TRUE.
                return
            end if
        end do
    else
        do iParam = nParamPhase(iDefinitionPhase-1)+1, nParamPhase(iDefinitionPhase)
            if (iRegularParam(iParam,1) /= 3) cycle
            iParamFace = iRegularParam(iParam,2:4)
            call SortInteger3(iParamFace)
            if (ALL(iParamFace == iFace)) then
                HasTernaryInteraction = .TRUE.
                return
            end if
        end do
    end if

    return
end function HasTernaryInteraction


!> \brief Map one low-discrepancy point to a simplex by sorted spacings.
subroutine SobolSimplexPoint(iPointIndex, nConstituent, dFraction)
    implicit none

    integer, intent(in) :: iPointIndex, nConstituent
    real(8), dimension(nConstituent), intent(out) :: dFraction
    real(8), dimension(MAX(1,nConstituent-1)) :: dCoordinate
    integer :: i

    if (nConstituent <= 1) then
        dFraction = 1D0
        return
    end if
    call SobolPoint(iPointIndex, nConstituent-1, dCoordinate(1:nConstituent-1))
    call SortReal(dCoordinate(1:nConstituent-1))
    dFraction(1) = dCoordinate(1)
    do i = 2, nConstituent-1
        dFraction(i) = dCoordinate(i)-dCoordinate(i-1)
    end do
    dFraction(nConstituent) = 1D0-dCoordinate(nConstituent-1)

    return
end subroutine SobolSimplexPoint


!> \brief Generate one Sobol point using fixed direction numbers through dimension eight.
subroutine SobolPoint(iPointIndex, nDimension, dPoint)
    implicit none

    integer, intent(in) :: iPointIndex, nDimension
    real(8), dimension(nDimension), intent(out) :: dPoint
    integer(8), dimension(8,nSobolBit), save :: iDirection = 0_8
    logical, save :: lInitialized = .FALSE.
    integer(8) :: iGray, iValue
    integer :: iDimension, iBit
    real(8), parameter :: dTwo32 = 4294967296D0

    if (.NOT.lInitialized) then
        call InitSobolDirections(iDirection)
        lInitialized = .TRUE.
    end if
    if (nDimension > 8) then
        call KroneckerPoint(iPointIndex, nDimension, dPoint)
        return
    end if

    iGray = IEOR(INT(iPointIndex,8),ISHFT(INT(iPointIndex,8),-1))
    do iDimension = 1, nDimension
        iValue = 0_8
        do iBit = 1, nSobolBit
            if (BTEST(iGray,iBit-1)) iValue = IEOR(iValue,iDirection(iDimension,iBit))
        end do
        dPoint(iDimension) = DBLE(iValue)/dTwo32
    end do

    return
end subroutine SobolPoint


!> \brief Initialize the fixed Sobol direction-number table.
subroutine InitSobolDirections(iDirection)
    implicit none

    integer(8), dimension(8,nSobolBit), intent(out) :: iDirection
    integer, dimension(8) :: iDegree, iCoefficient
    integer, dimension(8,5) :: iInitial
    integer :: iDimension, iBit, k, iDegreeLocal
    integer(8) :: iValue

    iDegree = (/0, 1, 2, 3, 3, 4, 4, 5/)
    iCoefficient = (/0, 0, 1, 1, 2, 1, 4, 2/)
    iInitial = 0
    iInitial(2,1:1) = (/1/)
    iInitial(3,1:2) = (/1,3/)
    iInitial(4,1:3) = (/1,3,1/)
    iInitial(5,1:3) = (/1,1,1/)
    iInitial(6,1:4) = (/1,3,5,13/)
    iInitial(7,1:4) = (/1,1,5,5/)
    iInitial(8,1:5) = (/1,3,3,9,7/)
    iDirection = 0_8

    do iBit = 1, nSobolBit
        iDirection(1,iBit) = ISHFT(1_8,nSobolBit-iBit)
    end do
    do iDimension = 2, 8
        iDegreeLocal = iDegree(iDimension)
        do iBit = 1, iDegreeLocal
            iDirection(iDimension,iBit) = &
                ISHFT(INT(iInitial(iDimension,iBit),8),nSobolBit-iBit)
        end do
        do iBit = iDegreeLocal+1, nSobolBit
            iValue = IEOR(iDirection(iDimension,iBit-iDegreeLocal), &
                ISHFT(iDirection(iDimension,iBit-iDegreeLocal),-iDegreeLocal))
            do k = 1, iDegreeLocal-1
                if (BTEST(iCoefficient(iDimension),iDegreeLocal-1-k)) then
                    iValue = IEOR(iValue,iDirection(iDimension,iBit-k))
                end if
            end do
            iDirection(iDimension,iBit) = iValue
        end do
    end do

    return
end subroutine InitSobolDirections


!> \brief Generate the deterministic Kronecker fallback above Sobol dimension eight.
subroutine KroneckerPoint(iPointIndex, nDimension, dPoint)
    implicit none

    integer, intent(in) :: iPointIndex, nDimension
    real(8), dimension(nDimension), intent(out) :: dPoint
    integer :: i
    real(8), parameter :: dGolden = 1.6180339887498948482D0

    do i = 1, nDimension
        dPoint(i) = MOD(DBLE(iPointIndex)*dGolden**(-DBLE(i)),1D0)
    end do

    return
end subroutine KroneckerPoint


!> \brief Return a small integer binomial coefficient for simplex face counts.
integer function BinomialCount(n, k)
    implicit none

    integer, intent(in) :: n, k
    integer :: i, kLocal

    if ((k < 0).OR.(k > n)) then
        BinomialCount = 0
        return
    end if
    kLocal = MIN(k,n-k)
    BinomialCount = 1
    do i = 1, kLocal
        BinomialCount = BinomialCount*(n-kLocal+i)/i
    end do

    return
end function BinomialCount


!> \brief Return the public disordered display phase for a sampled parent.
integer function GridDisplayPhase(iPhaseIndex)
    implicit none

    integer, intent(in) :: iPhaseIndex
    integer :: iCompanionPhase

    GridDisplayPhase = iPhaseIndex
    if (.NOT.allocated(iODCompanionPhase)) return
    if ((iPhaseIndex < 1).OR.(iPhaseIndex > SIZE(iODCompanionPhase))) return
    iCompanionPhase = iODCompanionPhase(iPhaseIndex)
    if (iCompanionPhase <= 0) return
    if (allocated(iODStandalonePhase)) then
        if (iPhaseIndex <= SIZE(iODStandalonePhase)) then
            if (iODStandalonePhase(iPhaseIndex) == iCompanionPhase) return
        end if
    end if
    GridDisplayPhase = iCompanionPhase

    return
end function GridDisplayPhase


!> \brief Recognize vacancy constituent spelling used by supported databases.
logical function IsVacancyName(cName)
    implicit none

    character(*), intent(in) :: cName
    character(8) :: cUpper

    cUpper = UpperGridName(cName)
    IsVacancyName = (TRIM(cUpper) == 'VA').OR.(TRIM(cUpper) == 'VACANCY')

    return
end function IsVacancyName


!> \brief Count substitutional constituents after excluding vacancy.
integer function CountNonVacancyConstituents(iSublatticePhase, iSublattice)
    implicit none

    integer, intent(in) :: iSublatticePhase, iSublattice
    integer :: iConstituent

    CountNonVacancyConstituents = 0
    do iConstituent = 1, nConstituentSublattice(iSublatticePhase,iSublattice)
        if (.NOT.IsVacancyName(cConstituentNameSUB(&
            iSublatticePhase,iSublattice,iConstituent))) then
            CountNonVacancyConstituents = CountNonVacancyConstituents + 1
        end if
    end do

    return
end function CountNonVacancyConstituents


!> \brief Count definition constituents that the sampled target can represent.
!!
!! \details Iterates the definition set and uses the target only as a membership
!! predicate.  A target-side superset therefore retains the complete definition
!! manifold; target-only constituents are not coordinates of that manifold.
integer function GridSimplexConstituentCount(iTargetPhase, iDefinitionPhase, &
    iSublatticePhase, iSublattice)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase
    integer, intent(in) :: iSublatticePhase, iSublattice
    integer :: iConstituent, nNonVacancy

    GridSimplexConstituentCount = 0
    nNonVacancy = CountNonVacancyConstituents(iSublatticePhase,iSublattice)
    do iConstituent = 1, nConstituentSublattice(iSublatticePhase,iSublattice)
        if ((nNonVacancy > 0).AND.IsVacancyName(cConstituentNameSUB(&
            iSublatticePhase,iSublattice,iConstituent))) cycle
        if (.NOT.TargetRepresentsDefinitionConstituent(iTargetPhase, iDefinitionPhase, &
            cConstituentNameSUB(iSublatticePhase,iSublattice,iConstituent))) cycle
        GridSimplexConstituentCount = GridSimplexConstituentCount + 1
    end do

    return
end function GridSimplexConstituentCount


!> \brief Map one filtered simplex coordinate to its definition constituent.
integer function GridSimplexConstituentID(iTargetPhase, iDefinitionPhase, &
    iSublatticePhase, iSublattice, iSimplexIndex)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase
    integer, intent(in) :: iSublatticePhase, iSublattice, iSimplexIndex
    integer :: iConstituent, iCount, nNonVacancy

    GridSimplexConstituentID = 0
    nNonVacancy = CountNonVacancyConstituents(iSublatticePhase,iSublattice)
    iCount = 0
    do iConstituent = 1, nConstituentSublattice(iSublatticePhase,iSublattice)
        if ((nNonVacancy > 0).AND.IsVacancyName(cConstituentNameSUB(&
            iSublatticePhase,iSublattice,iConstituent))) cycle
        if (.NOT.TargetRepresentsDefinitionConstituent(iTargetPhase, iDefinitionPhase, &
            cConstituentNameSUB(iSublatticePhase,iSublattice,iConstituent))) cycle
        iCount = iCount + 1
        if (iCount == iSimplexIndex) then
            GridSimplexConstituentID = iConstituent
            return
        end if
    end do

    return
end function GridSimplexConstituentID


!> \brief Test whether a companion constituent exists in its ordered target.
logical function TargetRepresentsDefinitionConstituent(iTargetPhase, iDefinitionPhase, cName)
    implicit none

    integer, intent(in) :: iTargetPhase, iDefinitionPhase
    character(*), intent(in) :: cName
    integer :: iTargetSublatticePhase, iSublattice, iConstituent
    character(8) :: cDefinitionName

    TargetRepresentsDefinitionConstituent = .FALSE.
    if (iTargetPhase == iDefinitionPhase) then
        TargetRepresentsDefinitionConstituent = .TRUE.
        return
    end if
    if ((iTargetPhase < 1).OR.(iTargetPhase > nSolnPhasesSys)) return
    if ((TRIM(cSolnPhaseType(iTargetPhase)) /= 'SUBL').AND.&
        (TRIM(cSolnPhaseType(iTargetPhase)) /= 'SUBLM').AND.&
        (TRIM(cSolnPhaseType(iTargetPhase)) /= 'SUBOM')) return
    iTargetSublatticePhase = iPhaseSublattice(iTargetPhase)
    if ((iTargetSublatticePhase < 1).OR.&
        (iTargetSublatticePhase > nCountSublattice)) return
    cDefinitionName = UpperGridName(cName)
    do iSublattice = 1, nSublatticePhase(iTargetSublatticePhase)
        do iConstituent = 1, nConstituentSublattice(iTargetSublatticePhase,iSublattice)
            if (IsVacancyName(cConstituentNameSUB(&
                iTargetSublatticePhase,iSublattice,iConstituent))) cycle
            if (UpperGridName(cConstituentNameSUB(iTargetSublatticePhase, &
                iSublattice,iConstituent)) == cDefinitionName) then
                TargetRepresentsDefinitionConstituent = .TRUE.
                return
            end if
        end do
    end do

    return
end function TargetRepresentsDefinitionConstituent


!> \brief Normalize a short constituent name to uppercase for structural checks.
character(8) function UpperGridName(cName)
    implicit none

    character(*), intent(in) :: cName
    integer :: i, iCode

    UpperGridName = ' '
    do i = 1, MIN(LEN_TRIM(cName),LEN(UpperGridName))
        iCode = IACHAR(cName(i:i))
        if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
            UpperGridName(i:i) = ACHAR(iCode-32)
        else
            UpperGridName(i:i) = cName(i:i)
        end if
    end do

    return
end function UpperGridName


!> \brief Sort a ternary constituent index set into canonical order.
subroutine SortInteger3(iValue)
    implicit none

    integer, dimension(3), intent(inout) :: iValue
    integer :: i, j, iKey

    do i = 2, 3
        iKey = iValue(i)
        j = i-1
        do while (j >= 1)
            if (iValue(j) <= iKey) exit
            iValue(j+1) = iValue(j)
            j = j-1
        end do
        iValue(j+1) = iKey
    end do

    return
end subroutine SortInteger3


!> \brief Sort low-discrepancy coordinates before the simplex spacing transform.
subroutine SortReal(dValue)
    implicit none

    real(8), dimension(:), intent(inout) :: dValue
    integer :: i, j
    real(8) :: dKey

    do i = 2, SIZE(dValue)
        dKey = dValue(i)
        j = i-1
        do while (j >= 1)
            if (dValue(j) <= dKey) exit
            dValue(j+1) = dValue(j)
            j = j-1
        end do
        dValue(j+1) = dKey
    end do

    return
end subroutine SortReal

end module GridDiscovery
