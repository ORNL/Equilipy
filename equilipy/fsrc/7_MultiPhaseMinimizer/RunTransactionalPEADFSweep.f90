!> \brief Own one atomic PEA driving-force candidate generation.
!!
!! \details Snapshots the accepted primal state once, evaluates all solution
!! candidates into staging storage, and either commits a canonically ordered
!! generation or restores the prior candidate state without partial evidence.
!-------------------------------------------------------------------------------------------------------------
!
!> \file    RunTransactionalPEADFSweep.f90
!> \brief   Atomic PEA driving-force sweep owner.
!> \author  S.Y. Kwon
!> \date    Jul. 19, 2026
!> \sa      CompMinSolnPoint.f90
!> \sa      CheckPhaseAssemblage.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Committed complete driving-force generations and restored the exact state after failed PEA sweeps.
    !
! Numerical assumptions:
! ======================
!
! - Candidate constitutions are packed phase-locally.
! - A failed or unknown record invalidates the complete sweep.
! - The static grid pool is read-only and is not part of transaction storage.
!-------------------------------------------------------------------------------------------------------------
subroutine InitPEADFSweepStorage
    USE ModuleThermo
    USE ModuleGEMSolver
    implicit none

    integer :: i,nLocal,nStart,nRecord,nFraction

    call ResetPEADFSweepStorage
    nRecord = 0
    nFraction = 0
    do i = 1,nSolnPhasesSys
        nLocal = nSpeciesPhase(i)-nSpeciesPhase(i-1)
        nStart = 1
        if (TRIM(cSolnPhaseType(i)) == 'SUBOM') nStart = nLocal
        nRecord = nRecord+nStart
        nFraction = nFraction+nStart*nLocal
    end do
    nRecord = MAX(1,nRecord)
    nFraction = MAX(1,nFraction)
    nPEADFSweepRecordCapacity = nRecord
    nPEADFSweepFractionCapacity = nFraction

    allocate(iPEADFSweepStagingParent(nRecord),iPEADFSweepStagingDisplay(nRecord))
    allocate(iPEADFSweepStagingRepresentation(nRecord),iPEADFSweepStagingOrdinal(nRecord))
    allocate(iPEADFSweepStagingStartRank(nRecord),iPEADFSweepStagingStatus(nRecord))
    allocate(iPEADFSweepStagingProof(nRecord),iPEADFSweepStagingBasis(nRecord))
    allocate(iPEADFSweepStagingFractionOffset(nRecord),iPEADFSweepStagingFractionSize(nRecord))
    allocate(iPEADFSweepCommittedParent(nRecord),iPEADFSweepCommittedDisplay(nRecord))
    allocate(iPEADFSweepCommittedRepresentation(nRecord),iPEADFSweepCommittedOrdinal(nRecord))
    allocate(iPEADFSweepCommittedStartRank(nRecord),iPEADFSweepCommittedStatus(nRecord))
    allocate(iPEADFSweepCommittedProof(nRecord),iPEADFSweepCommittedBasis(nRecord))
    allocate(iPEADFSweepCommittedFractionOffset(nRecord),iPEADFSweepCommittedFractionSize(nRecord))
    allocate(iPEADFSweepCommittedWitnessPending(nRecord))
    allocate(dPEADFSweepStagingDF(nRecord),dPEADFSweepStagingTiming(nRecord))
    allocate(dPEADFSweepStagingPlane(nRecord,MAX(1,nElements)))
    allocate(dPEADFSweepCommittedDF(nRecord),dPEADFSweepCommittedTiming(nRecord))
    allocate(dPEADFSweepCommittedPlane(nRecord,MAX(1,nElements)))
    allocate(dPEADFSweepStagingFraction(nFraction),dPEADFSweepCommittedFraction(nFraction))
    call ClearPEADFSweepBanks
    iPEAPlaneGeneration = 0
    iPEADFSweepGeneration = 0
    iPEADFSweepCommittedPlaneGeneration = 0
    iPEADFSweepOutcome = PEA_DF_SWEEP_NOT_RUN
    nPEADFSweepAttempted = 0
    nPEADFSweepCommitted = 0
    nPEADFSweepUnknown = 0
    nPEADFSweepFailed = 0
    nPEADFPendingWitness = 0
end subroutine InitPEADFSweepStorage


subroutine ClearPEADFSweepBanks
    USE ModuleGEMSolver
    implicit none

    iPEADFSweepStagingParent = 0
    iPEADFSweepStagingDisplay = 0
    iPEADFSweepStagingRepresentation = 0
    iPEADFSweepStagingOrdinal = 0
    iPEADFSweepStagingStartRank = 0
    iPEADFSweepStagingStatus = GEM_CERT_STATUS_UNKNOWN
    iPEADFSweepStagingProof = GEM_CERT_PROOF_UNKNOWN
    iPEADFSweepStagingBasis = 0
    iPEADFSweepStagingFractionOffset = 0
    iPEADFSweepStagingFractionSize = 0
    iPEADFSweepCommittedParent = 0
    iPEADFSweepCommittedDisplay = 0
    iPEADFSweepCommittedRepresentation = 0
    iPEADFSweepCommittedOrdinal = 0
    iPEADFSweepCommittedStartRank = 0
    iPEADFSweepCommittedStatus = GEM_CERT_STATUS_UNKNOWN
    iPEADFSweepCommittedProof = GEM_CERT_PROOF_UNKNOWN
    iPEADFSweepCommittedBasis = 0
    iPEADFSweepCommittedFractionOffset = 0
    iPEADFSweepCommittedFractionSize = 0
    iPEADFSweepCommittedWitnessPending = 0
    dPEADFSweepStagingDF = 0D0
    dPEADFSweepStagingTiming = 0D0
    dPEADFSweepStagingPlane = 0D0
    dPEADFSweepCommittedDF = 0D0
    dPEADFSweepCommittedTiming = 0D0
    dPEADFSweepCommittedPlane = 0D0
    dPEADFSweepStagingFraction = 0D0
    dPEADFSweepCommittedFraction = 0D0
end subroutine ClearPEADFSweepBanks


subroutine ResetPEADFSweepStorage
    USE ModuleGEMSolver
    implicit none

    if (allocated(iPEADFSweepStagingParent)) deallocate(iPEADFSweepStagingParent)
    if (allocated(iPEADFSweepStagingDisplay)) deallocate(iPEADFSweepStagingDisplay)
    if (allocated(iPEADFSweepStagingRepresentation)) deallocate(iPEADFSweepStagingRepresentation)
    if (allocated(iPEADFSweepStagingOrdinal)) deallocate(iPEADFSweepStagingOrdinal)
    if (allocated(iPEADFSweepStagingStartRank)) deallocate(iPEADFSweepStagingStartRank)
    if (allocated(iPEADFSweepStagingStatus)) deallocate(iPEADFSweepStagingStatus)
    if (allocated(iPEADFSweepStagingProof)) deallocate(iPEADFSweepStagingProof)
    if (allocated(iPEADFSweepStagingBasis)) deallocate(iPEADFSweepStagingBasis)
    if (allocated(iPEADFSweepStagingFractionOffset)) deallocate(iPEADFSweepStagingFractionOffset)
    if (allocated(iPEADFSweepStagingFractionSize)) deallocate(iPEADFSweepStagingFractionSize)
    if (allocated(iPEADFSweepCommittedParent)) deallocate(iPEADFSweepCommittedParent)
    if (allocated(iPEADFSweepCommittedDisplay)) deallocate(iPEADFSweepCommittedDisplay)
    if (allocated(iPEADFSweepCommittedRepresentation)) deallocate(iPEADFSweepCommittedRepresentation)
    if (allocated(iPEADFSweepCommittedOrdinal)) deallocate(iPEADFSweepCommittedOrdinal)
    if (allocated(iPEADFSweepCommittedStartRank)) deallocate(iPEADFSweepCommittedStartRank)
    if (allocated(iPEADFSweepCommittedStatus)) deallocate(iPEADFSweepCommittedStatus)
    if (allocated(iPEADFSweepCommittedProof)) deallocate(iPEADFSweepCommittedProof)
    if (allocated(iPEADFSweepCommittedBasis)) deallocate(iPEADFSweepCommittedBasis)
    if (allocated(iPEADFSweepCommittedFractionOffset)) deallocate(iPEADFSweepCommittedFractionOffset)
    if (allocated(iPEADFSweepCommittedFractionSize)) deallocate(iPEADFSweepCommittedFractionSize)
    if (allocated(iPEADFSweepCommittedWitnessPending)) deallocate(iPEADFSweepCommittedWitnessPending)
    if (allocated(dPEADFSweepStagingDF)) deallocate(dPEADFSweepStagingDF)
    if (allocated(dPEADFSweepStagingTiming)) deallocate(dPEADFSweepStagingTiming)
    if (allocated(dPEADFSweepStagingPlane)) deallocate(dPEADFSweepStagingPlane)
    if (allocated(dPEADFSweepStagingFraction)) deallocate(dPEADFSweepStagingFraction)
    if (allocated(dPEADFSweepCommittedDF)) deallocate(dPEADFSweepCommittedDF)
    if (allocated(dPEADFSweepCommittedTiming)) deallocate(dPEADFSweepCommittedTiming)
    if (allocated(dPEADFSweepCommittedPlane)) deallocate(dPEADFSweepCommittedPlane)
    if (allocated(dPEADFSweepCommittedFraction)) deallocate(dPEADFSweepCommittedFraction)
    nPEADFSweepRecordCapacity = 0
    nPEADFSweepFractionCapacity = 0
    nPEADFSweepStagingRecord = 0
    nPEADFSweepStagingFraction = 0
    nPEADFSweepCommittedRecord = 0
    nPEADFSweepCommittedFraction = 0
    iPEAPlaneGeneration = 0
    iPEADFSweepGeneration = 0
    iPEADFSweepCommittedPlaneGeneration = 0
    nPEADFSweepAttempted = 0
    nPEADFSweepCommitted = 0
    nPEADFSweepUnknown = 0
    nPEADFSweepFailed = 0
    nPEADFPendingWitness = 0
    iPEADFSweepOutcome = PEA_DF_SWEEP_NOT_RUN
    lPEADFSweepStagingActive = .FALSE.
end subroutine ResetPEADFSweepStorage


subroutine StagePEADFCandidate(iParent,iLevelRow,iBasis,iStatus,iProof,iRank,&
    dDrivingForce,dTiming,iFirstSpecies,iLastSpecies,dFraction)
    USE ModuleThermo
    USE ModuleGEMSolver
    implicit none

    integer, intent(in) :: iParent,iLevelRow,iBasis,iStatus,iProof,iRank
    integer, intent(in) :: iFirstSpecies,iLastSpecies
    real(8), intent(in) :: dDrivingForce,dTiming
    real(8), intent(in) :: dFraction(iLastSpecies-iFirstSpecies+1)
    integer :: iRecord,iOffset,nLocal,iCandidateLocal,iEffectiveCapacity

    if (.NOT.lPEADFSweepStagingActive) return
    if (iPEADFSweepOutcome == PEA_DF_SWEEP_FAILED) return
    iRecord = nPEADFSweepStagingRecord + 1
    if ((iPEADFSweepInjectFailureAfter > 0).AND.&
        (iRecord > iPEADFSweepInjectFailureAfter)) then
        iPEADFSweepOutcome = PEA_DF_SWEEP_FAILED
        return
    end if
    iEffectiveCapacity = nPEADFSweepRecordCapacity
    if (iPEADFSweepCapacityOverride > 0) iEffectiveCapacity = iPEADFSweepCapacityOverride
    nLocal = iLastSpecies-iFirstSpecies+1
    iOffset = nPEADFSweepStagingFraction + 1
    if ((iRecord > iEffectiveCapacity).OR.&
        (iOffset+nLocal-1 > nPEADFSweepFractionCapacity)) then
        iPEADFSweepOutcome = PEA_DF_SWEEP_FAILED
        return
    end if

    nPEADFSweepStagingRecord = iRecord
    nPEADFSweepStagingFraction = nPEADFSweepStagingFraction + nLocal
    iPEADFSweepStagingParent(iRecord) = iParent
    iPEADFSweepStagingDisplay(iRecord) = iParent
    iPEADFSweepStagingRepresentation(iRecord) = 0
    iPEADFSweepStagingOrdinal(iRecord) = 0
    if (allocated(iODCandidateClass)) then
        if ((iParent >= 1).AND.(iParent <= SIZE(iODCandidateClass))) &
            iPEADFSweepStagingRepresentation(iRecord) = iODCandidateClass(iParent)
    end if
    if (allocated(iLevelCandidateFromLevel)) then
        if ((iLevelRow >= 1).AND.(iLevelRow <= SIZE(iLevelCandidateFromLevel))) then
            iCandidateLocal = iLevelCandidateFromLevel(iLevelRow)
            if ((iCandidateLocal >= 1).AND.(iCandidateLocal <= nLevelCandidate)) then
                iPEADFSweepStagingDisplay(iRecord) = iLevelCandidateDisplayPhase(iCandidateLocal)
                iPEADFSweepStagingOrdinal(iRecord) = iLevelCandidateIdentityOrdinal(iCandidateLocal)
            end if
        end if
    end if
    iPEADFSweepStagingStartRank(iRecord) = iRank
    iPEADFSweepStagingStatus(iRecord) = iStatus
    if ((iPEADFSweepInjectUnknownAfter > 0).AND.&
        (iRecord == iPEADFSweepInjectUnknownAfter)) then
        iPEADFSweepStagingStatus(iRecord) = GEM_CERT_STATUS_UNKNOWN
    end if
    iPEADFSweepStagingProof(iRecord) = iProof
    iPEADFSweepStagingBasis(iRecord) = iBasis
    iPEADFSweepStagingFractionOffset(iRecord) = iOffset
    iPEADFSweepStagingFractionSize(iRecord) = nLocal
    dPEADFSweepStagingDF(iRecord) = dDrivingForce
    dPEADFSweepStagingTiming(iRecord) = dTiming
    dPEADFSweepStagingPlane(iRecord,:nElements) = dElementPotential
    dPEADFSweepStagingFraction(iOffset:iOffset+nLocal-1) = dFraction
end subroutine StagePEADFCandidate


subroutine RunTransactionalPEADFSweep
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
    implicit none

    integer :: nCertCountSave,nCertEmissionSave,nCertDroppedSave,nLevelCandidateSave,iOutcome
    real(8) :: dFunctionNormSave,dFunctionNormLastSave
    real(8), allocatable :: dMolFractionSave(:),dGibbsSolnPhaseSave(:),dChemicalPotentialSave(:)
    real(8), allocatable :: dLevelingChemicalPotentialSave(:),dPhasePotentialSave(:)
    real(8), allocatable :: dMolesSpeciesSave(:),dPartialExcessSave(:),dPartialHSave(:)
    real(8), allocatable :: dPartialSSave(:),dPartialCpSave(:),dPartialTotalHSave(:)
    real(8), allocatable :: dPartialTotalSSave(:),dPartialTotalCpSave(:)
    real(8), allocatable :: dMagGSave(:),dMagHSave(:),dMagSSave(:),dMagCpSave(:)
    real(8), allocatable :: dSiteFractionSave(:,:,:),dEffStoichSave(:,:),dSumFractionSave(:)
    real(8), allocatable :: dAtomFractionSave(:,:),dLevelCompositionSave(:,:),dStoichLevelSave(:,:)
    real(8), allocatable :: dLevelCandidateSave(:,:),dDrivingForceSave(:)
    integer, allocatable :: iPhaseLevelSave(:),iLevelCandidatePhaseSave(:),iLevelCandidateSourceSave(:)
    integer, allocatable :: iLevelCandidateStatusSave(:),iLevelCandidateParentSave(:)
    integer, allocatable :: iLevelCandidateDisplaySave(:),iLevelCandidateOrdinalSave(:)
    integer, allocatable :: iLevelCandidateStaticSave(:),iLevelCandidateFromLevelSave(:)
    integer, allocatable :: iSubMinStatusSave(:),iODClassSave(:),iODCompanionSave(:)
    real(8), allocatable :: dODCurrentSave(:),dODDisorderedSave(:),dODEigenSave(:)
    logical, allocatable :: lLevelingExcludedSave(:)
    integer, allocatable :: iCertPhaseSave(:),iCertCopySave(:),iCertLevelRowSave(:)
    integer, allocatable :: iCertBasisSave(:),iCertNormSave(:),iCertStatusSave(:),iCertProofSave(:)
    integer, allocatable :: iCertSubMinSave(:),iCertRankSave(:),iCertFirstSave(:),iCertLastSave(:)
    real(8), allocatable :: dCertDFSave(:),dCertTimingSave(:),dCertPlaneSave(:,:)

    nPEADFSweepAttempted = nPEADFSweepAttempted + 1
    iPEAPlaneGeneration = iPEAPlaneGeneration + 1
    iPEADFSweepOutcome = PEA_DF_SWEEP_NOT_RUN
    nPEADFSweepStagingRecord = 0
    nPEADFSweepStagingFraction = 0
    iPEADFSweepStagingParent = 0
    iPEADFSweepStagingStatus = GEM_CERT_STATUS_UNKNOWN
    iPEADFSweepStagingProof = GEM_CERT_PROOF_UNKNOWN
    dPEADFSweepStagingFraction = 0D0

    dMolFractionSave = dMolFraction
    dGibbsSolnPhaseSave = dGibbsSolnPhase
    dChemicalPotentialSave = dChemicalPotential
    dLevelingChemicalPotentialSave = dLevelingChemicalPotential
    dPhasePotentialSave = dPhasePotential
    dMolesSpeciesSave = dMolesSpecies
    dPartialExcessSave = dPartialExcessGibbs
    dPartialHSave = dPartialEnthalpyXS
    dPartialSSave = dPartialEntropyXS
    dPartialCpSave = dPartialHeatCapacityXS
    dPartialTotalHSave = dPartialEnthalpy
    dPartialTotalSSave = dPartialEntropy
    dPartialTotalCpSave = dPartialHeatCapacity
    dMagGSave = dMagGibbsEnergy
    dMagHSave = dMagEnthalpy
    dMagSSave = dMagEntropy
    dMagCpSave = dMagHeatCapacity
    if (allocated(dSiteFraction)) dSiteFractionSave = dSiteFraction
    dEffStoichSave = dEffStoichSolnPhase
    dSumFractionSave = dSumMolFractionSoln
    iSubMinStatusSave = iSubMinCandidateStatusSoln
    if (allocated(iODCandidateClass)) iODClassSave = iODCandidateClass
    if (allocated(iODCandidateCompanionPhase)) iODCompanionSave = iODCandidateCompanionPhase
    if (allocated(dODCandidateCurrentGibbs)) dODCurrentSave = dODCandidateCurrentGibbs
    if (allocated(dODCandidateDisorderedGibbs)) dODDisorderedSave = dODCandidateDisorderedGibbs
    if (allocated(dODCandidateOrderingEigenMin)) dODEigenSave = dODCandidateOrderingEigenMin
    dAtomFractionSave = dAtomFractionSpecies
    dLevelCompositionSave = dLevelingCompositionSpecies
    dStoichLevelSave = dStoichSpeciesLevel
    iPhaseLevelSave = iPhaseLevel
    dFunctionNormSave = dGEMFunctionNorm
    dFunctionNormLastSave = dGEMFunctionNormLast
    dDrivingForceSave = dDrivingForceSoln
    dLevelCandidateSave = dLevelCandidateMolFraction
    iLevelCandidatePhaseSave = iLevelCandidatePhase
    iLevelCandidateSourceSave = iLevelCandidateSource
    iLevelCandidateStatusSave = iLevelCandidateSubMinStatus
    iLevelCandidateParentSave = iLevelCandidateParentPhase
    iLevelCandidateDisplaySave = iLevelCandidateDisplayPhase
    iLevelCandidateOrdinalSave = iLevelCandidateIdentityOrdinal
    iLevelCandidateStaticSave = iLevelCandidateStaticRow
    iLevelCandidateFromLevelSave = iLevelCandidateFromLevel
    lLevelingExcludedSave = lLevelingRowExcluded
    iCertPhaseSave = iGEMCertPhase
    iCertCopySave = iGEMCertCopy
    iCertLevelRowSave = iGEMCertLevelRow
    iCertBasisSave = iGEMCertBasis
    iCertNormSave = iGEMCertNorm
    iCertStatusSave = iGEMCertStatus
    iCertProofSave = iGEMCertProof
    iCertSubMinSave = iGEMCertSubMinStatus
    iCertRankSave = iGEMCertRank
    iCertFirstSave = iGEMCertFirstSpecies
    iCertLastSave = iGEMCertLastSpecies
    dCertDFSave = dGEMCertDrivingForce
    dCertTimingSave = dGEMCertTiming
    dCertPlaneSave = dGEMCertPlane
    nCertCountSave = nGEMCertCount
    nCertEmissionSave = nGEMCertEmissionCount
    nCertDroppedSave = nGEMCertDropped
    nLevelCandidateSave = nLevelCandidate

    iGEMCertPhase = 0
    iGEMCertCopy = 0
    iGEMCertLevelRow = 0
    iGEMCertBasis = 0
    iGEMCertNorm = GEM_CERT_NORMALIZATION_PER_MOLE_ATOMS
    iGEMCertStatus = GEM_CERT_STATUS_UNKNOWN
    iGEMCertProof = GEM_CERT_PROOF_UNKNOWN
    iGEMCertSubMinStatus = SUBMIN_CANDIDATE_UNKNOWN
    iGEMCertRank = 0
    iGEMCertFirstSpecies = 0
    iGEMCertLastSpecies = 0
    dGEMCertDrivingForce = 9D5
    dGEMCertTiming = 0D0
    dGEMCertPlane = 0D0
    nGEMCertCount = 0
    nGEMCertEmissionCount = 0
    nGEMCertDropped = 0

    lPEADFSweepStagingActive = .TRUE.
    call CompMinSolnPoint
    lPEADFSweepStagingActive = .FALSE.

    iOutcome = iPEADFSweepOutcome
    if (iOutcome /= PEA_DF_SWEEP_FAILED) then
        if (INFOThermo /= 0) then
            iOutcome = PEA_DF_SWEEP_FAILED
        else if (nPEADFSweepStagingRecord <= 0) then
            iOutcome = PEA_DF_SWEEP_UNKNOWN
        else if (ANY(iPEADFSweepStagingStatus(:nPEADFSweepStagingRecord) == GEM_CERT_STATUS_UNKNOWN)) then
            iOutcome = PEA_DF_SWEEP_UNKNOWN
        else
            iOutcome = PEA_DF_SWEEP_COMMITTED
        end if
    end if

    if (iOutcome == PEA_DF_SWEEP_COMMITTED) then
        call CommitPEADFSweepGeneration
        iPEADFSweepGeneration = iPEADFSweepGeneration + 1
        iPEADFSweepCommittedPlaneGeneration = iPEAPlaneGeneration
        nPEADFSweepCommitted = nPEADFSweepCommitted + 1
    else
        dMolFraction = dMolFractionSave
        dGibbsSolnPhase = dGibbsSolnPhaseSave
        dChemicalPotential = dChemicalPotentialSave
        dLevelingChemicalPotential = dLevelingChemicalPotentialSave
        dPhasePotential = dPhasePotentialSave
        dMolesSpecies = dMolesSpeciesSave
        dPartialExcessGibbs = dPartialExcessSave
        dPartialEnthalpyXS = dPartialHSave
        dPartialEntropyXS = dPartialSSave
        dPartialHeatCapacityXS = dPartialCpSave
        dPartialEnthalpy = dPartialTotalHSave
        dPartialEntropy = dPartialTotalSSave
        dPartialHeatCapacity = dPartialTotalCpSave
        dMagGibbsEnergy = dMagGSave
        dMagEnthalpy = dMagHSave
        dMagEntropy = dMagSSave
        dMagHeatCapacity = dMagCpSave
        if (allocated(dSiteFractionSave)) dSiteFraction = dSiteFractionSave
        dEffStoichSolnPhase = dEffStoichSave
        dSumMolFractionSoln = dSumFractionSave
        iSubMinCandidateStatusSoln = iSubMinStatusSave
        if (allocated(iODClassSave)) iODCandidateClass = iODClassSave
        if (allocated(iODCompanionSave)) iODCandidateCompanionPhase = iODCompanionSave
        if (allocated(dODCurrentSave)) dODCandidateCurrentGibbs = dODCurrentSave
        if (allocated(dODDisorderedSave)) dODCandidateDisorderedGibbs = dODDisorderedSave
        if (allocated(dODEigenSave)) dODCandidateOrderingEigenMin = dODEigenSave
        dGEMFunctionNorm = dFunctionNormSave
        dGEMFunctionNormLast = dFunctionNormLastSave
        dDrivingForceSoln = dDrivingForceSave
        dAtomFractionSpecies = dAtomFractionSave
        dLevelingCompositionSpecies = dLevelCompositionSave
        dStoichSpeciesLevel = dStoichLevelSave
        iPhaseLevel = iPhaseLevelSave
        dLevelCandidateMolFraction = dLevelCandidateSave
        iLevelCandidatePhase = iLevelCandidatePhaseSave
        iLevelCandidateSource = iLevelCandidateSourceSave
        iLevelCandidateSubMinStatus = iLevelCandidateStatusSave
        iLevelCandidateParentPhase = iLevelCandidateParentSave
        iLevelCandidateDisplayPhase = iLevelCandidateDisplaySave
        iLevelCandidateIdentityOrdinal = iLevelCandidateOrdinalSave
        iLevelCandidateStaticRow = iLevelCandidateStaticSave
        iLevelCandidateFromLevel = iLevelCandidateFromLevelSave
        lLevelingRowExcluded = lLevelingExcludedSave
        iGEMCertPhase = iCertPhaseSave
        iGEMCertCopy = iCertCopySave
        iGEMCertLevelRow = iCertLevelRowSave
        iGEMCertBasis = iCertBasisSave
        iGEMCertNorm = iCertNormSave
        iGEMCertStatus = iCertStatusSave
        iGEMCertProof = iCertProofSave
        iGEMCertSubMinStatus = iCertSubMinSave
        iGEMCertRank = iCertRankSave
        iGEMCertFirstSpecies = iCertFirstSave
        iGEMCertLastSpecies = iCertLastSave
        dGEMCertDrivingForce = dCertDFSave
        dGEMCertTiming = dCertTimingSave
        dGEMCertPlane = dCertPlaneSave
        nGEMCertCount = nCertCountSave
        nGEMCertEmissionCount = nCertEmissionSave
        nGEMCertDropped = nCertDroppedSave
        nLevelCandidate = nLevelCandidateSave
        if (iOutcome == PEA_DF_SWEEP_UNKNOWN) then
            nPEADFSweepUnknown = nPEADFSweepUnknown + 1
        else
            nPEADFSweepFailed = nPEADFSweepFailed + 1
        end if
    end if
    iPEADFSweepOutcome = iOutcome
end subroutine RunTransactionalPEADFSweep


subroutine CommitPEADFSweepGeneration
    USE ModuleGEMSolver
    implicit none

    integer :: i,j,key,nLocal,iSourceOffset,iTargetOffset,nRecord
    integer, allocatable :: iOrder(:)

    nRecord = nPEADFSweepStagingRecord
    allocate(iOrder(nRecord))
    iOrder = [(i,i=1,nRecord)]
    do i = 2,nRecord
        key = iOrder(i)
        j = i-1
        do while (j >= 1)
            if (.NOT.PEADFCandidateKeyAfter(iOrder(j),key)) exit
            iOrder(j+1) = iOrder(j)
            j = j-1
        end do
        iOrder(j+1) = key
    end do

    iPEADFSweepCommittedParent = 0
    iPEADFSweepCommittedDisplay = 0
    iPEADFSweepCommittedRepresentation = 0
    iPEADFSweepCommittedOrdinal = 0
    iPEADFSweepCommittedStartRank = 0
    iPEADFSweepCommittedStatus = GEM_CERT_STATUS_UNKNOWN
    iPEADFSweepCommittedProof = GEM_CERT_PROOF_UNKNOWN
    iPEADFSweepCommittedBasis = 0
    iPEADFSweepCommittedFractionOffset = 0
    iPEADFSweepCommittedFractionSize = 0
    iPEADFSweepCommittedWitnessPending = 0
    dPEADFSweepCommittedDF = 0D0
    dPEADFSweepCommittedTiming = 0D0
    dPEADFSweepCommittedPlane = 0D0
    dPEADFSweepCommittedFraction = 0D0
    iTargetOffset = 1
    nPEADFPendingWitness = 0
    do i = 1,nRecord
        j = iOrder(i)
        iPEADFSweepCommittedParent(i) = iPEADFSweepStagingParent(j)
        iPEADFSweepCommittedDisplay(i) = iPEADFSweepStagingDisplay(j)
        iPEADFSweepCommittedRepresentation(i) = iPEADFSweepStagingRepresentation(j)
        iPEADFSweepCommittedOrdinal(i) = iPEADFSweepStagingOrdinal(j)
        iPEADFSweepCommittedStartRank(i) = iPEADFSweepStagingStartRank(j)
        iPEADFSweepCommittedStatus(i) = iPEADFSweepStagingStatus(j)
        iPEADFSweepCommittedProof(i) = iPEADFSweepStagingProof(j)
        iPEADFSweepCommittedBasis(i) = iPEADFSweepStagingBasis(j)
        iPEADFSweepCommittedFractionOffset(i) = iTargetOffset
        nLocal = iPEADFSweepStagingFractionSize(j)
        iPEADFSweepCommittedFractionSize(i) = nLocal
        dPEADFSweepCommittedDF(i) = dPEADFSweepStagingDF(j)
        dPEADFSweepCommittedTiming(i) = dPEADFSweepStagingTiming(j)
        dPEADFSweepCommittedPlane(i,:) = dPEADFSweepStagingPlane(j,:)
        iSourceOffset = iPEADFSweepStagingFractionOffset(j)
        dPEADFSweepCommittedFraction(iTargetOffset:iTargetOffset+nLocal-1) = &
            dPEADFSweepStagingFraction(iSourceOffset:iSourceOffset+nLocal-1)
        if (iPEADFSweepCommittedStatus(i) == GEM_CERT_STATUS_FAVORABLE) then
            iPEADFSweepCommittedWitnessPending(i) = 1
            nPEADFPendingWitness = nPEADFPendingWitness + 1
        end if
        iTargetOffset = iTargetOffset+nLocal
    end do
    nPEADFSweepCommittedRecord = nRecord
    nPEADFSweepCommittedFraction = iTargetOffset-1
    deallocate(iOrder)
contains
    logical function PEADFCandidateKeyAfter(iLeft,iRight)
        integer, intent(in) :: iLeft,iRight
        integer :: a(5),b(5),k
        a = [iPEADFSweepStagingParent(iLeft),iPEADFSweepStagingRepresentation(iLeft),&
            iPEADFSweepStagingDisplay(iLeft),iPEADFSweepStagingOrdinal(iLeft),&
            iPEADFSweepStagingStartRank(iLeft)]
        b = [iPEADFSweepStagingParent(iRight),iPEADFSweepStagingRepresentation(iRight),&
            iPEADFSweepStagingDisplay(iRight),iPEADFSweepStagingOrdinal(iRight),&
            iPEADFSweepStagingStartRank(iRight)]
        PEADFCandidateKeyAfter = .FALSE.
        do k = 1,5
            if (a(k) > b(k)) then
                PEADFCandidateKeyAfter = .TRUE.
                return
            else if (a(k) < b(k)) then
                return
            end if
        end do
    end function PEADFCandidateKeyAfter
end subroutine CommitPEADFSweepGeneration
