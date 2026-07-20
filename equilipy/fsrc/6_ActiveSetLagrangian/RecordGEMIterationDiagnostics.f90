!> \brief Reset the current Lagrangian GEM diagnostic scalars.
!!
!! \details Clears line-search and raw phase-amount diagnostic values at the
!! beginning of a Lagrangian iteration.  This prevents a Newton failure from
!! reporting stale line-search information from the previous iteration.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    RecordGEMIterationDiagnostics.f90
    !> \brief   Record per-iteration diagnostics for Lagrangian GEM.
    !> \author  S.Y. Kwon
    !> \date    Jun. 24, 2026
    !> \sa      RunLagrangianGEM.f90
    !> \sa      GEMNewton.f90
    !> \sa      GEMLineSearch.f90
    !> \sa      CompFunctionNorm.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Recorded phase, residual, line-search, complementarity, and active-set evidence for each GEM iteration.
    !
    ! Purpose:
    ! ========
    !
    !> \details These helper routines keep the Lagrangian diagnostics
    !! synchronized with the main solver loop without changing the numerical
    !! update.  Current-iteration scalar diagnostics are reset before Newton
    !! and copied into history arrays after the iteration state has been
    !! evaluated.
    !
    ! Required input variables:
    ! =========================
    !
    ! dMolesPhase             Current active-set phase amounts.
    ! iterGlobal              Current Lagrangian iteration index.
    ! dGEMFunctionNorm        Combined residual norm computed by CompFunctionNorm.
    ! dGEMMassBalanceNorm     Mass-balance residual norm.
    ! dGEMChemicalPotentialNorm
    !                         Chemical-potential residual norm.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! dGEM*History            Per-iteration diagnostic history arrays.
    ! dGEMLineSearch*         Current line-search diagnostic scalars and phase-mole snapshots.
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
    ! RunLagrangianGEM        Resets and records diagnostics around each fixed-active-set iteration.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Diagnostics are read-only with respect to the solver state that
    !   controls convergence.  They must not alter phase amounts, species
    !   amounts, elemental potentials, or phase-change decisions.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine ResetGEMIterationDiagnostics
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    iGEMLineSearchIterationCount = 0
    iGEMLineSearchNegativeFactorCount = 0
    iGEMLineSearchNegativePhaseCount = 0
    iGEMLineSearchFloorCount = 0
    iGEMLineSearchNoDescent = 0
    iGEMLineSearchNoDescentClass = 0
    iGEMPreLMNoDescent = 0
    iGEMPreLMNoDescentClass = 0
    iGEMPreLMNewtonInfo = 0
    iGEMPreLMKKTSize = 0
    iGEMPreLMNPrimal = 0
    iGEMPreLMNConstraint = 0
    iGEMPreLMPivot1x1 = 0
    iGEMPreLMPivot2x2 = 0
    iGEMPreLMPivotPositive = 0
    iGEMPreLMPivotNegative = 0
    iGEMPreLMPivotZero = 0
    iGEMNewtonKKTSize = 0
    iGEMNewtonPivot1x1Count = 0
    iGEMNewtonPivot2x2Count = 0
    iGEMNewtonPivotPositiveCount = 0
    iGEMNewtonPivotNegativeCount = 0
    iGEMNewtonPivotZeroCount = 0
    iGEMInertiaRegularizationAttemptedUsed = 0
    iGEMInertiaRegularizationAcceptedUsed = 0
    iGEMInertiaRegularizationFailedUsed = 0
    iGEMInertiaRegularizationStepLast = 0
    iGEMResidualLMUsed = 0
    iGEMResidualLMRawNegativeUsed = 0
    iGEMResidualLMNoDescentUsed = 0
    iGEMResidualLMNoDescentClass = 0
    iGEMInvalidCompBoundAttemptedUsed = 0
    iGEMInvalidCompBoundAcceptedUsed = 0
    iGEMInvalidCompBoundRejectedUsed = 0
    iGEMInvalidCompBoundVerdict = 0
    iGEMCoalesceUsed = 0
    iGEMStabilizeUsed = 0
    iGEMBoundaryRemovalUsed = 0
    iGEMRawNegativeRemovalUsed = 0
    iGEMBoundaryPinnedRemovalUsed = 0
    iGEMTinyBoundaryRemovalUsed = 0
    iGEMTraceRemoveUsed = 0
    iGEMTraceReinjectUsed = 0
    iGEMCEFRetryActivationUsed = 0
    iGEMSUBGQRideAlongUsed = 0
    iODPairCount = 0
    iODOrderingModeCount = 0
    iGEMRawNegComp = 0
    iGEMBoundaryRemovalSlot = 0
    iGEMBoundaryRemovalPhase = 0
    iGEMTinyBoundaryRemovalSlot = 0
    iGEMTinyBoundaryRemovalPhase = 0
    iGEMBoundaryRankGuardUsed = 0
    iGEMBoundaryRankGuardSlot = 0
    iGEMBoundaryRankGuardPhase = 0
    iGEMTraceRemoveSpecies = 0
    iGEMTraceRemovePhase = 0
    iGEMTraceRemoveCount = 0
    iGEMTraceReinjectSpecies = 0
    iGEMTraceReinjectPhase = 0
    iGEMTraceReinjectCount = 0

    dGEMLineSearchInitialNorm = dGEMFunctionNorm
    dGEMLineSearchBestNorm = dGEMFunctionNorm
    dGEMLineSearchFinalNorm = dGEMFunctionNorm
    dGEMLineSearchInitialStep = 0D0
    dGEMLineSearchBestStep = 0D0
    dGEMLineSearchFinalStep = 0D0
    dGEMLineSearchMinRawPhaseMoles = 0D0
    dGEMLineSearchMinFinalPhaseMoles = 0D0
    dGEMLineSearchInitialGibbs = 0D0
    dGEMLineSearchBestGibbs = 0D0
    dGEMLineSearchFinalGibbs = 0D0
    dGEMLineSearchInitialMerit = 0D0
    dGEMLineSearchBestMerit = 0D0
    dGEMLineSearchFinalMerit = 0D0
    dGEMLineSearchMeritCandNorm = 0D0
    dGEMLineSearchMeritCandMass = 0D0
    dGEMLineSearchMeritCandStep = 0D0
    dGEMLineSearchMeritCandGibbs = 0D0
    dGEMLineSearchMeritCandMerit = 0D0
    dGEMPreLMInitNorm = 0D0
    dGEMPreLMBestNorm = 0D0
    dGEMPreLMInitGibbs = 0D0
    dGEMPreLMBestGibbs = 0D0
    dGEMPreLMInitMerit = 0D0
    dGEMPreLMBestMerit = 0D0
    dGEMPreLMMeritCandNorm = 0D0
    dGEMPreLMMeritCandMass = 0D0
    dGEMPreLMMeritCandStep = 0D0
    dGEMPreLMMeritCandGibbs = 0D0
    dGEMPreLMMeritCandMerit = 0D0
    dGEMPreLMFinalNorm = 0D0
    dGEMPreLMFinalGibbs = 0D0
    dGEMPreLMFinalMerit = 0D0
    dGEMPreLMDirNormSlope = 0D0
    dGEMPreLMDirGibbsSlope = 0D0
    dGEMPreLMDirMeritSlope = 0D0
    dGEMPreLMDirectionNorm = 0D0
    dGEMPreLMMinPivotScale = 0D0
    dGEMPreLMMaxPivotScale = 0D0
    dGEMNewtonDirNormSlope = 0D0
    dGEMNewtonDirGibbsSlope = 0D0
    dGEMNewtonDirMeritSlope = 0D0
    dGEMNewtonMinPivotScale = 0D0
    dGEMNewtonMaxPivotScale = 0D0
    dGEMNewtonDirectionNorm = 0D0
    dODCompDist = 0D0
    dODOrderNorm = 0D0
    dODOrderingEigenMin = 0D0
    iGEMRawNegativePhaseSlot = 0
    iGEMRawNegativePhaseSoln = 0
    iGEMRawNegativePhaseSpecies = 0
    dGEMRawNegativePhaseAmount = 0D0
    dGEMRawNegativePhaseDirection = 0D0
    dGEMRawNegPhaseResidual = 0D0
    dGEMInvalidCompBoundPhi = 0D0

    if (allocated(dGEMLSRawPhaseMoles)) then
        dGEMLSRawPhaseMoles = dMolesPhase
    end if
    if (allocated(dGEMLSFinalPhaseMoles)) then
        dGEMLSFinalPhaseMoles = dMolesPhase
    end if

    return

end subroutine ResetGEMIterationDiagnostics


!> \brief Snapshot the analytical direction before residual-LM overwrites it.
!!
!! \details Copies the current Newton, line-search, and inertia diagnostics to
!! pre-LM scalar storage.  Residual-LM event rows use these values so C3 merit
!! design sees the failed analytical direction, not the fallback direction.
!!
subroutine SnapshotGEMPreLMDiagnostics
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    dGEMPreLMInitNorm = dGEMLineSearchInitialNorm
    dGEMPreLMBestNorm = dGEMLineSearchBestNorm
    dGEMPreLMFinalNorm = dGEMLineSearchFinalNorm
    dGEMPreLMInitGibbs = dGEMLineSearchInitialGibbs
    dGEMPreLMBestGibbs = dGEMLineSearchBestGibbs
    dGEMPreLMFinalGibbs = dGEMLineSearchFinalGibbs
    dGEMPreLMInitMerit = dGEMLineSearchInitialMerit
    dGEMPreLMBestMerit = dGEMLineSearchBestMerit
    dGEMPreLMFinalMerit = dGEMLineSearchFinalMerit
    dGEMPreLMMeritCandNorm = dGEMLineSearchMeritCandNorm
    dGEMPreLMMeritCandMass = dGEMLineSearchMeritCandMass
    dGEMPreLMMeritCandStep = dGEMLineSearchMeritCandStep
    dGEMPreLMMeritCandGibbs = dGEMLineSearchMeritCandGibbs
    dGEMPreLMMeritCandMerit = dGEMLineSearchMeritCandMerit
    dGEMPreLMDirNormSlope = dGEMNewtonDirNormSlope
    dGEMPreLMDirGibbsSlope = dGEMNewtonDirGibbsSlope
    dGEMPreLMDirMeritSlope = dGEMNewtonDirMeritSlope
    dGEMPreLMDirectionNorm = dGEMNewtonDirectionNorm
    dGEMPreLMMinPivotScale = dGEMNewtonMinPivotScale
    dGEMPreLMMaxPivotScale = dGEMNewtonMaxPivotScale
    iGEMPreLMNewtonInfo = iGEMNewtonDSYSVInfo
    iGEMPreLMKKTSize = iGEMNewtonKKTSize
    iGEMPreLMNPrimal = nGEMCEFPhaseVariables + nGEMCEFSiteVariables
    iGEMPreLMNConstraint = nElements
    iGEMPreLMPivot1x1 = iGEMNewtonPivot1x1Count
    iGEMPreLMPivot2x2 = iGEMNewtonPivot2x2Count
    iGEMPreLMPivotPositive = iGEMNewtonPivotPositiveCount
    iGEMPreLMPivotNegative = iGEMNewtonPivotNegativeCount
    iGEMPreLMPivotZero = iGEMNewtonPivotZeroCount
    call ComputePreLMODTransferAlignment

    return

end subroutine SnapshotGEMPreLMDiagnostics


!> \brief Measure rejected-direction alignment with active order/disorder transfer.
!!
!! \details C3-c1b records whether the Newton direction rejected into
!! residual-LM points mostly along a physical ordered/disordered phase-transfer
!! coordinate.  This is diagnostic only; it does not classify or accept steps.
!!
subroutine ComputePreLMODTransferAlignment
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: iOrderedPhase, iCompanionPhase, iPhaseVar, iSiteVar
    integer :: iOrderedVar, iCompanionVar
    integer :: nModeCapacity, nModeOut, iOrderingInfo, iEigenSet
    integer :: iPhaseID, iSub, jSub, iCon, iRef, iActive, iGroup
    integer :: nSiteOut, nGroup, nActiveCon, nTrialGroup, nTrialActive, nModeLocal
    integer :: iSiteA, iSiteB, iSiteC, iSiteD, iMode
    integer :: iSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
    integer :: iGroupSub(nMaxSublatticeSys), iTrialGroup(nMaxSublatticeSys)
    integer :: iActiveCon(nMaxConstituentSys), iTrialActive(nMaxConstituentSys)
    real(8) :: dDirectionNorm, dTransferComponent, dAlignment
    real(8) :: dOrderingEigenMinLocal
    real(8) :: dAverage, dSiteDirectionNorm, dProjection, dProjectionNorm
    logical :: lSameList
    real(8), allocatable :: dOrderingHessian(:,:), dOrderingGradient(:)
    real(8), allocatable :: dOrderingBasis(:,:), dSiteDirection(:), dTrialMode(:)

    dGEMPreLMODAlign = 0D0
    dGEMPreLMODEigen = 0D0
    iGEMPreLMODOrdPhase = 0
    iGEMPreLMODCompPhase = 0
    iEigenSet = 0

    if (nGEMCEFPhaseVariables <= 0) return
    if (.NOT.allocated(dGEMCEFPhaseDirection)) return
    if (.NOT.allocated(iGEMCEFPhaseSoln)) return
    if (.NOT.allocated(iDisorderedPhase)) return

    dDirectionNorm = 0D0
    dDirectionNorm = dDirectionNorm + &
        SUM(dGEMCEFPhaseDirection(1:nGEMCEFPhaseVariables)**2)
    if ((nGEMCEFSiteVariables > 0).AND.allocated(dGEMCEFSiteDirection)) then
        dDirectionNorm = dDirectionNorm + &
            SUM(dGEMCEFSiteDirection(1:nGEMCEFSiteVariables)**2)
    end if
    dDirectionNorm = DSQRT(DMAX1(dDirectionNorm, 0D0))
    if (dDirectionNorm <= 1D-300) return

    nModeCapacity = MAX(1, nMaxSublatticeSys*nMaxConstituentSys)
    allocate(dOrderingHessian(nModeCapacity,nModeCapacity))
    allocate(dOrderingGradient(nModeCapacity))
    allocate(dOrderingBasis(nModeCapacity,nModeCapacity))
    allocate(dSiteDirection(nModeCapacity))
    allocate(dTrialMode(nModeCapacity))

    do iOrderedPhase = 1, nSolnPhasesSys
        if (iOrderedPhase > SIZE(iDisorderedPhase)) exit
        if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') cycle
        iPhaseID = iPhaseSublattice(iOrderedPhase)
        if (iPhaseID <= 0) cycle

        iCompanionPhase = iDisorderedPhase(iOrderedPhase)
        if (allocated(iODCompanionPhase)) then
            if (iOrderedPhase <= SIZE(iODCompanionPhase)) then
                if ((iODCompanionPhase(iOrderedPhase) > 0).AND.&
                    (iODCompanionPhase(iOrderedPhase) <= nSolnPhasesSys)) then
                    iCompanionPhase = iODCompanionPhase(iOrderedPhase)
                end if
            end if
        end if
        if ((iCompanionPhase <= 0).OR.(iCompanionPhase > nSolnPhasesSys)) cycle

        iOrderedVar = 0
        iCompanionVar = 0
        do iPhaseVar = 1, nGEMCEFPhaseVariables
            if (iGEMCEFPhaseSoln(iPhaseVar) == iOrderedPhase) iOrderedVar = iPhaseVar
            if (iGEMCEFPhaseSoln(iPhaseVar) == iCompanionPhase) iCompanionVar = iPhaseVar
        end do

        if (iOrderedVar > 0) then
            call CompOrderingModeSUBOM(iOrderedPhase, nModeCapacity, dOrderingHessian, &
                dOrderingGradient, dOrderingEigenMinLocal, nModeOut, iOrderingInfo)
            if ((iOrderingInfo == 0).AND.(nModeOut > 0)) then
                if ((iEigenSet == 0).OR.&
                    (dOrderingEigenMinLocal < dGEMPreLMODEigen)) then
                    dGEMPreLMODEigen = dOrderingEigenMinLocal
                    iEigenSet = 1
                end if
            end if
        end if

        if (iOrderedVar > 0) then
            iSiteIndex = 0
            nSiteOut = 0
            do iSub = 1, nSublatticePhase(iPhaseID)
                do iCon = 1, nConstituentSublattice(iPhaseID,iSub)
                    nSiteOut = nSiteOut + 1
                    if (nSiteOut <= nModeCapacity) iSiteIndex(iSub,iCon) = nSiteOut
                end do
            end do

            iGroupSub = 0
            iActiveCon = 0
            nGroup = 0
            nActiveCon = 0
            do iSub = 1, nSublatticePhase(iPhaseID)
                if (nConstituentSublattice(iPhaseID,iSub) <= 1) cycle
                nTrialGroup = 1
                iTrialGroup = 0
                iTrialGroup(1) = iSub
                do jSub = iSub + 1, nSublatticePhase(iPhaseID)
                    if (nConstituentSublattice(iPhaseID,iSub) /= &
                        nConstituentSublattice(iPhaseID,jSub)) cycle
                    lSameList = .TRUE.
                    do iCon = 1, nConstituentSublattice(iPhaseID,iSub)
                        if (TRIM(cConstituentNameSUB(iPhaseID,iSub,iCon)) /= &
                            TRIM(cConstituentNameSUB(iPhaseID,jSub,iCon))) then
                            lSameList = .FALSE.
                            exit
                        end if
                    end do
                    if (lSameList) then
                        nTrialGroup = nTrialGroup + 1
                        iTrialGroup(nTrialGroup) = jSub
                    end if
                end do
                if (nTrialGroup < 2) cycle

                nTrialActive = 0
                iTrialActive = 0
                do iCon = 1, nConstituentSublattice(iPhaseID,iTrialGroup(1))
                    dAverage = 0D0
                    do iGroup = 1, nTrialGroup
                        dAverage = dAverage + &
                            dGEMCEFPhaseSiteLast(iOrderedVar,iTrialGroup(iGroup),iCon)
                    end do
                    dAverage = dAverage / DFLOAT(nTrialGroup)
                    if (dAverage > 1D-12) then
                        nTrialActive = nTrialActive + 1
                        iTrialActive(nTrialActive) = iCon
                    end if
                end do
                if (nTrialActive < 2) cycle
                nGroup = nTrialGroup
                nActiveCon = nTrialActive
                iGroupSub(1:nGroup) = iTrialGroup(1:nGroup)
                iActiveCon(1:nActiveCon) = iTrialActive(1:nActiveCon)
                exit
            end do

            if ((nGroup >= 2).AND.(nActiveCon >= 2).AND.(nSiteOut <= nModeCapacity)) then
                dOrderingBasis = 0D0
                nModeLocal = 0
                iRef = iActiveCon(1)
                do iGroup = 2, nGroup
                    iSub = iGroupSub(iGroup)
                    do iActive = 2, nActiveCon
                        iCon = iActiveCon(iActive)
                        iSiteA = iSiteIndex(iSub,iCon)
                        iSiteB = iSiteIndex(iSub,iRef)
                        iSiteC = iSiteIndex(iGroupSub(1),iCon)
                        iSiteD = iSiteIndex(iGroupSub(1),iRef)
                        if (MIN(iSiteA, iSiteB, iSiteC, iSiteD) <= 0) cycle
                        dTrialMode = 0D0
                        dTrialMode(iSiteA) = 1D0
                        dTrialMode(iSiteB) = -1D0
                        dTrialMode(iSiteC) = -1D0
                        dTrialMode(iSiteD) = 1D0
                        do iMode = 1, nModeLocal
                            dProjection = SUM(dTrialMode(1:nSiteOut) * &
                                dOrderingBasis(1:nSiteOut,iMode))
                            dTrialMode(1:nSiteOut) = dTrialMode(1:nSiteOut) - &
                                dProjection * dOrderingBasis(1:nSiteOut,iMode)
                        end do
                        dSiteDirectionNorm = DSQRT(SUM(dTrialMode(1:nSiteOut) * &
                            dTrialMode(1:nSiteOut)))
                        if (dSiteDirectionNorm <= 1D-12) cycle
                        if (nModeLocal >= nModeCapacity) exit
                        nModeLocal = nModeLocal + 1
                        dOrderingBasis(1:nSiteOut,nModeLocal) = &
                            dTrialMode(1:nSiteOut) / dSiteDirectionNorm
                    end do
                end do

                if (nModeLocal > 0) then
                    dSiteDirection = 0D0
                    do iSiteVar = 1, nGEMCEFSiteVariables
                        if (iGEMCEFVarSolnPhase(iSiteVar) /= iOrderedPhase) cycle
                        if (iGEMCEFVarPhaseID(iSiteVar) <= 0) cycle
                        iSub = iGEMCEFVarSub(iSiteVar)
                        iCon = iGEMCEFVarCon(iSiteVar)
                        iRef = iGEMCEFVarRef(iSiteVar)
                        if (MIN(iSub, iCon, iRef) <= 0) cycle
                        iSiteA = iSiteIndex(iSub,iCon)
                        iSiteB = iSiteIndex(iSub,iRef)
                        if (MIN(iSiteA, iSiteB) <= 0) cycle
                        dSiteDirection(iSiteA) = dSiteDirection(iSiteA) + &
                            dGEMCEFSiteDirection(iSiteVar)
                        dSiteDirection(iSiteB) = dSiteDirection(iSiteB) - &
                            dGEMCEFSiteDirection(iSiteVar)
                    end do
                    dProjectionNorm = 0D0
                    do iMode = 1, nModeLocal
                        dProjection = SUM(dSiteDirection(1:nSiteOut) * &
                            dOrderingBasis(1:nSiteOut,iMode))
                        dProjectionNorm = dProjectionNorm + dProjection*dProjection
                    end do
                    dAlignment = DSQRT(DMAX1(dProjectionNorm, 0D0)) / dDirectionNorm
                    if (dAlignment > dGEMPreLMODAlign) then
                        dGEMPreLMODAlign = dAlignment
                        iGEMPreLMODOrdPhase = iOrderedPhase
                        iGEMPreLMODCompPhase = iCompanionPhase
                    end if
                end if
            end if
        end if

        if ((iOrderedVar <= 0).OR.(iCompanionVar <= 0)) cycle

        dTransferComponent = &
            DABS(dGEMCEFPhaseDirection(iOrderedVar) - &
                 dGEMCEFPhaseDirection(iCompanionVar)) / DSQRT(2D0)
        dAlignment = dTransferComponent / dDirectionNorm
        if (dAlignment > dGEMPreLMODAlign) then
            dGEMPreLMODAlign = dAlignment
            iGEMPreLMODOrdPhase = iOrderedPhase
            iGEMPreLMODCompPhase = iCompanionPhase
        end if
    end do

    if (allocated(dOrderingHessian)) deallocate(dOrderingHessian)
    if (allocated(dOrderingGradient)) deallocate(dOrderingGradient)
    if (allocated(dOrderingBasis)) deallocate(dOrderingBasis)
    if (allocated(dSiteDirection)) deallocate(dSiteDirection)
    if (allocated(dTrialMode)) deallocate(dTrialMode)

    return

end subroutine ComputePreLMODTransferAlignment


!> \brief Record one successful residual-LM fallback activation.
!!
!! \details This call-site-stable buffer is independent of iterGlobal-indexed
!! histories, so it survives repeated PEA-polish and postprocess Lagrangian
!! calls within one equilibrium calculation.
!!
subroutine RecordGEMResidualLMEvent(iReason)
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer, intent(in) :: iReason
    integer :: iEvent, iCallSite

    if (.NOT.allocated(iGEMResidualLMEventCallSite)) then
        nGEMResidualLMEventCapacity = MAX(4096, 32*MAX(1, iterGlobalMax))
        allocate(iGEMResidualLMEventCallSite(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventIterPEA(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventIterGlobal(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventReason(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventNoDescentClass(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventNewtonInfo(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventKKTSize(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventNPrimal(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventNConstraint(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventPivot1x1(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventPivot2x2(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventPivotPositive(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventPivotNegative(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventPivotZero(nGEMResidualLMEventCapacity))
        allocate(iGEMLMEventODOrdPhase(nGEMResidualLMEventCapacity))
        allocate(iGEMLMEventODCompPhase(nGEMResidualLMEventCapacity))
        allocate(iGEMResidualLMEventAssemblage(nGEMResidualLMEventCapacity,nElements))
        allocate(dGEMResidualLMEventNormSlope(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventGibbsSlope(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventMeritSlope(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventInitialNorm(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventBestNorm(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventFinalNorm(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventInitialGibbs(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventBestGibbs(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventFinalGibbs(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventInitialMerit(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventBestMerit(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventFinalMerit(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventMeritCandNorm(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventMeritCandMass(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventMeritCandStep(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventMeritCandGibbs(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventMeritCandMerit(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventDirectionNorm(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventMinPivot(nGEMResidualLMEventCapacity))
        allocate(dGEMResidualLMEventMaxPivot(nGEMResidualLMEventCapacity))
        allocate(dGEMLMEventODAlign(nGEMResidualLMEventCapacity))
        allocate(dGEMLMEventODEigen(nGEMResidualLMEventCapacity))

        iGEMResidualLMEventCallSite = 0
        iGEMResidualLMEventIterPEA = 0
        iGEMResidualLMEventIterGlobal = 0
        iGEMResidualLMEventReason = 0
        iGEMResidualLMEventNoDescentClass = 0
        iGEMResidualLMEventNewtonInfo = 0
        iGEMResidualLMEventKKTSize = 0
        iGEMResidualLMEventNPrimal = 0
        iGEMResidualLMEventNConstraint = 0
        iGEMResidualLMEventPivot1x1 = 0
        iGEMResidualLMEventPivot2x2 = 0
        iGEMResidualLMEventPivotPositive = 0
        iGEMResidualLMEventPivotNegative = 0
        iGEMResidualLMEventPivotZero = 0
        iGEMLMEventODOrdPhase = 0
        iGEMLMEventODCompPhase = 0
        iGEMResidualLMEventAssemblage = 0
        dGEMResidualLMEventNormSlope = 0D0
        dGEMResidualLMEventGibbsSlope = 0D0
        dGEMResidualLMEventMeritSlope = 0D0
        dGEMResidualLMEventInitialNorm = 0D0
        dGEMResidualLMEventBestNorm = 0D0
        dGEMResidualLMEventFinalNorm = 0D0
        dGEMResidualLMEventInitialGibbs = 0D0
        dGEMResidualLMEventBestGibbs = 0D0
        dGEMResidualLMEventFinalGibbs = 0D0
        dGEMResidualLMEventInitialMerit = 0D0
        dGEMResidualLMEventBestMerit = 0D0
        dGEMResidualLMEventFinalMerit = 0D0
        dGEMResidualLMEventMeritCandNorm = 0D0
        dGEMResidualLMEventMeritCandMass = 0D0
        dGEMResidualLMEventMeritCandStep = 0D0
        dGEMResidualLMEventMeritCandGibbs = 0D0
        dGEMResidualLMEventMeritCandMerit = 0D0
        dGEMResidualLMEventDirectionNorm = 0D0
        dGEMResidualLMEventMinPivot = 0D0
        dGEMResidualLMEventMaxPivot = 0D0
        dGEMLMEventODAlign = 0D0
        dGEMLMEventODEigen = 0D0
    end if

    if (nGEMResidualLMEvent >= nGEMResidualLMEventCapacity) return

    iEvent = nGEMResidualLMEvent + 1
    nGEMResidualLMEvent = iEvent

    iCallSite = iGEMLagrangianCallSite
    if (lPEALagrangianPolishActive) iCallSite = 2
    if (lPostProcess.AND.((iCallSite < 3).OR.(iCallSite > 4))) iCallSite = 3
    if (iCallSite <= 0) iCallSite = 1

    iGEMResidualLMEventCallSite(iEvent) = iCallSite
    iGEMResidualLMEventIterPEA(iEvent) = iterPEA
    iGEMResidualLMEventIterGlobal(iEvent) = iterGlobal
    iGEMResidualLMEventReason(iEvent) = iReason
    iGEMResidualLMEventNoDescentClass(iEvent) = iGEMPreLMNoDescentClass
    iGEMResidualLMEventNewtonInfo(iEvent) = iGEMPreLMNewtonInfo
    iGEMResidualLMEventKKTSize(iEvent) = iGEMPreLMKKTSize
    iGEMResidualLMEventNPrimal(iEvent) = iGEMPreLMNPrimal
    iGEMResidualLMEventNConstraint(iEvent) = iGEMPreLMNConstraint
    iGEMResidualLMEventPivot1x1(iEvent) = iGEMPreLMPivot1x1
    iGEMResidualLMEventPivot2x2(iEvent) = iGEMPreLMPivot2x2
    iGEMResidualLMEventPivotPositive(iEvent) = iGEMPreLMPivotPositive
    iGEMResidualLMEventPivotNegative(iEvent) = iGEMPreLMPivotNegative
    iGEMResidualLMEventPivotZero(iEvent) = iGEMPreLMPivotZero
    iGEMLMEventODOrdPhase(iEvent) = iGEMPreLMODOrdPhase
    iGEMLMEventODCompPhase(iEvent) = iGEMPreLMODCompPhase
    iGEMResidualLMEventAssemblage(iEvent,1:nElements) = iAssemblage
    dGEMResidualLMEventNormSlope(iEvent) = dGEMPreLMDirNormSlope
    dGEMResidualLMEventGibbsSlope(iEvent) = dGEMPreLMDirGibbsSlope
    dGEMResidualLMEventMeritSlope(iEvent) = dGEMPreLMDirMeritSlope
    dGEMResidualLMEventInitialNorm(iEvent) = dGEMPreLMInitNorm
    dGEMResidualLMEventBestNorm(iEvent) = dGEMPreLMBestNorm
    dGEMResidualLMEventFinalNorm(iEvent) = dGEMPreLMFinalNorm
    dGEMResidualLMEventInitialGibbs(iEvent) = dGEMPreLMInitGibbs
    dGEMResidualLMEventBestGibbs(iEvent) = dGEMPreLMBestGibbs
    dGEMResidualLMEventFinalGibbs(iEvent) = dGEMPreLMFinalGibbs
    dGEMResidualLMEventInitialMerit(iEvent) = dGEMPreLMInitMerit
    dGEMResidualLMEventBestMerit(iEvent) = dGEMPreLMBestMerit
    dGEMResidualLMEventFinalMerit(iEvent) = dGEMPreLMFinalMerit
    dGEMResidualLMEventMeritCandNorm(iEvent) = dGEMPreLMMeritCandNorm
    dGEMResidualLMEventMeritCandMass(iEvent) = dGEMPreLMMeritCandMass
    dGEMResidualLMEventMeritCandStep(iEvent) = dGEMPreLMMeritCandStep
    dGEMResidualLMEventMeritCandGibbs(iEvent) = dGEMPreLMMeritCandGibbs
    dGEMResidualLMEventMeritCandMerit(iEvent) = dGEMPreLMMeritCandMerit
    dGEMResidualLMEventDirectionNorm(iEvent) = dGEMPreLMDirectionNorm
    dGEMResidualLMEventMinPivot(iEvent) = dGEMPreLMMinPivotScale
    dGEMResidualLMEventMaxPivot(iEvent) = dGEMPreLMMaxPivotScale
    dGEMLMEventODAlign(iEvent) = dGEMPreLMODAlign
    dGEMLMEventODEigen(iEvent) = dGEMPreLMODEigen

    return

end subroutine RecordGEMResidualLMEvent


!> \brief Store the latest Lagrangian GEM diagnostic values in history arrays.
!!
!! \details Copies residual families, Newton status, line-search status, raw
!! full-step phase amounts, accepted phase amounts, and trace-species counts
!! into iteration-indexed arrays exposed through f2py.
!!
subroutine RecordGEMIterationDiagnostics
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: iHistory, iSlot, iThermoPhase, iFirst, iLast, iOrdinal

    if (iterGlobal < 1) return
    if (iterGlobal > iterGlobalMax) return

    iHistory = iterGlobal
    call ComputeODDiagnostics

    dGEMLineSearchMinFinalPhaseMoles = MINVAL(dMolesPhase)
    if (allocated(dGEMLSFinalPhaseMoles)) then
        dGEMLSFinalPhaseMoles = dMolesPhase
    end if

    if (allocated(dGEMNormHist)) dGEMNormHist(iHistory) = dGEMFunctionNorm
    if (allocated(dGEMMassNormHist)) dGEMMassNormHist(iHistory) = dGEMMassBalanceNorm
    if (allocated(dGEMChemNormHist)) then
        dGEMChemNormHist(iHistory) = dGEMChemicalPotentialNorm
    end if
    if (allocated(dGEMSolnChemNormHist)) then
        dGEMSolnChemNormHist(iHistory) = dGEMSolutionChemicalPotentialNorm
    end if
    if (allocated(dGEMCondChemNormHist)) then
        dGEMCondChemNormHist(iHistory) = dGEMCondensedChemicalPotentialNorm
    end if
    if (allocated(dGEMSublExchangeNormHist)) then
        dGEMSublExchangeNormHist(iHistory) = dSublatticeExchangeNorm
    end if
    if (allocated(dGEMNewtonSymResidualHist)) then
        dGEMNewtonSymResidualHist(iHistory) = dGEMNewtonSymmetryResidual
    end if
    if (allocated(dGEMNewtonMinPivotScaleHist)) then
        dGEMNewtonMinPivotScaleHist(iHistory) = dGEMNewtonMinPivotScale
    end if
    if (allocated(dGEMNewtonMaxPivotScaleHist)) then
        dGEMNewtonMaxPivotScaleHist(iHistory) = dGEMNewtonMaxPivotScale
    end if
    if (allocated(dGEMNewtonDirectionNormHist)) then
        dGEMNewtonDirectionNormHist(iHistory) = dGEMNewtonDirectionNorm
    end if
    if (lSUBOMTwoSetCandidateEnabled.AND.allocated(iActiveSlotThermoPhase).AND.&
        allocated(iActiveSlotIdentityOrdinal).AND.allocated(dActiveSlotMolFraction)) then
        do iSlot = 1, MIN(nElements, SIZE(iActiveSlotThermoPhase))
            iThermoPhase = iActiveSlotThermoPhase(iSlot)
            if ((iThermoPhase < 1).OR.(iThermoPhase > nSolnPhasesSys)) cycle
            if (TRIM(cSolnPhaseType(iThermoPhase)) /= 'SUBOM') cycle
            iOrdinal = iActiveSlotIdentityOrdinal(iSlot)
            if (iOrdinal <= 0) cycle
            iFirst = nSpeciesPhase(iThermoPhase-1) + 1
            iLast  = nSpeciesPhase(iThermoPhase)
            if ((iFirst < 1).OR.(iLast > SIZE(dActiveSlotMolFraction,2))) cycle
            if (SUM(dActiveSlotMolFraction(iSlot,iFirst:iLast)) <= 1D-300) cycle
            call RecordSUBOMTwoSetTrace(SUBOM_TRACE_SOLVER_SLOT, iThermoPhase, iOrdinal, &
                iSlot, iLast-iFirst+1, dActiveSlotMolFraction(iSlot,iFirst:iLast))
        end do
    end if
    if (allocated(dGEMNewtonDirNormSlopeHist)) then
        dGEMNewtonDirNormSlopeHist(iHistory) = dGEMNewtonDirNormSlope
    end if
    if (allocated(dGEMNewtonDirGibbsSlopeHist)) then
        dGEMNewtonDirGibbsSlopeHist(iHistory) = dGEMNewtonDirGibbsSlope
    end if
    if (allocated(dGEMNewtonDirMeritSlopeHist)) then
        dGEMNewtonDirMeritSlopeHist(iHistory) = dGEMNewtonDirMeritSlope
    end if
    if (allocated(dGEMLSInitialNormHist)) then
        dGEMLSInitialNormHist(iHistory) = dGEMLineSearchInitialNorm
    end if
    if (allocated(dGEMLSBestNormHist)) then
        dGEMLSBestNormHist(iHistory) = dGEMLineSearchBestNorm
    end if
    if (allocated(dGEMLSFinalNormHist)) then
        dGEMLSFinalNormHist(iHistory) = dGEMLineSearchFinalNorm
    end if
    if (allocated(dGEMLSInitialStepHist)) then
        dGEMLSInitialStepHist(iHistory) = dGEMLineSearchInitialStep
    end if
    if (allocated(dGEMLSBestStepHist)) then
        dGEMLSBestStepHist(iHistory) = dGEMLineSearchBestStep
    end if
    if (allocated(dGEMLSFinalStepHist)) then
        dGEMLSFinalStepHist(iHistory) = dGEMLineSearchFinalStep
    end if
    if (allocated(dGEMLSMinRawPhaseMolesHist)) then
        dGEMLSMinRawPhaseMolesHist(iHistory) = dGEMLineSearchMinRawPhaseMoles
    end if
    if (allocated(dGEMLSMinFinalPhaseMolesHist)) then
        dGEMLSMinFinalPhaseMolesHist(iHistory) = dGEMLineSearchMinFinalPhaseMoles
    end if
    if (allocated(dGEMLSInitialGibbsHist)) then
        dGEMLSInitialGibbsHist(iHistory) = dGEMLineSearchInitialGibbs
    end if
    if (allocated(dGEMLSBestGibbsHist)) then
        dGEMLSBestGibbsHist(iHistory) = dGEMLineSearchBestGibbs
    end if
    if (allocated(dGEMLSFinalGibbsHist)) then
        dGEMLSFinalGibbsHist(iHistory) = dGEMLineSearchFinalGibbs
    end if
    if (allocated(dGEMLSInitialMeritHist)) then
        dGEMLSInitialMeritHist(iHistory) = dGEMLineSearchInitialMerit
    end if
    if (allocated(dGEMLSBestMeritHist)) then
        dGEMLSBestMeritHist(iHistory) = dGEMLineSearchBestMerit
    end if
    if (allocated(dGEMLSFinalMeritHist)) then
        dGEMLSFinalMeritHist(iHistory) = dGEMLineSearchFinalMerit
    end if
    if (allocated(dGEMPreLMInitNormHist)) then
        dGEMPreLMInitNormHist(iHistory) = dGEMPreLMInitNorm
    end if
    if (allocated(dGEMPreLMBestNormHist)) then
        dGEMPreLMBestNormHist(iHistory) = dGEMPreLMBestNorm
    end if
    if (allocated(dGEMPreLMInitGibbsHist)) then
        dGEMPreLMInitGibbsHist(iHistory) = dGEMPreLMInitGibbs
    end if
    if (allocated(dGEMPreLMBestGibbsHist)) then
        dGEMPreLMBestGibbsHist(iHistory) = dGEMPreLMBestGibbs
    end if
    if (allocated(dGEMPreLMInitMeritHist)) then
        dGEMPreLMInitMeritHist(iHistory) = dGEMPreLMInitMerit
    end if
    if (allocated(dGEMPreLMBestMeritHist)) then
        dGEMPreLMBestMeritHist(iHistory) = dGEMPreLMBestMerit
    end if
    if (allocated(dGEMPreLMMeritCandNormHist)) then
        dGEMPreLMMeritCandNormHist(iHistory) = dGEMPreLMMeritCandNorm
    end if
    if (allocated(dGEMPreLMMeritCandMassHist)) then
        dGEMPreLMMeritCandMassHist(iHistory) = dGEMPreLMMeritCandMass
    end if
    if (allocated(dGEMPreLMMeritCandStepHist)) then
        dGEMPreLMMeritCandStepHist(iHistory) = dGEMPreLMMeritCandStep
    end if
    if (allocated(dGEMPreLMMeritCandGibbsHist)) then
        dGEMPreLMMeritCandGibbsHist(iHistory) = dGEMPreLMMeritCandGibbs
    end if
    if (allocated(dGEMPreLMMeritCandMeritHist)) then
        dGEMPreLMMeritCandMeritHist(iHistory) = dGEMPreLMMeritCandMerit
    end if

    if (allocated(iGEMNewtonSolverHist)) iGEMNewtonSolverHist(iHistory) = iGEMNewtonSolver
    if (allocated(iGEMNewtonInfoHist)) iGEMNewtonInfoHist(iHistory) = iGEMNewtonDSYSVInfo
    if (allocated(iGEMNewtonKKTSizeHist)) iGEMNewtonKKTSizeHist(iHistory) = iGEMNewtonKKTSize
    if (allocated(iGEMNewtonPivot1x1Hist)) then
        iGEMNewtonPivot1x1Hist(iHistory) = iGEMNewtonPivot1x1Count
    end if
    if (allocated(iGEMNewtonPivot2x2Hist)) then
        iGEMNewtonPivot2x2Hist(iHistory) = iGEMNewtonPivot2x2Count
    end if
    if (allocated(iGEMNewtonPivotPositiveHist)) then
        iGEMNewtonPivotPositiveHist(iHistory) = iGEMNewtonPivotPositiveCount
    end if
    if (allocated(iGEMNewtonPivotNegativeHist)) then
        iGEMNewtonPivotNegativeHist(iHistory) = iGEMNewtonPivotNegativeCount
    end if
    if (allocated(iGEMNewtonPivotZeroHist)) then
        iGEMNewtonPivotZeroHist(iHistory) = iGEMNewtonPivotZeroCount
    end if
    if (allocated(iGEMLSIterHist)) then
        iGEMLSIterHist(iHistory) = iGEMLineSearchIterationCount
    end if
    if (allocated(iGEMLSNegFactorHist)) then
        iGEMLSNegFactorHist(iHistory) = iGEMLineSearchNegativeFactorCount
    end if
    if (allocated(iGEMLSNegPhaseHist)) then
        iGEMLSNegPhaseHist(iHistory) = iGEMLineSearchNegativePhaseCount
    end if
    if (allocated(iGEMLSFloorHist)) then
        iGEMLSFloorHist(iHistory) = iGEMLineSearchFloorCount
    end if
    if (allocated(iGEMLSNoDescentHist)) then
        iGEMLSNoDescentHist(iHistory) = iGEMLineSearchNoDescent
    end if
    if (allocated(iGEMLSNoDescentClassHist)) then
        iGEMLSNoDescentClassHist(iHistory) = iGEMLineSearchNoDescentClass
    end if
    if (allocated(iGEMPreLMNoDescentHist)) then
        iGEMPreLMNoDescentHist(iHistory) = iGEMPreLMNoDescent
    end if
    if (allocated(iGEMPreLMNoDescentClassHist)) then
        iGEMPreLMNoDescentClassHist(iHistory) = iGEMPreLMNoDescentClass
    end if
    if (allocated(iGEMTraceInactiveHist)) then
        if (allocated(lTraceSpeciesInactive)) then
            iGEMTraceInactiveHist(iHistory) = COUNT(lTraceSpeciesInactive)
        else
            iGEMTraceInactiveHist(iHistory) = 0
        end if
    end if
    if (allocated(iGEMTraceReinjectedHist)) then
        if (allocated(lTraceSpeciesReinjected)) then
            iGEMTraceReinjectedHist(iHistory) = COUNT(lTraceSpeciesReinjected)
        else
            iGEMTraceReinjectedHist(iHistory) = 0
        end if
    end if
    if (allocated(iGEMTraceRemoveHist)) iGEMTraceRemoveHist(iHistory) = iGEMTraceRemoveUsed
    if (allocated(iGEMTraceReinjectHist)) iGEMTraceReinjectHist(iHistory) = iGEMTraceReinjectUsed
    if (allocated(iGEMResidualLMHist)) iGEMResidualLMHist(iHistory) = iGEMResidualLMUsed
    if (allocated(iGEMResidualLMRawNegativeHist)) then
        iGEMResidualLMRawNegativeHist(iHistory) = iGEMResidualLMRawNegativeUsed
    end if
    if (allocated(iGEMResidualLMNoDescentHist)) then
        iGEMResidualLMNoDescentHist(iHistory) = iGEMResidualLMNoDescentUsed
    end if
    if (allocated(iGEMResidualLMNoDescentClassHist)) then
        iGEMResidualLMNoDescentClassHist(iHistory) = iGEMResidualLMNoDescentClass
    end if
    if (allocated(iGEMCoalesceHist)) iGEMCoalesceHist(iHistory) = iGEMCoalesceUsed
    if (allocated(iGEMStabilizeHist)) iGEMStabilizeHist(iHistory) = iGEMStabilizeUsed
    if (allocated(iGEMBoundaryRemovalHist)) then
        iGEMBoundaryRemovalHist(iHistory) = iGEMBoundaryRemovalUsed
    end if
    if (allocated(iGEMRawNegativeRemovalHist)) then
        iGEMRawNegativeRemovalHist(iHistory) = iGEMRawNegativeRemovalUsed
    end if
    if (allocated(iGEMRawNegativePhaseSlotHist)) then
        iGEMRawNegativePhaseSlotHist(iHistory) = iGEMRawNegativePhaseSlot
    end if
    if (allocated(iGEMRawNegativePhaseSolnHist)) then
        iGEMRawNegativePhaseSolnHist(iHistory) = iGEMRawNegativePhaseSoln
    end if
    if (allocated(iGEMRawNegativePhaseSpeciesHist)) then
        iGEMRawNegativePhaseSpeciesHist(iHistory) = iGEMRawNegativePhaseSpecies
    end if
    if (allocated(iGEMRawNegCompHist)) then
        iGEMRawNegCompHist(iHistory) = iGEMRawNegComp
    end if
    if (allocated(iGEMInvalidCompBoundAttemptedHist)) then
        iGEMInvalidCompBoundAttemptedHist(iHistory) = iGEMInvalidCompBoundAttemptedUsed
    end if
    if (allocated(iGEMInvalidCompBoundAcceptedHist)) then
        iGEMInvalidCompBoundAcceptedHist(iHistory) = iGEMInvalidCompBoundAcceptedUsed
    end if
    if (allocated(iGEMInvalidCompBoundRejectedHist)) then
        iGEMInvalidCompBoundRejectedHist(iHistory) = iGEMInvalidCompBoundRejectedUsed
    end if
    if (allocated(iGEMInvalidCompBoundVerdictHist)) then
        iGEMInvalidCompBoundVerdictHist(iHistory) = iGEMInvalidCompBoundVerdict
    end if
    if (allocated(dGEMRawNegativePhaseAmountHist)) then
        dGEMRawNegativePhaseAmountHist(iHistory) = dGEMRawNegativePhaseAmount
    end if
    if (allocated(dGEMRawNegativePhaseDirectionHist)) then
        dGEMRawNegativePhaseDirectionHist(iHistory) = dGEMRawNegativePhaseDirection
    end if
    if (allocated(dGEMRawNegPhaseResidualHist)) then
        dGEMRawNegPhaseResidualHist(iHistory) = dGEMRawNegPhaseResidual
    end if
    if (allocated(dGEMInvalidCompBoundPhiHist)) then
        dGEMInvalidCompBoundPhiHist(iHistory) = dGEMInvalidCompBoundPhi
    end if
    if (allocated(iGEMBoundaryPinnedRemovalHist)) then
        iGEMBoundaryPinnedRemovalHist(iHistory) = iGEMBoundaryPinnedRemovalUsed
    end if
    if (allocated(iGEMBoundaryRemovalSlotHist)) then
        iGEMBoundaryRemovalSlotHist(iHistory) = iGEMBoundaryRemovalSlot
    end if
    if (allocated(iGEMBoundaryRemovalPhaseHist)) then
        iGEMBoundaryRemovalPhaseHist(iHistory) = iGEMBoundaryRemovalPhase
    end if
    if (allocated(iGEMTinyBoundaryRemovalHist)) then
        iGEMTinyBoundaryRemovalHist(iHistory) = iGEMTinyBoundaryRemovalUsed
    end if
    if (allocated(iGEMTinyBoundaryRemovalSlotHist)) then
        iGEMTinyBoundaryRemovalSlotHist(iHistory) = iGEMTinyBoundaryRemovalSlot
    end if
    if (allocated(iGEMTinyBoundaryRemovalPhaseHist)) then
        iGEMTinyBoundaryRemovalPhaseHist(iHistory) = iGEMTinyBoundaryRemovalPhase
    end if
    if (allocated(iGEMBoundaryRankGuardHist)) then
        iGEMBoundaryRankGuardHist(iHistory) = iGEMBoundaryRankGuardUsed
    end if
    if (allocated(iGEMBoundaryRankGuardSlotHist)) then
        iGEMBoundaryRankGuardSlotHist(iHistory) = iGEMBoundaryRankGuardSlot
    end if
    if (allocated(iGEMBoundaryRankGuardPhaseHist)) then
        iGEMBoundaryRankGuardPhaseHist(iHistory) = iGEMBoundaryRankGuardPhase
    end if
    if (allocated(iGEMTraceRemoveSpeciesHist)) then
        iGEMTraceRemoveSpeciesHist(iHistory) = iGEMTraceRemoveSpecies
    end if
    if (allocated(iGEMTraceRemovePhaseHist)) then
        iGEMTraceRemovePhaseHist(iHistory) = iGEMTraceRemovePhase
    end if
    if (allocated(iGEMTraceRemoveCountHist)) then
        iGEMTraceRemoveCountHist(iHistory) = iGEMTraceRemoveCount
    end if
    if (allocated(iGEMTraceReinjectSpeciesHist)) then
        iGEMTraceReinjectSpeciesHist(iHistory) = iGEMTraceReinjectSpecies
    end if
    if (allocated(iGEMTraceReinjectPhaseHist)) then
        iGEMTraceReinjectPhaseHist(iHistory) = iGEMTraceReinjectPhase
    end if
    if (allocated(iGEMTraceReinjectCountHist)) then
        iGEMTraceReinjectCountHist(iHistory) = iGEMTraceReinjectCount
    end if
    if (allocated(iGEMCEFRetryActivationHist)) then
        iGEMCEFRetryActivationHist(iHistory) = iGEMCEFRetryActivationUsed
    end if
    if (allocated(iGEMSUBGQRideAlongHist)) then
        iGEMSUBGQRideAlongHist(iHistory) = iGEMSUBGQRideAlongUsed
    end if
    if (allocated(iODPairHist)) then
        iODPairHist(iHistory) = iODPairCount
    end if
    if (allocated(iODOrderingModeHist)) then
        iODOrderingModeHist(iHistory) = iODOrderingModeCount
    end if
    if (allocated(dODCompDistHist)) then
        dODCompDistHist(iHistory) = dODCompDist
    end if
    if (allocated(dODOrderNormHist)) then
        dODOrderNormHist(iHistory) = dODOrderNorm
    end if
    if (allocated(dODOrderingEigenMinHist)) then
        dODOrderingEigenMinHist(iHistory) = dODOrderingEigenMin
    end if

    if (allocated(dGEMLSRawPhaseMolesHist).AND.allocated(dGEMLSRawPhaseMoles)) then
        dGEMLSRawPhaseMolesHist(:,iHistory) = dGEMLSRawPhaseMoles
    end if
    if (allocated(dGEMLSFinalPhaseMolesHist).AND.allocated(dGEMLSFinalPhaseMoles)) then
        dGEMLSFinalPhaseMolesHist(:,iHistory) = dGEMLSFinalPhaseMoles
    end if

    return

end subroutine RecordGEMIterationDiagnostics


!> \brief Compute passive active order/disorder diagnostics for the current GEM state.
!!
!! \details Counts active ordered/DIS_PART helper pairs and records how far
!! active ordered phases are from their random manifold.  The values are
!! diagnostic only and must not control convergence or phase selection.
!!
subroutine ComputeODDiagnostics
    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: iOrderedPhase, iHelperPhase, iCompanionPhase, iOrderedSlot, iHelperSlot
    integer :: nModeCapacity, nModeOut, iOrderingInfo
    real(8) :: dCompositionDistance, dOrderingNorm, dOrderingEigenMinLocal
    real(8), dimension(nElements) :: dOrderedComposition, dHelperComposition
    real(8), allocatable :: dOrderingHessian(:,:), dOrderingGradient(:)

    iODPairCount = 0
    iODOrderingModeCount = 0
    dODCompDist = 0D0
    dODOrderNorm = 0D0
    dODOrderingEigenMin = 0D0

    if (.NOT.allocated(iDisorderedPhase)) return
    if (nSolnPhases <= 0) return

    nModeCapacity = MAX(1, nMaxSublatticeSys*nMaxConstituentSys)
    allocate(dOrderingHessian(nModeCapacity,nModeCapacity))
    allocate(dOrderingGradient(nModeCapacity))

    do iOrderedPhase = 1, nSolnPhasesSys
        if (iOrderedPhase > SIZE(iDisorderedPhase)) exit
        if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') cycle

        iHelperPhase = iDisorderedPhase(iOrderedPhase)
        if ((iHelperPhase <= 0).OR.(iHelperPhase > nSolnPhasesSys)) cycle
        if (allocated(iODCompanionPhase)) then
            if (iOrderedPhase <= SIZE(iODCompanionPhase)) then
                if ((iODCompanionPhase(iOrderedPhase) > 0).AND.&
                    (iODCompanionPhase(iOrderedPhase) <= nSolnPhasesSys)) then
                    iHelperPhase = iODCompanionPhase(iOrderedPhase)
                end if
            end if
        end if

        call FindActiveSolutionSlot(iOrderedPhase, iOrderedSlot)
        if (iOrderedSlot > 0) then
            dOrderingNorm = OrderedPhaseOrderingNorm(iOrderedPhase)
            dODOrderNorm = DMAX1(dODOrderNorm, dOrderingNorm)
            call CompOrderingModeSUBOM(iOrderedPhase, nModeCapacity, dOrderingHessian, &
                dOrderingGradient, dOrderingEigenMinLocal, nModeOut, iOrderingInfo)
            if ((iOrderingInfo == 0).AND.(nModeOut > 0)) then
                iODOrderingModeCount = iODOrderingModeCount + 1
                if (iODOrderingModeCount == 1) then
                    dODOrderingEigenMin = dOrderingEigenMinLocal
                else
                    dODOrderingEigenMin = DMIN1(dODOrderingEigenMin, dOrderingEigenMinLocal)
                end if
            end if
        end if

        iCompanionPhase = iHelperPhase
        call FindActiveSolutionSlot(iHelperPhase, iHelperSlot)
        if (iHelperSlot <= 0) then
            call FindActiveAliasCompanionSlot(iHelperPhase, iOrderedPhase, iCompanionPhase, iHelperSlot)
        end if
        if ((iOrderedSlot <= 0).OR.(iHelperSlot <= 0)) cycle

        iODPairCount = iODPairCount + 1
        call ComputeNormalizedSolutionComposition(iOrderedPhase, dOrderedComposition)
        call ComputeNormalizedSolutionComposition(iCompanionPhase, dHelperComposition)
        dCompositionDistance = SUM(DABS(dOrderedComposition - dHelperComposition)) / &
            DFLOAT(MAX(1,nElements))
        dODCompDist = &
            DMAX1(dODCompDist, dCompositionDistance)
    end do

    return

contains

    subroutine FindActiveSolutionSlot(iSolnPhaseIndex, iSlotOut)
        integer, intent(in)  :: iSolnPhaseIndex
        integer, intent(out) :: iSlotOut

        integer :: iSoln, iSlot

        iSlotOut = 0
        do iSoln = 1, nSolnPhases
            iSlot = nElements - iSoln + 1
            if (iAssemblage(iSlot) == -iSolnPhaseIndex) then
                iSlotOut = iSlot
                return
            end if
        end do

        return
    end subroutine FindActiveSolutionSlot


    subroutine FindActiveAliasCompanionSlot(iHelperPhaseIn, iOrderedPhaseIn, iCompanionPhaseOut, iSlotOut)
        integer, intent(in)  :: iHelperPhaseIn, iOrderedPhaseIn
        integer, intent(out) :: iCompanionPhaseOut, iSlotOut

        integer :: iSoln, iSlot, iCandidatePhase

        iCompanionPhaseOut = iHelperPhaseIn
        iSlotOut = 0

        do iSoln = 1, nSolnPhases
            iSlot = nElements - iSoln + 1
            iCandidatePhase = -iAssemblage(iSlot)
            if (iCandidatePhase <= 0) cycle
            if (iCandidatePhase == iOrderedPhaseIn) cycle
            if (OrderDisorderHelperAliasMatch(iHelperPhaseIn, iCandidatePhase)) then
                iCompanionPhaseOut = iCandidatePhase
                iSlotOut = iSlot
                return
            end if
        end do

        return
    end subroutine FindActiveAliasCompanionSlot


    logical function OrderDisorderHelperAliasMatch(iHelperPhaseIn, iCandidatePhaseIn)
        integer, intent(in) :: iHelperPhaseIn, iCandidatePhaseIn

        integer       :: iHelperKind, iCandidateKind, iOrderedPhase
        integer       :: OrderDisorderHelperAliasKind
        character(30) :: cHelperName, cCandidateName

        OrderDisorderHelperAliasMatch = .FALSE.
        if ((iHelperPhaseIn <= 0).OR.(iHelperPhaseIn > nSolnPhasesSys)) return
        if ((iCandidatePhaseIn <= 0).OR.(iCandidatePhaseIn > nSolnPhasesSys)) return

        cHelperName = UpperPhaseName(cSolnPhaseName(iHelperPhaseIn))
        cCandidateName = UpperPhaseName(cSolnPhaseName(iCandidatePhaseIn))

        if (TRIM(cHelperName) == TRIM(cCandidateName)) then
            OrderDisorderHelperAliasMatch = .TRUE.
            return
        end if

        if (allocated(iODTopologyClass).AND.allocated(iDisorderedPhase).AND.&
            allocated(iODStandalonePhase)) then
            do iOrderedPhase = 1, MIN(nSolnPhasesSys, SIZE(iODTopologyClass), &
                SIZE(iDisorderedPhase), SIZE(iODStandalonePhase))
                if (iODTopologyClass(iOrderedPhase) /= OD_TOPOLOGY_HELPER_STANDALONE) cycle
                if (((iDisorderedPhase(iOrderedPhase) == iHelperPhaseIn).AND.&
                    (iODStandalonePhase(iOrderedPhase) == iCandidatePhaseIn)).OR.&
                    ((iDisorderedPhase(iOrderedPhase) == iCandidatePhaseIn).AND.&
                    (iODStandalonePhase(iOrderedPhase) == iHelperPhaseIn))) then
                    OrderDisorderHelperAliasMatch = .TRUE.
                    return
                end if
            end do

            ! Typed TDB identity is authoritative.  The spelling fallback
            ! below exists only for metadata-free legacy databases.
            if (ANY((iODTopologyClass > 0).AND.&
                (iODTopologyClass /= OD_TOPOLOGY_LEGACY_NO_METADATA))) return
        end if

        iHelperKind = OrderDisorderHelperAliasKind(iHelperPhaseIn)
        iCandidateKind = OrderDisorderHelperAliasKind(iCandidatePhaseIn)

        if ((iHelperKind == 1).AND.(TRIM(cCandidateName) == 'BCC_A2')) then
            OrderDisorderHelperAliasMatch = .TRUE.
        else if ((iCandidateKind == 1).AND.(TRIM(cHelperName) == 'BCC_A2')) then
            OrderDisorderHelperAliasMatch = .TRUE.
        else if ((iHelperKind == 2).AND.(TRIM(cCandidateName) == 'FCC_A1')) then
            OrderDisorderHelperAliasMatch = .TRUE.
        else if ((iCandidateKind == 2).AND.(TRIM(cHelperName) == 'FCC_A1')) then
            OrderDisorderHelperAliasMatch = .TRUE.
        end if

        return
    end function OrderDisorderHelperAliasMatch


    subroutine ComputeNormalizedSolutionComposition(iSolnPhaseIndex, dCompositionOut)
        integer, intent(in) :: iSolnPhaseIndex
        real(8), dimension(nElements), intent(out) :: dCompositionOut

        integer :: iSpecies, iFirst, iLast, iElement
        real(8) :: dNorm

        dCompositionOut = 0D0
        iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
        iLast  = nSpeciesPhase(iSolnPhaseIndex)
        do iSpecies = iFirst, iLast
            do iElement = 1, nElements
                dCompositionOut(iElement) = dCompositionOut(iElement) + &
                    DMAX1(dMolFraction(iSpecies), 0D0) * dStoichSpecies(iSpecies,iElement) / &
                    DFLOAT(iParticlesPerMole(iSpecies))
            end do
        end do

        dNorm = SUM(DABS(dCompositionOut))
        if (dNorm > 1D-300) dCompositionOut = dCompositionOut / dNorm

        return
    end subroutine ComputeNormalizedSolutionComposition


    real(8) function OrderedPhaseOrderingNorm(iOrderedPhase)
        integer, intent(in) :: iOrderedPhase

        integer :: iSublPhase, iSub, iCon, iElement, nTerms
        real(8) :: dSite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), dimension(nElements) :: dComposition

        OrderedPhaseOrderingNorm = 0D0
        iSublPhase = iPhaseSublattice(iOrderedPhase)
        if (iSublPhase <= 0) return

        call ComputeNormalizedSolutionComposition(iOrderedPhase, dComposition)
        call BuildSolutionSiteFractions(iOrderedPhase, iSublPhase, dSite)

        nTerms = 0
        do iSub = 1, nSublatticePhase(iSublPhase)
            do iCon = 1, nConstituentSublattice(iSublPhase,iSub)
                if (IsVacancyConstituent(cConstituentNameSUB(iSublPhase,iSub,iCon))) cycle
                iElement = ConstituentElementIndex(cConstituentNameSUB(iSublPhase,iSub,iCon))
                if (iElement <= 0) cycle
                OrderedPhaseOrderingNorm = OrderedPhaseOrderingNorm + &
                    DABS(dSite(iSub,iCon) - dComposition(iElement))
                nTerms = nTerms + 1
            end do
        end do

        if (nTerms > 0) OrderedPhaseOrderingNorm = OrderedPhaseOrderingNorm / DFLOAT(nTerms)

        return
    end function OrderedPhaseOrderingNorm


    subroutine BuildSolutionSiteFractions(iSolnPhaseIndex, iSublPhase, dSiteOut)
        integer, intent(in) :: iSolnPhaseIndex, iSublPhase
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys), intent(out) :: dSiteOut

        integer :: iSpecies, iFirst, iLast, iLocal, iSub, iCon
        real(8) :: dSubSum

        dSiteOut = 0D0
        iFirst = nSpeciesPhase(iSolnPhaseIndex-1) + 1
        iLast  = nSpeciesPhase(iSolnPhaseIndex)
        do iSpecies = iFirst, iLast
            iLocal = iSpecies - iFirst + 1
            do iSub = 1, nSublatticePhase(iSublPhase)
                iCon = iConstituentSublattice(iSublPhase,iSub,iLocal)
                if (iCon > 0) dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                    DMAX1(dMolFraction(iSpecies), 0D0)
            end do
        end do

        do iSub = 1, nSublatticePhase(iSublPhase)
            dSubSum = SUM(dSiteOut(iSub,1:nConstituentSublattice(iSublPhase,iSub)))
            if (dSubSum > 1D-300) then
                dSiteOut(iSub,1:nConstituentSublattice(iSublPhase,iSub)) = &
                    dSiteOut(iSub,1:nConstituentSublattice(iSublPhase,iSub)) / dSubSum
            end if
        end do

        return
    end subroutine BuildSolutionSiteFractions


    integer function ConstituentElementIndex(cName)
        character(*), intent(in) :: cName

        integer :: iElement
        character(8) :: cTarget, cElement

        ConstituentElementIndex = 0
        cTarget = UpperName(cName)
        do iElement = 1, nElements
            cElement = UpperName(cElementName(iElement))
            if (TRIM(cTarget) == TRIM(cElement)) then
                ConstituentElementIndex = iElement
                return
            end if
        end do

        return
    end function ConstituentElementIndex


    logical function IsVacancyConstituent(cName)
        character(*), intent(in) :: cName

        IsVacancyConstituent = TRIM(UpperName(cName)) == 'VA'

        return
    end function IsVacancyConstituent


    character(30) function UpperPhaseName(cName)
        character(*), intent(in) :: cName

        integer :: iChar, iCode, nChar

        UpperPhaseName = ' '
        nChar = MIN(LEN_TRIM(cName), LEN(UpperPhaseName))
        do iChar = 1, nChar
            iCode = IACHAR(cName(iChar:iChar))
            if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                UpperPhaseName(iChar:iChar) = ACHAR(iCode - 32)
            else
                UpperPhaseName(iChar:iChar) = cName(iChar:iChar)
            end if
        end do

        return
    end function UpperPhaseName


    character(8) function UpperName(cName)
        character(*), intent(in) :: cName

        integer :: iChar, iCode, nChar

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

        return
    end function UpperName

end subroutine ComputeODDiagnostics
