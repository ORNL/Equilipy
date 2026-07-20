!> \brief Compute the passive projected ordering-mode curvature for a SUBOM phase.
!!
!! \details Evaluates an ordered SUBOM parent both at its current minimum and
!! on its exact disordered projection, then projects the site-fraction
!! gradient/Hessian onto ordering modes that change equivalent sublattices at
!! fixed averaged composition.  All thermodynamic state touched during
!! evaluation is restored before return.
!!
subroutine CompOrderingModeSUBOM(iSolnIndex, nModeDim, dOrderingHessian, &
    dOrderingGradient, dOrderingEigenMin, nModeOut, iInfo)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompOrderingModeSUBOM.f90
    !> \brief   Passive projected ordering-mode curvature diagnostic for SUBOM phases.
    !> \author  S.Y. Kwon
    !> \date    Jul. 03, 2026
    !> \sa      CompGradientSUBL.f90
    !> \sa      CompHessianSUBL.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Classified order/disorder candidates from current, projected, and site-Hessian thermodynamics without site-split tolerances.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details This routine measures whether a SUBOM phase has a local
    !! ordering instability at the random/disordered manifold without changing
    !! active-set behavior.  It first averages equivalent ordered
    !! substitutional sublattices, then builds ordering coordinates that
    !! preserve that averaged composition.  The returned Hessian is
    !! \f$Q^T K Q\f$, where \f$K\f$ is the redundant site Hessian from
    !! `CompHessianSUBL` and \f$Q\f$ is an orthonormal ordering-mode basis.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    integer, intent(in)  :: iSolnIndex, nModeDim
    real(8), intent(out) :: dOrderingHessian(nModeDim,nModeDim)
    real(8), intent(out) :: dOrderingGradient(nModeDim)
    real(8), intent(out) :: dOrderingEigenMin
    integer, intent(out) :: nModeOut, iInfo

    integer :: iPhaseID, nSublattice, nSiteCapacity, nSpeciesCapacity
    integer :: nSiteOut, nSpeciesOut, nGradientSiteOut, nModeLocal
    integer :: nGroup, nActiveCon, iGradientInfo, iHessianInfo
    integer :: iSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
    integer :: iGroupSub(nMaxSublatticeSys), iActiveCon(nMaxConstituentSys)
    integer, allocatable :: iSiteSub(:), iSiteCon(:)
    real(8) :: dScalarGibbs, dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity
    real(8) :: dCurrentScalarGibbs, dDisorderedScalarGibbs
    real(8) :: dCurrentSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8) :: dDisorderedSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8), allocatable :: dSiteGradient(:), dSiteGradientH(:), dSiteGradientS(:)
    real(8), allocatable :: dSiteGradientCp(:), dSiteHessian(:,:)
    real(8), allocatable :: dAmountHessian(:,:), dCompositionJacobian(:,:)
    real(8), allocatable :: dBasis(:,:)
    real(8), allocatable :: dMolFractionSave(:), dMolesSpeciesSave(:)
    real(8), allocatable :: dGibbsSolnPhaseSave(:), dSiteFractionSave(:,:,:)
    logical :: lStateSaved, lSiteOnlySave

    external DSYEV

    dOrderingHessian = 0D0
    dOrderingGradient = 0D0
    dOrderingEigenMin = 0D0
    nModeOut = 0
    iInfo = 0
    lStateSaved = .FALSE.
    dCurrentScalarGibbs = 0D0
    dDisorderedScalarGibbs = 0D0

    if (lODCandidateClassifyActive.AND.allocated(iODCandidateClass)) then
        if ((iSolnIndex >= 1).AND.(iSolnIndex <= SIZE(iODCandidateClass))) then
            iODCandidateClass(iSolnIndex) = OD_CANDIDATE_EVALUATION_FAILED
            dODCandidateCurrentGibbs(iSolnIndex) = 0D0
            dODCandidateDisorderedGibbs(iSolnIndex) = 0D0
            dODCandidateOrderingEigenMin(iSolnIndex) = 0D0
            dODCompanionEigenMin(iSolnIndex) = 0D0
        end if
    end if

    if (nModeDim <= 0) then
        iInfo = 3
        return
    end if
    if ((iSolnIndex <= 0).OR.(iSolnIndex > nSolnPhasesSys)) then
        iInfo = 1
        return
    end if
    if (TRIM(cSolnPhaseType(iSolnIndex)) /= 'SUBOM') then
        iInfo = 1
        return
    end if

    iPhaseID = iPhaseSublattice(iSolnIndex)
    if (iPhaseID <= 0) then
        iInfo = 2
        return
    end if

    nSublattice = nSublatticePhase(iPhaseID)
    nSiteCapacity = MAX(1, nMaxSublatticeSys*nMaxConstituentSys)
    nSpeciesCapacity = MAX(1, nSpeciesPhase(iSolnIndex) - nSpeciesPhase(iSolnIndex-1))
    allocate(iSiteSub(nSiteCapacity), iSiteCon(nSiteCapacity))
    allocate(dSiteGradient(nSiteCapacity), dSiteGradientH(nSiteCapacity))
    allocate(dSiteGradientS(nSiteCapacity), dSiteGradientCp(nSiteCapacity))
    allocate(dSiteHessian(nSiteCapacity,nSiteCapacity))
    allocate(dAmountHessian(nSpeciesCapacity,nSpeciesCapacity))
    allocate(dCompositionJacobian(nSpeciesCapacity,nSpeciesCapacity))
    allocate(dBasis(nSiteCapacity,nModeDim))

    call BuildSiteIndex(iPhaseID, nSiteCapacity, iSiteIndex, iSiteSub, iSiteCon, nSiteOut)
    if (nSiteOut <= 0) then
        iInfo = 2
        return
    end if

    call BuildNormalizedSiteFractions(iSolnIndex, iPhaseID, dCurrentSite)
    call FindEquivalentOrderingGroup(iPhaseID, dCurrentSite, iGroupSub, nGroup, &
        iActiveCon, nActiveCon, iInfo)
    if (iInfo /= 0) return

    call BuildDisorderedSite(iPhaseID, dCurrentSite, iGroupSub, nGroup, &
        iActiveCon, nActiveCon, dDisorderedSite)
    call BuildOrderingBasis(iSiteIndex, iGroupSub, nGroup, iActiveCon, nActiveCon, &
        nSiteCapacity, nModeDim, dBasis, nModeLocal, iInfo)
    if (iInfo /= 0) return
    if (nModeLocal <= 0) then
        iInfo = 4
        return
    end if

    call SaveThermoState(dMolFractionSave, dMolesSpeciesSave, dGibbsSolnPhaseSave, &
        dSiteFractionSave, lStateSaved)
    call CompGradientSUBL(iSolnIndex, nSiteCapacity, dSiteGradient, dCurrentScalarGibbs, &
        dSiteGradientH, dSiteGradientS, dSiteGradientCp, dScalarEnthalpy, dScalarEntropy, &
        dScalarHeatCapacity, nGradientSiteOut, iGradientInfo)
    if ((iGradientInfo /= 0).OR.(nGradientSiteOut /= nSiteOut)) then
        iInfo = 5
        call RestoreThermoState(dMolFractionSave, dMolesSpeciesSave, dGibbsSolnPhaseSave, &
            dSiteFractionSave, lStateSaved)
        return
    end if

    call SetPhaseFractionsFromSite(iSolnIndex, iPhaseID, dDisorderedSite, iInfo)
    if (iInfo /= 0) then
        call RestoreThermoState(dMolFractionSave, dMolesSpeciesSave, dGibbsSolnPhaseSave, &
            dSiteFractionSave, lStateSaved)
        return
    end if

    call CompGradientSUBL(iSolnIndex, nSiteCapacity, dSiteGradient, dScalarGibbs, &
        dSiteGradientH, dSiteGradientS, dSiteGradientCp, dScalarEnthalpy, dScalarEntropy, &
        dScalarHeatCapacity, nGradientSiteOut, iGradientInfo)
    if ((iGradientInfo /= 0).OR.(nGradientSiteOut /= nSiteOut)) then
        iInfo = 5
        call RestoreThermoState(dMolFractionSave, dMolesSpeciesSave, dGibbsSolnPhaseSave, &
            dSiteFractionSave, lStateSaved)
        return
    end if
    dDisorderedScalarGibbs = dScalarGibbs

    lSiteOnlySave = lSUBLHessianSiteOnlyActive
    lSUBLHessianSiteOnlyActive = .TRUE.
    call CompHessianSUBL(iSolnIndex, nSiteCapacity, nSpeciesCapacity, dSiteHessian, &
        dAmountHessian, dCompositionJacobian, nSiteOut, nSpeciesOut, iHessianInfo)
    lSUBLHessianSiteOnlyActive = lSiteOnlySave
    call RestoreThermoState(dMolFractionSave, dMolesSpeciesSave, dGibbsSolnPhaseSave, &
        dSiteFractionSave, lStateSaved)
    if (iHessianInfo /= 0) then
        iInfo = 6
        return
    end if

    call ProjectOrderingDerivatives(nSiteCapacity, nModeDim, nModeLocal, dBasis, &
        dSiteGradient, dSiteHessian, dOrderingGradient, dOrderingHessian)
    call ComputeMinimumEigenvalue(nModeDim, nModeLocal, dOrderingHessian, &
        dOrderingEigenMin, iInfo)
    if (iInfo == 0) then
        nModeOut = nModeLocal
        if (lODCandidateClassifyActive) call RecordCandidateIdentity
    end if

    return

contains

    subroutine BuildSiteIndex(iPhaseIDIn, nSiteCapacityIn, iSiteIndexOut, &
        iSiteSubOut, iSiteConOut, nSiteLocal)
        integer, intent(in) :: iPhaseIDIn, nSiteCapacityIn
        integer, intent(out) :: iSiteIndexOut(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(out) :: iSiteSubOut(nSiteCapacityIn), iSiteConOut(nSiteCapacityIn)
        integer, intent(out) :: nSiteLocal

        integer :: iSub, iCon

        iSiteIndexOut = 0
        iSiteSubOut = 0
        iSiteConOut = 0
        nSiteLocal = 0
        do iSub = 1, nSublatticePhase(iPhaseIDIn)
            do iCon = 1, nConstituentSublattice(iPhaseIDIn,iSub)
                nSiteLocal = nSiteLocal + 1
                if (nSiteLocal <= nSiteCapacityIn) then
                    iSiteIndexOut(iSub,iCon) = nSiteLocal
                    iSiteSubOut(nSiteLocal) = iSub
                    iSiteConOut(nSiteLocal) = iCon
                end if
            end do
        end do

        return
    end subroutine BuildSiteIndex


    subroutine BuildNormalizedSiteFractions(iSolnIndexIn, iPhaseIDIn, dSiteOut)
        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iFirstLocal, iLastLocal, iSpecies, iLocal, iSub, iCon
        real(8) :: dNorm

        dSiteOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnIndexIn)
        do iSpecies = iFirstLocal, iLastLocal
            iLocal = iSpecies - iFirstLocal + 1
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) then
                    dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                        DMAX1(dMolFraction(iSpecies), 0D0)
                end if
            end do
        end do

        do iSub = 1, nSublatticePhase(iPhaseIDIn)
            dNorm = SUM(dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)))
            if (dNorm > 1D-300) then
                dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)) = &
                    dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)) / dNorm
            end if
        end do

        return
    end subroutine BuildNormalizedSiteFractions


    subroutine FindEquivalentOrderingGroup(iPhaseIDIn, dSiteIn, iGroupSubOut, nGroupOut, &
        iActiveConOut, nActiveConOut, iInfoOut)
        integer, intent(in) :: iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(out) :: iGroupSubOut(nMaxSublatticeSys), iActiveConOut(nMaxConstituentSys)
        integer, intent(out) :: nGroupOut, nActiveConOut, iInfoOut

        integer :: iSub, jSub, nTrialGroup, nTrialActive
        integer :: iTrialGroup(nMaxSublatticeSys), iTrialActive(nMaxConstituentSys)

        iGroupSubOut = 0
        iActiveConOut = 0
        nGroupOut = 0
        nActiveConOut = 0
        iInfoOut = 4

        do iSub = 1, nSublatticePhase(iPhaseIDIn)
            if (nConstituentSublattice(iPhaseIDIn,iSub) <= 1) cycle
            nTrialGroup = 1
            iTrialGroup(1) = iSub
            do jSub = iSub + 1, nSublatticePhase(iPhaseIDIn)
                if (SameConstituentList(iPhaseIDIn, iSub, jSub)) then
                    nTrialGroup = nTrialGroup + 1
                    iTrialGroup(nTrialGroup) = jSub
                end if
            end do
            if (nTrialGroup < 2) cycle

            call StructuralOrderingConstituents(iPhaseIDIn, iTrialGroup, &
                nTrialGroup, iTrialActive, nTrialActive)
            if (nTrialActive < 2) cycle

            iGroupSubOut(1:nTrialGroup) = iTrialGroup(1:nTrialGroup)
            iActiveConOut(1:nTrialActive) = iTrialActive(1:nTrialActive)
            nGroupOut = nTrialGroup
            nActiveConOut = nTrialActive
            iInfoOut = 0
            return
        end do

        return
    end subroutine FindEquivalentOrderingGroup


    logical function SameConstituentList(iPhaseIDIn, iSubA, iSubB)
        integer, intent(in) :: iPhaseIDIn, iSubA, iSubB

        integer :: iCon

        SameConstituentList = .FALSE.
        if (nConstituentSublattice(iPhaseIDIn,iSubA) /= &
            nConstituentSublattice(iPhaseIDIn,iSubB)) return

        do iCon = 1, nConstituentSublattice(iPhaseIDIn,iSubA)
            if (TRIM(UpperName(cConstituentNameSUB(iPhaseIDIn,iSubA,iCon))) /= &
                TRIM(UpperName(cConstituentNameSUB(iPhaseIDIn,iSubB,iCon)))) return
        end do

        SameConstituentList = .TRUE.
        return
    end function SameConstituentList


    subroutine StructuralOrderingConstituents(iPhaseIDIn, iGroupSubIn, nGroupIn, &
        iActiveConOut, nActiveConOut)
        integer, intent(in) :: iPhaseIDIn, nGroupIn
        integer, intent(in) :: iGroupSubIn(nMaxSublatticeSys)
        integer, intent(out) :: iActiveConOut(nMaxConstituentSys), nActiveConOut

        integer :: iCon

        iActiveConOut = 0
        nActiveConOut = 0
        do iCon = 1, nConstituentSublattice(iPhaseIDIn,iGroupSubIn(1))
            nActiveConOut = nActiveConOut + 1
            iActiveConOut(nActiveConOut) = iCon
        end do

        return
    end subroutine StructuralOrderingConstituents


    subroutine BuildDisorderedSite(iPhaseIDIn, dSiteIn, iGroupSubIn, nGroupIn, &
        iActiveConIn, nActiveConIn, dSiteOut)
        integer, intent(in) :: iPhaseIDIn, nGroupIn, nActiveConIn
        integer, intent(in) :: iGroupSubIn(nMaxSublatticeSys), iActiveConIn(nMaxConstituentSys)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iGroup, iActive, iCon, iSub
        real(8) :: dAverage(nMaxConstituentSys), dNorm

        dSiteOut = dSiteIn
        dAverage = 0D0
        do iActive = 1, nActiveConIn
            iCon = iActiveConIn(iActive)
            do iGroup = 1, nGroupIn
                dAverage(iCon) = dAverage(iCon) + dSiteIn(iGroupSubIn(iGroup),iCon)
            end do
            dAverage(iCon) = dAverage(iCon) / DFLOAT(nGroupIn)
        end do
        dNorm = SUM(dAverage(iActiveConIn(1:nActiveConIn)))
        if (dNorm <= 1D-300) dNorm = 1D0

        do iGroup = 1, nGroupIn
            iSub = iGroupSubIn(iGroup)
            dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)) = 0D0
            do iActive = 1, nActiveConIn
                iCon = iActiveConIn(iActive)
                dSiteOut(iSub,iCon) = dAverage(iCon) / dNorm
            end do
        end do

        return
    end subroutine BuildDisorderedSite


    subroutine BuildOrderingBasis(iSiteIndexIn, iGroupSubIn, nGroupIn, iActiveConIn, &
        nActiveConIn, nSiteCapacityIn, nModeCapacityIn, dBasisOut, nModeLocal, iInfoOut)
        integer, intent(in) :: nGroupIn, nActiveConIn, nSiteCapacityIn, nModeCapacityIn
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iGroupSubIn(nMaxSublatticeSys), iActiveConIn(nMaxConstituentSys)
        real(8), intent(out) :: dBasisOut(nSiteCapacityIn,nModeCapacityIn)
        integer, intent(out) :: nModeLocal, iInfoOut

        integer :: iGroup, iActive, iSub, iRefSub, iCon, iRefCon
        integer :: iSiteA, iSiteB, iSiteC, iSiteD
        real(8) :: dTrial(nSiteCapacityIn)

        dBasisOut = 0D0
        nModeLocal = 0
        iInfoOut = 0
        iRefSub = iGroupSubIn(1)
        iRefCon = iActiveConIn(1)

        do iGroup = 2, nGroupIn
            iSub = iGroupSubIn(iGroup)
            do iActive = 2, nActiveConIn
                iCon = iActiveConIn(iActive)
                iSiteA = iSiteIndexIn(iSub,iCon)
                iSiteB = iSiteIndexIn(iSub,iRefCon)
                iSiteC = iSiteIndexIn(iRefSub,iCon)
                iSiteD = iSiteIndexIn(iRefSub,iRefCon)
                if (MIN(iSiteA, iSiteB, iSiteC, iSiteD) <= 0) cycle
                dTrial = 0D0
                dTrial(iSiteA) = 1D0
                dTrial(iSiteB) = -1D0
                dTrial(iSiteC) = -1D0
                dTrial(iSiteD) = 1D0
                call AddOrthonormalMode(nSiteCapacityIn, nModeCapacityIn, dTrial, &
                    dBasisOut, nModeLocal, iInfoOut)
                if (iInfoOut /= 0) return
            end do
        end do

        return
    end subroutine BuildOrderingBasis


    subroutine AddOrthonormalMode(nSiteCapacityIn, nModeCapacityIn, dTrialInOut, &
        dBasisInOut, nModeLocal, iInfoOut)
        integer, intent(in) :: nSiteCapacityIn, nModeCapacityIn
        real(8), intent(inout) :: dTrialInOut(nSiteCapacityIn)
        real(8), intent(inout) :: dBasisInOut(nSiteCapacityIn,nModeCapacityIn)
        integer, intent(inout) :: nModeLocal
        integer, intent(out) :: iInfoOut

        integer :: iMode
        real(8) :: dProjection, dNorm

        iInfoOut = 0
        do iMode = 1, nModeLocal
            dProjection = SUM(dTrialInOut * dBasisInOut(:,iMode))
            dTrialInOut = dTrialInOut - dProjection * dBasisInOut(:,iMode)
        end do

        dNorm = DSQRT(SUM(dTrialInOut*dTrialInOut))
        if (dNorm == 0D0) return
        if (nModeLocal >= nModeCapacityIn) then
            iInfoOut = 3
            return
        end if
        nModeLocal = nModeLocal + 1
        dBasisInOut(:,nModeLocal) = dTrialInOut / dNorm

        return
    end subroutine AddOrthonormalMode


    subroutine SaveThermoState(dMolFractionSaveOut, dMolesSpeciesSaveOut, &
        dGibbsSolnPhaseSaveOut, dSiteFractionSaveOut, lStateSavedOut)
        real(8), allocatable, intent(out) :: dMolFractionSaveOut(:), dMolesSpeciesSaveOut(:)
        real(8), allocatable, intent(out) :: dGibbsSolnPhaseSaveOut(:)
        real(8), allocatable, intent(out) :: dSiteFractionSaveOut(:,:,:)
        logical, intent(out) :: lStateSavedOut

        allocate(dMolFractionSaveOut(SIZE(dMolFraction)))
        allocate(dMolesSpeciesSaveOut(SIZE(dMolesSpecies)))
        allocate(dGibbsSolnPhaseSaveOut(SIZE(dGibbsSolnPhase)))
        allocate(dSiteFractionSaveOut(SIZE(dSiteFraction,1), SIZE(dSiteFraction,2), &
            SIZE(dSiteFraction,3)))
        dMolFractionSaveOut = dMolFraction
        dMolesSpeciesSaveOut = dMolesSpecies
        dGibbsSolnPhaseSaveOut = dGibbsSolnPhase
        dSiteFractionSaveOut = dSiteFraction
        lStateSavedOut = .TRUE.

        return
    end subroutine SaveThermoState


    subroutine RestoreThermoState(dMolFractionSaveIn, dMolesSpeciesSaveIn, &
        dGibbsSolnPhaseSaveIn, dSiteFractionSaveIn, lStateSavedIn)
        real(8), allocatable, intent(in) :: dMolFractionSaveIn(:), dMolesSpeciesSaveIn(:)
        real(8), allocatable, intent(in) :: dGibbsSolnPhaseSaveIn(:)
        real(8), allocatable, intent(in) :: dSiteFractionSaveIn(:,:,:)
        logical, intent(in) :: lStateSavedIn

        if (.NOT.lStateSavedIn) return
        dMolFraction = dMolFractionSaveIn
        dMolesSpecies = dMolesSpeciesSaveIn
        dGibbsSolnPhase = dGibbsSolnPhaseSaveIn
        dSiteFraction = dSiteFractionSaveIn

        return
    end subroutine RestoreThermoState


    subroutine SetPhaseFractionsFromSite(iSolnIndexIn, iPhaseIDIn, dSiteIn, iInfoOut)
        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(out) :: iInfoOut

        integer :: iFirstLocal, iLastLocal, iSpecies, iLocal, iSub, iCon
        real(8) :: dFractions(nSpeciesCapacity), dNorm

        iInfoOut = 0
        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnIndexIn)
        dFractions = 0D0
        do iSpecies = iFirstLocal, iLastLocal
            iLocal = iSpecies - iFirstLocal + 1
            dFractions(iLocal) = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) then
                    dFractions(iLocal) = dFractions(iLocal) * DMAX1(dSiteIn(iSub,iCon), 0D0)
                end if
            end do
        end do

        dNorm = SUM(dFractions(1:nSpeciesCapacity))
        if (dNorm <= 1D-300) then
            iInfoOut = 2
            return
        end if
        dFractions = dFractions / dNorm
        dMolFraction(iFirstLocal:iLastLocal) = dFractions(1:nSpeciesCapacity)
        dMolesSpecies(iFirstLocal:iLastLocal) = dFractions(1:nSpeciesCapacity)
        dGibbsSolnPhase(iSolnIndexIn) = 0D0
        dSiteFraction(iPhaseIDIn,:,:) = dSiteIn(:,:)

        return
    end subroutine SetPhaseFractionsFromSite


    subroutine ProjectOrderingDerivatives(nSiteCapacityIn, nModeCapacityIn, nModeLocal, &
        dBasisIn, dSiteGradientIn, dSiteHessianIn, dOrderingGradientOut, &
        dOrderingHessianOut)
        integer, intent(in) :: nSiteCapacityIn, nModeCapacityIn, nModeLocal
        real(8), intent(in) :: dBasisIn(nSiteCapacityIn,nModeCapacityIn)
        real(8), intent(in) :: dSiteGradientIn(nSiteCapacityIn)
        real(8), intent(in) :: dSiteHessianIn(nSiteCapacityIn,nSiteCapacityIn)
        real(8), intent(out) :: dOrderingGradientOut(nModeCapacityIn)
        real(8), intent(out) :: dOrderingHessianOut(nModeCapacityIn,nModeCapacityIn)

        integer :: iMode, jMode

        dOrderingGradientOut = 0D0
        dOrderingHessianOut = 0D0
        do iMode = 1, nModeLocal
            dOrderingGradientOut(iMode) = SUM(dBasisIn(:,iMode) * dSiteGradientIn(:))
            do jMode = 1, nModeLocal
                dOrderingHessianOut(iMode,jMode) = &
                    SUM(dBasisIn(:,iMode) * MATMUL(dSiteHessianIn, dBasisIn(:,jMode)))
            end do
        end do

        return
    end subroutine ProjectOrderingDerivatives


    subroutine ComputeMinimumEigenvalue(nModeCapacityIn, nModeLocal, dHessianIn, &
        dEigenMinOut, iInfoOut)
        integer, intent(in) :: nModeCapacityIn, nModeLocal
        real(8), intent(in) :: dHessianIn(nModeCapacityIn,nModeCapacityIn)
        real(8), intent(out) :: dEigenMinOut
        integer, intent(out) :: iInfoOut

        integer :: iLapackInfo, nWork
        real(8), allocatable :: dEigenMatrix(:,:), dEigen(:), dWork(:)

        iInfoOut = 0
        dEigenMinOut = 0D0
        if (nModeLocal <= 0) then
            iInfoOut = 4
            return
        end if

        allocate(dEigenMatrix(nModeLocal,nModeLocal))
        allocate(dEigen(nModeLocal))
        nWork = MAX(1, 3*nModeLocal - 1)
        allocate(dWork(nWork))
        dEigenMatrix = dHessianIn(1:nModeLocal,1:nModeLocal)
        call DSYEV('N', 'U', nModeLocal, dEigenMatrix, nModeLocal, dEigen, dWork, &
            nWork, iLapackInfo)
        if (iLapackInfo /= 0) then
            iInfoOut = 7
            return
        end if
        dEigenMinOut = dEigen(1)

        return
    end subroutine ComputeMinimumEigenvalue


    subroutine RecordCandidateIdentity

        if (.NOT.allocated(iODCandidateClass)) return
        if (iSolnIndex > SIZE(iODCandidateClass)) return

        dODCandidateCurrentGibbs(iSolnIndex) = dCurrentScalarGibbs
        dODCandidateDisorderedGibbs(iSolnIndex) = dDisorderedScalarGibbs
        dODCandidateOrderingEigenMin(iSolnIndex) = dOrderingEigenMin

        if ((.NOT.ieee_is_finite(dCurrentScalarGibbs)).OR.&
            (.NOT.ieee_is_finite(dDisorderedScalarGibbs)).OR.&
            (.NOT.ieee_is_finite(dOrderingEigenMin))) then
            iODCandidateClass(iSolnIndex) = OD_CANDIDATE_EVALUATION_FAILED
        else if (dOrderingEigenMin == 0D0) then
            ! Exact rank loss is a physical node fact.  It is neither an ordered
            ! nor a regular disordered minimum and must remain typed.
            iODCandidateClass(iSolnIndex) = OD_CANDIDATE_AMBIGUOUS_NODE
        else if (ALL(dCurrentSite == dDisorderedSite)) then
            if (dOrderingEigenMin < 0D0) then
                ! Exact symmetry with negative ordering curvature means the
                ! selected subminimum is a saddle, not a disordered minimum.
                iODCandidateClass(iSolnIndex) = OD_CANDIDATE_AMBIGUOUS_UNSTABLE
            else
                iODCandidateClass(iSolnIndex) = OD_CANDIDATE_DISORDERED
            end if
        else if (dCurrentScalarGibbs == dDisorderedScalarGibbs) then
            if (dOrderingEigenMin < 0D0) then
                iODCandidateClass(iSolnIndex) = OD_CANDIDATE_AMBIGUOUS_UNSTABLE
            else
                ! An asymmetric floating representation with exactly the same
                ! parent-model energy as its structural projection is the
                ! disordered branch.  This is an exact projection identity,
                ! not an ordering-degree cutoff.
                iODCandidateClass(iSolnIndex) = OD_CANDIDATE_DISORDERED
            end if
        else if (dCurrentScalarGibbs < dDisorderedScalarGibbs) then
            if ((dOrderingEigenMin > 0D0).AND.&
                WithinOneRepresentationStepBelow(dCurrentScalarGibbs, &
                    dDisorderedScalarGibbs)) then
                ! The two values overlap at floating-point representation
                ! precision.  NEAREST defines that arithmetic fact exactly;
                ! no thermodynamic or composition-distance cutoff is used.
                iODCandidateClass(iSolnIndex) = OD_CANDIDATE_AMBIGUOUS_ROUNDOFF
            else
                iODCandidateClass(iSolnIndex) = OD_CANDIDATE_ORDERED
            end if
        else
            if (dOrderingEigenMin > 0D0) then
                ! The exact symmetric projection is both lower and locally
                ! stable, so the asymmetric submin result cannot be the phase
                ! minimum.  Preserve that rejected-branch fact while selecting
                ! the structural disordered projection.
                iODCandidateClass(iSolnIndex) = OD_CANDIDATE_DISORDERED_PROJECTED
            else
                iODCandidateClass(iSolnIndex) = OD_CANDIDATE_AMBIGUOUS_UNSTABLE
            end if
        end if

        return
    end subroutine RecordCandidateIdentity


    logical function WithinOneRepresentationStepBelow(dValue, dReference)
        real(8), intent(in) :: dValue, dReference

        WithinOneRepresentationStepBelow = .FALSE.
        if (dValue >= dReference) return
        WithinOneRepresentationStepBelow = NEAREST(dValue, 1D0) >= dReference

        return
    end function WithinOneRepresentationStepBelow


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

    pure logical function ieee_is_finite(dValue)
        ! Portable finiteness test: the manylinux2014 aarch64 gfortran lacks
        ! the ieee_arithmetic intrinsic module. NaN fails dValue == dValue;
        ! the infinities exceed HUGE; HUGE itself is finite, hence <=.
        real(8), intent(in) :: dValue
        ieee_is_finite = (dValue == dValue) .AND. (ABS(dValue) <= HUGE(dValue))
    end function ieee_is_finite
!
end subroutine CompOrderingModeSUBOM
