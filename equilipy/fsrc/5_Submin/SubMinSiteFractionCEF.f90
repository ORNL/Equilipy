!-------------------------------------------------------------------------------------------------------------
!
!> \file    SubMinSiteFractionCEF.f90
!> \brief   Minimize CEF solution-phase driving force in site-fraction space.
!> \author  S.Y. Kwon
!> \date    Jun. 24, 2026
!> \sa      Subminimization.f90
!> \sa      CompGradientSUBL.f90
!> \sa      CompHessianSUBL.f90
!> \sa      SubMinDrivingForce.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   06/24/2026      S.Y. Kwon           Added CEF site-fraction Subminimization path
!   06/24/2026      S.Y. Kwon           Added charge-neutrality and reduced active-set solves
!   06/24/2026      S.Y. Kwon           Matched the site gradient to the line-search objective
!   06/24/2026      S.Y. Kwon           Switched CEF Subminimization to direct scalar site-gradient
!   06/24/2026      S.Y. Kwon           Moved the routine prologue above the subroutine declaration
!   06/24/2026      S.Y. Kwon           Recorded CEF path diagnostics for classic fallback max-out cases
!   06/24/2026      S.Y. Kwon           Clarified CEF handled-state diagnostics for Subminimization
!   06/28/2026      S.Y. Kwon           Added a residual-stagnation exit for the CEF site-fraction solve.
!   07/01/2026      S.Y. Kwon           Retried a projected site-gradient direction when a descent Newton
!                                       direction fails the CEF Subminimization line search.
!
!
! Purpose:
! ========
!
!> \details The endmember-fraction Subminimization Newton system becomes very
!! large for CEF phases such as B2 and FCC_4SL because every endmember
!! combination is an optimization variable.  This routine instead works in
!! independent site fractions on each sublattice.  The residual is the phase
!! grand-potential derivative
!! \f$\partial G/\partial y_{s,c} - \lambda \cdot \partial M/\partial y_{s,c}\f$.
!! For a one-sublattice CEF phase this becomes the ordinary exchange residual
!! between endmember mole fractions.  The scalar objective and first derivative
!! are evaluated directly in site-fraction coordinates by CompGradientSUBL.
!! The Hessian used for the Newton proposal is the analytical fixed-state CEF
!! Hessian from CompHessianSUBL projected to the independent site-fraction
!! basis.  Ionic CEF phases add a linear charge-neutrality row after
!! per-sublattice normalization has been eliminated by the reference variables.
!
!
! Required input variables:
! =========================
!
!> \param[in] iSolnPhaseIndex  Absolute solution phase index.
!
! cSolnPhaseType        Solution model type.  SUBL, SUBLM, and SUBOM are handled.
! iPhaseSublattice      Runtime CEF phase topology index.
! dMolFraction          Initial endmember fractions from Leveling or a previous PEA estimate.
! dElementPotential     Current elemental potentials used by dChemicalPotentialStar.
! dChemicalPotentialStar
!                       Element-potential projection for active endmembers, initialized by SubMinInit.
! dSublatticeCharge     Constituent charge used to construct the ionic CEF charge row.
!
!
! Output/updated variables:
! =========================
!
!> \param[out] lHandled  True when this routine converged tightly enough to pass the CEF submin result.
!
! dMolFraction          Updated endmember fractions corresponding to the accepted site fractions.
! dChemicalPotential    Updated by SubMinChemicalPotential at the accepted site state.
! dDrivingForce         Updated phase driving force at the accepted site state.
! dSubMinFunctionNorm   Site-normalization/charge residual guard for this path.
! dMaxPotentialVector   Site-gradient infinity norm for the CEF solve.
! lSubMinConverged      True when the site-gradient and step criteria are satisfied.
!
!
! Called subroutines/functions:
! =============================
!
! SubMinChemicalPotential      Evaluates phase chemical potentials from current endmember fractions.
! ComputeCEFSubMinGradient     Computes the active site-gradient from the direct CEF scalar gradient.
! EvaluateCEFSubMinObjective   Evaluates the phase grand potential and reporting driving force.
! CompGradientSUBL             Provides the direct fixed-state CEF scalar and site-gradient.
! CompHessianSUBL              Provides analytical CEF site Hessian for Newton proposals.
! DSYSV                        Solves the small symmetric site-fraction Newton system.
!
!
! Primary callers:
! ================
!
! Subminimization              Calls this before the legacy endmember Newton loop for CEF phases.
!
!
! Numerical assumptions:
! ======================
!
! - Per-sublattice normalization is eliminated by using one reference
!   constituent per sublattice.  Ionic CEF phases add one charge-neutrality
!   row/column to the projected Newton system.
! - The Newton merit is the unnormalized phase grand potential.  The legacy
!   per-real-atom dDrivingForce is updated only for PEA phase-entry reporting.
! - The active-set gradient is evaluated from the same fixed-state CEF scalar
!   used by the line search.
! - Constituent element projections assume ordinary CEF constituents whose
!   names map to active element names; VA/blank constituents do not contribute
!   to conserved element mass.
! - A small site-fraction floor is used only as a numerical interior-point
!   guard for log terms and line-search feasibility.
! - Lower-bound site constituents satisfy complementarity in the CEF
!   Subminimization solve.
! - A CEF submin start exits early when the active site-gradient residual no
!   longer changes within 1D-8.  This mirrors the classic Subminimization stall
!   exit and avoids spending the full CEF iteration budget on a flat path.
! - The analytical Hessian supplies the primary Newton direction.  If that
!   direction is a formal descent direction but fails the actual Armijo line
!   search, the same active site variables are retried with a projected
!   steepest-descent direction before the start is rejected.
!
!-------------------------------------------------------------------------------------------------------------



subroutine SubMinSiteFractionCEF(iSolnPhaseIndex, lHandled)
    USE ModuleThermo
    USE ModuleSubMin
    USE ModuleGEMSolver
    USE ModuleThermoIO

    implicit none

    integer, intent(in)  :: iSolnPhaseIndex
    logical, intent(out) :: lHandled

    integer, parameter :: iSiteIterMax = 80
    real(8), parameter :: dSiteFloor = 1D-12
    real(8), parameter :: dArmijo = 1D-4
    real(8), parameter :: dSiteGradientTolerance = 1D-3
    real(8), parameter :: dChargeMeritPenalty = 1D3
    real(8), parameter :: dSiteActiveSetFloor = 1D-6
    real(8), parameter :: dCEFGradientStagnationTol = 1D-8

    integer :: iPhaseID, nSiteCapacity, nSpeciesCapacity
    integer :: nIndependent, nActiveIndependent, iIter, iInfo
    real(8) :: dObjective, dObjectiveLast, dGradientNorm, dStepNorm, dConvergenceGradientTol
    real(8) :: dGradientNormLast
    real(8) :: dDirectionalDerivative, dMeritDirectionalDerivative
    real(8) :: dChargeResidual, dChargeTolerance
    logical :: lAccepted, lChargeConstraint
    logical :: lGradientStagnated
    logical :: lIsCEFPhase
    logical :: lUsedNewtonDirection

    integer, dimension(:), allocatable :: iVarSub, iVarCon, iVarRef
    integer, dimension(:), allocatable :: iVarSite, iVarRefSite
    integer, dimension(:), allocatable :: iActiveVar
    real(8), dimension(:), allocatable :: dGradient, dDirection, dChargeJacobian
    real(8), dimension(:,:), allocatable :: dSite, dSiteTrial, dIndependentHessian
    logical, dimension(:), allocatable :: lVarActive

    lHandled = .FALSE.

    lIsCEFPhase = (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBL').OR.&
        (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBLM').OR.&
        (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBOM')
    if (.NOT.lIsCEFPhase) return

    iPhaseID = iPhaseSublattice(iSolnPhaseIndex)
    if (iPhaseID <= 0) return
    lSubMinCEFAttempted = .TRUE.
    iSubMinCEFPhaseIndex = iSolnPhaseIndex
    lChargeConstraint = iPhaseElectronID(iSolnPhaseIndex) /= 0

    nSiteCapacity = nMaxSublatticeSys * nMaxConstituentSys
    nSpeciesCapacity = MAX(1, nVar)

    allocate(dSite(nMaxSublatticeSys,nMaxConstituentSys))
    allocate(dSiteTrial(nMaxSublatticeSys,nMaxConstituentSys))
    allocate(iVarSub(nSiteCapacity), iVarCon(nSiteCapacity), iVarRef(nSiteCapacity))
    allocate(iVarSite(nSiteCapacity), iVarRefSite(nSiteCapacity))
    allocate(iActiveVar(nSiteCapacity))
    allocate(dGradient(nSiteCapacity), dDirection(nSiteCapacity), dChargeJacobian(nSiteCapacity))
    allocate(dIndependentHessian(nSiteCapacity,nSiteCapacity))
    allocate(lVarActive(nSiteCapacity))

    call BuildCEFSubMinSite(iSolnPhaseIndex, iPhaseID, dSite)
    call NormalizeCEFSubMinSite(iPhaseID, dSite, dSiteFloor)
    call SetCEFSubMinMoleFractions(iSolnPhaseIndex, iPhaseID, dSite)
    call EvaluateCEFSubMinObjective(iSolnPhaseIndex, iPhaseID, dSite, dObjective)

    dObjectiveLast = dObjective
    dGradientNorm = HUGE(1D0)
    dStepNorm = HUGE(1D0)
    dConvergenceGradientTol = DMIN1(dMaxPotentialTol, dSiteGradientTolerance)
    dChargeTolerance = DMAX1(dTolerance(1), 1D-12)
    dChargeResidual = 0D0
    dChargeJacobian = 0D0
    dGradientNormLast = HUGE(1D0)
    lGradientStagnated = .FALSE.
    lSubMinConverged = .FALSE.

    LOOP_SiteIter: do iIter = 1, iSiteIterMax
        iterSub = iIter
        call BuildCEFSubMinVariables(iPhaseID, dSite, nSiteCapacity, nIndependent, &
            iVarSub, iVarCon, iVarRef, iVarSite, iVarRefSite)
        iSubMinCEFIter = iIter
        nSubMinCEFIndependent = nIndependent
        if (lChargeConstraint) then
            call ComputeCEFSubMinChargeConstraint(iPhaseID, dSite, nIndependent, &
                iVarSub, iVarCon, iVarRef, dChargeResidual, dChargeJacobian)
        else
            dChargeResidual = 0D0
            dChargeJacobian = 0D0
        end if
        dSubMinCEFChargeResidual = dChargeResidual

        if (nIndependent <= 0) then
            lSubMinConverged = (.NOT.lChargeConstraint).OR.(DABS(dChargeResidual) < dChargeTolerance)
            exit LOOP_SiteIter
        end if

        call ComputeCEFSubMinGradient(iSolnPhaseIndex, iPhaseID, dSite, nIndependent, &
            iVarSub, iVarCon, iVarRef, iVarSite, iVarRefSite, dGradient)
        call BuildCEFSubMinActiveSet(iPhaseID, dSite, nIndependent, iVarSub, iVarCon, &
            dGradient, dSiteActiveSetFloor, lVarActive, iActiveVar, &
            nActiveIndependent, dGradientNorm)
        nSubMinCEFActiveIndependent = nActiveIndependent
        if (lChargeConstraint) dGradientNorm = DMAX1(dGradientNorm, DABS(dChargeResidual))
        dMaxPotentialVector = dGradientNorm
        dSubMinCEFGradientNorm = dGradientNorm

        if ((dGradientNorm < dConvergenceGradientTol).AND.&
            ((.NOT.lChargeConstraint).OR.(DABS(dChargeResidual) < dChargeTolerance))) then
            lSubMinConverged = .TRUE.
            exit LOOP_SiteIter
        end if

        lGradientStagnated = (DABS(dGradientNorm - dGradientNormLast) < &
            dCEFGradientStagnationTol)
        dGradientNormLast = dGradientNorm
        if (lGradientStagnated) exit LOOP_SiteIter

        call ComputeCEFSubMinHessian(iSolnPhaseIndex, nSiteCapacity, nSpeciesCapacity, &
            nIndependent, iVarSite, iVarRefSite, dIndependentHessian, iInfo)

        lUsedNewtonDirection = .FALSE.
        if (iInfo == 0) then
            call SolveCEFSubMinActiveDirection(nIndependent, nActiveIndependent, iActiveVar, &
                dIndependentHessian, dGradient, lChargeConstraint, dChargeResidual, &
                dChargeJacobian, dDirection)
            lUsedNewtonDirection = .TRUE.
        else
            dDirection = 0D0
            call SetCEFSubMinActiveGradientDirection(nIndependent, lVarActive, dGradient, dDirection)
        end if

        call LimitCEFSubMinDirection(nIndependent, dDirection)
        dDirectionalDerivative = SUM(dGradient(1:nIndependent) * dDirection(1:nIndependent))
        dMeritDirectionalDerivative = dDirectionalDerivative
        if (lChargeConstraint) then
            dMeritDirectionalDerivative = dMeritDirectionalDerivative + &
                2D0 * dChargeMeritPenalty * dChargeResidual * &
                SUM(dChargeJacobian(1:nIndependent) * dDirection(1:nIndependent))
        end if

        if (dMeritDirectionalDerivative >= 0D0) then
            dDirection = 0D0
            call SetCEFSubMinActiveGradientDirection(nIndependent, lVarActive, dGradient, dDirection)
            call LimitCEFSubMinDirection(nIndependent, dDirection)
            lUsedNewtonDirection = .FALSE.
        end if

        call LineSearchCEFSubMin(iSolnPhaseIndex, iPhaseID, dSite, nIndependent, &
            iVarSub, iVarCon, iVarRef, dGradient, dDirection, dObjective, &
            dSiteFloor, dArmijo, lChargeConstraint, dChargeResidual, dChargeJacobian, &
            dChargeMeritPenalty, dSiteTrial, lAccepted, dStepNorm)

        if ((.NOT.lAccepted).AND.lUsedNewtonDirection) then
            dDirection = 0D0
            call SetCEFSubMinActiveGradientDirection(nIndependent, lVarActive, dGradient, dDirection)
            call LimitCEFSubMinDirection(nIndependent, dDirection)
            call LineSearchCEFSubMin(iSolnPhaseIndex, iPhaseID, dSite, nIndependent, &
                iVarSub, iVarCon, iVarRef, dGradient, dDirection, dObjective, &
                dSiteFloor, dArmijo, lChargeConstraint, dChargeResidual, dChargeJacobian, &
                dChargeMeritPenalty, dSiteTrial, lAccepted, dStepNorm)
        end if

        lSubMinCEFLineSearchAccepted = lAccepted
        dSubMinCEFStepNorm = dStepNorm

        if (.NOT.lAccepted) exit LOOP_SiteIter

        dSite = dSiteTrial
        dObjectiveLast = dObjective
        dSubMinCEFObjective = dObjective

        if ((dStepNorm < dSubMinTolerance).AND.(dGradientNorm < dConvergenceGradientTol)) then
            lSubMinConverged = .TRUE.
            exit LOOP_SiteIter
        end if
    end do LOOP_SiteIter

    call SetCEFSubMinMoleFractions(iSolnPhaseIndex, iPhaseID, dSite)
    call EvaluateCEFSubMinObjective(iSolnPhaseIndex, iPhaseID, dSite, dObjective)
    dSubMinFunctionNorm = 0D0
    if (lChargeConstraint) dSubMinFunctionNorm = DABS(dChargeResidual)
    dDrivingForceLast = dObjectiveLast
    dMaxPotentialVector = dGradientNorm
    dSubMinCEFObjective = dObjective
    dSubMinCEFGradientNorm = dGradientNorm
    dSubMinCEFStepNorm = dStepNorm
    dSubMinCEFChargeResidual = dChargeResidual

    lHandled = lSubMinConverged.AND.(dSubMinFunctionNorm < dTolerance(1))
    lSubMinCEFHandled = lHandled

    if (allocated(dSite)) deallocate(dSite)
    if (allocated(dSiteTrial)) deallocate(dSiteTrial)
    if (allocated(iVarSub)) deallocate(iVarSub)
    if (allocated(iVarCon)) deallocate(iVarCon)
    if (allocated(iVarRef)) deallocate(iVarRef)
    if (allocated(iVarSite)) deallocate(iVarSite)
    if (allocated(iVarRefSite)) deallocate(iVarRefSite)
    if (allocated(iActiveVar)) deallocate(iActiveVar)
    if (allocated(dGradient)) deallocate(dGradient)
    if (allocated(dDirection)) deallocate(dDirection)
    if (allocated(dChargeJacobian)) deallocate(dChargeJacobian)
    if (allocated(dIndependentHessian)) deallocate(dIndependentHessian)
    if (allocated(lVarActive)) deallocate(lVarActive)

    return

contains

    subroutine BuildCEFSubMinSite(iSolnPhaseIndexIn, iPhaseIDIn, dSiteOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: i, m, s, c, iFirstLocal, iLastLocal

        dSiteOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)

        do i = iFirstLocal, iLastLocal
            m = i - iFirstLocal + 1
            do s = 1, nSublatticePhase(iPhaseIDIn)
                c = iConstituentSublattice(iPhaseIDIn,s,m)
                if (c > 0) dSiteOut(s,c) = dSiteOut(s,c) + DMAX1(dMolFraction(i), 0D0)
            end do
        end do

        return

    end subroutine BuildCEFSubMinSite

    subroutine NormalizeCEFSubMinSite(iPhaseIDIn, dSiteInOut, dFloor)

        implicit none

        integer, intent(in) :: iPhaseIDIn
        real(8), intent(in) :: dFloor
        real(8), intent(inout) :: dSiteInOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: s, c
        real(8) :: dSum

        do s = 1, nSublatticePhase(iPhaseIDIn)
            dSum = 0D0
            do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                dSiteInOut(s,c) = DMAX1(dSiteInOut(s,c), dFloor)
                dSum = dSum + dSiteInOut(s,c)
            end do
            if (dSum > 0D0) then
                do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                    dSiteInOut(s,c) = dSiteInOut(s,c) / dSum
                end do
            end if
        end do

        return

    end subroutine NormalizeCEFSubMinSite

    subroutine SetCEFSubMinMoleFractions(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: i, m, s, c, iFirstLocal, iLastLocal
        real(8) :: dSum, dProduct

        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        dSum = 0D0

        do i = iFirstLocal, iLastLocal
            m = i - iFirstLocal + 1
            dProduct = 1D0
            do s = 1, nSublatticePhase(iPhaseIDIn)
                c = iConstituentSublattice(iPhaseIDIn,s,m)
                if (c > 0) dProduct = dProduct * DMAX1(dSiteIn(s,c), dMinMoleFraction)
            end do
            dMolFraction(i) = dProduct
            dSum = dSum + dProduct
        end do

        if (dSum > 0D0) then
            do i = iFirstLocal, iLastLocal
                dMolFraction(i) = dMolFraction(i) / dSum
            end do
        end if

        return

    end subroutine SetCEFSubMinMoleFractions

    subroutine EvaluateCEFSubMinObjective(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, dObjectiveOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dObjectiveOut

        integer :: i, s, c, iFirstLocal, iLastLocal, nSiteCapacityLocal, nSiteOut, iInfoLocal
        real(8) :: dAtomDenom, dPlane, dScalarGibbs
        real(8) :: dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity
        real(8), allocatable :: dSiteGradient(:), dSiteGradientH(:)
        real(8), allocatable :: dSiteGradientS(:), dSiteGradientCp(:)

        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        nSiteCapacityLocal = nMaxSublatticeSys * nMaxConstituentSys

        call SubMinChemicalPotential(iSolnPhaseIndexIn)
        allocate(dSiteGradient(nSiteCapacityLocal), dSiteGradientH(nSiteCapacityLocal), &
            dSiteGradientS(nSiteCapacityLocal), dSiteGradientCp(nSiteCapacityLocal))
        call CompGradientSUBL(iSolnPhaseIndexIn, nSiteCapacityLocal, dSiteGradient, &
            dScalarGibbs, dSiteGradientH, dSiteGradientS, dSiteGradientCp, &
            dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity, nSiteOut, iInfoLocal)
        if (iInfoLocal /= 0) then
            dObjectiveOut = 9D5
            dDrivingForce = 9D5
            if (allocated(dSiteGradient)) deallocate(dSiteGradient)
            if (allocated(dSiteGradientH)) deallocate(dSiteGradientH)
            if (allocated(dSiteGradientS)) deallocate(dSiteGradientS)
            if (allocated(dSiteGradientCp)) deallocate(dSiteGradientCp)
            return
        end if

        dObjectiveOut = dScalarGibbs
        dPlane = 0D0
        dAtomDenom = 0D0

        do i = iFirstLocal, iLastLocal
            dAtomDenom = dAtomDenom + dMolFraction(i) * &
                SUM(dStoichSpecies(i,1:nElements)) / DBLE(iParticlesPerMole(i))
        end do

        do s = 1, nSublatticePhase(iPhaseIDIn)
            do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                dPlane = dPlane + dSiteIn(s,c) * ComputeCEFSubMinElementProjection(iPhaseIDIn, s, c)
            end do
        end do

        dObjectiveOut = dObjectiveOut - dPlane
        if (DABS(dAtomDenom) > 1D-300) then
            dDrivingForce = dObjectiveOut / dAtomDenom
        else
            dDrivingForce = 9D5
        end if

        if (allocated(dSiteGradient)) deallocate(dSiteGradient)
        if (allocated(dSiteGradientH)) deallocate(dSiteGradientH)
        if (allocated(dSiteGradientS)) deallocate(dSiteGradientS)
        if (allocated(dSiteGradientCp)) deallocate(dSiteGradientCp)

        return

    end subroutine EvaluateCEFSubMinObjective

    subroutine BuildCEFSubMinVariables(iPhaseIDIn, dSiteIn, nSiteCapacityIn, nIndependentOut, &
        iVarSubOut, iVarConOut, iVarRefOut, iVarSiteOut, iVarRefSiteOut)

        implicit none

        integer, intent(in) :: iPhaseIDIn, nSiteCapacityIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(out) :: nIndependentOut
        integer, intent(out) :: iVarSubOut(nSiteCapacityIn), iVarConOut(nSiteCapacityIn)
        integer, intent(out) :: iVarRefOut(nSiteCapacityIn)
        integer, intent(out) :: iVarSiteOut(nSiteCapacityIn), iVarRefSiteOut(nSiteCapacityIn)

        integer :: s, c, iRefCon, iCursor
        real(8) :: dRefValue

        nIndependentOut = 0
        iVarSubOut = 0
        iVarConOut = 0
        iVarRefOut = 0
        iVarSiteOut = 0
        iVarRefSiteOut = 0
        iCursor = 0

        do s = 1, nSublatticePhase(iPhaseIDIn)
            iRefCon = 1
            dRefValue = dSiteIn(s,1)
            do c = 2, nConstituentSublattice(iPhaseIDIn,s)
                if (dSiteIn(s,c) > dRefValue) then
                    iRefCon = c
                    dRefValue = dSiteIn(s,c)
                end if
            end do

            do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                if (c /= iRefCon) then
                    nIndependentOut = nIndependentOut + 1
                    iVarSubOut(nIndependentOut) = s
                    iVarConOut(nIndependentOut) = c
                    iVarRefOut(nIndependentOut) = iRefCon
                    iVarSiteOut(nIndependentOut) = iCursor + c
                    iVarRefSiteOut(nIndependentOut) = iCursor + iRefCon
                end if
            end do

            iCursor = iCursor + nConstituentSublattice(iPhaseIDIn,s)
        end do

        return

    end subroutine BuildCEFSubMinVariables

    subroutine ComputeCEFSubMinChargeConstraint(iPhaseIDIn, dSiteIn, nIndependentIn, &
        iVarSubIn, iVarConIn, iVarRefIn, dChargeResidualOut, dChargeJacobianOut)

        implicit none

        integer, intent(in) :: iPhaseIDIn, nIndependentIn
        integer, intent(in) :: iVarSubIn(:), iVarConIn(:), iVarRefIn(:)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dChargeResidualOut
        real(8), intent(out) :: dChargeJacobianOut(:)

        integer :: k, s, c, r

        dChargeResidualOut = 0D0
        dChargeJacobianOut = 0D0
        if (.NOT.allocated(dSublatticeCharge)) return

        do s = 1, nSublatticePhase(iPhaseIDIn)
            do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                dChargeResidualOut = dChargeResidualOut + dStoichSublattice(iPhaseIDIn,s) * &
                    dSublatticeCharge(iPhaseIDIn,s,c) * dSiteIn(s,c)
            end do
        end do

        do k = 1, nIndependentIn
            s = iVarSubIn(k)
            c = iVarConIn(k)
            r = iVarRefIn(k)
            dChargeJacobianOut(k) = dStoichSublattice(iPhaseIDIn,s) * &
                (dSublatticeCharge(iPhaseIDIn,s,c) - dSublatticeCharge(iPhaseIDIn,s,r))
        end do

        return

    end subroutine ComputeCEFSubMinChargeConstraint

    subroutine ComputeCEFSubMinGradient(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, nIndependentIn, &
        iVarSubIn, iVarConIn, iVarRefIn, iVarSiteIn, iVarRefSiteIn, dGradientOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn, nIndependentIn
        integer, intent(in) :: iVarSubIn(:), iVarConIn(:), iVarRefIn(:)
        integer, intent(in) :: iVarSiteIn(:), iVarRefSiteIn(:)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dGradientOut(:)

        integer :: k, s, c, r, nSiteCapacityLocal, nSiteOut, iInfoLocal
        real(8) :: dScalarGibbs, dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity
        real(8) :: dProjectionDerivative
        real(8), allocatable :: dSiteGradient(:), dSiteGradientH(:)
        real(8), allocatable :: dSiteGradientS(:), dSiteGradientCp(:)

        dGradientOut = 0D0
        nSiteCapacityLocal = nMaxSublatticeSys * nMaxConstituentSys

        call SetCEFSubMinMoleFractions(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn)
        allocate(dSiteGradient(nSiteCapacityLocal), dSiteGradientH(nSiteCapacityLocal), &
            dSiteGradientS(nSiteCapacityLocal), dSiteGradientCp(nSiteCapacityLocal))
        call CompGradientSUBL(iSolnPhaseIndexIn, nSiteCapacityLocal, dSiteGradient, &
            dScalarGibbs, dSiteGradientH, dSiteGradientS, dSiteGradientCp, &
            dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity, nSiteOut, iInfoLocal)
        if (iInfoLocal /= 0) then
            dGradientOut(1:nIndependentIn) = 9D5
            if (allocated(dSiteGradient)) deallocate(dSiteGradient)
            if (allocated(dSiteGradientH)) deallocate(dSiteGradientH)
            if (allocated(dSiteGradientS)) deallocate(dSiteGradientS)
            if (allocated(dSiteGradientCp)) deallocate(dSiteGradientCp)
            return
        end if

        do k = 1, nIndependentIn
            s = iVarSubIn(k)
            c = iVarConIn(k)
            r = iVarRefIn(k)

            dProjectionDerivative = ComputeCEFSubMinElementProjection(iPhaseIDIn, s, c) - &
                ComputeCEFSubMinElementProjection(iPhaseIDIn, s, r)
            dGradientOut(k) = dSiteGradient(iVarSiteIn(k)) - dSiteGradient(iVarRefSiteIn(k)) - &
                dProjectionDerivative
        end do

        if (allocated(dSiteGradient)) deallocate(dSiteGradient)
        if (allocated(dSiteGradientH)) deallocate(dSiteGradientH)
        if (allocated(dSiteGradientS)) deallocate(dSiteGradientS)
        if (allocated(dSiteGradientCp)) deallocate(dSiteGradientCp)

        return

    end subroutine ComputeCEFSubMinGradient

    subroutine BuildCEFSubMinActiveSet(iPhaseIDIn, dSiteIn, nIndependentIn, &
        iVarSubIn, iVarConIn, dGradientIn, dActiveFloorIn, lVarActiveOut, &
        iActiveVarOut, nActiveOut, dActiveGradientNormOut)

        implicit none

        integer, intent(in) :: iPhaseIDIn, nIndependentIn
        integer, intent(in) :: iVarSubIn(:), iVarConIn(:)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(in) :: dGradientIn(:), dActiveFloorIn
        logical, intent(out) :: lVarActiveOut(:)
        integer, intent(out) :: iActiveVarOut(:), nActiveOut
        real(8), intent(out) :: dActiveGradientNormOut

        integer :: k, s, c
        logical :: lLowerBoundSatisfied

        lVarActiveOut = .FALSE.
        iActiveVarOut = 0
        nActiveOut = 0
        dActiveGradientNormOut = 0D0

        do k = 1, nIndependentIn
            s = iVarSubIn(k)
            c = iVarConIn(k)
            lLowerBoundSatisfied = (dSiteIn(s,c) <= dActiveFloorIn).AND.(dGradientIn(k) >= 0D0)

            if (.NOT.lLowerBoundSatisfied) then
                nActiveOut = nActiveOut + 1
                lVarActiveOut(k) = .TRUE.
                iActiveVarOut(nActiveOut) = k
                dActiveGradientNormOut = DMAX1(dActiveGradientNormOut, DABS(dGradientIn(k)))
            end if
        end do

        return

    end subroutine BuildCEFSubMinActiveSet

    real(8) function ComputeCEFSubMinElementProjection(iPhaseIDIn, iSublatticeIn, iConstituentIn)

        implicit none

        integer, intent(in) :: iPhaseIDIn, iSublatticeIn, iConstituentIn
        integer :: iElement

        ComputeCEFSubMinElementProjection = 0D0
        iElement = FindCEFSubMinElementIndex(cConstituentNameSUB(iPhaseIDIn,iSublatticeIn,iConstituentIn))
        if (iElement > 0) then
            ComputeCEFSubMinElementProjection = dStoichSublattice(iPhaseIDIn,iSublatticeIn) * &
                dElementPotential(iElement)
        end if

        return

    end function ComputeCEFSubMinElementProjection

    integer function FindCEFSubMinElementIndex(cNameIn)

        implicit none

        character(*), intent(in) :: cNameIn
        character(16) :: cName, cElement, cBase
        integer :: iElement, iPlus, iMinus, iCut

        FindCEFSubMinElementIndex = 0
        cName = NormalizeCEFSubMinName(cNameIn)
        if ((LEN_TRIM(cName) == 0).OR.(TRIM(cName) == 'VA')) return

        cBase = cName
        iPlus = INDEX(cBase, '+')
        iMinus = INDEX(cBase, '-')
        iCut = 0
        if (iPlus > 1) iCut = iPlus
        if ((iMinus > 1).AND.((iCut == 0).OR.(iMinus < iCut))) iCut = iMinus
        if (iCut > 1) cBase = cBase(1:iCut-1)

        do iElement = 1, nElements
            cElement = NormalizeCEFSubMinName(cElementName(iElement))
            if ((TRIM(cName) == TRIM(cElement)).OR.(TRIM(cBase) == TRIM(cElement))) then
                FindCEFSubMinElementIndex = iElement
                return
            end if
        end do

        return

    end function FindCEFSubMinElementIndex

    character(16) function NormalizeCEFSubMinName(cNameIn)

        implicit none

        character(*), intent(in) :: cNameIn
        integer :: iChar, iCode

        NormalizeCEFSubMinName = ADJUSTL(cNameIn)
        do iChar = 1, LEN(NormalizeCEFSubMinName)
            iCode = IACHAR(NormalizeCEFSubMinName(iChar:iChar))
            if ((iCode >= IACHAR('a')).AND.(iCode <= IACHAR('z'))) then
                NormalizeCEFSubMinName(iChar:iChar) = ACHAR(iCode - 32)
            end if
        end do

        return

    end function NormalizeCEFSubMinName

    subroutine ComputeCEFSubMinHessian(iSolnPhaseIndexIn, nSiteCapacityIn, nSpeciesCapacityIn, &
        nIndependentIn, iVarSiteIn, iVarRefSiteIn, dIndependentHessianOut, iInfoOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, nSiteCapacityIn, nSpeciesCapacityIn, nIndependentIn
        integer, intent(in) :: iVarSiteIn(:), iVarRefSiteIn(:)
        real(8), intent(out) :: dIndependentHessianOut(:,:)
        integer, intent(out) :: iInfoOut

        integer :: i, j, nSiteOut, nSpeciesOut
        real(8), allocatable :: dSiteHessian(:,:), dAmountHessian(:,:), dJacobian(:,:)

        dIndependentHessianOut = 0D0
        allocate(dSiteHessian(nSiteCapacityIn,nSiteCapacityIn))
        allocate(dAmountHessian(nSpeciesCapacityIn,nSpeciesCapacityIn))
        allocate(dJacobian(nSpeciesCapacityIn,nSpeciesCapacityIn))

        call CompHessianSUBL(iSolnPhaseIndexIn, nSiteCapacityIn, nSpeciesCapacityIn, &
            dSiteHessian, dAmountHessian, dJacobian, nSiteOut, nSpeciesOut, iInfoOut)

        if (iInfoOut == 0) then
            do i = 1, nIndependentIn
                do j = 1, nIndependentIn
                    dIndependentHessianOut(i,j) = &
                        dSiteHessian(iVarSiteIn(i),iVarSiteIn(j)) - &
                        dSiteHessian(iVarSiteIn(i),iVarRefSiteIn(j)) - &
                        dSiteHessian(iVarRefSiteIn(i),iVarSiteIn(j)) + &
                        dSiteHessian(iVarRefSiteIn(i),iVarRefSiteIn(j))
                end do
            end do
        end if

        if (allocated(dSiteHessian)) deallocate(dSiteHessian)
        if (allocated(dAmountHessian)) deallocate(dAmountHessian)
        if (allocated(dJacobian)) deallocate(dJacobian)

        return

    end subroutine ComputeCEFSubMinHessian

    subroutine SolveCEFSubMinDirection(nIndependentIn, dHessianIn, dGradientIn, &
        lChargeConstraintIn, dChargeResidualIn, dChargeJacobianIn, dDirectionOut)

        implicit none

        integer, intent(in) :: nIndependentIn
        real(8), intent(in) :: dHessianIn(:,:), dGradientIn(:)
        logical, intent(in) :: lChargeConstraintIn
        real(8), intent(in) :: dChargeResidualIn, dChargeJacobianIn(:)
        real(8), intent(out) :: dDirectionOut(:)

        integer :: i, iTry, iInfo, lWork, nKKT
        integer, allocatable :: iPivot(:)
        real(8) :: dShift, dWorkQuery(1), dChargeNorm, dChargeCorrection
        real(8), allocatable :: dMatrix(:,:), dRHSWork(:), dWork(:)

        dDirectionOut = 0D0
        if (nIndependentIn <= 0) return

        nKKT = nIndependentIn
        if (lChargeConstraintIn) nKKT = nKKT + 1

        allocate(iPivot(nKKT))
        allocate(dMatrix(nKKT,nKKT))
        allocate(dRHSWork(nKKT))

        do iTry = 1, 8
            dMatrix = 0D0
            dRHSWork = 0D0
            dMatrix(1:nIndependentIn,1:nIndependentIn) = dHessianIn(1:nIndependentIn,1:nIndependentIn)
            if (iTry == 1) then
                dShift = 0D0
            else
                dShift = 10D0**(iTry-9)
            end if
            do i = 1, nIndependentIn
                dMatrix(i,i) = dMatrix(i,i) + dShift
                dRHSWork(i) = -dGradientIn(i)
            end do

            if (lChargeConstraintIn) then
                do i = 1, nIndependentIn
                    dMatrix(i,nKKT) = dChargeJacobianIn(i)
                    dMatrix(nKKT,i) = dChargeJacobianIn(i)
                end do
                dRHSWork(nKKT) = -dChargeResidualIn
            end if

            lWork = -1
            dWorkQuery = 0D0
            call DSYSV('U', nKKT, 1, dMatrix, nKKT, iPivot, &
                dRHSWork, nKKT, dWorkQuery, lWork, iInfo)

            if (iInfo == 0) then
                lWork = MAX(1, INT(dWorkQuery(1)))
                allocate(dWork(lWork))
                call DSYSV('U', nKKT, 1, dMatrix, nKKT, iPivot, &
                    dRHSWork, nKKT, dWork, lWork, iInfo)
                if (allocated(dWork)) deallocate(dWork)
            end if

            if (iInfo == 0) then
                dDirectionOut(1:nIndependentIn) = dRHSWork(1:nIndependentIn)
                if (lChargeConstraintIn) then
                    exit
                else if (SUM(dGradientIn(1:nIndependentIn) * dDirectionOut(1:nIndependentIn)) < 0D0) then
                    exit
                end if
            end if
        end do

        if ((iInfo /= 0).OR.&
            ((.NOT.lChargeConstraintIn).AND.&
            (SUM(dGradientIn(1:nIndependentIn) * dDirectionOut(1:nIndependentIn)) >= 0D0))) then
            dDirectionOut(1:nIndependentIn) = -dGradientIn(1:nIndependentIn)
        end if

        if (lChargeConstraintIn) then
            dChargeNorm = SUM(dChargeJacobianIn(1:nIndependentIn)**2)
            if (dChargeNorm > 1D-300) then
                dChargeCorrection = (SUM(dChargeJacobianIn(1:nIndependentIn) * &
                    dDirectionOut(1:nIndependentIn)) + dChargeResidualIn) / dChargeNorm
                dDirectionOut(1:nIndependentIn) = dDirectionOut(1:nIndependentIn) - &
                    dChargeCorrection * dChargeJacobianIn(1:nIndependentIn)
            end if
        end if

        if (allocated(iPivot)) deallocate(iPivot)
        if (allocated(dMatrix)) deallocate(dMatrix)
        if (allocated(dRHSWork)) deallocate(dRHSWork)
        if (allocated(dWork)) deallocate(dWork)

        return

    end subroutine SolveCEFSubMinDirection

    subroutine SolveCEFSubMinActiveDirection(nIndependentIn, nActiveIn, iActiveVarIn, &
        dHessianIn, dGradientIn, lChargeConstraintIn, dChargeResidualIn, dChargeJacobianIn, &
        dDirectionOut)

        implicit none

        integer, intent(in) :: nIndependentIn, nActiveIn
        integer, intent(in) :: iActiveVarIn(:)
        real(8), intent(in) :: dHessianIn(:,:), dGradientIn(:)
        logical, intent(in) :: lChargeConstraintIn
        real(8), intent(in) :: dChargeResidualIn, dChargeJacobianIn(:)
        real(8), intent(out) :: dDirectionOut(:)

        integer :: i, j, iFull, jFull
        real(8), allocatable :: dReducedHessian(:,:), dReducedGradient(:)
        real(8), allocatable :: dReducedChargeJacobian(:), dReducedDirection(:)

        dDirectionOut = 0D0
        if (nActiveIn <= 0) return

        allocate(dReducedHessian(nActiveIn,nActiveIn))
        allocate(dReducedGradient(nActiveIn), dReducedChargeJacobian(nActiveIn))
        allocate(dReducedDirection(nActiveIn))

        dReducedHessian = 0D0
        dReducedGradient = 0D0
        dReducedChargeJacobian = 0D0
        dReducedDirection = 0D0

        do i = 1, nActiveIn
            iFull = iActiveVarIn(i)
            dReducedGradient(i) = dGradientIn(iFull)
            dReducedChargeJacobian(i) = dChargeJacobianIn(iFull)
            do j = 1, nActiveIn
                jFull = iActiveVarIn(j)
                dReducedHessian(i,j) = dHessianIn(iFull,jFull)
            end do
        end do

        call SolveCEFSubMinDirection(nActiveIn, dReducedHessian, dReducedGradient, &
            lChargeConstraintIn, dChargeResidualIn, dReducedChargeJacobian, dReducedDirection)

        do i = 1, nActiveIn
            dDirectionOut(iActiveVarIn(i)) = dReducedDirection(i)
        end do

        if (allocated(dReducedHessian)) deallocate(dReducedHessian)
        if (allocated(dReducedGradient)) deallocate(dReducedGradient)
        if (allocated(dReducedChargeJacobian)) deallocate(dReducedChargeJacobian)
        if (allocated(dReducedDirection)) deallocate(dReducedDirection)

        return

    end subroutine SolveCEFSubMinActiveDirection

    subroutine SetCEFSubMinActiveGradientDirection(nIndependentIn, lVarActiveIn, dGradientIn, dDirectionOut)

        implicit none

        integer, intent(in) :: nIndependentIn
        logical, intent(in) :: lVarActiveIn(:)
        real(8), intent(in) :: dGradientIn(:)
        real(8), intent(inout) :: dDirectionOut(:)

        integer :: k

        do k = 1, nIndependentIn
            if (lVarActiveIn(k)) dDirectionOut(k) = -dGradientIn(k)
        end do

        return

    end subroutine SetCEFSubMinActiveGradientDirection

    subroutine LimitCEFSubMinDirection(nIndependentIn, dDirectionInOut)

        implicit none

        integer, intent(in) :: nIndependentIn
        real(8), intent(inout) :: dDirectionInOut(:)

        real(8) :: dMaxDirection

        if (nIndependentIn <= 0) return
        dMaxDirection = MAXVAL(DABS(dDirectionInOut(1:nIndependentIn)))
        if (dMaxDirection > 0.2D0) then
            dDirectionInOut(1:nIndependentIn) = 0.2D0 * dDirectionInOut(1:nIndependentIn) / dMaxDirection
        end if

        return

    end subroutine LimitCEFSubMinDirection

    subroutine LineSearchCEFSubMin(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, nIndependentIn, &
        iVarSubIn, iVarConIn, iVarRefIn, dGradientIn, dDirectionIn, dObjectiveInOut, &
        dFloor, dArmijoIn, lChargeConstraintIn, dChargeResidualIn, dChargeJacobianIn, &
        dChargePenaltyIn, dSiteAcceptedOut, lAcceptedOut, dStepNormOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn, nIndependentIn
        integer, intent(in) :: iVarSubIn(:), iVarConIn(:), iVarRefIn(:)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(in) :: dGradientIn(:), dDirectionIn(:), dFloor, dArmijoIn
        logical, intent(in) :: lChargeConstraintIn
        real(8), intent(in) :: dChargeResidualIn, dChargeJacobianIn(:), dChargePenaltyIn
        real(8), intent(inout) :: dObjectiveInOut
        real(8), intent(out) :: dSiteAcceptedOut(nMaxSublatticeSys,nMaxConstituentSys)
        logical, intent(out) :: lAcceptedOut
        real(8), intent(out) :: dStepNormOut

        integer :: k, iLine, s, c, r
        real(8) :: dAlpha, dTrialObjective, dDirectionalDerivative
        real(8) :: dMeritBase, dMeritTrial, dMeritDirectionalDerivative
        real(8) :: dChargeTrial, dChargeDirectionalDerivative
        logical :: lFeasible
        real(8) :: dSiteWork(nMaxSublatticeSys,nMaxConstituentSys)
        real(8) :: dChargeJacobianTrial(nIndependentIn)

        lAcceptedOut = .FALSE.
        dSiteAcceptedOut = dSiteIn
        dStepNormOut = 0D0
        dDirectionalDerivative = SUM(dGradientIn(1:nIndependentIn) * dDirectionIn(1:nIndependentIn))
        dChargeDirectionalDerivative = 0D0
        dMeritBase = dObjectiveInOut
        dMeritDirectionalDerivative = dDirectionalDerivative

        if (lChargeConstraintIn) then
            dChargeDirectionalDerivative = SUM(dChargeJacobianIn(1:nIndependentIn) * &
                dDirectionIn(1:nIndependentIn))
            dMeritBase = dMeritBase + dChargePenaltyIn * dChargeResidualIn**2
            dMeritDirectionalDerivative = dMeritDirectionalDerivative + &
                2D0 * dChargePenaltyIn * dChargeResidualIn * dChargeDirectionalDerivative
        end if

        dAlpha = 1D0
        do iLine = 1, 30
            dSiteWork = dSiteIn
            lFeasible = .TRUE.

            do k = 1, nIndependentIn
                s = iVarSubIn(k)
                c = iVarConIn(k)
                r = iVarRefIn(k)
                dSiteWork(s,c) = dSiteWork(s,c) + dAlpha * dDirectionIn(k)
                dSiteWork(s,r) = dSiteWork(s,r) - dAlpha * dDirectionIn(k)
            end do

            do s = 1, nSublatticePhase(iPhaseIDIn)
                do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                    if (dSiteWork(s,c) <= dFloor) lFeasible = .FALSE.
                end do
            end do

            if (lFeasible) then
                call SetCEFSubMinMoleFractions(iSolnPhaseIndexIn, iPhaseIDIn, dSiteWork)
                call EvaluateCEFSubMinObjective(iSolnPhaseIndexIn, iPhaseIDIn, dSiteWork, dTrialObjective)
                dMeritTrial = dTrialObjective

                if (lChargeConstraintIn) then
                    call ComputeCEFSubMinChargeConstraint(iPhaseIDIn, dSiteWork, nIndependentIn, &
                        iVarSubIn, iVarConIn, iVarRefIn, dChargeTrial, dChargeJacobianTrial)
                    dMeritTrial = dMeritTrial + dChargePenaltyIn * dChargeTrial**2
                end if

                if (dMeritTrial <= dMeritBase + dArmijoIn*dAlpha*dMeritDirectionalDerivative) then
                    lAcceptedOut = .TRUE.
                    dSiteAcceptedOut = dSiteWork
                    dStepNormOut = dAlpha * MAXVAL(DABS(dDirectionIn(1:nIndependentIn)))
                    dObjectiveInOut = dTrialObjective
                    exit
                end if
            end if

            dAlpha = 0.5D0 * dAlpha
        end do

        call SetCEFSubMinMoleFractions(iSolnPhaseIndexIn, iPhaseIDIn, dSiteAcceptedOut)
        call EvaluateCEFSubMinObjective(iSolnPhaseIndexIn, iPhaseIDIn, dSiteAcceptedOut, dTrialObjective)

        return

    end subroutine LineSearchCEFSubMin

end subroutine SubMinSiteFractionCEF
