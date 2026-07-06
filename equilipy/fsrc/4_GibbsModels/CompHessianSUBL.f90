!> \brief Compute analytical fixed-state CEF Hessians for SUBL/SUBOM phases.
!!
!! \details This diagnostic routine returns the site-fraction Hessian and its
!! projection to endmember space for one CEF phase.  Magnetic phases include
!! chain-rule curvature of the same scalar magnetic contribution used by the
!! Gibbs-energy evaluator.  It does not alter the minimizer state and it does
!! not replace the current GEM Newton system.
!!
!! \param[in]  iSolnIndex                      Absolute solution phase index.
!! \param[in]  nSiteDim                        Caller-provided site Hessian capacity.
!! \param[in]  nSpeciesDim                     Caller-provided species Hessian capacity.
!! \param[out] dSiteHessian                    Site-space Hessian in internal G/(RT) units.
!! \param[out] dEndmemberAmountHessian         Projected amount Hessian for unit phase amount.
!! \param[out] dEndmemberCompositionJacobian   Redundant endmember-composition Jacobian.
!! \param[out] nSiteOut                        Number of active site variables written.
!! \param[out] nSpeciesOut                     Number of active phase endmembers written.
!! \param[out] iInfo                           Zero on success; nonzero for unsupported or undersized calls.



subroutine CompHessianSUBL(iSolnIndex, nSiteDim, nSpeciesDim, dSiteHessian, &
    dEndmemberAmountHessian, dEndmemberCompositionJacobian, nSiteOut, nSpeciesOut, iInfo)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompHessianSUBL.f90
    !> \brief   Compute analytical fixed-state CEF Hessians for SUBL/SUBOM phases.
    !> \author  S.Y. Kwon
    !> \date    Jun. 24, 2026
    !> \sa      CompExcessGibbsEnergySUBL.f90
    !> \sa      solution_models.md
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/24/2026      S.Y. Kwon           Added diagnostic analytical CEF fixed-state Hessian.
    !   06/24/2026      S.Y. Kwon           Added magnetic scalar Hessian contribution.
    !   06/24/2026      S.Y. Kwon           Removed singular RK linear-factor Hessian evaluation.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details This routine evaluates the site-fraction Hessian
    !! \f$K_{AB} = \partial^2 g / \partial y_A \partial y_B\f$ for the
    !! reference, ideal configurational, nonmagnetic excess, and magnetic
    !! scalar terms of SUBL phases.  For SUBOM phases it applies the fixed ordered-state
    !! order/disorder correction
    !! \f$K^{fix}=K^{ord}+C^T K^{dis} C-M^T K^{ord,r} M\f$.
    !
    ! Required input variables:
    ! =========================
    !
    !> \param[in]  iSolnIndex   Absolute solution phase index.
    !> \param[in]  nSiteDim     First/second dimension of dSiteHessian.
    !> \param[in]  nSpeciesDim  First/second dimension of endmember Hessian/Jacobian arrays.
    !
    ! cSolnPhaseType             Solution model type.  SUBL, SUBLM, and SUBOM are accepted.
    ! dMolFraction               Current endmember fractions used to rebuild site fractions.
    ! dStdGibbsEnergy            Dimensionless standard endmember Gibbs energies.
    ! dExcessGibbsParam          Dimensionless nonmagnetic excess parameters.
    ! iRegularParam/iSUBLParamData
    !                            Runtime CEF parameter topology.
    ! iDisorderedPhase           Ordered-to-disordered phase mapping for SUBOM.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] dSiteHessian                   Site-space Hessian in G/(RT) units.
    !> \param[out] dEndmemberAmountHessian        Symmetric projected amount Hessian for unit phase amount.
    !> \param[out] dEndmemberCompositionJacobian  Redundant endmember-composition Jacobian.
    !> \param[out] nSiteOut                       Number of active site variables written.
    !> \param[out] nSpeciesOut                    Number of active phase endmembers written.
    !> \param[out] iInfo                          Zero on success.  1=unsupported phase, 2=bad topology,
    !!                                             3=caller arrays too small.
    !
    ! dSiteFraction              Rebuilt for iSolnIndex from dMolFraction.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! BuildNormalizedSiteFractions       Builds normalized site fractions from endmember fractions.
    ! AccumulateRawSUBLSiteHessian       Adds reference, ideal, excess, and magnetic CEF curvature.
    ! AccumulateMagneticSiteHessian      Adds scalar magnetic curvature matching CompGibbsMagneticSoln.
    ! AccumulateMagneticVariableDerivatives
    !                                    Computes site derivatives of magnetic TCrit and B variables.
    ! EvaluateMagneticScalarDerivatives  Computes scalar derivatives with respect to TCrit and B.
    ! AddProductLinearTermHessian        Differentiates product/linear-factor CEF interaction terms.
    ! AddSUBOMFixedCorrection            Applies fixed ordered-state order/disorder Hessian mapping.
    ! ProjectSiteHessianToEndmembers     Projects site curvature to endmember amount/composition coordinates.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! Python diagnostics/devop tests      Validate analytical curvature against finite differences.
    ! Future site-fraction minimizer      Candidate source for a Hessian-based local minimization step.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Returned Hessians use the internal dimensionless Gibbs-energy convention, G/(RT).
    ! - Site fractions are treated as fixed ordered-state coordinates.
    ! - Magnetic Hessian terms use analytical site derivatives of TCrit and B with numerical derivatives only for
    !   the two-variable Hillert-Jarl scalar.  This keeps the diagnostic Hessian consistent with
    !   CompExcessGibbsEnergy for magnetic SUBL/SUBOM phases while avoiding site-coordinate finite differences.
    ! - The full site-space matrix is redundant; callers should project to independent site coordinates before
    !   inverting or using it in a Newton/trust-region step.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer, intent(in)  :: iSolnIndex, nSiteDim, nSpeciesDim
    real(8), intent(out) :: dSiteHessian(nSiteDim,nSiteDim)
    real(8), intent(out) :: dEndmemberAmountHessian(nSpeciesDim,nSpeciesDim)
    real(8), intent(out) :: dEndmemberCompositionJacobian(nSpeciesDim,nSpeciesDim)
    integer, intent(out) :: nSiteOut, nSpeciesOut, iInfo

    integer :: iPhaseID, nSublattice, iFirst, iLast
    integer :: iSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
    integer :: iSiteSub(nSiteDim), iSiteCon(nSiteDim)
    real(8) :: dWorkSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8) :: dWorkHessian(nSiteDim,nSiteDim)

    dSiteHessian = 0D0
    dEndmemberAmountHessian = 0D0
    dEndmemberCompositionJacobian = 0D0
    nSiteOut = 0
    nSpeciesOut = 0
    iInfo = 0

    if ((iSolnIndex <= 0).OR.(iSolnIndex > nSolnPhasesSys)) then
        iInfo = 1
        return
    end if
    if ((TRIM(cSolnPhaseType(iSolnIndex)) /= 'SUBL').AND. &
        (TRIM(cSolnPhaseType(iSolnIndex)) /= 'SUBLM').AND. &
        (TRIM(cSolnPhaseType(iSolnIndex)) /= 'SUBOM')) then
        iInfo = 1
        return
    end if

    iPhaseID = iPhaseSublattice(iSolnIndex)
    if (iPhaseID <= 0) then
        iInfo = 2
        return
    end if

    nSublattice = nSublatticePhase(iPhaseID)
    iFirst = nSpeciesPhase(iSolnIndex-1) + 1
    iLast = nSpeciesPhase(iSolnIndex)
    nSpeciesOut = iLast - iFirst + 1

    call BuildSiteIndex(iPhaseID, nSiteDim, iSiteIndex, iSiteSub, iSiteCon, nSiteOut)
    if ((nSiteOut <= 0).OR.(nSpeciesOut <= 0)) then
        iInfo = 2
        return
    end if
    if ((nSiteOut > nSiteDim).OR.(nSpeciesOut > nSpeciesDim)) then
        iInfo = 3
        return
    end if

    call BuildNormalizedSiteFractions(iSolnIndex, iPhaseID, dWorkSite)
    dSiteFraction(iPhaseID,:,:) = dWorkSite(:,:)

    call AccumulateRawSUBLSiteHessian(iSolnIndex, iPhaseID, nSiteDim, iSiteIndex, &
        iSiteSub, iSiteCon, nSiteOut, dWorkSite, dWorkHessian)

    if (TRIM(cSolnPhaseType(iSolnIndex)) == 'SUBOM') then
        call AddSUBOMFixedCorrection(iSolnIndex, nSiteDim, iSiteIndex, iSiteSub, iSiteCon, &
            nSiteOut, dWorkSite, dWorkHessian)
    end if

    dSiteHessian(1:nSiteOut,1:nSiteOut) = dWorkHessian(1:nSiteOut,1:nSiteOut)
    call ProjectSiteHessianToEndmembers(iSolnIndex, iPhaseID, nSiteDim, nSpeciesDim, &
        iSiteIndex, iSiteSub, iSiteCon, nSiteOut, nSpeciesOut, dWorkSite, dSiteHessian, &
        dEndmemberAmountHessian, dEndmemberCompositionJacobian)

    return

contains

    subroutine BuildSiteIndex(iPhaseIDIn, nSiteCapacity, iSiteIndexOut, iSiteSubOut, iSiteConOut, nSiteLocal)

        implicit none

        integer, intent(in) :: iPhaseIDIn, nSiteCapacity
        integer, intent(out) :: iSiteIndexOut(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(out) :: iSiteSubOut(nSiteCapacity), iSiteConOut(nSiteCapacity)
        integer, intent(out) :: nSiteLocal

        integer :: s, c

        iSiteIndexOut = 0
        iSiteSubOut = 0
        iSiteConOut = 0
        nSiteLocal = 0

        do s = 1, nSublatticePhase(iPhaseIDIn)
            do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                nSiteLocal = nSiteLocal + 1
                if (nSiteLocal <= nSiteCapacity) then
                    iSiteIndexOut(s,c) = nSiteLocal
                    iSiteSubOut(nSiteLocal) = s
                    iSiteConOut(nSiteLocal) = c
                end if
            end do
        end do

        return

    end subroutine BuildSiteIndex

    subroutine BuildNormalizedSiteFractions(iSolnIndexIn, iPhaseIDIn, dSiteOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iFirstLocal, iLastLocal, i, m, s, c
        real(8) :: dNorm

        dSiteOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnIndexIn)

        do i = iFirstLocal, iLastLocal
            m = i - iFirstLocal + 1
            do s = 1, nSublatticePhase(iPhaseIDIn)
                c = iConstituentSublattice(iPhaseIDIn,s,m)
                dSiteOut(s,c) = dSiteOut(s,c) + DMAX1(dMolFraction(i), 0D0)
            end do
        end do

        do s = 1, nSublatticePhase(iPhaseIDIn)
            dNorm = 0D0
            do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                dNorm = dNorm + dSiteOut(s,c)
            end do
            if (dNorm > 0D0) then
                do c = 1, nConstituentSublattice(iPhaseIDIn,s)
                    dSiteOut(s,c) = dSiteOut(s,c) / dNorm
                end do
            end if
        end do

        return

    end subroutine BuildNormalizedSiteFractions

    subroutine AccumulateRawSUBLSiteHessian(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
        iSiteSubIn, iSiteConIn, nSiteLocal, dSiteIn, dHessianOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn, nSiteCapacity, nSiteLocal
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dHessianOut(nSiteCapacity,nSiteCapacity)

        integer :: iFirstLocal, iLastLocal, i, m, s, t, a, b, c
        real(8) :: dProduct, dSiteA, dSiteB

        dHessianOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnIndexIn)

        do i = iFirstLocal, iLastLocal
            m = i - iFirstLocal + 1
            dProduct = 1D0
            do s = 1, nSublatticePhase(iPhaseIDIn)
                c = iConstituentSublattice(iPhaseIDIn,s,m)
                dProduct = dProduct * DMAX1(dSiteIn(s,c), 1D-75)
            end do

            do s = 1, nSublatticePhase(iPhaseIDIn)
                a = iSiteIndexIn(s,iConstituentSublattice(iPhaseIDIn,s,m))
                dSiteA = DMAX1(dSiteIn(s,iConstituentSublattice(iPhaseIDIn,s,m)), 1D-75)
                do t = 1, nSublatticePhase(iPhaseIDIn)
                    if (s == t) cycle
                    b = iSiteIndexIn(t,iConstituentSublattice(iPhaseIDIn,t,m))
                    dSiteB = DMAX1(dSiteIn(t,iConstituentSublattice(iPhaseIDIn,t,m)), 1D-75)
                    dHessianOut(a,b) = dHessianOut(a,b) + dStdGibbsEnergy(i) * dProduct / (dSiteA * dSiteB)
                end do
            end do
        end do

        do a = 1, nSiteLocal
            s = iSiteSubIn(a)
            c = iSiteConIn(a)
            dHessianOut(a,a) = dHessianOut(a,a) + dStoichSublattice(iPhaseIDIn,s) / &
                DMAX1(dSiteIn(s,c), 1D-75)
        end do

        do i = nParamPhase(iSolnIndexIn-1) + 1, nParamPhase(iSolnIndexIn)
            call AccumulateOneSUBLParameter(i, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
                iSiteSubIn, iSiteConIn, dSiteIn, dHessianOut)
        end do

        if (HasMagneticTerms(iSolnIndexIn)) then
            call AccumulateMagneticSiteHessian(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
                iSiteSubIn, iSiteConIn, nSiteLocal, dSiteIn, dHessianOut)
        end if

        return

    end subroutine AccumulateRawSUBLSiteHessian

    subroutine AccumulateMagneticSiteHessian(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
        iSiteSubIn, iSiteConIn, nSiteLocal, dSiteIn, dHessianOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn, nSiteCapacity, nSiteLocal
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(inout) :: dHessianOut(nSiteCapacity,nSiteCapacity)

        integer :: a, b
        real(8) :: dTCrit, dB, dMagG, dMagTCrit, dMagB, dMagTT, dMagTB, dMagBB
        real(8) :: dTCritGrad(nSiteCapacity), dBGrad(nSiteCapacity)
        real(8) :: dTCritHessian(nSiteCapacity,nSiteCapacity), dBHessian(nSiteCapacity,nSiteCapacity)

        if (.NOT.HasMagneticTerms(iSolnIndexIn)) return
        if ((TRIM(cSolnPhaseType(iSolnIndexIn)) /= 'SUBLM').AND. &
            (TRIM(cSolnPhaseType(iSolnIndexIn)) /= 'SUBOM')) return

        call AccumulateMagneticVariableDerivatives(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, &
            iSiteIndexIn, iSiteSubIn, iSiteConIn, dSiteIn, dTCrit, dB, dTCritGrad, dBGrad, &
            dTCritHessian, dBHessian)
        call EvaluateMagneticScalarDerivatives(iSolnIndexIn, dTCrit, dB, dMagG, &
            dMagTCrit, dMagB, dMagTT, dMagTB, dMagBB)

        do a = 1, nSiteLocal
            do b = 1, nSiteLocal
                dHessianOut(a,b) = dHessianOut(a,b) + &
                    dMagTT * dTCritGrad(a) * dTCritGrad(b) + &
                    dMagTB * (dTCritGrad(a) * dBGrad(b) + dBGrad(a) * dTCritGrad(b)) + &
                    dMagBB * dBGrad(a) * dBGrad(b) + &
                    dMagTCrit * dTCritHessian(a,b) + dMagB * dBHessian(a,b)
            end do
        end do

        return

    end subroutine AccumulateMagneticSiteHessian

    subroutine AccumulateMagneticVariableDerivatives(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, &
        iSiteIndexIn, iSiteSubIn, iSiteConIn, dSiteIn, dTCrit, dB, dTCritGrad, dBGrad, &
        dTCritHessian, dBHessian)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn, nSiteCapacity
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dTCrit, dB
        real(8), intent(out) :: dTCritGrad(nSiteCapacity), dBGrad(nSiteCapacity)
        real(8), intent(out) :: dTCritHessian(nSiteCapacity,nSiteCapacity), dBHessian(nSiteCapacity,nSiteCapacity)

        integer :: iFirstLocal, iLastLocal, i, m, s, c, a, k, iParam, nParamCon, iExponent
        real(8) :: dNu(nSiteCapacity)
        real(8) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)

        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnIndexIn)
        dTCrit = 0D0
        dB = 0D0
        dTCritGrad = 0D0
        dBGrad = 0D0
        dTCritHessian = 0D0
        dBHessian = 0D0

        do i = iFirstLocal, iLastLocal
            m = i - iFirstLocal + 1
            dNu = 0D0
            do s = 1, nSublatticePhase(iPhaseIDIn)
                c = iConstituentSublattice(iPhaseIDIn,s,m)
                a = iSiteIndexIn(s,c)
                if (a > 0) dNu(a) = dNu(a) + 1D0
            end do
            call AccumulateProductLinearTermDerivatives(nSiteCapacity, dCoeffGibbsMagnetic(i,1), dNu, &
                0, dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
                dTCrit, dTCritGrad, dTCritHessian)
            call AccumulateProductLinearTermDerivatives(nSiteCapacity, dCoeffGibbsMagnetic(i,2), dNu, &
                0, dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
                dB, dBGrad, dBHessian)
        end do

        do iParam = nMagParamPhase(iSolnIndexIn-1) + 1, nMagParamPhase(iSolnIndexIn)
            nParamCon = iMagneticParam(iParam,1)
            if (nParamCon <= 0) cycle

            dNu = 0D0
            dLinearValue = 0D0
            dLinearPower = 0D0
            dLinearCoeff = 0D0
            do k = 2, nParamCon + 1
                c = MOD(iMagneticParam(iParam,k), 10000)
                s = (iMagneticParam(iParam,k) - c) / 10000
                a = iSiteIndexIn(s,c)
                if (a > 0) dNu(a) = dNu(a) + 1D0
                if (k == 2) then
                    dLinearCoeff(1,a) = 1D0
                    dLinearValue(1) = dLinearValue(1) + dSiteIn(s,c)
                else if (k == 3) then
                    dLinearCoeff(1,a) = dLinearCoeff(1,a) - 1D0
                    dLinearValue(1) = dLinearValue(1) - dSiteIn(s,c)
                end if
            end do

            iExponent = iMagneticParam(iParam,nParamCon+2)
            if (iExponent > 0) dLinearPower(1) = DBLE(iExponent)
            call AccumulateProductLinearTermDerivatives(nSiteCapacity, dMagneticParam(iParam,1), dNu, &
                MERGE(1, 0, iExponent > 0), dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, &
                iSiteSubIn, iSiteConIn, dTCrit, dTCritGrad, dTCritHessian)
            call AccumulateProductLinearTermDerivatives(nSiteCapacity, dMagneticParam(iParam,2), dNu, &
                MERGE(1, 0, iExponent > 0), dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, &
                iSiteSubIn, iSiteConIn, dB, dBGrad, dBHessian)
        end do

        return

    end subroutine AccumulateMagneticVariableDerivatives

    subroutine AccumulateProductLinearTermDerivatives(nSiteCapacity, dParameter, dNu, nLinear, &
        dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
        dValueOut, dGradientOut, dHessianOut)

        implicit none

        integer, intent(in) :: nSiteCapacity, nLinear
        real(8), intent(in) :: dParameter, dNu(nSiteCapacity)
        real(8), intent(in) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(inout) :: dValueOut, dGradientOut(nSiteCapacity)
        real(8), intent(inout) :: dHessianOut(nSiteCapacity,nSiteCapacity)

        real(8) :: dValueTerm
        real(8) :: dGradientTerm(nSiteCapacity)
        real(8) :: dHessianTerm(nSiteCapacity,nSiteCapacity)

        call EvaluateProductLinearDerivatives(nSiteCapacity, dParameter, dNu, nLinear, &
            dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
            dValueTerm, dGradientTerm, dHessianTerm)

        dValueOut = dValueOut + dValueTerm
        dGradientOut = dGradientOut + dGradientTerm
        dHessianOut = dHessianOut + dHessianTerm

        return

    end subroutine AccumulateProductLinearTermDerivatives

    subroutine EvaluateMagneticScalarDerivatives(iSolnIndexIn, dTCritIn, dBIn, dMagG, &
        dMagTCrit, dMagB, dMagTT, dMagTB, dMagBB)

        implicit none

        integer, intent(in) :: iSolnIndexIn
        real(8), intent(in) :: dTCritIn, dBIn
        real(8), intent(out) :: dMagG, dMagTCrit, dMagB, dMagTT, dMagTB, dMagBB

        real(8) :: dTCritStep, dBStep, dPlusG, dMinusG, dPPG, dPMG, dMPG, dMMG

        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn, dBIn, dMagG)
        dMagTCrit = 0D0
        dMagB = 0D0
        dMagTT = 0D0
        dMagTB = 0D0
        dMagBB = 0D0
        if (dTCritIn == 0D0) return

        dTCritStep = DMAX1(1D-4, 1D-6 * DABS(dTCritIn))
        dBStep = DMAX1(1D-7, 1D-6 * DMAX1(1D0, DABS(dBIn)))

        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn + dTCritStep, dBIn, dPlusG)
        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn - dTCritStep, dBIn, dMinusG)
        dMagTCrit = (dPlusG - dMinusG) / (2D0 * dTCritStep)
        dMagTT = (dPlusG - 2D0*dMagG + dMinusG) / (dTCritStep * dTCritStep)

        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn, dBIn + dBStep, dPlusG)
        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn, dBIn - dBStep, dMinusG)
        dMagB = (dPlusG - dMinusG) / (2D0 * dBStep)
        dMagBB = (dPlusG - 2D0*dMagG + dMinusG) / (dBStep * dBStep)

        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn + dTCritStep, dBIn + dBStep, dPPG)
        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn + dTCritStep, dBIn - dBStep, dPMG)
        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn - dTCritStep, dBIn + dBStep, dMPG)
        call EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn - dTCritStep, dBIn - dBStep, dMMG)
        dMagTB = (dPPG - dPMG - dMPG + dMMG) / (4D0 * dTCritStep * dBStep)

        return

    end subroutine EvaluateMagneticScalarDerivatives

    subroutine EvaluateMagneticScalarFromTCritB(iSolnIndexIn, dTCritIn, dBIn, dMagG)

        implicit none

        integer, intent(in) :: iSolnIndexIn
        real(8), intent(in) :: dTCritIn, dBIn
        real(8), intent(out) :: dMagG

        integer :: iFirstLocal
        real(8) :: dTCritLocal, dBLocal, dStructureFactor, dP, dInvPMinusOne
        real(8) :: dTau, dD, dA, dBTemp, dC, dF

        dMagG = 0D0
        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        dStructureFactor = dCoeffGibbsMagnetic(iFirstLocal,3)
        dP = dCoeffGibbsMagnetic(iFirstLocal,4)
        if (dP == 0D0) return

        dTCritLocal = dTCritIn
        dBLocal = dBIn
        if (dBLocal < 0D0) dBLocal = -dBLocal * dStructureFactor
        if (dTCritLocal < 0D0) dTCritLocal = -dTCritLocal * dStructureFactor
        if ((dTCritLocal == 0D0).OR.(dBLocal <= -1D0)) return

        dTau = dTemperature / dTCritLocal
        dInvPMinusOne = 1D0 / dP - 1D0
        dD = (518D0/1125D0) + (11692D0/15975D0) * dInvPMinusOne

        if (dTau > 1D0) then
            dA = dTau**(-5)
            dBTemp = dA**3
            dC = dA * dA * dBTemp
            dF = -(dA/10D0 + dBTemp/315D0 + dC/1500D0) / dD
        else
            dA = dTau**3
            dBTemp = dA**3
            dC = dA * dA * dBTemp
            dF = 1D0 - (79D0/(140D0*dP*dTau) + &
                (474D0/497D0)*dInvPMinusOne*(dA/6D0 + dBTemp/135D0 + dC/600D0)) / dD
        end if

        dMagG = DLOG(1D0 + dBLocal) * dF

        return

    end subroutine EvaluateMagneticScalarFromTCritB

    subroutine AccumulateOneSUBLParameter(iParam, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
        iSiteSubIn, iSiteConIn, dSiteIn, dHessianOut)

        implicit none

        integer, intent(in) :: iParam, iPhaseIDIn, nSiteCapacity
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(inout) :: dHessianOut(nSiteCapacity,nSiteCapacity)

        integer :: k, nParamCon, iMixType, iExponent, iTernaryCon
        integer :: iFirstParam, iSecondParam, iThirdParam, iSubParam
        integer :: iFirstParam2, iSecondParam2, iSubParam2, iTempParam
        integer :: s, c, a, iLinearCount
        real(8) :: dNu(nSiteCapacity)
        real(8) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)
        real(8) :: dTempValue, dTernaryWeight
        logical :: lTernaryWeighted

        nParamCon = iRegularParam(iParam,1)
        if (nParamCon <= 0) return

        dNu = 0D0
        dLinearValue = 0D0
        dLinearPower = 0D0
        dLinearCoeff = 0D0
        iLinearCount = 0
        iMixType = 0
        iFirstParam = 0
        iSecondParam = 0
        iThirdParam = 0
        iSubParam = 0
        iFirstParam2 = 0
        iSecondParam2 = 0
        iSubParam2 = 0
        lTernaryWeighted = .FALSE.

        do k = 2, nParamCon + 1
            c = MOD(iRegularParam(iParam,k), 10000)
            s = (iRegularParam(iParam,k) - c) / 10000
            if ((s <= 0).OR.(c <= 0)) return
            a = iSiteIndexIn(s,c)
            if ((a <= 0).OR.(a > nSiteCapacity)) return
            dNu(a) = dNu(a) + 1D0
        end do

        if ((iSUBLParamData(iParam,1) == 1).AND.(iSUBLParamData(iParam,3) == 2)) then
            iMixType = 2
            do k = 2, nParamCon + 1
                c = MOD(iRegularParam(iParam,k), 10000)
                s = (iRegularParam(iParam,k) - c) / 10000
                if (k == iSUBLParamData(iParam,2)) then
                    iFirstParam = c
                    iSubParam = s
                else if (k == iSUBLParamData(iParam,2) + 1) then
                    iSecondParam = c
                end if
            end do

            iExponent = iRegularParam(iParam,nParamCon+2)
            if (iExponent > 0) then
                iLinearCount = 1
                dLinearPower(iLinearCount) = DBLE(iExponent)
                a = iSiteIndexIn(iSubParam,iFirstParam)
                dLinearCoeff(iLinearCount,a) = 1D0
                a = iSiteIndexIn(iSubParam,iSecondParam)
                dLinearCoeff(iLinearCount,a) = -1D0
                dLinearValue(iLinearCount) = dSiteIn(iSubParam,iFirstParam) - dSiteIn(iSubParam,iSecondParam)
            end if
        else if ((iSUBLParamData(iParam,1) == 1).AND.(iSUBLParamData(iParam,3) == 3)) then
            iMixType = 3
            iTernaryCon = iRegularParam(iParam,nParamCon+2)
            lTernaryWeighted = (iRegularParam(iParam,nParamCon+3) == 1)

            do k = 2, nParamCon + 1
                c = MOD(iRegularParam(iParam,k), 10000)
                s = (iRegularParam(iParam,k) - c) / 10000
                if (k == iSUBLParamData(iParam,2) + iTernaryCon) then
                    iFirstParam = c
                    iSubParam = s
                else if (k == iSUBLParamData(iParam,2) + MOD(iTernaryCon + 1, 3)) then
                    iSecondParam = c
                else if (k == iSUBLParamData(iParam,2) + MOD(iTernaryCon + 2, 3)) then
                    iThirdParam = c
                end if
            end do

            if (lTernaryWeighted) then
                iLinearCount = 1
                dLinearPower(iLinearCount) = 1D0
                dTernaryWeight = dSiteIn(iSubParam,iFirstParam)
                dTempValue = 0D0
                do c = 1, nConstituentSublattice(iPhaseIDIn,iSubParam)
                    dTempValue = dTempValue + dSiteIn(iSubParam,c)
                end do
                dTernaryWeight = dTernaryWeight + &
                    (dTempValue - dSiteIn(iSubParam,iFirstParam) - &
                    dSiteIn(iSubParam,iSecondParam) - dSiteIn(iSubParam,iThirdParam)) / 3D0
                dLinearValue(iLinearCount) = dTernaryWeight

                do c = 1, nConstituentSublattice(iPhaseIDIn,iSubParam)
                    a = iSiteIndexIn(iSubParam,c)
                    dLinearCoeff(iLinearCount,a) = 1D0 / 3D0
                end do
                a = iSiteIndexIn(iSubParam,iFirstParam)
                dLinearCoeff(iLinearCount,a) = 1D0
                a = iSiteIndexIn(iSubParam,iSecondParam)
                dLinearCoeff(iLinearCount,a) = 0D0
                a = iSiteIndexIn(iSubParam,iThirdParam)
                dLinearCoeff(iLinearCount,a) = 0D0
            end if
        else if ((iSUBLParamData(iParam,1) == 2).AND.(iSUBLParamData(iParam,3) == 2).AND. &
            (iSUBLParamData(iParam,5) == 2)) then
            iMixType = 4
            do k = 2, nParamCon + 1
                c = MOD(iRegularParam(iParam,k), 10000)
                s = (iRegularParam(iParam,k) - c) / 10000
                if (k == iSUBLParamData(iParam,2)) then
                    iFirstParam = c
                    iSubParam = s
                else if (k == iSUBLParamData(iParam,2) + 1) then
                    iSecondParam = c
                else if (k == iSUBLParamData(iParam,4)) then
                    iFirstParam2 = c
                    iSubParam2 = s
                else if (k == iSUBLParamData(iParam,4) + 1) then
                    iSecondParam2 = c
                end if
            end do

            if (MOD(iRegularParam(iParam,nParamCon+2),2) == 0) then
                iExponent = iRegularParam(iParam,nParamCon+2) / 2
            else
                iExponent = (iRegularParam(iParam,nParamCon+2) - 1) / 2
                iTempParam = iFirstParam
                iFirstParam = iFirstParam2
                iFirstParam2 = iTempParam
                iTempParam = iSecondParam
                iSecondParam = iSecondParam2
                iSecondParam2 = iTempParam
                iTempParam = iSubParam
                iSubParam = iSubParam2
                iSubParam2 = iTempParam
            end if

            if (iExponent > 0) then
                iLinearCount = 1
                dLinearPower(iLinearCount) = DBLE(iExponent)
                a = iSiteIndexIn(iSubParam,iFirstParam)
                dLinearCoeff(iLinearCount,a) = 1D0
                a = iSiteIndexIn(iSubParam,iSecondParam)
                dLinearCoeff(iLinearCount,a) = -1D0
                dLinearValue(iLinearCount) = dSiteIn(iSubParam,iFirstParam) - dSiteIn(iSubParam,iSecondParam)
            end if
        else
            return
        end if

        if (iMixType > 0) then
            call AddProductLinearTermHessian(nSiteCapacity, dExcessGibbsParam(iParam), dNu, &
                iLinearCount, dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
                dHessianOut)
        end if

        return

    end subroutine AccumulateOneSUBLParameter

    subroutine AddProductLinearTermHessian(nSiteCapacity, dParameter, dNu, nLinear, dLinearValue, &
        dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, dHessianOut)

        implicit none

        integer, intent(in) :: nSiteCapacity, nLinear
        real(8), intent(in) :: dParameter, dNu(nSiteCapacity)
        real(8), intent(in) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(inout) :: dHessianOut(nSiteCapacity,nSiteCapacity)

        real(8) :: dValueTerm
        real(8) :: dGradientTerm(nSiteCapacity)
        real(8) :: dHessianTerm(nSiteCapacity,nSiteCapacity)

        call EvaluateProductLinearDerivatives(nSiteCapacity, dParameter, dNu, nLinear, &
            dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
            dValueTerm, dGradientTerm, dHessianTerm)

        dHessianOut = dHessianOut + dHessianTerm

        return

    end subroutine AddProductLinearTermHessian

    subroutine EvaluateProductLinearDerivatives(nSiteCapacity, dParameter, dNu, nLinear, &
        dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
        dValueTerm, dGradientTerm, dHessianTerm)

        implicit none

        integer, intent(in) :: nSiteCapacity, nLinear
        real(8), intent(in) :: dParameter, dNu(nSiteCapacity)
        real(8), intent(in) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(out) :: dValueTerm, dGradientTerm(nSiteCapacity)
        real(8), intent(out) :: dHessianTerm(nSiteCapacity,nSiteCapacity)

        integer :: a, b, q, s, c, iPower
        real(8) :: dProduct, dLinearProduct, dFactor, dSiteA, dSiteB
        real(8) :: dProductGradient(nSiteCapacity), dLinearGradient(nSiteCapacity)
        real(8) :: dProductHessian(nSiteCapacity,nSiteCapacity)
        real(8) :: dLinearHessian(nSiteCapacity,nSiteCapacity)
        real(8) :: dFactorGradient(nSiteCapacity), dFactorHessian(nSiteCapacity,nSiteCapacity)
        real(8) :: dNewLinearProduct, dNewLinearGradient(nSiteCapacity)
        real(8) :: dNewLinearHessian(nSiteCapacity,nSiteCapacity)

        dValueTerm = 0D0
        dGradientTerm = 0D0
        dHessianTerm = 0D0
        if (dParameter == 0D0) return

        dProduct = dParameter
        dProductGradient = 0D0
        dProductHessian = 0D0

        do a = 1, nSiteCapacity
            if (iSiteSubIn(a) <= 0) cycle
            if (dNu(a) == 0D0) cycle
            s = iSiteSubIn(a)
            c = iSiteConIn(a)
            dProduct = dProduct * DMAX1(dSiteIn(s,c), 1D-75)**dNu(a)
        end do

        do a = 1, nSiteCapacity
            if (iSiteSubIn(a) <= 0) cycle
            s = iSiteSubIn(a)
            c = iSiteConIn(a)
            dSiteA = DMAX1(dSiteIn(s,c), 1D-75)
            dProductGradient(a) = dProduct * dNu(a) / dSiteA

            do b = 1, nSiteCapacity
                if (iSiteSubIn(b) <= 0) cycle
                dSiteB = DMAX1(dSiteIn(iSiteSubIn(b),iSiteConIn(b)), 1D-75)
                dProductHessian(a,b) = dProduct * dNu(a) * dNu(b) / (dSiteA * dSiteB)
                if (a == b) dProductHessian(a,b) = dProductHessian(a,b) - &
                    dProduct * dNu(a) / (dSiteA * dSiteA)
            end do
        end do

        dLinearProduct = 1D0
        dLinearGradient = 0D0
        dLinearHessian = 0D0

        do q = 1, nLinear
            iPower = NINT(dLinearPower(q))
            if (iPower <= 0) cycle

            dFactor = dLinearValue(q)**iPower
            dFactorGradient = 0D0
            dFactorHessian = 0D0

            do a = 1, nSiteCapacity
                if (iSiteSubIn(a) <= 0) cycle
                dFactorGradient(a) = DBLE(iPower) * dLinearValue(q)**(iPower - 1) * dLinearCoeff(q,a)

                if (iPower > 1) then
                    do b = 1, nSiteCapacity
                        if (iSiteSubIn(b) <= 0) cycle
                        dFactorHessian(a,b) = DBLE(iPower * (iPower - 1)) * &
                            dLinearValue(q)**(iPower - 2) * dLinearCoeff(q,a) * dLinearCoeff(q,b)
                    end do
                end if
            end do

            dNewLinearProduct = dLinearProduct * dFactor
            dNewLinearGradient = dLinearGradient * dFactor + dLinearProduct * dFactorGradient
            dNewLinearHessian = dLinearHessian * dFactor

            do a = 1, nSiteCapacity
                if (iSiteSubIn(a) <= 0) cycle
                do b = 1, nSiteCapacity
                    if (iSiteSubIn(b) <= 0) cycle
                    dNewLinearHessian(a,b) = dNewLinearHessian(a,b) + &
                        dLinearGradient(a) * dFactorGradient(b) + &
                        dFactorGradient(a) * dLinearGradient(b) + &
                        dLinearProduct * dFactorHessian(a,b)
                end do
            end do

            dLinearProduct = dNewLinearProduct
            dLinearGradient = dNewLinearGradient
            dLinearHessian = dNewLinearHessian
        end do

        dValueTerm = dProduct * dLinearProduct
        dGradientTerm = dProductGradient * dLinearProduct + dProduct * dLinearGradient
        dHessianTerm = dProductHessian * dLinearProduct

        do a = 1, nSiteCapacity
            if (iSiteSubIn(a) <= 0) cycle
            do b = 1, nSiteCapacity
                if (iSiteSubIn(b) <= 0) cycle
                dHessianTerm(a,b) = dHessianTerm(a,b) + &
                    dProductGradient(a) * dLinearGradient(b) + &
                    dLinearGradient(a) * dProductGradient(b) + &
                    dProduct * dLinearHessian(a,b)
            end do
        end do

        return

    end subroutine EvaluateProductLinearDerivatives

    subroutine AddSUBOMFixedCorrection(iSolnIndexIn, nSiteCapacity, iOrdSiteIndex, iOrdSiteSub, iOrdSiteCon, &
        nOrdSite, dOrdSite, dOrdHessian)

        implicit none

        integer, intent(in) :: iSolnIndexIn, nSiteCapacity, nOrdSite
        integer, intent(in) :: iOrdSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iOrdSiteSub(nSiteCapacity), iOrdSiteCon(nSiteCapacity)
        real(8), intent(in) :: dOrdSite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(inout) :: dOrdHessian(nSiteCapacity,nSiteCapacity)

        integer :: iDisIndex, iOrdPhaseID, iDisPhaseID, nDisSite, nRandomSite
        integer :: so, sd, c, d, a, b, p, q, iBestDisSub, iBestCount, iCountDis
        integer :: iDisSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
        integer :: iDisSiteSub(nSiteCapacity), iDisSiteCon(nSiteCapacity)
        integer :: iRandomSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
        integer :: iRandomSiteSub(nSiteCapacity), iRandomSiteCon(nSiteCapacity)
        integer :: iOrdToDis(nMaxSublatticeSys)
        integer :: iOrdConToDisCon(nMaxSublatticeSys,nMaxConstituentSys)
        real(8) :: dGroupStoich(nMaxSublatticeSys), dCollapsedSite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8) :: dRandomSite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8) :: dDisHessian(nSiteCapacity,nSiteCapacity)
        real(8) :: dRandomHessian(nSiteCapacity,nSiteCapacity)
        real(8) :: dCollapseMap(nSiteCapacity,nSiteCapacity)
        real(8) :: dRandomMap(nSiteCapacity,nSiteCapacity)
        real(8) :: dWeight
        logical :: lSubset, lFound

        iDisIndex = 0
        if (allocated(iDisorderedPhase)) iDisIndex = iDisorderedPhase(iSolnIndexIn)
        if ((iDisIndex <= 0).OR.(iDisIndex > nSolnPhasesSys)) return

        iOrdPhaseID = iPhaseSublattice(iSolnIndexIn)
        iDisPhaseID = iPhaseSublattice(iDisIndex)
        if ((iOrdPhaseID <= 0).OR.(iDisPhaseID <= 0)) return

        iOrdToDis = 0
        iOrdConToDisCon = 0
        dGroupStoich = 0D0

        do so = 1, nSublatticePhase(iOrdPhaseID)
            iBestDisSub = 0
            iBestCount = nMaxConstituentSys + 1
            do sd = 1, nSublatticePhase(iDisPhaseID)
                lSubset = .TRUE.
                do c = 1, nConstituentSublattice(iOrdPhaseID,so)
                    lFound = .FALSE.
                    do d = 1, nConstituentSublattice(iDisPhaseID,sd)
                        if (TRIM(cConstituentNameSUB(iOrdPhaseID,so,c)) == &
                            TRIM(cConstituentNameSUB(iDisPhaseID,sd,d))) then
                            lFound = .TRUE.
                            exit
                        end if
                    end do
                    if (.NOT.lFound) lSubset = .FALSE.
                end do

                iCountDis = nConstituentSublattice(iDisPhaseID,sd)
                if (lSubset.AND.(iCountDis < iBestCount)) then
                    iBestDisSub = sd
                    iBestCount = iCountDis
                end if
            end do
            if (iBestDisSub == 0) return

            iOrdToDis(so) = iBestDisSub
            dGroupStoich(iBestDisSub) = dGroupStoich(iBestDisSub) + dStoichSublattice(iOrdPhaseID,so)

            do c = 1, nConstituentSublattice(iOrdPhaseID,so)
                do d = 1, nConstituentSublattice(iDisPhaseID,iBestDisSub)
                    if (TRIM(cConstituentNameSUB(iOrdPhaseID,so,c)) == &
                        TRIM(cConstituentNameSUB(iDisPhaseID,iBestDisSub,d))) then
                        iOrdConToDisCon(so,c) = d
                        exit
                    end if
                end do
                if (iOrdConToDisCon(so,c) == 0) return
            end do
        end do

        do sd = 1, nSublatticePhase(iDisPhaseID)
            if (dGroupStoich(sd) <= 0D0) return
        end do

        call BuildSiteIndex(iDisPhaseID, nSiteCapacity, iDisSiteIndex, iDisSiteSub, iDisSiteCon, nDisSite)
        call BuildSiteIndex(iOrdPhaseID, nSiteCapacity, iRandomSiteIndex, iRandomSiteSub, iRandomSiteCon, nRandomSite)
        if ((nDisSite > nSiteCapacity).OR.(nRandomSite > nSiteCapacity)) return

        dCollapsedSite = 0D0
        dCollapseMap = 0D0
        do so = 1, nSublatticePhase(iOrdPhaseID)
            sd = iOrdToDis(so)
            dWeight = dStoichSublattice(iOrdPhaseID,so) / dGroupStoich(sd)
            do c = 1, nConstituentSublattice(iOrdPhaseID,so)
                d = iOrdConToDisCon(so,c)
                dCollapsedSite(sd,d) = dCollapsedSite(sd,d) + dWeight * dOrdSite(so,c)
                a = iDisSiteIndex(sd,d)
                b = iOrdSiteIndex(so,c)
                dCollapseMap(a,b) = dCollapseMap(a,b) + dWeight
            end do
        end do

        dRandomSite = 0D0
        dRandomMap = 0D0
        do so = 1, nSublatticePhase(iOrdPhaseID)
            sd = iOrdToDis(so)
            do c = 1, nConstituentSublattice(iOrdPhaseID,so)
                d = iOrdConToDisCon(so,c)
                dRandomSite(so,c) = dCollapsedSite(sd,d)
                a = iRandomSiteIndex(so,c)
                do b = 1, nOrdSite
                    dRandomMap(a,b) = dCollapseMap(iDisSiteIndex(sd,d),b)
                end do
            end do
        end do

        call AccumulateRawSUBLSiteHessian(iDisIndex, iDisPhaseID, nSiteCapacity, iDisSiteIndex, &
            iDisSiteSub, iDisSiteCon, nDisSite, dCollapsedSite, dDisHessian)
        call AccumulateRawSUBLSiteHessian(iSolnIndexIn, iOrdPhaseID, nSiteCapacity, iRandomSiteIndex, &
            iRandomSiteSub, iRandomSiteCon, nRandomSite, dRandomSite, dRandomHessian)

        do a = 1, nOrdSite
            do b = 1, nOrdSite
                do p = 1, nDisSite
                    do q = 1, nDisSite
                        dOrdHessian(a,b) = dOrdHessian(a,b) + dCollapseMap(p,a) * dDisHessian(p,q) * &
                            dCollapseMap(q,b)
                    end do
                end do
                do p = 1, nRandomSite
                    do q = 1, nRandomSite
                        dOrdHessian(a,b) = dOrdHessian(a,b) - dRandomMap(p,a) * dRandomHessian(p,q) * &
                            dRandomMap(q,b)
                    end do
                end do
            end do
        end do

        return

    end subroutine AddSUBOMFixedCorrection

    subroutine ProjectSiteHessianToEndmembers(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, nSpeciesCapacity, &
        iSiteIndexIn, iSiteSubIn, iSiteConIn, nSiteLocal, nSpeciesLocal, dSiteIn, dSiteHessianIn, &
        dAmountHessianOut, dCompositionJacobianOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn, nSiteCapacity, nSpeciesCapacity
        integer, intent(in) :: nSiteLocal, nSpeciesLocal
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(in) :: dSiteHessianIn(nSiteCapacity,nSiteCapacity)
        real(8), intent(out) :: dAmountHessianOut(nSpeciesCapacity,nSpeciesCapacity)
        real(8), intent(out) :: dCompositionJacobianOut(nSpeciesCapacity,nSpeciesCapacity)

        integer :: i, j, a, b, s, c, iFirstLocal
        real(8) :: dBi, dBj, dChiJ

        dAmountHessianOut = 0D0
        dCompositionJacobianOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1

        do i = 1, nSpeciesLocal
            do j = 1, nSpeciesLocal
                do a = 1, nSiteLocal
                    s = iSiteSubIn(a)
                    c = iSiteConIn(a)
                    dBi = -dSiteIn(s,c)
                    if (iConstituentSublattice(iPhaseIDIn,s,i) == c) dBi = dBi + 1D0

                    do b = 1, nSiteLocal
                        dBj = -dSiteIn(iSiteSubIn(b),iSiteConIn(b))
                        if (iConstituentSublattice(iPhaseIDIn,iSiteSubIn(b),j) == iSiteConIn(b)) dBj = dBj + 1D0
                        dAmountHessianOut(i,j) = dAmountHessianOut(i,j) + dBi * dSiteHessianIn(a,b) * dBj

                        dChiJ = 0D0
                        if (iConstituentSublattice(iPhaseIDIn,iSiteSubIn(b),j) == iSiteConIn(b)) dChiJ = 1D0
                        dCompositionJacobianOut(i,j) = dCompositionJacobianOut(i,j) + &
                            dBi * dSiteHessianIn(a,b) * dChiJ
                    end do
                end do
            end do
        end do

        return

    end subroutine ProjectSiteHessianToEndmembers

    logical function HasMagneticTerms(iSolnIndexIn)

        implicit none

        integer, intent(in) :: iSolnIndexIn
        integer :: iFirstLocal, iLastLocal

        HasMagneticTerms = .FALSE.
        if (allocated(nMagParamPhase)) then
            if (nMagParamPhase(iSolnIndexIn) > nMagParamPhase(iSolnIndexIn-1)) HasMagneticTerms = .TRUE.
        end if
        if (allocated(dCoeffGibbsMagnetic)) then
            iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
            iLastLocal = nSpeciesPhase(iSolnIndexIn)
            if (MAXVAL(DABS(dCoeffGibbsMagnetic(iFirstLocal:iLastLocal,1:2))) > 0D0) HasMagneticTerms = .TRUE.
        end if

        return

    end function HasMagneticTerms

end subroutine CompHessianSUBL
