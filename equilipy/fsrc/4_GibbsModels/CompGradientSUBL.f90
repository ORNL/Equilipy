!> \brief Compute analytical fixed-state CEF scalar thermodynamics and site-gradients.
!!
!! \details This diagnostic/minimizer routine returns the same fixed-state
!! CEF scalars used by the SUBL Hessian diagnostic together with their first
!! derivatives in full redundant site-fraction coordinates.
!!
!! \param[in]  iSolnIndex      Absolute solution phase index.
!! \param[in]  nSiteDim        Caller-provided site-gradient capacity.
!! \param[out] dSiteGradient   Site-space gradient in internal G/(RT) units.
!! \param[out] dScalarGibbs    Scalar CEF Gibbs energy in internal G/(RT) units.
!! \param[out] dSiteGradientH  Site-space gradient in internal H/(RT) units.
!! \param[out] dSiteGradientS  Site-space gradient in internal S/R units.
!! \param[out] dSiteGradientCp Site-space gradient in internal Cp/R units.
!! \param[out] dScalarEnthalpy Scalar CEF enthalpy in internal H/(RT) units.
!! \param[out] dScalarEntropy  Scalar CEF entropy in internal S/R units.
!! \param[out] dScalarHeatCapacity Scalar CEF heat capacity in internal Cp/R units.
!! \param[out] nSiteOut        Number of active site variables written.
!! \param[out] iInfo           Zero on success; nonzero for unsupported or undersized calls.



subroutine CompGradientSUBL(iSolnIndex, nSiteDim, dSiteGradient, dScalarGibbs, &
    dSiteGradientH, dSiteGradientS, dSiteGradientCp, dScalarEnthalpy, dScalarEntropy, &
    dScalarHeatCapacity, nSiteOut, iInfo)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompGradientSUBL.f90
    !> \brief   Compute analytical fixed-state CEF scalar thermodynamics and site-gradients.
    !> \author  S.Y. Kwon
    !> \date    Jun. 24, 2026
    !> \sa      CompHessianSUBL.f90
    !> \sa      CompExcessGibbsEnergySUBL.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Corrected CEF chemical-potential and magnetic-temperature derivatives for consistent sublattice gradients.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details This routine evaluates the fixed-state CEF scalar Gibbs energy,
    !! enthalpy, entropy, heat capacity, and their first derivatives with respect
    !! to site fractions for SUBL, SUBLM, and SUBOM phases.  The scalars include
    !! reference endmember products, ideal configurational terms where applicable,
    !! nonmagnetic excess terms, and magnetic scalar terms.  For SUBOM phases,
    !! the same fixed ordered-state correction used by CompHessianSUBL is applied:
    !! \f$q^{fix}=q^{ord}+q^{dis}(Cy)-q^{ord,r}(My)\f$.
    !
    !
    ! Required input variables:
    ! =========================
    !
    !> \param[in] iSolnIndex  Absolute solution phase index.
    !> \param[in] nSiteDim    Caller-provided site-gradient capacity.
    !
    ! cSolnPhaseType             Solution model type.  SUBL, SUBLM, and SUBOM are accepted.
    ! dMolFraction               Current endmember fractions used to rebuild site fractions.
    ! dStdGibbsEnergy            Dimensionless standard endmember Gibbs energies.
    ! dStdEnthalpy/Entropy/HeatCapacity
    !                            Dimensionless standard endmember thermodynamic properties.
    ! dExcessGibbsParam/HParam/SParam/CpParam
    !                            Dimensionless nonmagnetic excess parameters.
    ! iRegularParam/iSUBLParamData
    !                            Runtime CEF parameter topology.
    ! iDisorderedPhase           Ordered-to-disordered phase mapping for SUBOM.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] dSiteGradient  Redundant site-space gradient in G/(RT) units.
    !> \param[out] dScalarGibbs   Scalar CEF Gibbs energy in G/(RT) units.
    !> \param[out] dSiteGradientH Redundant site-space gradient in H/(RT) units.
    !> \param[out] dSiteGradientS Redundant site-space gradient in S/R units.
    !> \param[out] dSiteGradientCp
    !!                            Redundant site-space gradient in Cp/R units.
    !> \param[out] dScalarEnthalpy
    !!                            Scalar CEF enthalpy in H/(RT) units.
    !> \param[out] dScalarEntropy Scalar CEF entropy in S/R units.
    !> \param[out] dScalarHeatCapacity
    !!                            Scalar CEF heat capacity in Cp/R units.
    !> \param[out] nSiteOut       Number of active site variables written.
    !> \param[out] iInfo          Zero on success.  1=unsupported phase, 2=bad topology,
    !!                             3=caller arrays too small, 4=unsupported CEF parameter topology.
    !
    ! dSiteFraction              Rebuilt for iSolnIndex from dMolFraction.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! BuildSiteIndex                   Builds redundant site-coordinate indexing.
    ! BuildNormalizedSiteFractions     Builds normalized site fractions from endmember fractions.
    ! ValidateSUBLParameterTopology    Refuses CEF parameter topologies this analytical path cannot evaluate.
    ! AccumulateRawSUBLGradient        Adds reference, ideal, excess, and magnetic CEF scalars/gradients.
    ! AccumulateMagneticSiteGradient   Adds scalar magnetic first derivatives.
    ! AddSUBOMFixedGradient            Applies fixed ordered-state order/disorder scalar/gradient mapping.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! SubMinSiteFractionCEF            Uses this scalar and gradient for CEF site-fraction minimization.
    ! Python diagnostics/devop tests   Validate analytical first derivatives against finite differences.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Returned values use internal dimensionless conventions: G/(RT), H/(RT), S/R, and Cp/R.
    ! - The site-space gradient is redundant.  Callers that eliminate per-sublattice normalization should
    !   subtract the reference-constituent gradient on each sublattice.
    ! - SUBOM uses the fixed ordered-state scalar consistent with CompHessianSUBL.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO

    implicit none

    integer, intent(in)  :: iSolnIndex, nSiteDim
    real(8), intent(out) :: dSiteGradient(nSiteDim)
    real(8), intent(out) :: dScalarGibbs
    real(8), intent(out) :: dSiteGradientH(nSiteDim), dSiteGradientS(nSiteDim)
    real(8), intent(out) :: dSiteGradientCp(nSiteDim)
    real(8), intent(out) :: dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity
    integer, intent(out) :: nSiteOut, iInfo

    integer :: iPhaseID
    integer :: iSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
    integer :: iSiteSub(nSiteDim), iSiteCon(nSiteDim)
    real(8) :: dWorkSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8) :: dWorkGradientG(nSiteDim), dWorkGradientH(nSiteDim)
    real(8) :: dWorkGradientS(nSiteDim), dWorkGradientCp(nSiteDim)

    dSiteGradient = 0D0
    dSiteGradientH = 0D0
    dSiteGradientS = 0D0
    dSiteGradientCp = 0D0
    dScalarGibbs = 0D0
    dScalarEnthalpy = 0D0
    dScalarEntropy = 0D0
    dScalarHeatCapacity = 0D0
    nSiteOut = 0
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

    call BuildSiteIndex(iPhaseID, nSiteDim, iSiteIndex, iSiteSub, iSiteCon, nSiteOut)
    if (nSiteOut <= 0) then
        iInfo = 2
        return
    end if
    if (nSiteOut > nSiteDim) then
        iInfo = 3
        return
    end if

    call ValidateSUBLParameterTopology(iSolnIndex, iInfo)
    if (iInfo /= 0) return
    if ((TRIM(cSolnPhaseType(iSolnIndex)) == 'SUBOM').AND.allocated(iDisorderedPhase)) then
        if ((iDisorderedPhase(iSolnIndex) > 0).AND.(iDisorderedPhase(iSolnIndex) <= nSolnPhasesSys)) then
            call ValidateSUBLParameterTopology(iDisorderedPhase(iSolnIndex), iInfo)
            if (iInfo /= 0) return
        end if
    end if

    call BuildNormalizedSiteFractions(iSolnIndex, iPhaseID, dWorkSite)
    dSiteFraction(iPhaseID,:,:) = dWorkSite(:,:)

    call AccumulateRawSUBLGradient(iSolnIndex, iPhaseID, nSiteDim, iSiteIndex, &
        iSiteSub, iSiteCon, nSiteOut, dWorkSite, dScalarGibbs, dScalarEnthalpy, &
        dScalarEntropy, dScalarHeatCapacity, dWorkGradientG, dWorkGradientH, &
        dWorkGradientS, dWorkGradientCp)

    if (TRIM(cSolnPhaseType(iSolnIndex)) == 'SUBOM') then
        call AddSUBOMFixedGradient(iSolnIndex, nSiteDim, iSiteIndex, iSiteSub, iSiteCon, &
            nSiteOut, dWorkSite, dScalarGibbs, dScalarEnthalpy, dScalarEntropy, &
            dScalarHeatCapacity, dWorkGradientG, dWorkGradientH, dWorkGradientS, dWorkGradientCp)
    end if

    dSiteGradient(1:nSiteOut) = dWorkGradientG(1:nSiteOut)
    dSiteGradientH(1:nSiteOut) = dWorkGradientH(1:nSiteOut)
    dSiteGradientS(1:nSiteOut) = dWorkGradientS(1:nSiteOut)
    dSiteGradientCp(1:nSiteOut) = dWorkGradientCp(1:nSiteOut)

    return

contains

    subroutine ValidateSUBLParameterTopology(iSolnIndexIn, iInfoOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn
        integer, intent(out) :: iInfoOut

        integer :: iParam

        iInfoOut = 0
        do iParam = nParamPhase(iSolnIndexIn-1) + 1, nParamPhase(iSolnIndexIn)
            if (iRegularParam(iParam,1) <= 0) cycle
            if ((iSUBLParamData(iParam,1) == 1).AND.(iSUBLParamData(iParam,3) == 2)) cycle
            if ((iSUBLParamData(iParam,1) == 1).AND.(iSUBLParamData(iParam,3) == 3)) cycle
            if ((iSUBLParamData(iParam,1) == 2).AND.(iSUBLParamData(iParam,3) == 2).AND. &
                (iSUBLParamData(iParam,5) == 2)) cycle
            iInfoOut = 4
            return
        end do

        return

    end subroutine ValidateSUBLParameterTopology

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

    subroutine AccumulateRawSUBLGradient(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
        iSiteSubIn, iSiteConIn, nSiteLocal, dSiteIn, dValueGOut, dValueHOut, dValueSOut, dValueCpOut, &
        dGradientGOut, dGradientHOut, dGradientSOut, dGradientCpOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn, nSiteCapacity, nSiteLocal
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dValueGOut, dValueHOut, dValueSOut, dValueCpOut
        real(8), intent(out) :: dGradientGOut(nSiteCapacity), dGradientHOut(nSiteCapacity)
        real(8), intent(out) :: dGradientSOut(nSiteCapacity), dGradientCpOut(nSiteCapacity)

        integer :: iFirstLocal, iLastLocal, i, m, s, a, c
        real(8) :: dProduct, dSiteA

        dValueGOut = 0D0
        dValueHOut = 0D0
        dValueSOut = 0D0
        dValueCpOut = 0D0
        dGradientGOut = 0D0
        dGradientHOut = 0D0
        dGradientSOut = 0D0
        dGradientCpOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnIndexIn)

        do i = iFirstLocal, iLastLocal
            m = i - iFirstLocal + 1
            dProduct = 1D0
            do s = 1, nSublatticePhase(iPhaseIDIn)
                c = iConstituentSublattice(iPhaseIDIn,s,m)
                dProduct = dProduct * DMAX1(dSiteIn(s,c), 1D-75)
            end do

            dValueGOut = dValueGOut + dStdGibbsEnergy(i) * dProduct
            dValueHOut = dValueHOut + dStdEnthalpy(i) * dProduct
            dValueSOut = dValueSOut + dStdEntropy(i) * dProduct
            dValueCpOut = dValueCpOut + dStdHeatCapacity(i) * dProduct

            do s = 1, nSublatticePhase(iPhaseIDIn)
                c = iConstituentSublattice(iPhaseIDIn,s,m)
                a = iSiteIndexIn(s,c)
                dSiteA = DMAX1(dSiteIn(s,c), 1D-75)
                dGradientGOut(a) = dGradientGOut(a) + dStdGibbsEnergy(i) * dProduct / dSiteA
                dGradientHOut(a) = dGradientHOut(a) + dStdEnthalpy(i) * dProduct / dSiteA
                dGradientSOut(a) = dGradientSOut(a) + dStdEntropy(i) * dProduct / dSiteA
                dGradientCpOut(a) = dGradientCpOut(a) + dStdHeatCapacity(i) * dProduct / dSiteA
            end do
        end do

        do a = 1, nSiteLocal
            s = iSiteSubIn(a)
            c = iSiteConIn(a)
            dSiteA = DMAX1(dSiteIn(s,c), 1D-75)
            dValueGOut = dValueGOut + dStoichSublattice(iPhaseIDIn,s) * dSiteA * DLOG(dSiteA)
            dGradientGOut(a) = dGradientGOut(a) + dStoichSublattice(iPhaseIDIn,s) * (DLOG(dSiteA) + 1D0)
            dValueSOut = dValueSOut - dStoichSublattice(iPhaseIDIn,s) * dSiteA * DLOG(dSiteA)
            dGradientSOut(a) = dGradientSOut(a) - dStoichSublattice(iPhaseIDIn,s) * (DLOG(dSiteA) + 1D0)
        end do

        do i = nParamPhase(iSolnIndexIn-1) + 1, nParamPhase(iSolnIndexIn)
            call AccumulateOneSUBLParameterGradients(i, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
                iSiteSubIn, iSiteConIn, dSiteIn, dValueGOut, dValueHOut, dValueSOut, dValueCpOut, &
                dGradientGOut, dGradientHOut, dGradientSOut, dGradientCpOut)
        end do

        if (HasMagneticTerms(iSolnIndexIn)) then
            call AccumulateMagneticSiteGradient(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
                iSiteSubIn, iSiteConIn, nSiteLocal, dSiteIn, dValueGOut, dValueHOut, dValueSOut, &
                dValueCpOut, dGradientGOut, dGradientHOut, dGradientSOut, dGradientCpOut)
        end if

        return

    end subroutine AccumulateRawSUBLGradient

    subroutine AccumulateMagneticSiteGradient(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
        iSiteSubIn, iSiteConIn, nSiteLocal, dSiteIn, dValueGOut, dValueHOut, dValueSOut, dValueCpOut, &
        dGradientGOut, dGradientHOut, dGradientSOut, dGradientCpOut)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn, nSiteCapacity, nSiteLocal
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(inout) :: dValueGOut, dValueHOut, dValueSOut, dValueCpOut
        real(8), intent(inout) :: dGradientGOut(nSiteCapacity), dGradientHOut(nSiteCapacity)
        real(8), intent(inout) :: dGradientSOut(nSiteCapacity), dGradientCpOut(nSiteCapacity)

        integer :: a
        real(8) :: dTCrit, dB
        real(8) :: dMagValue(4), dMagTCrit(4), dMagB(4)
        real(8) :: dTCritGrad(nSiteCapacity), dBGrad(nSiteCapacity)

        if (.NOT.HasMagneticTerms(iSolnIndexIn)) return
        if ((TRIM(cSolnPhaseType(iSolnIndexIn)) /= 'SUBLM').AND. &
            (TRIM(cSolnPhaseType(iSolnIndexIn)) /= 'SUBOM')) return

        call AccumulateMagneticVariableGradient(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, &
            iSiteIndexIn, iSiteSubIn, iSiteConIn, dSiteIn, dTCrit, dB, dTCritGrad, dBGrad)
        call EvaluateMagneticScalarFirstDerivatives(iSolnIndexIn, dTCrit, dB, dMagValue, dMagTCrit, dMagB)

        dValueGOut = dValueGOut + dMagValue(1)
        dValueHOut = dValueHOut + dMagValue(2)
        dValueSOut = dValueSOut + dMagValue(3)
        dValueCpOut = dValueCpOut + dMagValue(4)
        do a = 1, nSiteLocal
            dGradientGOut(a) = dGradientGOut(a) + dMagTCrit(1) * dTCritGrad(a) + dMagB(1) * dBGrad(a)
            dGradientHOut(a) = dGradientHOut(a) + dMagTCrit(2) * dTCritGrad(a) + dMagB(2) * dBGrad(a)
            dGradientSOut(a) = dGradientSOut(a) + dMagTCrit(3) * dTCritGrad(a) + dMagB(3) * dBGrad(a)
            dGradientCpOut(a) = dGradientCpOut(a) + dMagTCrit(4) * dTCritGrad(a) + dMagB(4) * dBGrad(a)
        end do

        return

    end subroutine AccumulateMagneticSiteGradient

    subroutine AccumulateMagneticVariableGradient(iSolnIndexIn, iPhaseIDIn, nSiteCapacity, &
        iSiteIndexIn, iSiteSubIn, iSiteConIn, dSiteIn, dTCrit, dB, dTCritGrad, dBGrad)

        implicit none

        integer, intent(in) :: iSolnIndexIn, iPhaseIDIn, nSiteCapacity
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dTCrit, dB, dTCritGrad(nSiteCapacity), dBGrad(nSiteCapacity)

        integer :: iFirstLocal, iLastLocal, i, m, s, c, a, k, iParam, nParamCon, iExponent
        real(8) :: dNu(nSiteCapacity)
        real(8) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)

        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnIndexIn)
        dTCrit = 0D0
        dB = 0D0
        dTCritGrad = 0D0
        dBGrad = 0D0

        do i = iFirstLocal, iLastLocal
            m = i - iFirstLocal + 1
            dNu = 0D0
            do s = 1, nSublatticePhase(iPhaseIDIn)
                c = iConstituentSublattice(iPhaseIDIn,s,m)
                a = iSiteIndexIn(s,c)
                if (a > 0) dNu(a) = dNu(a) + 1D0
            end do
            call AccumulateProductLinearTermGradient(nSiteCapacity, dCoeffGibbsMagnetic(i,1), dNu, &
                0, dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
                dTCrit, dTCritGrad)
            call AccumulateProductLinearTermGradient(nSiteCapacity, dCoeffGibbsMagnetic(i,2), dNu, &
                0, dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
                dB, dBGrad)
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
            call AccumulateProductLinearTermGradient(nSiteCapacity, dMagneticParam(iParam,1), dNu, &
                MERGE(1, 0, iExponent > 0), dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, &
                iSiteSubIn, iSiteConIn, dTCrit, dTCritGrad)
            call AccumulateProductLinearTermGradient(nSiteCapacity, dMagneticParam(iParam,2), dNu, &
                MERGE(1, 0, iExponent > 0), dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, &
                iSiteSubIn, iSiteConIn, dB, dBGrad)
        end do

        return

    end subroutine AccumulateMagneticVariableGradient

    subroutine EvaluateMagneticScalarFirstDerivatives(iSolnIndexIn, dTCritIn, dBIn, dMagValue, dMagTCrit, dMagB)

        implicit none

        integer, intent(in) :: iSolnIndexIn
        real(8), intent(in) :: dTCritIn, dBIn
        real(8), intent(out) :: dMagValue(4), dMagTCrit(4), dMagB(4)

        integer :: iFirstLocal
        real(8) :: dTCritLocal, dBLocal, dTCritSign, dBSign
        real(8) :: dStructureFactor, dP, dInvPMinusOne
        real(8) :: dTau, dD, dA, dBTemp, dC, dF, dFp, dFdp, dFtp
        real(8) :: dQp, dQdp, dQtp
        real(8) :: dLogMoment

        dMagValue = 0D0
        dMagTCrit = 0D0
        dMagB = 0D0

        if (iSolnIndexIn <= 0) return

        iFirstLocal = nSpeciesPhase(iSolnIndexIn-1) + 1
        dStructureFactor = dCoeffGibbsMagnetic(iFirstLocal,3)
        dP = dCoeffGibbsMagnetic(iFirstLocal,4)
        if (dP == 0D0) return

        dTCritLocal = dTCritIn
        dBLocal = dBIn
        dTCritSign = 1D0
        dBSign = 1D0
        if (dBLocal < 0D0) then
            dBLocal = -dBLocal * dStructureFactor
            dBSign = -dStructureFactor
        end if
        if (dTCritLocal < 0D0) then
            dTCritLocal = -dTCritLocal * dStructureFactor
            dTCritSign = -dStructureFactor
        end if
        if ((dTCritLocal == 0D0).OR.(dBLocal <= -1D0)) return

        dTau = dTemperature / dTCritLocal
        dInvPMinusOne = 1D0 / dP - 1D0
        dD = (518D0/1125D0) + (11692D0/15975D0) * dInvPMinusOne

        if (dTau > 1D0) then
            dA = dTau**(-5)
            dBTemp = dA**3
            dC = dA * dA * dBTemp

            dF = -(dA/10D0 + dBTemp/315D0 + dC/1500D0) / dD
            dFp = (1D0 / (dD * dTCritLocal)) * (dA/2D0 + dBTemp/21D0 + dC/60D0)
            dFdp = -(1D0 / (dD * dTCritLocal)) * &
                (3D0*dA + 16D0*dBTemp/21D0 + 13D0*dC/30D0)
            dFtp = +(1D0 / (dD * dTCritLocal)) * &
                (21D0*dA + 272D0*dBTemp/21D0 + 117D0*dC/10D0)
            dQp = -dFp
            dQdp = -dFdp
            dQtp = -dFtp
        else
            dA = dTau**3
            dBTemp = dA**3
            dC = dA * dA * dBTemp

            dF = 1D0 - (79D0/(140D0*dP*dTau) + &
                (474D0/497D0)*dInvPMinusOne*(dA/6D0 + dBTemp/135D0 + dC/600D0)) / dD
            dFp = (79D0/(140D0*dP*dTemperature) - &
                (474D0/(497D0*dTCritLocal))*dInvPMinusOne*(dA/2D0 + dBTemp/15D0 + dC/40D0)) / dD
            dFdp = -(79D0/(70D0*dP*dTemperature) + &
                (474D0/(497D0*dTCritLocal))*dInvPMinusOne*(dA + 8D0*dBTemp/15D0 + 7D0*dC/20D0)) / dD
            dFtp = (237D0/(70D0*dP*dTemperature) - &
                (474D0/(497D0*dTCritLocal))*dInvPMinusOne*(dA + 56D0*dBTemp/15D0 + &
                91D0*dC/20D0)) / dD
            dQp = -dFp
            dQdp = -dFdp
            dQtp = -dFtp
        end if

        dLogMoment = DLOG(1D0 + dBLocal)
        dMagValue(1) = dLogMoment * dF
        dMagValue(2) = -dLogMoment * dTCritLocal * dFp
        dMagValue(3) = -dLogMoment * (dF + dTCritLocal * dFp)
        dMagValue(4) = -dLogMoment * dTCritLocal * (2D0*dFp + dFdp)

        dMagTCrit(1) = dLogMoment * dQp * dTCritSign
        dMagTCrit(2) = -dLogMoment * (dQp + dQdp) * dTCritSign
        dMagTCrit(3) = -dLogMoment * (2D0*dQp + dQdp) * dTCritSign
        dMagTCrit(4) = -dLogMoment * (2D0*dQp + 4D0*dQdp + dQtp) * dTCritSign

        dMagB(1) = dF * dBSign / (1D0 + dBLocal)
        dMagB(2) = -dTCritLocal * dFp * dBSign / (1D0 + dBLocal)
        dMagB(3) = -(dF + dTCritLocal * dFp) * dBSign / (1D0 + dBLocal)
        dMagB(4) = -dTCritLocal * (2D0*dFp + dFdp) * dBSign / (1D0 + dBLocal)

        return

    end subroutine EvaluateMagneticScalarFirstDerivatives

    subroutine AccumulateOneSUBLParameterGradients(iParam, iPhaseIDIn, nSiteCapacity, iSiteIndexIn, &
        iSiteSubIn, iSiteConIn, dSiteIn, dValueGOut, dValueHOut, dValueSOut, dValueCpOut, &
        dGradientGOut, dGradientHOut, dGradientSOut, dGradientCpOut)

        implicit none

        integer, intent(in) :: iParam, iPhaseIDIn, nSiteCapacity
        integer, intent(in) :: iSiteIndexIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(inout) :: dValueGOut, dValueHOut, dValueSOut, dValueCpOut
        real(8), intent(inout) :: dGradientGOut(nSiteCapacity), dGradientHOut(nSiteCapacity)
        real(8), intent(inout) :: dGradientSOut(nSiteCapacity), dGradientCpOut(nSiteCapacity)

        integer :: k, nParamCon, iMixType, iExponent, iTernaryCon
        integer :: iFirstParam, iSecondParam, iThirdParam, iSubParam
        integer :: iFirstParam2, iSecondParam2, iSubParam2, iTempParam
        integer :: s, c, a, iLinearCount
        real(8) :: dNu(nSiteCapacity)
        real(8) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)
        real(8) :: dBaseValue, dBaseGradient(nSiteCapacity)
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
            call EvaluateProductLinearFirstDerivatives(nSiteCapacity, 1D0, dNu, iLinearCount, &
                dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
                dBaseValue, dBaseGradient)

            dValueGOut = dValueGOut + dExcessGibbsParam(iParam) * dBaseValue
            dValueHOut = dValueHOut + dExcessHParam(iParam) * dBaseValue
            dValueSOut = dValueSOut + dExcessSParam(iParam) * dBaseValue
            dValueCpOut = dValueCpOut + dExcessCpParam(iParam) * dBaseValue

            dGradientGOut = dGradientGOut + dExcessGibbsParam(iParam) * dBaseGradient
            dGradientHOut = dGradientHOut + dExcessHParam(iParam) * dBaseGradient
            dGradientSOut = dGradientSOut + dExcessSParam(iParam) * dBaseGradient
            dGradientCpOut = dGradientCpOut + dExcessCpParam(iParam) * dBaseGradient
        end if

        return

    end subroutine AccumulateOneSUBLParameterGradients

    subroutine AccumulateProductLinearTermGradient(nSiteCapacity, dParameter, dNu, nLinear, &
        dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
        dValueOut, dGradientOut)

        implicit none

        integer, intent(in) :: nSiteCapacity, nLinear
        real(8), intent(in) :: dParameter, dNu(nSiteCapacity)
        real(8), intent(in) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(inout) :: dValueOut, dGradientOut(nSiteCapacity)

        real(8) :: dValueTerm
        real(8) :: dGradientTerm(nSiteCapacity)

        call EvaluateProductLinearFirstDerivatives(nSiteCapacity, dParameter, dNu, nLinear, &
            dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
            dValueTerm, dGradientTerm)

        dValueOut = dValueOut + dValueTerm
        dGradientOut = dGradientOut + dGradientTerm

        return

    end subroutine AccumulateProductLinearTermGradient

    subroutine EvaluateProductLinearFirstDerivatives(nSiteCapacity, dParameter, dNu, nLinear, &
        dLinearValue, dLinearPower, dLinearCoeff, dSiteIn, iSiteSubIn, iSiteConIn, &
        dValueTerm, dGradientTerm)

        implicit none

        integer, intent(in) :: nSiteCapacity, nLinear
        real(8), intent(in) :: dParameter, dNu(nSiteCapacity)
        real(8), intent(in) :: dLinearValue(2), dLinearPower(2), dLinearCoeff(2,nSiteCapacity)
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iSiteSubIn(nSiteCapacity), iSiteConIn(nSiteCapacity)
        real(8), intent(out) :: dValueTerm, dGradientTerm(nSiteCapacity)

        integer :: a, q, s, c, iPower
        real(8) :: dProduct, dLinearProduct, dFactor, dSiteA
        real(8) :: dProductGradient(nSiteCapacity), dLinearGradient(nSiteCapacity)
        real(8) :: dFactorGradient(nSiteCapacity)
        real(8) :: dNewLinearProduct, dNewLinearGradient(nSiteCapacity)

        dValueTerm = 0D0
        dGradientTerm = 0D0
        if (dParameter == 0D0) return

        dProduct = dParameter
        dProductGradient = 0D0

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
        end do

        dLinearProduct = 1D0
        dLinearGradient = 0D0

        do q = 1, nLinear
            iPower = NINT(dLinearPower(q))
            if (iPower <= 0) cycle

            dFactor = dLinearValue(q)**iPower
            dFactorGradient = 0D0

            do a = 1, nSiteCapacity
                if (iSiteSubIn(a) <= 0) cycle
                dFactorGradient(a) = DBLE(iPower) * dLinearValue(q)**(iPower - 1) * dLinearCoeff(q,a)
            end do

            dNewLinearProduct = dLinearProduct * dFactor
            dNewLinearGradient = dLinearGradient * dFactor + dLinearProduct * dFactorGradient

            dLinearProduct = dNewLinearProduct
            dLinearGradient = dNewLinearGradient
        end do

        dValueTerm = dProduct * dLinearProduct
        dGradientTerm = dProductGradient * dLinearProduct + dProduct * dLinearGradient

        return

    end subroutine EvaluateProductLinearFirstDerivatives

    subroutine AddSUBOMFixedGradient(iSolnIndexIn, nSiteCapacity, iOrdSiteIndex, iOrdSiteSub, iOrdSiteCon, &
        nOrdSite, dOrdSite, dOrdValueG, dOrdValueH, dOrdValueS, dOrdValueCp, &
        dOrdGradientG, dOrdGradientH, dOrdGradientS, dOrdGradientCp)

        implicit none

        integer, intent(in) :: iSolnIndexIn, nSiteCapacity, nOrdSite
        integer, intent(in) :: iOrdSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
        integer, intent(in) :: iOrdSiteSub(nSiteCapacity), iOrdSiteCon(nSiteCapacity)
        real(8), intent(in) :: dOrdSite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(inout) :: dOrdValueG, dOrdValueH, dOrdValueS, dOrdValueCp
        real(8), intent(inout) :: dOrdGradientG(nSiteCapacity), dOrdGradientH(nSiteCapacity)
        real(8), intent(inout) :: dOrdGradientS(nSiteCapacity), dOrdGradientCp(nSiteCapacity)

        integer :: iDisIndex, iOrdPhaseID, iDisPhaseID, nDisSite, nRandomSite
        integer :: so, sd, c, d, a, b, p, iBestDisSub, iBestCount, iCountDis
        integer :: iDisSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
        integer :: iDisSiteSub(nSiteCapacity), iDisSiteCon(nSiteCapacity)
        integer :: iRandomSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
        integer :: iRandomSiteSub(nSiteCapacity), iRandomSiteCon(nSiteCapacity)
        integer :: iOrdToDis(nMaxSublatticeSys)
        integer :: iOrdConToDisCon(nMaxSublatticeSys,nMaxConstituentSys)
        real(8) :: dGroupStoich(nMaxSublatticeSys), dCollapsedSite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8) :: dRandomSite(nMaxSublatticeSys,nMaxConstituentSys)
        real(8) :: dDisValueG, dDisValueH, dDisValueS, dDisValueCp
        real(8) :: dRandomValueG, dRandomValueH, dRandomValueS, dRandomValueCp
        real(8) :: dDisGradientG(nSiteCapacity), dDisGradientH(nSiteCapacity)
        real(8) :: dDisGradientS(nSiteCapacity), dDisGradientCp(nSiteCapacity)
        real(8) :: dRandomGradientG(nSiteCapacity), dRandomGradientH(nSiteCapacity)
        real(8) :: dRandomGradientS(nSiteCapacity), dRandomGradientCp(nSiteCapacity)
        real(8) :: dOrdTCrit, dOrdB, dDisTCrit, dDisB, dRandomTCrit, dRandomB
        real(8) :: dPartTCrit, dPartB, dPartTCritGrad, dPartBGrad
        real(8) :: dOrdMagValue(4), dOrdMagTCrit(4), dOrdMagB(4)
        real(8) :: dDisMagValue(4), dDisMagTCrit(4), dDisMagB(4)
        real(8) :: dRandomMagValue(4), dRandomMagTCrit(4), dRandomMagB(4)
        real(8) :: dPartMagValue(4), dPartMagTCrit(4), dPartMagB(4)
        real(8) :: dOrdTCritGrad(nSiteCapacity), dOrdBGrad(nSiteCapacity)
        real(8) :: dDisTCritGrad(nSiteCapacity), dDisBGrad(nSiteCapacity)
        real(8) :: dRandomTCritGrad(nSiteCapacity), dRandomBGrad(nSiteCapacity)
        real(8) :: dCollapseMap(nSiteCapacity,nSiteCapacity)
        real(8) :: dRandomMap(nSiteCapacity,nSiteCapacity)
        real(8) :: dWeight
        logical :: lSubset, lFound
        logical :: lMagneticCorrection

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

        call AccumulateRawSUBLGradient(iDisIndex, iDisPhaseID, nSiteCapacity, iDisSiteIndex, &
            iDisSiteSub, iDisSiteCon, nDisSite, dCollapsedSite, dDisValueG, dDisValueH, dDisValueS, &
            dDisValueCp, dDisGradientG, dDisGradientH, dDisGradientS, dDisGradientCp)
        call AccumulateRawSUBLGradient(iSolnIndexIn, iOrdPhaseID, nSiteCapacity, iRandomSiteIndex, &
            iRandomSiteSub, iRandomSiteCon, nRandomSite, dRandomSite, dRandomValueG, dRandomValueH, &
            dRandomValueS, dRandomValueCp, dRandomGradientG, dRandomGradientH, dRandomGradientS, &
            dRandomGradientCp)

        dOrdValueG = dOrdValueG + dDisValueG - dRandomValueG
        dOrdValueH = dOrdValueH + dDisValueH - dRandomValueH
        dOrdValueS = dOrdValueS + dDisValueS - dRandomValueS
        dOrdValueCp = dOrdValueCp + dDisValueCp - dRandomValueCp

        do a = 1, nOrdSite
            do p = 1, nDisSite
                dOrdGradientG(a) = dOrdGradientG(a) + dCollapseMap(p,a) * dDisGradientG(p)
                dOrdGradientH(a) = dOrdGradientH(a) + dCollapseMap(p,a) * dDisGradientH(p)
                dOrdGradientS(a) = dOrdGradientS(a) + dCollapseMap(p,a) * dDisGradientS(p)
                dOrdGradientCp(a) = dOrdGradientCp(a) + dCollapseMap(p,a) * dDisGradientCp(p)
            end do
            do p = 1, nRandomSite
                dOrdGradientG(a) = dOrdGradientG(a) - dRandomMap(p,a) * dRandomGradientG(p)
                dOrdGradientH(a) = dOrdGradientH(a) - dRandomMap(p,a) * dRandomGradientH(p)
                dOrdGradientS(a) = dOrdGradientS(a) - dRandomMap(p,a) * dRandomGradientS(p)
                dOrdGradientCp(a) = dOrdGradientCp(a) - dRandomMap(p,a) * dRandomGradientCp(p)
            end do
        end do

        lMagneticCorrection = HasMagneticTerms(iSolnIndexIn).OR.HasMagneticTerms(iDisIndex)
        if (lMagneticCorrection) then
            call AccumulateMagneticVariableGradient(iSolnIndexIn, iOrdPhaseID, nSiteCapacity, &
                iOrdSiteIndex, iOrdSiteSub, iOrdSiteCon, dOrdSite, dOrdTCrit, dOrdB, &
                dOrdTCritGrad, dOrdBGrad)
            call AccumulateMagneticVariableGradient(iDisIndex, iDisPhaseID, nSiteCapacity, &
                iDisSiteIndex, iDisSiteSub, iDisSiteCon, dCollapsedSite, dDisTCrit, dDisB, &
                dDisTCritGrad, dDisBGrad)
            call AccumulateMagneticVariableGradient(iSolnIndexIn, iOrdPhaseID, nSiteCapacity, &
                iRandomSiteIndex, iRandomSiteSub, iRandomSiteCon, dRandomSite, dRandomTCrit, dRandomB, &
                dRandomTCritGrad, dRandomBGrad)

            dPartTCrit = dOrdTCrit + dDisTCrit - dRandomTCrit
            dPartB = dOrdB + dDisB - dRandomB
            call EvaluateMagneticScalarFirstDerivatives(iSolnIndexIn, dOrdTCrit, dOrdB, &
                dOrdMagValue, dOrdMagTCrit, dOrdMagB)
            call EvaluateMagneticScalarFirstDerivatives(iDisIndex, dDisTCrit, dDisB, &
                dDisMagValue, dDisMagTCrit, dDisMagB)
            call EvaluateMagneticScalarFirstDerivatives(iSolnIndexIn, dRandomTCrit, dRandomB, &
                dRandomMagValue, dRandomMagTCrit, dRandomMagB)
            call EvaluateMagneticScalarFirstDerivatives(iSolnIndexIn, dPartTCrit, dPartB, &
                dPartMagValue, dPartMagTCrit, dPartMagB)

            dOrdValueG = dOrdValueG + dPartMagValue(1) - dOrdMagValue(1) - dDisMagValue(1) + dRandomMagValue(1)
            dOrdValueH = dOrdValueH + dPartMagValue(2) - dOrdMagValue(2) - dDisMagValue(2) + dRandomMagValue(2)
            dOrdValueS = dOrdValueS + dPartMagValue(3) - dOrdMagValue(3) - dDisMagValue(3) + dRandomMagValue(3)
            dOrdValueCp = dOrdValueCp + dPartMagValue(4) - dOrdMagValue(4) - dDisMagValue(4) + dRandomMagValue(4)

            do a = 1, nOrdSite
                dPartTCritGrad = dOrdTCritGrad(a)
                dPartBGrad = dOrdBGrad(a)
                do p = 1, nDisSite
                    dPartTCritGrad = dPartTCritGrad + dCollapseMap(p,a) * dDisTCritGrad(p)
                    dPartBGrad = dPartBGrad + dCollapseMap(p,a) * dDisBGrad(p)
                end do
                do p = 1, nRandomSite
                    dPartTCritGrad = dPartTCritGrad - dRandomMap(p,a) * dRandomTCritGrad(p)
                    dPartBGrad = dPartBGrad - dRandomMap(p,a) * dRandomBGrad(p)
                end do

                dOrdGradientG(a) = dOrdGradientG(a) + &
                    dPartMagTCrit(1) * dPartTCritGrad + dPartMagB(1) * dPartBGrad - &
                    dOrdMagTCrit(1) * dOrdTCritGrad(a) - dOrdMagB(1) * dOrdBGrad(a)
                dOrdGradientH(a) = dOrdGradientH(a) + &
                    dPartMagTCrit(2) * dPartTCritGrad + dPartMagB(2) * dPartBGrad - &
                    dOrdMagTCrit(2) * dOrdTCritGrad(a) - dOrdMagB(2) * dOrdBGrad(a)
                dOrdGradientS(a) = dOrdGradientS(a) + &
                    dPartMagTCrit(3) * dPartTCritGrad + dPartMagB(3) * dPartBGrad - &
                    dOrdMagTCrit(3) * dOrdTCritGrad(a) - dOrdMagB(3) * dOrdBGrad(a)
                dOrdGradientCp(a) = dOrdGradientCp(a) + &
                    dPartMagTCrit(4) * dPartTCritGrad + dPartMagB(4) * dPartBGrad - &
                    dOrdMagTCrit(4) * dOrdTCritGrad(a) - dOrdMagB(4) * dOrdBGrad(a)

                do p = 1, nDisSite
                    dOrdGradientG(a) = dOrdGradientG(a) - dCollapseMap(p,a) * &
                        (dDisMagTCrit(1) * dDisTCritGrad(p) + dDisMagB(1) * dDisBGrad(p))
                    dOrdGradientH(a) = dOrdGradientH(a) - dCollapseMap(p,a) * &
                        (dDisMagTCrit(2) * dDisTCritGrad(p) + dDisMagB(2) * dDisBGrad(p))
                    dOrdGradientS(a) = dOrdGradientS(a) - dCollapseMap(p,a) * &
                        (dDisMagTCrit(3) * dDisTCritGrad(p) + dDisMagB(3) * dDisBGrad(p))
                    dOrdGradientCp(a) = dOrdGradientCp(a) - dCollapseMap(p,a) * &
                        (dDisMagTCrit(4) * dDisTCritGrad(p) + dDisMagB(4) * dDisBGrad(p))
                end do

                do p = 1, nRandomSite
                    dOrdGradientG(a) = dOrdGradientG(a) + dRandomMap(p,a) * &
                        (dRandomMagTCrit(1) * dRandomTCritGrad(p) + dRandomMagB(1) * dRandomBGrad(p))
                    dOrdGradientH(a) = dOrdGradientH(a) + dRandomMap(p,a) * &
                        (dRandomMagTCrit(2) * dRandomTCritGrad(p) + dRandomMagB(2) * dRandomBGrad(p))
                    dOrdGradientS(a) = dOrdGradientS(a) + dRandomMap(p,a) * &
                        (dRandomMagTCrit(3) * dRandomTCritGrad(p) + dRandomMagB(3) * dRandomBGrad(p))
                    dOrdGradientCp(a) = dOrdGradientCp(a) + dRandomMap(p,a) * &
                        (dRandomMagTCrit(4) * dRandomTCritGrad(p) + dRandomMagB(4) * dRandomBGrad(p))
                end do
            end do
        end if

        return

    end subroutine AddSUBOMFixedGradient

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

end subroutine CompGradientSUBL
