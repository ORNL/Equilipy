!> \brief Evaluate one static-grid phase block with single-core G-only kernels.
!!
!! \details Loads phase topology once, then evaluates every constitution in
!! canonical row order without changing shared thermodynamic work arrays.
module GridBatchKernels
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GridBatchKernels.f90
    !> \brief   Single-core phase-block Gibbs kernels for static-grid Leveling rows.
    !> \author  S.Y. Kwon
    !> \date    Jul. 19, 2026
    !> \sa      GridDiscovery.f90
    !> \sa      CompGradientSUBL.f90
    !> \sa      CompExcessGibbsEnergyRKMP.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Added batched Gibbs-energy evaluation for CEF, SUBOM, and RKMP phase blocks.
    !
    ! Purpose:
    ! ========
    !
    !> \details Static-grid rows require only the scalar Gibbs energy at fixed
    !! constitution. These kernels use direct scalar-G arithmetic while
    !! omitting partial properties, gradients, Hessians, and per-row
    !! module-state setup.
    !
    ! Required input variables:
    ! =========================
    !
    !> \param[in] iPhaseIndex   Active-system solution-phase index.
    !> \param[in] dFractionBlock Endmember constitutions in immutable row order.
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] dGibbsBlock  Scalar Gibbs energy in internal G/(RT) units.
    !> \param[out] iMode        Typed kernel-family result.
    !> \param[out] iInfo        Zero on success; nonzero requests phase-level scalar fallback.
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - The caller groups rows by phase and preserves their original order.
    ! - No thread or parallel runtime participates; evaluation is single-core.
    ! - CEF parameter topology is accepted only when the analytical production
    !   CEF path already supports it.  Other topology falls back as one phase.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: dTemperature
    USE ModuleGEMSolver, ONLY: GRID_EVALUATION_BATCH_CEF, GRID_EVALUATION_BATCH_RKMP

    implicit none

    private
    public :: EvaluateGridPhaseGibbsBatch
    public :: EvaluateGridPhaseGibbsScalar

contains

subroutine EvaluateGridPhaseGibbsBatch(iPhaseIndex, dFractionBlock, dGibbsBlock, iMode, iInfo)
    implicit none

    integer, intent(in) :: iPhaseIndex
    real(8), dimension(:,:), intent(in) :: dFractionBlock
    real(8), dimension(:), intent(out) :: dGibbsBlock
    integer, intent(out) :: iMode, iInfo

    iMode = 0
    iInfo = 0
    dGibbsBlock = 0D0

    select case (TRIM(cSolnPhaseType(iPhaseIndex)))
        case ('SUBL', 'SUBLM', 'SUBOM')
            call EvaluateCEFGibbsBlock(iPhaseIndex, dFractionBlock, dGibbsBlock, iInfo)
            if (iInfo == 0) iMode = GRID_EVALUATION_BATCH_CEF
        case ('RKMP')
            call EvaluateRKMPGibbsBlock(iPhaseIndex, dFractionBlock, dGibbsBlock, iInfo)
            if (iInfo == 0) iMode = GRID_EVALUATION_BATCH_RKMP
        case default
            iInfo = 1
    end select

    return
end subroutine EvaluateGridPhaseGibbsBatch


!> \brief Evaluate one phase constitution with the direct scalar-G reference path.
subroutine EvaluateGridPhaseGibbsScalar(iPhaseIndex, dFraction, dGibbs, iInfo)
    implicit none

    integer, intent(in) :: iPhaseIndex
    real(8), dimension(:), intent(in) :: dFraction
    real(8), intent(out) :: dGibbs
    integer, intent(out) :: iInfo
    integer :: iMode
    real(8) :: dFractionBlock(1,SIZE(dFraction)), dGibbsBlock(1)

    dFractionBlock(1,:) = dFraction
    call EvaluateGridPhaseGibbsBatch(iPhaseIndex, dFractionBlock, dGibbsBlock, iMode, iInfo)
    dGibbs = dGibbsBlock(1)

    return
end subroutine EvaluateGridPhaseGibbsScalar


subroutine EvaluateCEFGibbsBlock(iPhaseIndex, dFractionBlock, dGibbsBlock, iInfo)
    implicit none

    integer, intent(in) :: iPhaseIndex
    real(8), dimension(:,:), intent(in) :: dFractionBlock
    real(8), dimension(:), intent(out) :: dGibbsBlock
    integer, intent(out) :: iInfo

    integer :: iPhaseID, iDisPhase, iDisPhaseID, iRow
    integer :: iOrdToDis(nMaxSublatticeSys)
    integer :: iOrdConToDis(nMaxSublatticeSys,nMaxConstituentSys)
    real(8) :: dGroupStoich(nMaxSublatticeSys)
    real(8) :: dOrdSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8) :: dDisSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8) :: dRandomSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8) :: dOrdG, dDisG, dRandomG
    real(8) :: dOrdTCrit, dOrdB, dDisTCrit, dDisB, dRandomTCrit, dRandomB
    real(8) :: dOrdMag, dDisMag, dRandomMag, dPartMag

    iInfo = 0
    iPhaseID = iPhaseSublattice(iPhaseIndex)
    if (iPhaseID <= 0) then
        iInfo = 2
        return
    end if
    call ValidateCEFParameterPlan(iPhaseIndex, iInfo)
    if (iInfo /= 0) return

    iDisPhase = 0
    iDisPhaseID = 0
    if (TRIM(cSolnPhaseType(iPhaseIndex)) == 'SUBOM') then
        if (allocated(iDisorderedPhase)) iDisPhase = iDisorderedPhase(iPhaseIndex)
        if ((iDisPhase <= 0).OR.(iDisPhase > nSolnPhasesSys)) then
            iInfo = 3
            return
        end if
        iDisPhaseID = iPhaseSublattice(iDisPhase)
        call ValidateCEFParameterPlan(iDisPhase, iInfo)
        if (iInfo /= 0) return
        call BuildPartitionPlan(iPhaseID, iDisPhaseID, iOrdToDis, iOrdConToDis, dGroupStoich, iInfo)
        if (iInfo /= 0) return
    end if

    do iRow = 1, SIZE(dFractionBlock,1)
        call BuildCEFSiteFractions(iPhaseIndex, iPhaseID, dFractionBlock(iRow,:), dOrdSite)
        call EvaluateRawCEFGibbs(iPhaseIndex, iPhaseID, dOrdSite, dOrdG, dOrdTCrit, dOrdB, dOrdMag)
        if (TRIM(cSolnPhaseType(iPhaseIndex)) /= 'SUBOM') then
            dGibbsBlock(iRow) = dOrdG
            cycle
        end if

        call MapPartitionSites(iPhaseID, iDisPhaseID, iOrdToDis, iOrdConToDis, dGroupStoich, &
            dOrdSite, dDisSite, dRandomSite)
        call EvaluateRawCEFGibbs(iDisPhase, iDisPhaseID, dDisSite, dDisG, dDisTCrit, dDisB, dDisMag)
        call EvaluateRawCEFGibbs(iPhaseIndex, iPhaseID, dRandomSite, dRandomG, &
            dRandomTCrit, dRandomB, dRandomMag)
        call EvaluateMagneticGibbs(iPhaseIndex, dOrdTCrit+dDisTCrit-dRandomTCrit, &
            dOrdB+dDisB-dRandomB, dPartMag)

        dGibbsBlock(iRow) = dOrdG + dDisG - dRandomG
        dGibbsBlock(iRow) = dGibbsBlock(iRow) + dPartMag - dOrdMag - dDisMag + dRandomMag
    end do

    return
end subroutine EvaluateCEFGibbsBlock


subroutine ValidateCEFParameterPlan(iPhaseIndex, iInfo)
    implicit none

    integer, intent(in) :: iPhaseIndex
    integer, intent(out) :: iInfo
    integer :: iParam

    iInfo = 0
    do iParam = nParamPhase(iPhaseIndex-1)+1, nParamPhase(iPhaseIndex)
        if (iRegularParam(iParam,1) <= 0) cycle
        if ((iSUBLParamData(iParam,1) == 1).AND.(iSUBLParamData(iParam,3) == 2)) cycle
        if ((iSUBLParamData(iParam,1) == 1).AND.(iSUBLParamData(iParam,3) == 3)) cycle
        if ((iSUBLParamData(iParam,1) == 2).AND.(iSUBLParamData(iParam,3) == 2).AND. &
            (iSUBLParamData(iParam,5) == 2)) cycle
        iInfo = 4
        return
    end do

    return
end subroutine ValidateCEFParameterPlan


subroutine BuildCEFSiteFractions(iPhaseIndex, iPhaseID, dFraction, dSite)
    implicit none

    integer, intent(in) :: iPhaseIndex, iPhaseID
    real(8), dimension(:), intent(in) :: dFraction
    real(8), intent(out) :: dSite(nMaxSublatticeSys,nMaxConstituentSys)
    integer :: i, m, s, c
    real(8) :: dNorm

    dSite = 0D0
    do i = 1, SIZE(dFraction)
        m = i
        do s = 1, nSublatticePhase(iPhaseID)
            c = iConstituentSublattice(iPhaseID,s,m)
            dSite(s,c) = dSite(s,c) + DMAX1(dFraction(i),0D0)
        end do
    end do
    do s = 1, nSublatticePhase(iPhaseID)
        dNorm = 0D0
        do c = 1, nConstituentSublattice(iPhaseID,s)
            dNorm = dNorm + dSite(s,c)
        end do
        if (dNorm <= 0D0) cycle
        do c = 1, nConstituentSublattice(iPhaseID,s)
            dSite(s,c) = dSite(s,c) / dNorm
        end do
    end do

    return
end subroutine BuildCEFSiteFractions


subroutine EvaluateRawCEFGibbs(iPhaseIndex, iPhaseID, dSite, dGibbs, dTCrit, dB, dMagnetic)
    implicit none

    integer, intent(in) :: iPhaseIndex, iPhaseID
    real(8), intent(in) :: dSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8), intent(out) :: dGibbs, dTCrit, dB, dMagnetic
    integer :: iFirst, iLast, i, m, s, c, iParam
    real(8) :: dProduct, dSiteValue, dFactor

    dGibbs = 0D0
    iFirst = nSpeciesPhase(iPhaseIndex-1) + 1
    iLast = nSpeciesPhase(iPhaseIndex)
    do i = iFirst, iLast
        m = i - iFirst + 1
        dProduct = 1D0
        do s = 1, nSublatticePhase(iPhaseID)
            c = iConstituentSublattice(iPhaseID,s,m)
            dProduct = dProduct * DMAX1(dSite(s,c),1D-75)
        end do
        dGibbs = dGibbs + dStdGibbsEnergy(i)*dProduct
    end do
    do s = 1, nSublatticePhase(iPhaseID)
        do c = 1, nConstituentSublattice(iPhaseID,s)
            dSiteValue = DMAX1(dSite(s,c),1D-75)
            dGibbs = dGibbs + dStoichSublattice(iPhaseID,s)*dSiteValue*DLOG(dSiteValue)
        end do
    end do
    do iParam = nParamPhase(iPhaseIndex-1)+1, nParamPhase(iPhaseIndex)
        call EvaluateCEFParameterFactor(iParam, iPhaseID, dSite, dFactor)
        dGibbs = dGibbs + dExcessGibbsParam(iParam)*dFactor
    end do

    call EvaluateMagneticVariables(iPhaseIndex, iPhaseID, dSite, dTCrit, dB)
    call EvaluateMagneticGibbs(iPhaseIndex, dTCrit, dB, dMagnetic)
    dGibbs = dGibbs + dMagnetic

    return
end subroutine EvaluateRawCEFGibbs


subroutine EvaluateCEFParameterFactor(iParam, iPhaseID, dSite, dFactor)
    implicit none

    integer, intent(in) :: iParam, iPhaseID
    real(8), intent(in) :: dSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8), intent(out) :: dFactor
    integer :: nParamCon, k, s, c, q, iExponent, iTernaryCon
    integer :: iFirst, iSecond, iThird, iSub, iFirst2, iSecond2, iSub2, iTemp
    real(8) :: dLinearValue(2), dLinearPower(2), dTempValue
    logical :: lWeighted

    dLinearValue = 0D0
    dLinearPower = 0D0
    dFactor = 0D0
    nParamCon = iRegularParam(iParam,1)
    if (nParamCon <= 0) return
    iFirst = 0; iSecond = 0; iThird = 0; iSub = 0
    iFirst2 = 0; iSecond2 = 0; iSub2 = 0
    if ((iSUBLParamData(iParam,1) == 1).AND.(iSUBLParamData(iParam,3) == 2)) then
        do k = 2, nParamCon+1
            c = MOD(iRegularParam(iParam,k),10000)
            s = (iRegularParam(iParam,k)-c)/10000
            if (k == iSUBLParamData(iParam,2)) then
                iFirst = c; iSub = s
            else if (k == iSUBLParamData(iParam,2)+1) then
                iSecond = c
            end if
        end do
        iExponent = iRegularParam(iParam,nParamCon+2)
        if (iExponent > 0) then
            dLinearPower(1) = DBLE(iExponent)
            dLinearValue(1) = dSite(iSub,iFirst)-dSite(iSub,iSecond)
        end if
    else if ((iSUBLParamData(iParam,1) == 1).AND.(iSUBLParamData(iParam,3) == 3)) then
        iTernaryCon = iRegularParam(iParam,nParamCon+2)
        lWeighted = iRegularParam(iParam,nParamCon+3) == 1
        do k = 2, nParamCon+1
            c = MOD(iRegularParam(iParam,k),10000)
            s = (iRegularParam(iParam,k)-c)/10000
            if (k == iSUBLParamData(iParam,2)+iTernaryCon) then
                iFirst = c; iSub = s
            else if (k == iSUBLParamData(iParam,2)+MOD(iTernaryCon+1,3)) then
                iSecond = c
            else if (k == iSUBLParamData(iParam,2)+MOD(iTernaryCon+2,3)) then
                iThird = c
            end if
        end do
        if (lWeighted) then
            dTempValue = 0D0
            do c = 1, nConstituentSublattice(iPhaseID,iSub)
                dTempValue = dTempValue + dSite(iSub,c)
            end do
            dLinearPower(1) = 1D0
            dLinearValue(1) = dSite(iSub,iFirst) + &
                (dTempValue-dSite(iSub,iFirst)-dSite(iSub,iSecond)-dSite(iSub,iThird))/3D0
        end if
    else
        do k = 2, nParamCon+1
            c = MOD(iRegularParam(iParam,k),10000)
            s = (iRegularParam(iParam,k)-c)/10000
            if (k == iSUBLParamData(iParam,2)) then
                iFirst = c; iSub = s
            else if (k == iSUBLParamData(iParam,2)+1) then
                iSecond = c
            else if (k == iSUBLParamData(iParam,4)) then
                iFirst2 = c; iSub2 = s
            else if (k == iSUBLParamData(iParam,4)+1) then
                iSecond2 = c
            end if
        end do
        if (MOD(iRegularParam(iParam,nParamCon+2),2) == 0) then
            iExponent = iRegularParam(iParam,nParamCon+2)/2
        else
            iExponent = (iRegularParam(iParam,nParamCon+2)-1)/2
            iTemp=iFirst; iFirst=iFirst2; iFirst2=iTemp
            iTemp=iSecond; iSecond=iSecond2; iSecond2=iTemp
            iTemp=iSub; iSub=iSub2; iSub2=iTemp
        end if
        if (iExponent > 0) then
            dLinearPower(1) = DBLE(iExponent)
            dLinearValue(1) = dSite(iSub,iFirst)-dSite(iSub,iSecond)
        end if
    end if

    dFactor = 1D0
    do k = 2, nParamCon+1
        c = MOD(iRegularParam(iParam,k),10000)
        s = (iRegularParam(iParam,k)-c)/10000
        dFactor = dFactor*DMAX1(dSite(s,c),1D-75)
    end do
    do q = 1, 2
        iExponent = NINT(dLinearPower(q))
        if (iExponent > 0) dFactor = dFactor*dLinearValue(q)**iExponent
    end do

    return
end subroutine EvaluateCEFParameterFactor


subroutine EvaluateMagneticVariables(iPhaseIndex, iPhaseID, dSite, dTCrit, dB)
    implicit none

    integer, intent(in) :: iPhaseIndex, iPhaseID
    real(8), intent(in) :: dSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8), intent(out) :: dTCrit, dB
    integer :: iFirst, iLast, i, m, s, c, iParam, k, nParamCon, iExponent
    real(8) :: dProduct, dLinear

    dTCrit = 0D0
    dB = 0D0
    iFirst = nSpeciesPhase(iPhaseIndex-1)+1
    iLast = nSpeciesPhase(iPhaseIndex)
    do i = iFirst, iLast
        m = i-iFirst+1
        dProduct = 1D0
        do s = 1, nSublatticePhase(iPhaseID)
            c = iConstituentSublattice(iPhaseID,s,m)
            dProduct = dProduct*DMAX1(dSite(s,c),1D-75)
        end do
        dTCrit = dTCrit + dCoeffGibbsMagnetic(i,1)*dProduct
        dB = dB + dCoeffGibbsMagnetic(i,2)*dProduct
    end do
    do iParam = nMagParamPhase(iPhaseIndex-1)+1, nMagParamPhase(iPhaseIndex)
        nParamCon = iMagneticParam(iParam,1)
        if (nParamCon <= 0) cycle
        dProduct = 1D0
        dLinear = 0D0
        do k = 2, nParamCon+1
            c = MOD(iMagneticParam(iParam,k),10000)
            s = (iMagneticParam(iParam,k)-c)/10000
            dProduct = dProduct*DMAX1(dSite(s,c),1D-75)
            if (k == 2) dLinear = dLinear + dSite(s,c)
            if (k == 3) dLinear = dLinear - dSite(s,c)
        end do
        iExponent = iMagneticParam(iParam,nParamCon+2)
        if (iExponent > 0) dProduct = dProduct*dLinear**iExponent
        dTCrit = dTCrit + dMagneticParam(iParam,1)*dProduct
        dB = dB + dMagneticParam(iParam,2)*dProduct
    end do

    return
end subroutine EvaluateMagneticVariables


subroutine EvaluateMagneticGibbs(iPhaseIndex, dTCritInput, dBInput, dMagnetic)
    implicit none

    integer, intent(in) :: iPhaseIndex
    real(8), intent(in) :: dTCritInput, dBInput
    real(8), intent(out) :: dMagnetic
    integer :: iFirst
    real(8) :: dTCrit, dB, dStructure, dP, dTau, dD, dA, dTemp, dC, dF

    dMagnetic = 0D0
    iFirst = nSpeciesPhase(iPhaseIndex-1)+1
    dStructure = dCoeffGibbsMagnetic(iFirst,3)
    dP = dCoeffGibbsMagnetic(iFirst,4)
    if (dP == 0D0) return
    dTCrit = dTCritInput
    dB = dBInput
    if (dB < 0D0) dB = -dB*dStructure
    if (dTCrit < 0D0) dTCrit = -dTCrit*dStructure
    if ((dTCrit == 0D0).OR.(dB <= -1D0)) return

    dTau = dTemperature/dTCrit
    dD = 518D0/1125D0 + 11692D0/15975D0*(1D0/dP-1D0)
    if (dTau > 1D0) then
        dA = dTau**(-5)
        dTemp = dA**3
        dC = dA*dA*dTemp
        dF = -(dA/10D0+dTemp/315D0+dC/1500D0)/dD
    else
        dA = dTau**3
        dTemp = dA**3
        dC = dA*dA*dTemp
        dF = 1D0-(79D0/(140D0*dP*dTau)+ &
            474D0/497D0*(1D0/dP-1D0)*(dA/6D0+dTemp/135D0+dC/600D0))/dD
    end if
    dMagnetic = DLOG(1D0+dB)*dF

    return
end subroutine EvaluateMagneticGibbs


subroutine BuildPartitionPlan(iOrdPhaseID, iDisPhaseID, iOrdToDis, iOrdConToDis, dGroupStoich, iInfo)
    implicit none

    integer, intent(in) :: iOrdPhaseID, iDisPhaseID
    integer, intent(out) :: iOrdToDis(nMaxSublatticeSys)
    integer, intent(out) :: iOrdConToDis(nMaxSublatticeSys,nMaxConstituentSys)
    real(8), intent(out) :: dGroupStoich(nMaxSublatticeSys)
    integer, intent(out) :: iInfo
    integer :: so, sd, c, d, iBest, iBestCount
    logical :: lSubset, lFound

    iOrdToDis = 0
    iOrdConToDis = 0
    dGroupStoich = 0D0
    iInfo = 0
    do so = 1, nSublatticePhase(iOrdPhaseID)
        iBest = 0
        iBestCount = nMaxConstituentSys+1
        do sd = 1, nSublatticePhase(iDisPhaseID)
            lSubset = .TRUE.
            do c = 1, nConstituentSublattice(iOrdPhaseID,so)
                lFound = .FALSE.
                do d = 1, nConstituentSublattice(iDisPhaseID,sd)
                    if (TRIM(cConstituentNameSUB(iOrdPhaseID,so,c)) == &
                        TRIM(cConstituentNameSUB(iDisPhaseID,sd,d))) lFound = .TRUE.
                end do
                if (.NOT.lFound) lSubset = .FALSE.
            end do
            if (lSubset.AND.(nConstituentSublattice(iDisPhaseID,sd) < iBestCount)) then
                iBest = sd
                iBestCount = nConstituentSublattice(iDisPhaseID,sd)
            end if
        end do
        if (iBest == 0) then
            iInfo = 5
            return
        end if
        iOrdToDis(so) = iBest
        dGroupStoich(iBest) = dGroupStoich(iBest)+dStoichSublattice(iOrdPhaseID,so)
        do c = 1, nConstituentSublattice(iOrdPhaseID,so)
            do d = 1, nConstituentSublattice(iDisPhaseID,iBest)
                if (TRIM(cConstituentNameSUB(iOrdPhaseID,so,c)) == &
                    TRIM(cConstituentNameSUB(iDisPhaseID,iBest,d))) iOrdConToDis(so,c) = d
            end do
            if (iOrdConToDis(so,c) == 0) then
                iInfo = 5
                return
            end if
        end do
    end do

    return
end subroutine BuildPartitionPlan


subroutine MapPartitionSites(iOrdPhaseID, iDisPhaseID, iOrdToDis, iOrdConToDis, dGroupStoich, &
    dOrdSite, dDisSite, dRandomSite)
    implicit none

    integer, intent(in) :: iOrdPhaseID, iDisPhaseID
    integer, intent(in) :: iOrdToDis(nMaxSublatticeSys)
    integer, intent(in) :: iOrdConToDis(nMaxSublatticeSys,nMaxConstituentSys)
    real(8), intent(in) :: dGroupStoich(nMaxSublatticeSys)
    real(8), intent(in) :: dOrdSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8), intent(out) :: dDisSite(nMaxSublatticeSys,nMaxConstituentSys)
    real(8), intent(out) :: dRandomSite(nMaxSublatticeSys,nMaxConstituentSys)
    integer :: so, sd, c, d
    real(8) :: dWeight

    dDisSite = 0D0
    do so = 1, nSublatticePhase(iOrdPhaseID)
        sd = iOrdToDis(so)
        dWeight = dStoichSublattice(iOrdPhaseID,so)/dGroupStoich(sd)
        do c = 1, nConstituentSublattice(iOrdPhaseID,so)
            d = iOrdConToDis(so,c)
            dDisSite(sd,d) = dDisSite(sd,d)+dWeight*dOrdSite(so,c)
        end do
    end do
    dRandomSite = 0D0
    do so = 1, nSublatticePhase(iOrdPhaseID)
        sd = iOrdToDis(so)
        do c = 1, nConstituentSublattice(iOrdPhaseID,so)
            d = iOrdConToDis(so,c)
            dRandomSite(so,c) = dDisSite(sd,d)
        end do
    end do

    return
end subroutine MapPartitionSites


subroutine EvaluateRKMPGibbsBlock(iPhaseIndex, dFractionBlock, dGibbsBlock, iInfo)
    implicit none

    integer, intent(in) :: iPhaseIndex
    real(8), dimension(:,:), intent(in) :: dFractionBlock
    real(8), dimension(:), intent(out) :: dGibbsBlock
    integer, intent(out) :: iInfo
    integer :: iRow, iParam, i, j, nLocal, iExponent
    integer :: i1, i2, i3, i4, iSelector
    real(8) :: x1, x2, x3, x4, dx, xprod, xj, dIdeal, dExcess

    iInfo = 0
    nLocal = SIZE(dFractionBlock,2)
    do iParam = nParamPhase(iPhaseIndex-1)+1, nParamPhase(iPhaseIndex)
        if ((iRegularParam(iParam,1) < 2).OR.(iRegularParam(iParam,1) > 4)) then
            iInfo = 6
            return
        end if
    end do

    do iRow = 1, SIZE(dFractionBlock,1)
        dIdeal = 0D0
        do i = 1, nLocal
            j = nSpeciesPhase(iPhaseIndex-1)+i
            dIdeal = dIdeal+dFractionBlock(iRow,i)*dStdGibbsEnergy(j)
            if ((dFractionBlock(iRow,i) >= 1D-75).AND. &
                (dFractionBlock(iRow,i) < 1D0-1D-75)) then
                dIdeal = dIdeal+dFractionBlock(iRow,i)*DLOG(dFractionBlock(iRow,i))
            end if
        end do

        dExcess = 0D0
        do iParam = nParamPhase(iPhaseIndex-1)+1, nParamPhase(iPhaseIndex)
            i1 = iRegularParam(iParam,2)
            i2 = iRegularParam(iParam,3)
            x1 = dFractionBlock(iRow,i1)
            x2 = dFractionBlock(iRow,i2)
            xprod = x1*x2
            dx = x1-x2
            if (iRegularParam(iParam,1) == 2) then
                iExponent = iRegularParam(iParam,4)
                dExcess = dExcess+dExcessGibbsParam(iParam)*xprod*dx**iExponent
            else if (iRegularParam(iParam,1) == 3) then
                i3 = iRegularParam(iParam,4)
                iSelector = iRegularParam(iParam,5)
                x3 = dFractionBlock(iRow,i3)
                xj = dFractionBlock(iRow,iSelector)
                xprod = xprod*x3
                dExcess = dExcess+dExcessGibbsParam(iParam)*xprod* &
                    ((1D0-x1-x2-x3)/3D0+xj)
            else
                i3 = iRegularParam(iParam,4)
                i4 = iRegularParam(iParam,5)
                x3 = dFractionBlock(iRow,i3)
                x4 = dFractionBlock(iRow,i4)
                xprod = xprod*x3*x4
                dExcess = dExcess+dExcessGibbsParam(iParam)*xprod
            end if
        end do
        dGibbsBlock(iRow) = dIdeal+dExcess
    end do

    return
end subroutine EvaluateRKMPGibbsBlock

end module GridBatchKernels
