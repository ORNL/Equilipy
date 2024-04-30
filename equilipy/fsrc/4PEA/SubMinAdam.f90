!
subroutine SubMinAdam(iSolnPhaseIndex)
    !---------------------------------------------------------------------------
    !
    ! Purpose:
    ! ========
    !
    ! The purpose of this subroutine is to compute mole fraction of next iteration
    ! based on Lagrangian method. An enhanced algorithm (PEA+Adam+MDMM) is used here
    !
    ! Pertinent variables:
    ! ====================
    !
    ! INFO          An integer scalar used by LAPACK to identify a successful
    !                exit (INFO = 0) or an error (INFO /= 0).
    ! nEqn          The number of equations in the Hessian matrix.
    ! dRHS          A double real vector representing the right hand side (RHS)
    !                of the Hessian matrix and returns the direction vector.
    !
    !---------------------------------------------------------------------------
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin
!
    implicit none
!
    integer :: i, j, k, iSolnPhaseIndex, INFO, nEqn, m
    real(8),dimension(nVar):: dChemicalPotentialSubmin, dMolFractionNew, dUpdateDelta, dMolFractionLast
    real(8) :: dMaxDiff, dConstraintFactor
!
!
    ! Initialize Adam parameters
    dConstraintFactor= dMIN1(float(iterSub)/2,9D+10)
    dMolFractionNew =0D0
    dMaxDiff = 0D0
    dMolFractionLast = dMolFraction(iFirstSUB:iLastSUB)

    dChemicalPotentialSubmin = dChemicalPotential(iFirstSUB:iLastSUB)-dChemicalPotentialStar
    dSubminGibbsEst = sum(dChemicalPotentialSubmin)/nVar
    dPotentialVector = dSubminGibbsEst-dChemicalPotentialSubmin
    dMaxPotentialVector = MAXVAL(DABS(dPotentialVector))
    
    dPiDot(:nVar) = dMolFractionLast - DABS(dMolFractionLast)
    dPiDot(nVar+1) = 1-sum(dMolFractionLast)   

    ! Constraints for charge neutrality (We can implement charge here in the future)
    if (iPhaseElectronID(iSolnPhaseIndex) /= 0) then
        k    = iPhaseElectronID(iSolnPhaseIndex)
        dPiDot(nVar+2) = sum(dStoichSpecies(iFirstSUB:iLastSUB,k) * dMolFractionLast)
    end if

    ! Update lambda array
    dLambda = beta1*dLambdaLast+(1-beta1)*dPiDot

    ! Update pi array
    dPi = dPi + alpha*dLambda

    ! Update gradient of Lagrange with respect to mole fraction
    do j = 1, nVar
        ! Add the term for gradient of objective function and for inequality constraints x>0
        if(dMolFractionLast(j)<0) then
            dxDot(j) = dPotentialVector(j) - dPi(j)*2 - dConstraintFactor*dPiDot(j)
        else
            dxDot(j) = dPotentialVector(j) - dConstraintFactor*dPiDot(j)
        end if

        ! Add constraint term for sum of all mole fraction to be one
        dxDot(j) = dxDot(j)  + dPi(nVar+1) - dConstraintFactor*dPiDot(nVar+1)
       
        if (iPhaseElectronID(iSolnPhaseIndex) /= 0) then
        ! Add term for charge neutrality 
            dxDot(j) = dxDot(j)+dStoichSpecies(j,k)*(dPi(nVar+2)+dConstraintFactor*(dPiDot(nVar+2)))
        end if
    end do

    ! We don't need exponential decaying function the dxDot decays
    dvAdam = (beta1*dvAdamLast+(1-beta1)*dxDot)
    dsAdam = (beta2*dsAdamLast+(1-beta2)*dxDot**2)
!
    ! Update mole fraction
    ! Note that MDMM can results in small negative moles. In that case, refine the calculation

    do j = 1, nVar
        dUpdateDelta(j)=alpha*dvAdam(j)/(sqrt(dsAdam(j))+epsilon)
        dMolFractionNew(j) = dMolFractionLast(j) + dUpdateDelta(j)
        if(dMolFractionNew(j)<0) then
            ! Approximate mole fraction from chemical potential assuming Henrian
            ! ln(ri) = mui-ln(Xi) where ln(ri) is assumed constant
            ! dPotentialVector(j) = 0 = mu_bar-ln(ri)-ln(Xi)
            dMolFractionNew(j)= dMolFractionLast(j)*exp(dPotentialVector(j))
            iAdamNeg = iAdamNeg +1            
        end if
    end do

    ! Normalize mole fraction
    dMolFractionNew = dMolFractionNew/sum(dMolFractionNew)
    dMaxDiff = MAXVAL(DABS(dMolFractionLast-dMolFractionNew))
    dMolFraction(iFirstSUB:iLastSUB) = dMolFractionNew
    
    ! When Adam doesn't progress well move on to Lagrangian
    if((MAXVAL(ABS(dUpdateDelta))<dSubMinTolerance)&
    .or. (dMaxDiff<1D-2)&
    .or. (iAdamNeg>100)&
    .or.dMaxPotentialVector<3D0) then
        lNegativeFraction = .TRUE.
    end if


    ! Compute the chemical potentials of solution phase constituents:
    call SubMinChemicalPotential(iSolnPhaseIndex)
    !
    ! Compute the driving force of this solution phase:
    call SubMinDrivingForce
    
    ! Update 
    dvAdamLast = dvAdam
    dsAdamLast = dsAdam
    dLambdaLast= dLambda
    dDrivingForceLast = dDrivingForce
!
end subroutine SubMinAdam
!
!
