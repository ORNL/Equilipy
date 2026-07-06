!-------------------------------------------------------------------------------------------------------------
!
!> \file    SubMinLineSearch.f90
!> \brief   Backtrack the Subminimization Newton step for one solution phase.
!> \author  M.H.A. Piro
!> \date    Aug. 21, 2012
!> \sa      Subminimization.f90
!> \sa      SubMinNewton.f90
!> \sa      SubMinTraceSpeciesControl.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   08/21/2012      M.H.A. Piro         Original Subminimization line-search routine
!   06/24/2026      S.Y. Kwon           Backtrack rejected Newton trials from the original base state
!   06/24/2026      S.Y. Kwon           Documented the classic-only scope of this line search
!   06/27/2026      S.Y. Kwon           Accepted classic Newton trials with a driving-force and
!                                        mole-fraction-weighted residual merit function.
!   06/28/2026      S.Y. Kwon           Required SUBG/SUBQ line-search trials to improve strict pair
!                                        stationarity instead of only lowering driving force.
!   06/28/2026      S.Y. Kwon           Removed the driving-force guard from SUBG/SUBQ trial acceptance so
!                                        pair phases converge as stationarity root solves.
!   06/28/2026      S.Y. Kwon           Rejected non-finite line-search trials before testing pair
!                                        stationarity improvement.
!   06/28/2026      S.Y. Kwon           Removed the generic multiplicative residual correction from
!                                        SUBG/SUBQ pair phases.
!   06/29/2026      S.Y. Kwon           Switched SUBG/SUBQ trial acceptance to Thermochimica-style
!                                        driving-force decrease while retaining strict residual diagnostics.
!
!
! Purpose:
! ========
!
!> \details This routine performs a line search on the direction vector
!! computed by SubMinNewton.  Candidate mole fractions are generated from a
!! fixed base composition and a damped Newton update.  Trial acceptance uses a
!! merit function that combines the scalar phase driving force with a
!! mole-fraction-weighted active residual for ordinary non-CEF phases.  SUBG
!! and SUBQ use pair/quadruplet variables; following Thermochimica's MQM
!! subminimization, their line search accepts a damped trial when the phase
!! driving force decreases.  The strict pair stationarity residual is still
!! updated and reported, but it does not block a useful downhill driving-force
!! step during PEA phase-entry checks.  This routine is not used by the CEF
!! site-fraction path.
!
!
! Required input variables:
! =========================
!
!> \param[in] iSolnPhaseIndex Absolute solution phase index being subminimized.
!
! dMolFraction          Current active endmember fractions before the Newton trial.
! dRHS                  Newton update direction from SubMinNewton.
! dDrivingForceLast     Driving force from the previous accepted submin state.
! dMaxPotentialVector   Current active endmember chemical-potential residual.
! lSubMinTraceInactive  Phase-local trace species excluded from the active Newton solve.
!
!
! Output/updated variables:
! =========================
!
! dMolFraction          Accepted active endmember fractions after damping.
! dChemicalPotential    Updated by SubMinChemicalPotential for the accepted state.
! dDrivingForce         Updated by SubMinDrivingForce for the accepted state.
! dMaxPotentialVector   Updated by UpdateSubMinPotentialVector for the accepted state.
! lSubMinConverged      Set true when composition change and potential residual are small.
!
!
! Called subroutines/functions:
! =============================
!
! SubMinChemicalPotential      Evaluates phase chemical potentials for a trial composition.
! SubMinDrivingForce           Evaluates the scalar phase driving force.
! UpdateSubMinPotentialVector  Updates the active stationarity residual.
!
!
! Primary callers:
! ================
!
! Subminimization              Calls after SubMinNewton during PEA candidate minimization.
!
!
! Numerical assumptions:
! ======================
!
! - All backtracking trials are generated from the original base composition
!   for the current Newton step.
! - Trace-inactive endmembers remain at zero and do not participate in the
!   active residual or step-size limit.
! - dMaxPotentialVector remains the strict convergence residual.  SUBG/SUBQ
!   acceptance uses driving-force decrease, matching Thermochimica's MQM
!   subminimization behavior.  Other non-CEF phases retain the weighted-residual
!   merit so trace endmembers do not block every downhill driving-force step.
! - If the damped Newton update stalls on an ordinary non-pair, non-ionic phase
!   with mixed formula-site scaling, a bounded multiplicative residual
!   correction is tried and accepted only through the same merit function.
!
!-------------------------------------------------------------------------------------------------------------



subroutine SubMinLineSearch(iSolnPhaseIndex)
    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin
!
    implicit none
!
    integer :: i, j, k, l, iAlpha, iSolnPhaseIndex
    real(8), parameter :: dMeritResidualWeight = 0.5D0
    real(8) :: dStepLength, dTemp, dMaxChange
    real(8) :: dBaseDrivingForce, dBasePotentialNorm, dBaseWeightedNorm, dBaseMerit
    real(8) :: dBestMerit, dBestPotentialNorm, dBestDrivingForce
    real(8) :: dTrialWeightedNorm, dTrialMerit, dDrivingForceGuard
    real(8) :: dExpAlpha, dExpArg, dExpSign, dSum, dWeight
    real(8) :: dAtomFactor, dAtomFactorMin, dAtomFactorMax
    real(8), dimension(nVar) :: dMolFractionBase, dMolFractionBest, dMolFractionTrial
    logical :: lAcceptedTrial, lTrialAcceptable, lTrialFinite, lPairPhase

    dStepLength = 1D0
    dMaxChange  = 0.2D0
    dMolFractionBase = dMolFraction(iFirstSUB:iLastSUB)
    dBaseDrivingForce = dDrivingForce
    dBasePotentialNorm = dMaxPotentialVector
    dBestPotentialNorm = dBasePotentialNorm
    dBestDrivingForce = dBaseDrivingForce
    dDrivingForceGuard = 1D-4 * DMAX1(1D0, DABS(dBaseDrivingForce))
    dBaseWeightedNorm = 0D0

    do j = 1, nVar
        if (allocated(lSubMinTraceInactive)) then
            if (lSubMinTraceInactive(j)) cycle
        end if
        dWeight = DMAX1(dMolFractionBase(j), dMinMoleFraction)
        dBaseWeightedNorm = dBaseWeightedNorm + dWeight * dPotentialVector(j)**2
    end do
    dBaseWeightedNorm = DSQRT(dBaseWeightedNorm)
    dBaseMerit = dBaseDrivingForce + dMeritResidualWeight * dBaseWeightedNorm**2
    dBestMerit = dBaseMerit
    dMolFractionBest = dMolFractionBase
    lAcceptedTrial = .FALSE.
    lPairPhase = (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBG').OR.&
        (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBQ')

    ! Initialize steplength (determine a steplength that prevents negative mole fractions
    ! and constrains maximum change):
    do j = 1, nVar
!
        ! Absolute species index:
        i = iFirstSUB + j - 1 
        if (allocated(lSubMinTraceInactive)) then
            if (lSubMinTraceInactive(j)) then
                dRHS(j) = 0D0
                dMolFraction(i) = 0D0
                cycle
            end if
        end if
!
        
        ! Check if the mole fraction of this constituent is driven to be negative:
        if (dMolFraction(i) + dRHS(j) <= 0D0) then
!
            ! Determine a step length that reduces the mole fraction by a factor of 100:
            dTemp = -0.99D0 * dMolFraction(i) / dRHS(j)
!
            ! Update the step length:
            dStepLength = DMIN1(dStepLength, dTemp)
!
        end if
!
        ! Compute a step length that constrains the maximum change by dMaxChange:
        dTemp=dStepLength
        if(dRHS(j) /= 0D0)then
           dTemp       = DABS(dMaxChange / dRHS(j))
        end if
        dStepLength = DMIN1(dStepLength, dTemp)
!
    end do
!
    ! Iterate to satisfy the merit decrease conditions:
    LOOP_WOLFE: do k = 1, 8
        dSum = 0D0
        dTemp = 0D0
        do j = 1, nVar
            i = iFirstSUB + j - 1 
            if (allocated(lSubMinTraceInactive)) then
                if (lSubMinTraceInactive(j)) then
                    dMolFraction(i) = 0D0
                    cycle
                end if
            end if
            dMolFraction(i) = DMAX1(dMolFractionBase(j) + dStepLength*dRHS(j), dMinMoleFraction)
            dSum = dSum + dMolFraction(i)
        end do
        if (dSum > 0D0) then
            do j = 1, nVar
                i = iFirstSUB + j - 1
                if (allocated(lSubMinTraceInactive)) then
                    if (lSubMinTraceInactive(j)) cycle
                end if
                dMolFraction(i) = dMolFraction(i) / dSum
                dTemp = DMAX1(dTemp, DABS(dMolFraction(i) - dMolFractionBase(j)))
            end do
        end if

        ! Compute the chemical potentials of solution phase constituents:
        call SubMinChemicalPotential(iSolnPhaseIndex)
!
        ! Compute the driving force of this solution phase:
        call SubMinDrivingForce
        
        call UpdateSubMinPotentialVector
        lTrialFinite = (dDrivingForce == dDrivingForce).AND.&
            (dMaxPotentialVector == dMaxPotentialVector).AND.&
            ALL(dMolFraction(iFirstSUB:iLastSUB) == dMolFraction(iFirstSUB:iLastSUB))
!
        dTrialWeightedNorm = 0D0
        do j = 1, nVar
            i = iFirstSUB + j - 1
            if (allocated(lSubMinTraceInactive)) then
                if (lSubMinTraceInactive(j)) cycle
            end if
            dWeight = DMAX1(dMolFraction(i), dMinMoleFraction)
            dTrialWeightedNorm = dTrialWeightedNorm + dWeight * dPotentialVector(j)**2
        end do
        dTrialWeightedNorm = DSQRT(dTrialWeightedNorm)
        dTrialMerit = dDrivingForce + dMeritResidualWeight * dTrialWeightedNorm**2

        if (.NOT.lTrialFinite) then
            lTrialAcceptable = .FALSE.
        else if (lPairPhase) then
            lTrialAcceptable = (dDrivingForce < dBestDrivingForce)
        else
            lTrialAcceptable = (dDrivingForce <= dBaseDrivingForce + dDrivingForceGuard).AND.&
                (dTrialMerit < dBestMerit)
            if ((dMaxPotentialVector < dBasePotentialNorm).AND.&
                (dDrivingForce <= dBaseDrivingForce + dDrivingForceGuard)) lTrialAcceptable = .TRUE.
        end if

        if (lTrialAcceptable) then
            dBestMerit = dTrialMerit
            dBestPotentialNorm = dMaxPotentialVector
            dBestDrivingForce = dDrivingForce
            dMolFractionBest = dMolFraction(iFirstSUB:iLastSUB)
            lAcceptedTrial = .TRUE.
            exit LOOP_WOLFE
!
        else
            ! Divergence has been detected.  Damping the sub-system:
            !SYMF
            dStepLength = dStepLength * 0.5D0
            dTemp       = 0D0
!
            cycle LOOP_WOLFE
!
        end if
!
    end do LOOP_WOLFE

    if (.NOT.lAcceptedTrial) then
        dMolFraction(iFirstSUB:iLastSUB) = dMolFractionBase
        dTemp = 0D0
    else
        dMolFraction(iFirstSUB:iLastSUB) = dMolFractionBest
        dTemp = MAXVAL(DABS(dMolFractionBest - dMolFractionBase))
    end if
    call SubMinChemicalPotential(iSolnPhaseIndex)
    call SubMinDrivingForce
    call UpdateSubMinPotentialVector

    ! If an ordinary non-pair Newton step is being over-damped, try bounded
    ! multiplicative residual corrections and accept only a lower merit state.
    dAtomFactorMin = 1D100
    dAtomFactorMax = -1D100
    do j = 1, nVar
        i = iFirstSUB + j - 1
        if (allocated(lSubMinTraceInactive)) then
            if (lSubMinTraceInactive(j)) cycle
        end if
        dAtomFactor = SUM(dStoichSpecies(i,1:nElements)) / DFLOAT(iParticlesPerMole(i))
        dAtomFactorMin = DMIN1(dAtomFactorMin, dAtomFactor)
        dAtomFactorMax = DMAX1(dAtomFactorMax, dAtomFactor)
    end do
!
    if ((.NOT.lPairPhase).AND.&
        (iPhaseElectronID(iSolnPhaseIndex) == 0).AND.&
        (dMaxPotentialVector > dMaxPotentialTol).AND.&
        (dAtomFactorMin > 1D-8).AND.&
        ((dAtomFactorMax-dAtomFactorMin) > 1D-8)) then
!
        dMolFractionBase  = dMolFraction(iFirstSUB:iLastSUB)
        dMolFractionBest  = dMolFractionBase
        dBaseDrivingForce = dDrivingForce
        dBasePotentialNorm = dMaxPotentialVector
        dBestPotentialNorm = dBasePotentialNorm
        dBestDrivingForce = dBaseDrivingForce
        dBaseWeightedNorm = 0D0
        do j = 1, nVar
            i = iFirstSUB + j - 1
            if (allocated(lSubMinTraceInactive)) then
                if (lSubMinTraceInactive(j)) cycle
            end if
            dWeight = DMAX1(dMolFractionBase(j), dMinMoleFraction)
            dBaseWeightedNorm = dBaseWeightedNorm + dWeight * dPotentialVector(j)**2
        end do
        dBaseWeightedNorm = DSQRT(dBaseWeightedNorm)
        dBaseMerit = dBaseDrivingForce + dMeritResidualWeight * dBaseWeightedNorm**2
        dDrivingForceGuard = 1D-4 * DMAX1(1D0, DABS(dBaseDrivingForce))
        dBestMerit = dBaseMerit
!
        do l = 1, 10
!
            iAlpha = (l + 1) / 2
            if (iAlpha == 1) then
                dExpAlpha = 0.25D0
            else if (iAlpha == 2) then
                dExpAlpha = 0.5D0
            else if (iAlpha == 3) then
                dExpAlpha = 1D0
            else if (iAlpha == 4) then
                dExpAlpha = 2D0
            else
                dExpAlpha = 5D0
            end if
            if (MOD(l,2) == 1) then
                dExpSign = 1D0
            else
                dExpSign = -1D0
            end if
!
            dSum = 0D0
            do j = 1, nVar
                if (allocated(lSubMinTraceInactive)) then
                    if (lSubMinTraceInactive(j)) then
                        dMolFractionTrial(j) = 0D0
                        cycle
                    end if
                end if
                dExpArg = DMAX1(-60D0, DMIN1(60D0, dExpSign*dExpAlpha*dPotentialVector(j)))
                dMolFractionTrial(j) = DMAX1(dMolFractionBase(j)*EXP(dExpArg), dMinMoleFraction)
                dSum = dSum + dMolFractionTrial(j)
            end do
            if (dSum > 0D0) dMolFractionTrial = dMolFractionTrial / dSum
            dMolFraction(iFirstSUB:iLastSUB) = dMolFractionTrial
!
            call SubMinChemicalPotential(iSolnPhaseIndex)
            call SubMinDrivingForce
!
            call UpdateSubMinPotentialVector
            lTrialFinite = (dDrivingForce == dDrivingForce).AND.&
                (dMaxPotentialVector == dMaxPotentialVector).AND.&
                ALL(dMolFraction(iFirstSUB:iLastSUB) == dMolFraction(iFirstSUB:iLastSUB))

            dTrialWeightedNorm = 0D0
            do j = 1, nVar
                i = iFirstSUB + j - 1
                if (allocated(lSubMinTraceInactive)) then
                    if (lSubMinTraceInactive(j)) cycle
                end if
                dWeight = DMAX1(dMolFraction(i), dMinMoleFraction)
                dTrialWeightedNorm = dTrialWeightedNorm + dWeight * dPotentialVector(j)**2
            end do
            dTrialWeightedNorm = DSQRT(dTrialWeightedNorm)
            dTrialMerit = dDrivingForce + dMeritResidualWeight * dTrialWeightedNorm**2
!
            if (.NOT.lTrialFinite) then
                lTrialAcceptable = .FALSE.
            else if (lPairPhase) then
                lTrialAcceptable = (dDrivingForce < dBestDrivingForce)
            else
                lTrialAcceptable = (dDrivingForce <= dBaseDrivingForce + dDrivingForceGuard).AND.&
                    (dTrialMerit < dBestMerit)
            end if
!
            if (lTrialAcceptable) then
                dBestMerit = dTrialMerit
                dBestPotentialNorm = dMaxPotentialVector
                dBestDrivingForce  = dDrivingForce
                dMolFractionBest   = dMolFractionTrial
            end if
!
        end do
!
        dMolFraction(iFirstSUB:iLastSUB) = dMolFractionBest
        call SubMinChemicalPotential(iSolnPhaseIndex)
        call SubMinDrivingForce
!
        call UpdateSubMinPotentialVector
        dTemp = DMAX1(dTemp, MAXVAL(DABS(dMolFractionBest-dMolFractionBase)))
!
    end if
!
    ! Check convergence (maximum change to any mole fraction):
    if ((dTemp <= dSubMinTolerance).AND.&
        (dMaxPotentialVector < dMaxPotentialTol)) lSubMinConverged = .TRUE.
!
    ! Store the driving force from the last iteration:
    dDrivingForceLast = dDrivingForce
!
end subroutine SubMinLineSearch
!
!
