!> \brief   Renormalize active endmembers after phase-local trace removal.
!!
!! \details Trace-inactive endmembers are held at zero and the remaining
!! active endmember fractions are normalized to unity for the current
!! Subminimization phase only.



subroutine NormalizeSubMinTraceActiveSpecies

    USE ModuleThermo
    USE ModuleSubMin

    implicit none

    integer :: i, j
    real(8) :: dActiveSum

    if (.NOT.allocated(lSubMinTraceInactive)) return

    dActiveSum = 0D0
    do j = 1, nVar
        i = iFirstSUB + j - 1
        if (lSubMinTraceInactive(j)) then
            dMolFraction(i) = 0D0
        else
            dActiveSum = dActiveSum + dMolFraction(i)
        end if
    end do

    if (dActiveSum <= 0D0) return

    do j = 1, nVar
        i = iFirstSUB + j - 1
        if (.NOT.lSubMinTraceInactive(j)) dMolFraction(i) = dMolFraction(i) / dActiveSum
    end do

    return

end subroutine NormalizeSubMinTraceActiveSpecies


!> \brief   Compute the active Subminimization chemical-potential residual.
!!
!! \details Builds the phase-local residual vector used by classic
!! SubMinNewton/SubMinLineSearch.  Trace-inactive endmembers are skipped so
!! dMaxPotentialVector measures only the reduced active set.



subroutine UpdateSubMinPotentialVector

    USE ModuleThermo
    USE ModuleSubMin

    implicit none

    integer :: i, j
    real(8) :: dResidual

    dSubminGibbsEst = 0D0
    dPotentialVector = 0D0
    dMaxPotentialVector = 0D0

    do j = 1, nVar
        if (allocated(lSubMinTraceInactive)) then
            if (lSubMinTraceInactive(j)) cycle
        end if
        i = iFirstSUB + j - 1
        dSubminGibbsEst = dSubminGibbsEst + &
            (dChemicalPotential(i) - dChemicalPotentialStar(j)) * dMolFraction(i)
    end do

    do j = 1, nVar
        if (allocated(lSubMinTraceInactive)) then
            if (lSubMinTraceInactive(j)) cycle
        end if
        i = iFirstSUB + j - 1
        dResidual = dSubminGibbsEst - (dChemicalPotential(i) - dChemicalPotentialStar(j))
        dPotentialVector(j) = dResidual
        dMaxPotentialVector = DMAX1(dMaxPotentialVector, DABS(dResidual))
    end do

    return

end subroutine UpdateSubMinPotentialVector


!> \brief   Approximate active trace endmembers with a Henrian-region update.
!!
!! \details For an active trace endmember, use the previous evaluated
!! residual level as the local equilibrium constant and hold the activity
!! coefficient fixed:
!!
!!    ln(x_i) = mu_e - (g_i^0 - mu_i*) - ln(gamma_i^H)
!!
!! where ``mu_e`` is represented by ``dSubminGibbsEst`` in the classic submin
!! residual equations, and ``ln(gamma_i^H)`` is estimated from the current
!! chemical potential.  The estimate is accepted only when the old and
!! predicted mole fractions are both in the trace region.  CEF phases are
!! skipped because their reference term is a site-fraction/order-disorder
!! construction rather than a simple endmember ``g_i^0``.  SUBG/SUBQ pair
!! phases are skipped because their variables are MQM pair/quadruplet
!! populations, not ordinary endmember mole fractions; their trace relation
!! must include coordination-number scaling.  Ionic phases are skipped until
!! the same update is derived with the charge-neutrality row.
!!
!! \param[in]  iSolnPhaseIndex  Absolute solution phase index being subminimized.
!! \param[out] lChanged         True when one or more trace estimates were applied.



subroutine ApplySubMinHenrianTraceEstimate(iSolnPhaseIndex, lChanged)

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer, intent(in) :: iSolnPhaseIndex
    logical, intent(out) :: lChanged
    integer :: i, j, nActiveSpecies, nEstimatedSpecies
    real(8) :: dOldFraction, dLogGamma, dLogEstimate, dEstimate
    real(8) :: dTraceSum, dNonTraceSum, dScale, dMaxChange
    real(8) :: dTraceLogLimit, dLogFloor
    real(8), dimension(nVar) :: dTrialFraction
    logical, dimension(nVar) :: lEstimated

    lChanged = .FALSE.

    if (.NOT.lTraceSpeciesControlEnabled) return
    if (lPostProcess) return
    if (.NOT.allocated(lSubMinTraceInactive)) return
    if ((TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBL').OR.&
        (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBLM').OR.&
        (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBOM').OR.&
        (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBG').OR.&
        (TRIM(cSolnPhaseType(iSolnPhaseIndex)) == 'SUBQ')) return
    if (iPhaseElectronID(iSolnPhaseIndex) /= 0) return
    if (nVar <= 1) return
    if (dSubminGibbsEst /= dSubminGibbsEst) return

    dTrialFraction = dMolFraction(iFirstSUB:iLastSUB)
    lEstimated = .FALSE.
    dTraceSum = 0D0
    dNonTraceSum = 0D0
    nActiveSpecies = 0
    nEstimatedSpecies = 0
    dTraceLogLimit = DLOG(dSubMinHenrianTraceThreshold)
    dLogFloor = DLOG(dMinMoleFraction)

    do j = 1, nVar
        i = iFirstSUB + j - 1
        if (lSubMinTraceInactive(j)) then
            dTrialFraction(j) = 0D0
            cycle
        end if

        nActiveSpecies = nActiveSpecies + 1
        dOldFraction = dMolFraction(i)

        if ((dOldFraction > 0D0).AND.(dOldFraction <= dSubMinHenrianTraceThreshold)) then
            dLogGamma = dChemicalPotential(i) - dStdGibbsEnergy(i) - &
                DLOG(DMAX1(dOldFraction, dMinMoleFraction))
            if (dLogGamma /= dLogGamma) then
                dNonTraceSum = dNonTraceSum + dOldFraction
                cycle
            end if

            dLogEstimate = dSubminGibbsEst - &
                (dStdGibbsEnergy(i) - dChemicalPotentialStar(j)) - dLogGamma

            if (dLogEstimate <= dTraceLogLimit) then
                dEstimate = DEXP(DMAX1(dLogFloor, dLogEstimate))
                dEstimate = DMAX1(dEstimate, dMinMoleFraction)
                dTrialFraction(j) = dEstimate
                lEstimated(j) = .TRUE.
                dTraceSum = dTraceSum + dEstimate
                nEstimatedSpecies = nEstimatedSpecies + 1
            else
                dNonTraceSum = dNonTraceSum + dOldFraction
            end if
        else
            dNonTraceSum = dNonTraceSum + dOldFraction
        end if
    end do

    if (nEstimatedSpecies <= 0) return
    if (nActiveSpecies <= nEstimatedSpecies) return
    if (dTraceSum >= 0.5D0) return
    if (dNonTraceSum <= 0D0) return

    dScale = (1D0 - dTraceSum) / dNonTraceSum
    do j = 1, nVar
        i = iFirstSUB + j - 1
        if (lSubMinTraceInactive(j)) then
            dTrialFraction(j) = 0D0
        else if (.NOT.lEstimated(j)) then
            dTrialFraction(j) = dMolFraction(i) * dScale
        end if
    end do

    dMaxChange = MAXVAL(DABS(dTrialFraction - dMolFraction(iFirstSUB:iLastSUB)))
    if (dMaxChange <= 1D-14) return

    dMolFraction(iFirstSUB:iLastSUB) = dTrialFraction
    call SubMinChemicalPotential(iSolnPhaseIndex)
    call SubMinDrivingForce
    call UpdateSubMinPotentialVector
    call SubMinFunctionNorm(iSolnPhaseIndex)

    iSubMinHenrianTraceEstimateCount = iSubMinHenrianTraceEstimateCount + nEstimatedSpecies
    lChanged = .TRUE.

    return

end subroutine ApplySubMinHenrianTraceEstimate


!> \brief   Remove tiny active endmembers from the classic submin Newton solve.
!!
!! \details Marks phase-local trace endmembers inactive when their mole
!! fraction falls below dTraceSpeciesRemoveFraction.  This is a PEA candidate
!! repair used to solve the reduced phase constitution more cleanly; it is not
!! a global mass-conserving Lagrangian trace operation.
!!
!! \param[in]  iSolnPhaseIndex  Absolute solution phase index being subminimized.
!! \param[out] lChanged         True when one or more endmembers were removed.



subroutine RemoveSubMinTraceSpeciesBelowThreshold(iSolnPhaseIndex, lChanged)

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer, intent(in) :: iSolnPhaseIndex
    logical, intent(out) :: lChanged
    integer :: i, j, nActiveSpecies

    lChanged = .FALSE.

    if (.NOT.lTraceSpeciesControlEnabled) return
    if (lPostProcess) return
    if (.NOT.allocated(lSubMinTraceInactive)) return

    nActiveSpecies = 0
    do j = 1, nVar
        if (.NOT.lSubMinTraceInactive(j)) nActiveSpecies = nActiveSpecies + 1
    end do
    if (nActiveSpecies <= 1) return

    do j = 1, nVar
        if (lSubMinTraceInactive(j)) cycle
        if (allocated(lSubMinTraceReinjected)) then
            if (lSubMinTraceReinjected(j)) cycle
        end if

        i = iFirstSUB + j - 1
        if ((dMolFraction(i) > 0D0).AND.&
            (dMolFraction(i) <= dTraceSpeciesRemoveFraction).AND.&
            (nActiveSpecies > 1)) then
            lSubMinTraceInactive(j) = .TRUE.
            dMolFraction(i) = 0D0
            nActiveSpecies = nActiveSpecies - 1
            lChanged = .TRUE.
        end if
    end do

    if (lChanged) then
        call NormalizeSubMinTraceActiveSpecies
        call SubMinChemicalPotential(iSolnPhaseIndex)
        call SubMinDrivingForce
        call UpdateSubMinPotentialVector
        call SubMinFunctionNorm(iSolnPhaseIndex)
        iSubMinTraceSlowProgressCount = 0
        iterSubMinTraceLastRemoval = iterSub
        dSubMinTraceReducedNormLast = dMaxPotentialVector
    end if

    return

end subroutine RemoveSubMinTraceSpeciesBelowThreshold


!> \brief   Reinsert phase-local trace endmembers into classic submin.
!!
!! \details Reintroduces trace-inactive endmembers at
!! dTraceSpeciesReinjectFraction and rescales the remaining active fractions
!! inside the current candidate phase.  Reinserted endmembers are marked so
!! they are not immediately removed again.
!!
!! \param[in]  iSolnPhaseIndex  Absolute solution phase index being subminimized.
!! \param[out] lChanged         True when one or more endmembers were reinserted.



subroutine ReinjectSubMinTraceSpecies(iSolnPhaseIndex, lChanged)

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    integer, intent(in) :: iSolnPhaseIndex
    logical, intent(out) :: lChanged
    integer :: i, j, nInactiveSpecies
    real(8) :: dActiveSum, dTotalInject, dScale

    lChanged = .FALSE.

    if (.NOT.lTraceSpeciesControlEnabled) return
    if (lPostProcess) return
    if (.NOT.allocated(lSubMinTraceInactive)) return
    if (.NOT.ANY(lSubMinTraceInactive)) return

    nInactiveSpecies = COUNT(lSubMinTraceInactive)
    dTotalInject = DFLOAT(nInactiveSpecies) * dTraceSpeciesReinjectFraction
    if (dTotalInject <= 0D0) return
    if (dTotalInject >= 0.5D0) return

    dActiveSum = 0D0
    do j = 1, nVar
        i = iFirstSUB + j - 1
        if (.NOT.lSubMinTraceInactive(j)) dActiveSum = dActiveSum + dMolFraction(i)
    end do
    if (dActiveSum <= 0D0) return

    dScale = (1D0 - dTotalInject) / dActiveSum
    do j = 1, nVar
        i = iFirstSUB + j - 1
        if (lSubMinTraceInactive(j)) then
            dMolFraction(i) = dTraceSpeciesReinjectFraction
            lSubMinTraceInactive(j) = .FALSE.
            if (allocated(lSubMinTraceReinjected)) lSubMinTraceReinjected(j) = .TRUE.
            lChanged = .TRUE.
        else
            dMolFraction(i) = dMolFraction(i) * dScale
        end if
    end do

    if (lChanged) then
        call SubMinChemicalPotential(iSolnPhaseIndex)
        call SubMinDrivingForce
        call UpdateSubMinPotentialVector
        call SubMinFunctionNorm(iSolnPhaseIndex)
        iSubMinTraceSlowProgressCount = 0
        iterSubMinTraceLastReinject = iterSub
        dSubMinTraceReducedNormLast = HUGE(1D0)
    end if

    return

end subroutine ReinjectSubMinTraceSpecies


!> \brief   Decide whether reduced classic submin progress has stalled.
!!
!! \details Requests trace reinjection once the reduced active-set residual is
!! already small or its relative improvement has been below the slow-progress
!! threshold for the configured window.
!!
!! \param[out] lShouldReinject  True when trace endmembers should be reinserted.



subroutine ShouldReinjectSubMinTraceSpeciesForSlowProgress(lShouldReinject)

    USE ModuleGEMSolver
    USE ModuleSubMin

    implicit none

    logical, intent(out) :: lShouldReinject
    real(8) :: dDenom, dRelativeImprovement

    lShouldReinject = .FALSE.

    if (.NOT.lTraceSpeciesControlEnabled) return
    if (lPostProcess) return
    if (.NOT.allocated(lSubMinTraceInactive)) return
    if (.NOT.ANY(lSubMinTraceInactive)) return
    if ((iterSub - iterSubMinTraceLastRemoval) < iTraceSpeciesMinIterBeforeReinject) return
    if (iterSubMinTraceLastReinject > 0) then
        if ((iterSub - iterSubMinTraceLastReinject) < iTraceSpeciesMinIterBeforeReinject) return
    end if

    if (dMaxPotentialVector < dTraceSpeciesReducedSolveTolerance) then
        lShouldReinject = .TRUE.
        return
    end if

    if ((dSubMinTraceReducedNormLast <= 0D0).OR.&
        (dSubMinTraceReducedNormLast >= 0.1D0*HUGE(1D0))) then
        dSubMinTraceReducedNormLast = dMaxPotentialVector
        iSubMinTraceSlowProgressCount = 0
        return
    end if

    dDenom = DMAX1(DABS(dSubMinTraceReducedNormLast), 1D-300)
    dRelativeImprovement = (dSubMinTraceReducedNormLast - dMaxPotentialVector) / dDenom

    if (dRelativeImprovement < dTraceSpeciesSlowProgressTolerance) then
        iSubMinTraceSlowProgressCount = iSubMinTraceSlowProgressCount + 1
    else
        iSubMinTraceSlowProgressCount = 0
    end if

    dSubMinTraceReducedNormLast = dMaxPotentialVector

    if (iSubMinTraceSlowProgressCount >= iTraceSpeciesSlowProgressWindow) then
        lShouldReinject = .TRUE.
    end if

    return

end subroutine ShouldReinjectSubMinTraceSpeciesForSlowProgress
