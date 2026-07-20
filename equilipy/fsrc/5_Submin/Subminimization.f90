!-------------------------------------------------------------------------------------------------------------
!
!> \file    Subminimization.f90
!> \brief   Minimize one inactive solution phase against the current Gibbs plane.
!> \author  M.H.A. Piro
!> \date    Aug. 21, 2012
!> \sa      CheckPhaseAssemblage.f90
!> \sa      SubMinSiteFractionCEF.f90
!> \sa      SubMinNewton.f90
!> \sa      SubMinLineSearch.f90
!> \sa      SubMinTraceSpeciesControl.f90
!
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   08/21/2012      M.H.A. Piro         Original code
    !   08/30/2012      M.H.A. Piro         Enforced line-search Wolfe conditions and miscibility duplicate checks
    !   02/28/2013      M.H.A. Piro         Added ionic phase charge-neutrality constraints
    !   11/26/2021      S.Y. Kwon           Changed driving-force basis to Joule per mole of atoms
    !   07/20/2026      S.Y. Kwon           Integrated CEF, analytical-curvature, trace, and miscibility-gap paths into robust solution-phase subminimization.
    !
    !
! Purpose:
! ========
!
!> \details The purpose of this subroutine is to minimize the phase driving
!! force for one inactive solution phase in CheckPhaseAssemblage/PEA.  A
!! positive phase-entry result means the phase has a constitution below the
!! current Gibbs plane and should be considered for the active assemblage.  CEF
!! phases such as SUBL/SUBLM/SUBOM use the site-fraction path because their
!! endmember space grows combinatorially.  A feasible negative-driving-force
!! point from any solution model is already enough evidence that the phase lies
!! below the current Gibbs plane; it is accepted as a PEA entry witness even if
!! the phase-local submin stationarity residual is not fully minimized.
!! SUBG/SUBQ/RKMP/QKTO phases use the analytical solution-Hessian branch inside
!! SubMinNewton with the same line search.  Other non-CEF solution models
!! remain on the classic endmember Newton branch.
!
!
! Required input variables:
! =========================
!
!> \param[in] iSolnPhaseIndex  Absolute solution phase index being tested.
!
! dElementPotential     Element-potential plane from the current GEM/PEA state.
! dMolFraction          Initial endmember fractions from Leveling or previous PEA estimates.
! cSolnPhaseType        Used by SubMinSiteFractionCEF to decide whether the CEF path applies.
! lMiscibility          Enables duplicate-minimum checks for miscibility-gap phases.
!
!
! Output/updated variables:
! =========================
!
!> \param[out] lPhasePass  True when the phase should be considered for addition.
!
! dDrivingForce          Subminimized phase driving force.
! dDrivingForceSoln      Stores dDrivingForce for this solution phase.
! dPhasePotential        Stores the same phase potential for endmembers in this phase.
! iSubMinCandidateStatusSoln
!                        Stores whether the latest subminimized candidate row is usable by PEA Leveling.
! lSubMinConverged       True when either CEF or classic submin path satisfies its residuals.
! lSubMinMaxHit          True when the classic path reaches the iteration guard.
!
!
! Called subroutines/functions:
! =============================
!
! SubMinInit                    Initializes phase-local arrays and residuals.
! SubMinSiteFractionCEF         Handles SUBL/SUBLM/SUBOM phases in site-fraction variables.
! SubMinNewton                  Builds classic or analytical-Hessian Newton directions.
! SubMinLineSearch              Accepts/damps classic Newton endmember updates.
! ApplySubMinHenrianTraceEstimate
!                               Updates active trace endmembers from the local Henrian approximation.
! SubMinTraceSpeciesControl     Removes/reinjects trace endmembers for classic submin.
! SubMinFunctionNorm            Updates classic submin residual norms.
! SubMinCheckDuplicate          Rejects duplicate miscibility-gap minima.
!
!
! Primary callers:
! ================
!
! CheckPhaseAssemblage          Uses subminimized driving forces during PEA active-set repair.
! CompMolFraction               Uses this routine to refine solution-phase candidates.
! CompInitMinSolnPoint          Uses this routine while discovering initial solution minima.
! CompMinSolnPoint              Uses this routine while refining solution minima.
!
!
! Numerical assumptions:
! ======================
!
! - CEF site-fraction submin owns SUBL/SUBLM/SUBOM phases.  A converged CEF
!   solve is preferred, but any unconverged feasible point with negative
!   driving force is still a valid phase-entry witness because it proves the
!   current Gibbs plane is too high for that phase.  Positive unconverged
!   candidates are rejected.
! - SUBG/SUBQ submin works in formal pair/quadruplet fractions.  RKMP/QKTO
!   submin works in one-sublattice species fractions.  Both branches project
!   the analytical Hessian to the simplex and use the same phase-local line
!   search as the classic path.
! - Classic submin works in endmember mole fractions and keeps phase-local
!   trace handling inside PEA.  This trace handling is not a global
!   mass-conserving Lagrangian trace algorithm.
! - A classic submin start exits early when the max potential residual no
!   longer changes within 1D-8.  This lets CompInitMinSolnPoint/CompMinSolnPoint
!   move on to the next start instead of spending the full iteration budget on
!   a line-search-stalled local path.
! - Active classic trace endmembers below the Henrian trace threshold can be
!   re-estimated by holding the activity coefficient fixed from the previous
!   evaluated state.  The update is accepted only while the predicted fraction
!   remains in the trace region.
! - dDrivingForce is trusted for exact phase entry when the chosen path
!   converges.  For CEF only, a negative feasible point can also pass as a
!   phase-entry witness even if the site-gradient is not fully minimized.
! - A duplicate miscibility-gap minimum is neutralized by setting the driving
!   force to a large positive value.
!
!-------------------------------------------------------------------------------------------------------------



subroutine Subminimization(iSolnPhaseIndex,lPhasePass)
    USE ModuleThermo
    USE ModuleSubMin
    USE ModuleGEMSolver
    USE ModuleThermoIO

    implicit none

    integer :: iSolnPhaseIndex, iterSubMax, j
    logical :: lPhasePass, lDuplicate, lHitSubMinMax, lHandledByCEFSubmin
    logical :: lNegativeDrivingForceWitness
    logical :: lPrintSubMinMaxWarning

    lPhasePass = .FALSE.
    lDuplicate = .FALSE.
    lNegativeDrivingForceWitness = .FALSE.

    iFirstSUB = nSpeciesPhase(iSolnPhaseIndex-1) + 1
    iLastSUB  = nSpeciesPhase(iSolnPhaseIndex)
    nVar      = iLastSUB - iFirstSUB + 1

    call SubMinInit(iSolnPhaseIndex,iterSubMax)

    call SubMinSiteFractionCEF(iSolnPhaseIndex, lHandledByCEFSubmin)
    if (.NOT.lSubMinCEFAttempted) then
        call RunClassicSubMinNewtonPath(iSolnPhaseIndex, iterSubMax, lDuplicate)
    end if

    lHitSubMinMax = ((.NOT.lSubMinConverged).AND.(INFOThermo == 0).AND.(iterSub > iterSubMax))
    lPrintSubMinMaxWarning = lDebugMode .AND. lHitSubMinMax
    if (lHitSubMinMax) then
        lSubMinMaxHit             = .TRUE.
        iSubMinMaxPhaseIndex      = iSolnPhaseIndex
        iSubMinMaxIter            = iterSubMax
        dSubMinMaxFunctionNorm    = dSubMinFunctionNorm
        dSubMinMaxPotentialVector = dMaxPotentialVector
        dSubMinMaxMinFraction     = MINVAL(dMolFraction(iFirstSUB:iLastSUB))
        dSubMinMaxMaxFraction     = MAXVAL(dMolFraction(iFirstSUB:iLastSUB))
        if (allocated(dSubMinMaxElementPotential)) deallocate(dSubMinMaxElementPotential)
        if (allocated(dSubMinMaxInitialMolFraction)) deallocate(dSubMinMaxInitialMolFraction)
        if (allocated(dSubMinMaxFinalMolFraction)) deallocate(dSubMinMaxFinalMolFraction)
        allocate(dSubMinMaxElementPotential(nElements))
        allocate(dSubMinMaxInitialMolFraction(nVar), dSubMinMaxFinalMolFraction(nVar))
        dSubMinMaxElementPotential = dElementPotential(1:nElements)
        dSubMinMaxInitialMolFraction = dSubMinInitialMolFraction
        dSubMinMaxFinalMolFraction = dMolFraction(iFirstSUB:iLastSUB)
    end if

    if (lPrintSubMinMaxWarning) then
        print *, 'WARNING: Subminimization hit max iterations ', &
            ' phase=', iSolnPhaseIndex, &
            ' phaseName=', trim(cSolnPhaseName(iSolnPhaseIndex)), &
            ' iterGlobal=', iterGlobal, &
            ' iterPEA=', iterPEA, &
            ' iterSubMax=', iterSubMax, &
            ' T=', dTemperature, &
            ' P=', dPressure, &
            ' dDrivingForce=', dDrivingForce, &
            ' dSubMinFunctionNorm=', dSubMinFunctionNorm, &
            ' dMaxPotentialVector=', dMaxPotentialVector, &
            ' minX=', MINVAL(dMolFraction(iFirstSUB:iLastSUB)), &
            ' maxX=', MAXVAL(dMolFraction(iFirstSUB:iLastSUB))
        do j = 1, nElements
            print *, 'WARNING: Subminimization condition element ', &
                ' name=', trim(cElementName(j)), &
                ' moles=', dMolesElement(j), &
                ' potential=', dElementPotential(j)
        end do
        print *, 'WARNING: Subminimization CEF path ', &
            ' attempted=', lSubMinCEFAttempted, &
            ' handled=', lSubMinCEFHandled, &
            ' acceptedLineSearch=', lSubMinCEFLineSearchAccepted, &
            ' iter=', iSubMinCEFIter, &
            ' nIndependent=', nSubMinCEFIndependent, &
            ' nActive=', nSubMinCEFActiveIndependent, &
            ' objective=', dSubMinCEFObjective, &
            ' gradientNorm=', dSubMinCEFGradientNorm, &
            ' stepNorm=', dSubMinCEFStepNorm, &
            ' chargeResidual=', dSubMinCEFChargeResidual
        print *, 'WARNING: Subminimization Newton analytical Hessian ', &
            ' attempted=', lSubMinNewtonAnalyticalHessianAttempted, &
            ' accepted=', lSubMinNewtonAnalyticalHessianAccepted, &
            ' dsysvInfo=', iSubMinNewtonDSYSVInfo, &
            ' symmetryResidual=', dSubMinNewtonSymmetryResidual
    end if

    lNegativeDrivingForceWitness = (dDrivingForce < dTolerance(4)).AND.&
        (dSubMinFunctionNorm <= 10D0*dTolerance(1))

    if (allocated(iSubMinCandidateStatusSoln)) then
        if ((iSolnPhaseIndex >= 1).AND.(iSolnPhaseIndex <= SIZE(iSubMinCandidateStatusSoln))) then
            if (lDuplicate) then
                iSubMinCandidateStatusSoln(iSolnPhaseIndex) = SUBMIN_CANDIDATE_DUPLICATE
            else if (lSubMinConverged) then
                iSubMinCandidateStatusSoln(iSolnPhaseIndex) = SUBMIN_CANDIDATE_CONVERGED
            else if (lNegativeDrivingForceWitness) then
                iSubMinCandidateStatusSoln(iSolnPhaseIndex) = SUBMIN_CANDIDATE_NEGATIVE_WITNESS
            else if (lHitSubMinMax) then
                iSubMinCandidateStatusSoln(iSolnPhaseIndex) = SUBMIN_CANDIDATE_MAX_ITER
            else
                iSubMinCandidateStatusSoln(iSolnPhaseIndex) = SUBMIN_CANDIDATE_REJECTED
            end if
        end if
    end if

    if (lDuplicate.OR.(.NOT.(lSubMinConverged.OR.lNegativeDrivingForceWitness))) dDrivingForce = 9D5
    if ((dDrivingForce < dTolerance(4)).AND.(lSubMinConverged.OR.lNegativeDrivingForceWitness)) lPhasePass = .TRUE.

    dDrivingForceSoln(iSolnPhaseIndex)  = dDrivingForce
    dPhasePotential(iFirstSUB:iLastSUB) = dDrivingForce

    if (allocated(dChemicalPotentialStar)) deallocate(dChemicalPotentialStar)
    if (allocated(dRHS)) deallocate(dRHS)
    if (allocated(dHessian)) deallocate(dHessian)
    if (allocated(iHessian)) deallocate(iHessian)
    if (allocated(dPotentialVector)) deallocate(dPotentialVector)
    if (allocated(lSubMinTraceInactive)) deallocate(lSubMinTraceInactive)
    if (allocated(lSubMinTraceReinjected)) deallocate(lSubMinTraceReinjected)
    if (allocated(dSubMinInitialMolFraction)) deallocate(dSubMinInitialMolFraction)

    return

contains

    subroutine RunClassicSubMinNewtonPath(iSolnPhaseIndexIn, iterSubMaxIn, lDuplicateOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iterSubMaxIn
        logical, intent(inout) :: lDuplicateOut

        real(8), parameter :: dMaxPotentialStagnationTol = 1D-8

        real(8) :: dMaxPotentialVectorLast

        logical :: lPairPhaseForTrace, lTraceChanged, lTraceEstimated, lTraceSlowReinject
        logical :: lPotentialStagnated

        lPairPhaseForTrace = (TRIM(cSolnPhaseType(iSolnPhaseIndexIn)) == 'SUBG').OR.&
            (TRIM(cSolnPhaseType(iSolnPhaseIndexIn)) == 'SUBQ')
        dMaxPotentialVectorLast = HUGE(1D0)
        lPotentialStagnated     = .FALSE.

        LOOP_ClassicSubMin: do iterSub = 1, iterSubMaxIn
            call SubMinNewton(iSolnPhaseIndexIn)
            call SubMinLineSearch(iSolnPhaseIndexIn)

            if (.NOT.lPairPhaseForTrace) then
                call ApplySubMinHenrianTraceEstimate(iSolnPhaseIndexIn, lTraceEstimated)
                if (lTraceEstimated) lSubMinConverged = .FALSE.
            end if

            call SubMinFunctionNorm(iSolnPhaseIndexIn)

            if (.NOT.lPairPhaseForTrace) then
                call RemoveSubMinTraceSpeciesBelowThreshold(iSolnPhaseIndexIn, lTraceChanged)
                if (lTraceChanged) cycle LOOP_ClassicSubMin
            end if

            if (MinimumActiveSubMinFraction() <= 0D0) exit LOOP_ClassicSubMin

            lPotentialStagnated = (ABS(dMaxPotentialVector - dMaxPotentialVectorLast) < &
                dMaxPotentialStagnationTol)
            dMaxPotentialVectorLast = dMaxPotentialVector

            if ((dSubMinFunctionNorm < dTolerance(1)).AND.&
                (dMaxPotentialVector < dMaxPotentialTol)) lSubMinConverged = .TRUE.

            if (lMiscibility(iSolnPhaseIndexIn)) call SubMinCheckDuplicate(lDuplicateOut)
            if (lDuplicateOut) exit LOOP_ClassicSubMin

            lTraceSlowReinject = .FALSE.
            if ((.NOT.lPairPhaseForTrace).AND.(.NOT.lSubMinConverged)) then
                call ShouldReinjectSubMinTraceSpeciesForSlowProgress(lTraceSlowReinject)
            end if

            if ((.NOT.lPairPhaseForTrace).AND.(lSubMinConverged .OR. lTraceSlowReinject)) then
                call ReinjectSubMinTraceSpecies(iSolnPhaseIndexIn, lTraceChanged)
                if (lTraceChanged) then
                    lSubMinConverged = .FALSE.
                    cycle LOOP_ClassicSubMin
                end if
            end if

            if (lSubMinConverged .OR. (INFOThermo /= 0)) exit LOOP_ClassicSubMin
            if (lPotentialStagnated) exit LOOP_ClassicSubMin
        end do LOOP_ClassicSubMin

        return

    end subroutine RunClassicSubMinNewtonPath

    real(8) function MinimumActiveSubMinFraction()

        implicit none

        integer :: i, j

        MinimumActiveSubMinFraction = HUGE(1D0)
        do j = 1, nVar
            if (allocated(lSubMinTraceInactive)) then
                if (lSubMinTraceInactive(j)) cycle
            end if
            i = iFirstSUB + j - 1
            MinimumActiveSubMinFraction = DMIN1(MinimumActiveSubMinFraction, dMolFraction(i))
        end do

        return

    end function MinimumActiveSubMinFraction

end subroutine Subminimization
