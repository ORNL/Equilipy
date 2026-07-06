!-------------------------------------------------------------------------------------------------------------
!
!> \file    SubMinInit.f90
!> \brief   Initialize shared state for one solution-phase Subminimization call.
!> \author  M.H.A. Piro
!> \date    Aug. 21, 2012
!> \sa      Subminimization.f90
!> \sa      SubMinSiteFractionCEF.f90
!> \sa      SubMinNewton.f90
!> \sa      SubMinTraceSpeciesControl.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   08/21/2012      M.H.A. Piro         Original Subminimization initialization routine
!   06/24/2026      S.Y. Kwon           Stored starting fractions for max-iteration diagnostics
!   06/24/2026      S.Y. Kwon           Documented shared CEF/classic Subminimization setup
!   06/24/2026      S.Y. Kwon           Reset CEF path diagnostics for each Subminimization call
!   06/26/2026      S.Y. Kwon           Reset the per-phase PEA candidate status before each submin call.
!   06/27/2026      S.Y. Kwon           Reset Henrian trace-estimate diagnostics for each submin call.
!   06/28/2026      S.Y. Kwon           Reset analytical-Hessian Newton diagnostics for each submin call.
!
!
! Purpose:
! ========
!
!> \details This routine prepares the module state for one solution-phase
!! Subminimization call.  It is shared by the CEF site-fraction path and the
!! classic endmember Newton path.  It clamps starting endmember fractions to the
!! phase-local trace floor, normalizes the phase, computes dChemicalPotentialStar
!! from the current elemental potentials, evaluates the initial chemical
!! potentials/driving force, and resets trace and Newton diagnostics.
!
!
! Required input variables:
! =========================
!
!> \param[in] iSolnPhaseIndex  Absolute solution phase index being subminimized.
!
! iFirstSUB            Absolute index of first species in phase.
! iLastSUB             Absolute index of last species in phase.
! nVar                 Number of endmembers/species in the phase.
! dElementPotential    Current element-potential plane.
! dMolFraction         Starting endmember fractions for this phase.
! iPhaseElectronID     Determines whether an extra ionic charge row is needed.
!
!
! Output/updated variables:
! =========================
!
!> \param[out] iterSubMax  Maximum classic Subminimization iterations.
!
! dChemicalPotentialStar    Element-potential projection for each endmember.
! dRHS                     Newton direction array.
! dHessian, iHessian        Classic Newton KKT matrix and DSYSV pivots.
! dPotentialVector          Active chemical-potential residual vector.
! dSubMinInitialMolFraction Normalized starting fractions for diagnostics.
! lSubMinTraceInactive      Phase-local trace-active mask.
! lSubMinTraceReinjected    Reinject-cycle guard for trace endmembers.
!
!
! Called subroutines/functions:
! =============================
!
! SubMinChemicalPotential      Computes initial phase chemical potentials.
! SubMinDrivingForce           Computes initial phase driving force.
! UpdateSubMinPotentialVector  Computes initial residual vector.
!
!
! Primary callers:
! ================
!
! Subminimization              Calls this before selecting CEF or classic submin path.
!
!
! Numerical assumptions:
! ======================
!
! - dChemicalPotentialStar uses the same per-particle scaling as the classic
!   endmember Newton equations.
! - All starting endmember fractions are made positive before normalization so
!   log-based Gibbs models can be evaluated.
! - The trace masks are allocated for both paths because CEF may decline to
!   handle the phase, after which the classic trace-aware path runs.
!
!-------------------------------------------------------------------------------------------------------------



subroutine SubMinInit(iSolnPhaseIndex,iterSubMax)
    USE ModuleThermo
    USE ModuleSubMin
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, iSolnPhaseIndex, iterSubMax

    dDrivingForceLast = 10D0
    lSubMinConverged  = .FALSE.
    dDrivingForceSoln(iSolnPhaseIndex) = 0D0
    if (allocated(iSubMinCandidateStatusSoln)) then
        if ((iSolnPhaseIndex >= 1).AND.(iSolnPhaseIndex <= SIZE(iSubMinCandidateStatusSoln))) then
            iSubMinCandidateStatusSoln(iSolnPhaseIndex) = SUBMIN_CANDIDATE_UNKNOWN
        end if
    end if

    if (allocated(dChemicalPotentialStar)) deallocate(dChemicalPotentialStar)
    if (allocated(dRHS))                   deallocate(dRHS)
    if (allocated(dHessian))               deallocate(dHessian)
    if (allocated(iHessian))               deallocate(iHessian)
    if (allocated(dPotentialVector))        deallocate(dPotentialVector)
    if (allocated(dSubMinInitialMolFraction)) deallocate(dSubMinInitialMolFraction)
    if (allocated(lSubMinTraceInactive))    deallocate(lSubMinTraceInactive)
    if (allocated(lSubMinTraceReinjected)) deallocate(lSubMinTraceReinjected)

    if (iPhaseElectronID(iSolnPhaseIndex) == 0) then
        i = 1
    else
        i = 2
    end if

    allocate(dChemicalPotentialStar(nVar), dPotentialVector(nVar))
    allocate(dSubMinInitialMolFraction(nVar))
    allocate(dRHS(nVar+i), iHessian(nVar+i), dHessian(nVar+i,nVar+i))
    allocate(lSubMinTraceInactive(nVar), lSubMinTraceReinjected(nVar))

    dRHS                   = 0D0
    dChemicalPotentialStar = 0D0
    dSubminGibbsEst        = 0D0
    dPotentialVector       = 0D0
    lSubMinTraceInactive   = .FALSE.
    lSubMinTraceReinjected = .FALSE.
    iSubMinTraceSlowProgressCount  = 0
    iSubMinHenrianTraceEstimateCount = 0
    iSubMinNewtonDSYSVInfo         = 0
    iterSubMinTraceLastRemoval = 0
    iterSubMinTraceLastReinject = 0
    dSubMinTraceReducedNormLast = HUGE(1D0)
    dSubMinNewtonSymmetryResidual = 0D0
    dSubMinInitialMolFraction = 0D0
    iSubMinCEFPhaseIndex = 0
    iSubMinCEFIter = 0
    nSubMinCEFIndependent = 0
    nSubMinCEFActiveIndependent = 0
    dSubMinCEFObjective = 0D0
    dSubMinCEFGradientNorm = 0D0
    dSubMinCEFStepNorm = 0D0
    dSubMinCEFChargeResidual = 0D0
    lSubMinCEFAttempted = .FALSE.
    lSubMinCEFHandled = .FALSE.
    lSubMinCEFLineSearchAccepted = .FALSE.
    lSubMinNewtonAnalyticalHessianAttempted = .FALSE.
    lSubMinNewtonAnalyticalHessianAccepted = .FALSE.

    iterSub = 1
    iterSubMax = 1000
    dMaxPotentialTol = 1D-3

    do k = 1, nVar
        i = iFirstSUB + k - 1
        dChemicalPotentialStar(k) = 0D0
        do j = 1, nElements
            dChemicalPotentialStar(k) = dChemicalPotentialStar(k) + dElementPotential(j) &
                * dStoichSpecies(i,j)
        end do
        dChemicalPotentialStar(k) = dChemicalPotentialStar(k) / DFLOAT(iParticlesPerMole(i))

        dMolFraction(i) = DMIN1(dMolFraction(i),1D0)
        dMolFraction(i) = DMAX1(dMolFraction(i), dMinMoleFraction)
    end do

    dMolFraction(iFirstSUB:iLastSUB) = &
        dMolFraction(iFirstSUB:iLastSUB) / SUM(dMolFraction(iFirstSUB:iLastSUB))
    dSubMinInitialMolFraction = dMolFraction(iFirstSUB:iLastSUB)

    call SubMinChemicalPotential(iSolnPhaseIndex)
    call SubMinDrivingForce
    call UpdateSubMinPotentialVector

    if (lMiscibility(iSolnPhaseIndex)) then
        if (cSolnPhaseName(iSolnPhaseIndex) == cSolnPhaseName(iSolnPhaseIndex-1)) then
            iSolnPhaseIndexOther = iSolnPhaseIndex - 1
        else
            LOOP_SolnSys: do i = 1, nSolnPhasesSys
                if (i == iSolnPhaseIndex) cycle LOOP_SolnSys
                if (cSolnPhaseName(iSolnPhaseIndex) == cSolnPhaseName(i)) then
                    iSolnPhaseIndexOther = i
                    exit LOOP_SolnSys
                end if
            end do LOOP_SolnSys
        end if
    end if

    return

end subroutine SubMinInit
