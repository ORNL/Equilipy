subroutine InitTraceSpeciesMask

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    lTraceSpeciesControlEnabled = .TRUE.
    iTraceSpeciesSlowProgressCount = 0
    iterTraceSpeciesLastRemoval = 0
    iterTraceSpeciesLastReinject = 0
    dTraceSpeciesReducedNormLast = HUGE(1D0)

    if (allocated(lTraceSpeciesInactive)) deallocate(lTraceSpeciesInactive)
    allocate(lTraceSpeciesInactive(nSpecies))
    lTraceSpeciesInactive = .FALSE.

    if (allocated(lTraceSpeciesReinjected)) deallocate(lTraceSpeciesReinjected)
    allocate(lTraceSpeciesReinjected(nSpecies))
    lTraceSpeciesReinjected = .FALSE.

    return

end subroutine InitTraceSpeciesMask


subroutine NormalizeTraceActiveSolutionPhases

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, m, n, iPhaseSlot
    real(8) :: dPhaseMoles

    if (.NOT.allocated(lTraceSpeciesInactive)) return

    do j = 1, nSolnPhases
        iPhaseSlot = nElements - j + 1
        k = -iAssemblage(iPhaseSlot)
        if (k <= 0) cycle

        m = nSpeciesPhase(k-1) + 1
        n = nSpeciesPhase(k)
        dPhaseMoles = SUM(dMolesSpecies(m:n))
        dMolesPhase(iPhaseSlot) = dPhaseMoles

        if (dPhaseMoles > 0D0) then
            do i = m, n
                if (lTraceSpeciesInactive(i)) then
                    dMolFraction(i) = 0D0
                else
                    dMolFraction(i) = dMolesSpecies(i) / dPhaseMoles
                end if
            end do
        end if
    end do

    return

end subroutine NormalizeTraceActiveSolutionPhases


!> \brief Preserve explicitly reinjected trace fractions after model projection.
!!
!! \details CEF property evaluation can rebuild local endmember fractions from
!! site fractions.  This helper reapplies the explicit trace-reinjection
!! fraction and rescales the remaining active endmembers phase-locally so the
!! next reduced Lagrangian iteration starts from the intended reinjected state.
!!
subroutine NormalizeTraceReinjectedSolutionPhases

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: i, j, k, m, n, iPhaseSlot, nReinjected
    real(8) :: dPhaseMoles, dTraceTotal, dActiveSum, dScale

    if (.NOT.allocated(lTraceSpeciesReinjected)) return
    if (.NOT.ANY(lTraceSpeciesReinjected)) return

    do j = 1, nSolnPhases
        iPhaseSlot = nElements - j + 1
        k = -iAssemblage(iPhaseSlot)
        if (k <= 0) cycle

        m = nSpeciesPhase(k-1) + 1
        n = nSpeciesPhase(k)
        nReinjected = COUNT(lTraceSpeciesReinjected(m:n))
        if (nReinjected <= 0) cycle

        dTraceTotal = DBLE(nReinjected) * dTraceSpeciesReinjectFraction
        if (dTraceTotal >= 1D0) cycle

        dActiveSum = 0D0
        do i = m, n
            if (.NOT.lTraceSpeciesReinjected(i)) then
                dActiveSum = dActiveSum + DMAX1(dMolFraction(i), 0D0)
            end if
        end do
        if (dActiveSum <= 0D0) cycle

        dScale = (1D0 - dTraceTotal) / dActiveSum
        do i = m, n
            if (lTraceSpeciesReinjected(i)) then
                dMolFraction(i) = dTraceSpeciesReinjectFraction
            else
                dMolFraction(i) = DMAX1(dMolFraction(i), 0D0) * dScale
            end if
        end do

        dPhaseMoles = dMolesPhase(iPhaseSlot)
        if (dPhaseMoles <= 0D0) dPhaseMoles = SUM(dMolesSpecies(m:n))
        if (dPhaseMoles > 0D0) then
            do i = m, n
                dMolesSpecies(i) = dPhaseMoles * dMolFraction(i)
            end do
        end if
    end do

    return

end subroutine NormalizeTraceReinjectedSolutionPhases


!> \brief Remove tiny active solution endmembers from the legacy GEM reduced solve.
!!
!! \details Trace removal accelerates the classic active-set Lagrangian solve,
!! but SUBG/SUBQ pair fractions are skipped because low-fraction pairs can
!! carry the MQM pair-exchange residual.
!!
subroutine RemoveTraceSpeciesBelowThreshold(lChanged)

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TraceSpeciesControl.f90
    !> \brief   Remove trace species from non-MQM active solution phases during Lagrangian GEM.
    !> \author  S.Y. Kwon
    !> \date    Jun. 24, 2026
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Removed trace species without breaking coupled SUBG and SUBQ phase constitutions.
    !
    ! Purpose:
    ! ========
    !
    !> \details Temporarily removes tiny endmembers from the legacy active-set
    !! Lagrangian solve when the solution model can tolerate a reduced
    !! endmember set.  SUBG/SUBQ phases are skipped because their pair fractions
    !! are constrained by internal exchange equilibria rather than only by
    !! elemental mass balance.
    !
    ! Required input variables:
    ! =========================
    !
    ! iAssemblage                    Current active phase assemblage.
    ! dMolesSpecies                  Active species amounts.
    ! dMolesPhase                    Active phase amounts.
    ! cSolnPhaseType                 Solution-model type for each active solution phase.
    ! lTraceSpeciesInactive          Current trace-active mask.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] lChanged          True when at least one endmember is removed.
    !
    ! lTraceSpeciesInactive          Marks removed species.
    ! dMolesSpecies, dMolFraction    Removed species are set to zero for the reduced solve.
    ! dGEMFunctionNorm               Recomputed when a removal occurs.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! NormalizeTraceActiveSolutionPhases  Renormalizes active species after removal.
    ! CompChemicalPotential               Recomputes model potentials.
    ! CompFunctionNorm                    Recomputes the GEM residual.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! RunLagrangianGEM                    Calls after the legacy non-CEF line search.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - CEF site-fraction Lagrangian does not use this trace deletion path.
    ! - SUBG/SUBQ pair fractions are not removed here because low-fraction
    !   pairs can be required by the MQM exchange stationarity condition.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lChanged
    integer :: i, j, k, m, n, iPhaseSlot, nActiveSpecies
    real(8) :: dPhaseMoles, dSpeciesFraction
    logical :: lCompEverything

    lChanged = .FALSE.
    lCompEverything = .FALSE.

    if (.NOT.lTraceSpeciesControlEnabled) return
    if (lPostProcess) return
    if (.NOT.allocated(lTraceSpeciesInactive)) return

    do j = 1, nSolnPhases
        iPhaseSlot = nElements - j + 1
        k = -iAssemblage(iPhaseSlot)
        if (k <= 0) cycle
        if ((TRIM(cSolnPhaseType(k)) == 'SUBG').OR.&
            (TRIM(cSolnPhaseType(k)) == 'SUBQ')) cycle

        m = nSpeciesPhase(k-1) + 1
        n = nSpeciesPhase(k)
        dPhaseMoles = dMolesPhase(iPhaseSlot)
        if (dPhaseMoles <= 0D0) cycle

        nActiveSpecies = 0
        do i = m, n
            if (.NOT.lTraceSpeciesInactive(i)) nActiveSpecies = nActiveSpecies + 1
        end do
        if (nActiveSpecies <= 1) cycle

        do i = m, n
            if (lTraceSpeciesInactive(i)) cycle
            if (allocated(lTraceSpeciesReinjected)) then
                if (lTraceSpeciesReinjected(i)) cycle
            end if
            dSpeciesFraction = dMolesSpecies(i) / dPhaseMoles
            if ((dSpeciesFraction > 0D0).AND.&
                (dSpeciesFraction < dTraceSpeciesRemoveFraction).AND.&
                (nActiveSpecies > 1)) then
                lTraceSpeciesInactive(i) = .TRUE.
                dMolesSpecies(i) = 0D0
                dMolFraction(i) = 0D0
                iGEMTraceRemoveSpecies = i
                iGEMTraceRemovePhase = k
                iGEMTraceRemoveCount = iGEMTraceRemoveCount + 1
                nActiveSpecies = nActiveSpecies - 1
                lChanged = .TRUE.
            end if
        end do
    end do

    if (lChanged) then
        call NormalizeTraceActiveSolutionPhases
        call CompChemicalPotential(lCompEverything)
        call NormalizeTraceActiveSolutionPhases
        call CompFunctionNorm
        do i = 1, nSpecies
            if (lTraceSpeciesInactive(i)) then
                dMolesSpecies(i) = 0D0
                dMolFraction(i) = 0D0
            end if
        end do
        iTraceSpeciesSlowProgressCount = 0
        iterTraceSpeciesLastRemoval = iterGlobal
        dTraceSpeciesReducedNormLast = dGEMFunctionNorm
    end if

    return

end subroutine RemoveTraceSpeciesBelowThreshold


!> \brief Reinsert previously removed trace endmembers into active solution phases.
!!
!! \details Reinjection returns the legacy reduced Lagrangian solve to the full
!! endmember set after the non-trace variables have progressed.  Element-bearing
!! endmembers borrow a trace amount from an active endmember with identical
!! elemental stoichiometry; all-vacancy endmembers are added directly because
!! they carry no conserved elements.
!!
subroutine ReinjectTraceSpecies(lChanged)

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    TraceSpeciesControl.f90
    !> \brief   Reinsert trace species removed from the legacy Lagrangian GEM solve.
    !> \author  S.Y. Kwon
    !> \date    Jul. 02, 2026
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Reinjected trace species while preserving explicit solution fractions and diagnostic identity.
    !
    ! Purpose:
    ! ========
    !
    !> \details Restore previously removed trace endmembers before checking the
    !! final active-set Lagrangian solution.  CEF solution phases can rebuild
    !! endmember fractions during property evaluation, so this routine reapplies
    !! the explicit reinjection fraction after recomputing chemical potentials.
    !
    ! Required input variables:
    ! =========================
    !
    ! iAssemblage                    Current active phase assemblage.
    ! dMolesSpecies, dMolesPhase     Current species and phase amounts.
    ! lTraceSpeciesInactive          Species removed by the trace-control pass.
    ! dTraceSpeciesReinjectFraction  Fraction used to seed removed trace species.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] lChanged          True when at least one trace species is reinserted.
    !
    ! lTraceSpeciesInactive          Reinserted species are restored to the active set.
    ! lTraceSpeciesReinjected        Marks species seeded by this routine.
    ! dMolesSpecies, dMolFraction    Updated phase-local endmember amounts/fractions.
    ! dGEMFunctionNorm               Recomputed when reinjection occurs.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! CompChemicalPotential                 Recomputes model potentials.
    ! NormalizeTraceReinjectedSolutionPhases Reapplies the explicit reinjected fraction.
    ! CompFunctionNorm                      Recomputes the GEM residual.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! RunLagrangianGEM                      Calls after reduced-solve convergence or stagnation.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Element-bearing traces are only reinserted by mass-neutral transfer
    !   from an active endmember with identical elemental stoichiometry.
    ! - All-vacancy endmembers can be seeded directly because they do not
    !   change conserved element balances.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lChanged
    integer :: i, ii, j, k, m, n, iPhaseSlot, nInactiveSpecies, iCompensate
    real(8) :: dPhaseMoles, dBasePhaseMoles, dInjectMoles, dStoichMax, dStoichDiff
    logical :: lCompEverything
    logical :: lPhaseChanged

    lChanged = .FALSE.
    lCompEverything = .FALSE.

    if (.NOT.lTraceSpeciesControlEnabled) return
    if (lPostProcess) return
    if (.NOT.allocated(lTraceSpeciesInactive)) return
    if (.NOT.ANY(lTraceSpeciesInactive)) return

    do j = 1, nSolnPhases
        iPhaseSlot = nElements - j + 1
        k = -iAssemblage(iPhaseSlot)
        if (k <= 0) cycle

        m = nSpeciesPhase(k-1) + 1
        n = nSpeciesPhase(k)
        nInactiveSpecies = COUNT(lTraceSpeciesInactive(m:n))
        if (nInactiveSpecies <= 0) cycle

        dBasePhaseMoles = dMolesPhase(iPhaseSlot)
        if (dBasePhaseMoles <= 0D0) dBasePhaseMoles = SUM(dMolesSpecies(m:n))
        if (dBasePhaseMoles <= 0D0) cycle

        lPhaseChanged = .FALSE.
        do i = m, n
            if (.NOT.lTraceSpeciesInactive(i)) cycle

            dInjectMoles = dBasePhaseMoles * dTraceSpeciesReinjectFraction
            if (dInjectMoles <= 0D0) cycle

            ! All-vacancy endmembers do not carry elements.  Add their trace
            ! moles directly, then let the phase amount grow on renormalization.
            dStoichMax = MAXVAL(DABS(dStoichSpecies(i,1:nElements)))
            if (dStoichMax < 1D-30) then
                dMolesSpecies(i) = dInjectMoles
                lTraceSpeciesInactive(i) = .FALSE.
                if (allocated(lTraceSpeciesReinjected)) lTraceSpeciesReinjected(i) = .TRUE.
                iGEMTraceReinjectSpecies = i
                iGEMTraceReinjectPhase = k
                iGEMTraceReinjectCount = iGEMTraceReinjectCount + 1
                lPhaseChanged = .TRUE.
                lChanged = .TRUE.
                cycle
            end if

            ! For element-bearing trace species, only reinject when an active
            ! species with identical elemental stoichiometry can donate the
            ! trace amount.  This keeps element mass exactly unchanged.
            iCompensate = 0
            do ii = m, n
                if (ii == i) cycle
                if (lTraceSpeciesInactive(ii)) cycle
                if (dMolesSpecies(ii) <= dInjectMoles + dTolerance(8)) cycle

                dStoichDiff = MAXVAL(DABS(dStoichSpecies(ii,1:nElements) - dStoichSpecies(i,1:nElements)))
                if (dStoichDiff < 1D-12) then
                    iCompensate = ii
                    exit
                end if
            end do

            if (iCompensate > 0) then
                dMolesSpecies(iCompensate) = dMolesSpecies(iCompensate) - dInjectMoles
                dMolesSpecies(i) = dInjectMoles
                lTraceSpeciesInactive(i) = .FALSE.
                if (allocated(lTraceSpeciesReinjected)) lTraceSpeciesReinjected(i) = .TRUE.
                iGEMTraceReinjectSpecies = i
                iGEMTraceReinjectPhase = k
                iGEMTraceReinjectCount = iGEMTraceReinjectCount + 1
                lPhaseChanged = .TRUE.
                lChanged = .TRUE.
            end if
        end do

        if (lPhaseChanged) then
            dPhaseMoles = SUM(dMolesSpecies(m:n))
            dMolesPhase(iPhaseSlot) = dPhaseMoles
            if (dPhaseMoles > 0D0) then
                do i = m, n
                    dMolFraction(i) = dMolesSpecies(i) / dPhaseMoles
                end do
            end if
        end if
    end do

    if (lChanged) then
        iTraceSpeciesSlowProgressCount = 0
        iterTraceSpeciesLastReinject = iterGlobal
        dTraceSpeciesReducedNormLast = HUGE(1D0)
        call CompChemicalPotential(lCompEverything)
        call NormalizeTraceReinjectedSolutionPhases
        call CompFunctionNorm
        call NormalizeTraceReinjectedSolutionPhases
    end if

    return

end subroutine ReinjectTraceSpecies


subroutine ShouldReinjectTraceSpeciesForSlowProgress(lShouldReinject)

    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lShouldReinject
    real(8) :: dDenom, dRelativeImprovement

    lShouldReinject = .FALSE.

    if (.NOT.lTraceSpeciesControlEnabled) return
    if (lPostProcess) return
    if (.NOT.allocated(lTraceSpeciesInactive)) return
    if (.NOT.ANY(lTraceSpeciesInactive)) return
    if ((iterGlobal - iterTraceSpeciesLastRemoval) < iTraceSpeciesMinIterBeforeReinject) return
    if (iterTraceSpeciesLastReinject > 0) then
        if ((iterGlobal - iterTraceSpeciesLastReinject) < iTraceSpeciesMinIterBeforeReinject) return
    end if

    if (dGEMFunctionNorm < dTraceSpeciesReducedSolveTolerance) then
        lShouldReinject = .TRUE.
        return
    end if

    if ((dTraceSpeciesReducedNormLast <= 0D0).OR.&
        (dTraceSpeciesReducedNormLast >= 0.1D0*HUGE(1D0))) then
        dTraceSpeciesReducedNormLast = dGEMFunctionNorm
        iTraceSpeciesSlowProgressCount = 0
        return
    end if

    dDenom = DMAX1(DABS(dTraceSpeciesReducedNormLast), 1D-300)
    dRelativeImprovement = (dTraceSpeciesReducedNormLast - dGEMFunctionNorm) / dDenom

    if (dRelativeImprovement < dTraceSpeciesSlowProgressTolerance) then
        iTraceSpeciesSlowProgressCount = iTraceSpeciesSlowProgressCount + 1
    else
        iTraceSpeciesSlowProgressCount = 0
    end if

    dTraceSpeciesReducedNormLast = dGEMFunctionNorm

    if (iTraceSpeciesSlowProgressCount >= iTraceSpeciesSlowProgressWindow) then
        lShouldReinject = .TRUE.
    end if

    return

end subroutine ShouldReinjectTraceSpeciesForSlowProgress
