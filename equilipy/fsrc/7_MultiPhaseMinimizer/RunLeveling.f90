!> \brief Initialize and run classical Leveling before active-set refinement.
!!
!! \details Builds the initial Leveling active-set candidate and marks true
!! compound-only systems and provisional compounds-only Leveling results.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    RunLeveling.f90
    !> \brief   Run the initial classical Leveling stage for multiphase minimization.
    !> \author  S.Y. Kwon
    !> \date    Jun. 26, 2026
    !> \sa      InitGEMSolver.f90
    !> \sa      LevelingSolver.f90
    !> \sa      MultiPhaseMinimizer.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Added certified compounds-only screening and eligibility-gated static-grid seeding.
    !
    ! Purpose:
    ! ========
    !
    !> \details This routine initializes GEM state, runs classical Leveling over
    !! stoichiometric compounds and solution endmembers, and sets the initial
    !! phase-change reason consumed by MultiPhaseMinimizer.
    !
    ! Required input variables:
    ! =========================
    !
    ! nSolnPhasesSys     Number of solution phases available in the current system.
    ! nElements          Number of active conserved components after preprocessing.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! lCompbdOnly        True for systems without solution phases or for a
    !                    provisional compounds-only Leveling result that must
    !                    be certified by MultiPhaseMinimizer.
    ! lPhaseChange       Marks that initial Leveling still needs a downstream
    !                    PEA or direct Lagrangian handoff.
    ! iPhaseChangeReason Initial Leveling phase-change reason.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! InitGEMSolver      Allocates and initializes Leveling/Lagrangian state.
    ! LevelingSolver     Solves the classical Leveling problem.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! MultiPhaseMinimizer Starts the production minimization workflow here.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - A one-component system is not automatically compound-only.  Oxide
    !   pseudo-component inputs such as SiO2 may reduce to one active component
    !   while still containing solution phases.  MultiPhaseMinimizer decides
    !   whether that Leveling result needs PEA or direct Lagrangian conversion.
    ! - When solution phases exist, lCompbdOnly is only a cheap endmember-level
    !   screen.  MultiPhaseMinimizer must certify the full solution minima
    !   before treating it as a terminal result.
    ! - The provisional screen uses a local Leveling assemblage residual.  It
    !   must not overwrite the staged solver's unconverged dGEMFunctionNorm
    !   sentinel before CheckPhaseAssemblage initializes its own state.
    ! - A Leveling solution endmember is not a compound-only terminal phase,
    !   even when it spans a unary active-component system by itself.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine RunLeveling

    USE ModuleThermo
    USE ModuleThermoIO, ONLY: INFOThermo
    USE ModuleGEMSolver
    USE ModuleParseCS, ONLY: iPhaseCS
    USE GridDiscovery, ONLY: GenerateStaticGrid

    implicit none

    integer :: i, j
    real(8) :: dLevelingFunctionNorm, dNormComponent
    real(8) :: dTimerStart, dTimerStop
    logical :: lGridFrontEndRequested, lGridSeedsLevelingEligible
    logical :: lGridSolutionCoverageAvailable, lLevelingAssemblageCompoundsOnly

    ! Step 1. Initialize the physical endmember and compound Leveling problem.
    lGridFrontEndRequested = lGridFrontEndActive
    lGridSolutionCoverageAvailable = allocated(cSolnPhaseType)
    if (lGridSolutionCoverageAvailable) then
        lGridSolutionCoverageAvailable = .NOT.ANY(cSolnPhaseType == 'SUBQ')
    end if

    ! Grid seeds Leveling only when the full, unselected multiphase universe
    ! has v0.3.3 model coverage; restricted solves stay on the classical path.
    lGridSeedsLevelingEligible = lGridFrontEndRequested .AND. &
        (.NOT.ANY(iPhaseCS < 0)) .AND. (nSolnPhasesSys >= 2) .AND. &
        (.NOT.lCompbdOnly) .AND. (.NOT.lPostProcess) .AND. &
        (nChargedConstraints == 0) .AND. lGridSolutionCoverageAvailable
    lGridFrontEndActive = lGridSeedsLevelingEligible

    lCompbdOnly  = .False.
    call InitGEMSolver
    lPhaseChange       = .True.
    iPhaseChangeReason = PHASE_CHANGE_REASON_LEVELING_INITIAL

    ! Step 2. Append immutable fixed-constitution rows when static discovery is
    ! enabled.  With the flag off, Leveling sees the original endmember set.
    if (lGridFrontEndActive) then
        call cpu_time(dTimerStart)
        call GenerateStaticGrid
        call cpu_time(dTimerStop)
        dGEMTimingGridGeneration = MAX(0D0,dTimerStop-dTimerStart)
        if (INFOThermo /= 0) return
    end if

    ! Step 3. Solve the finite Leveling exchange over the prepared row pool.
    call LevelingSolver

    ! Evaluate the residual of the final classical Leveling assemblage.  The
    ! InitGEMSolver value is an unconverged sentinel and cannot gate a public
    ! compounds-only certificate path.
    dLevelingFunctionNorm = 0D0
    do j = 1, nElements
        dNormComponent = dMolesElement(j)
        do i = 1, nElements
            dNormComponent = dNormComponent - dMolesPhase(i)*dStoichSpeciesGEM(i,j)
        end do
        dLevelingFunctionNorm = dLevelingFunctionNorm + dNormComponent**2
    end do
    do i = 1, nElements
        dNormComponent = DOT_PRODUCT(dElementPotential, dAtomFractionSpeciesGEM(i,:)) - &
            dChemicalPotentialGEM(i)
        dLevelingFunctionNorm = dLevelingFunctionNorm + dNormComponent**2
    end do
    dLevelingFunctionNorm = DSQRT(dLevelingFunctionNorm)

    lLevelingAssemblageCompoundsOnly = .TRUE.
    do i = 1, nElements
        if ((iAssemblage(i) <= 0).OR.(iAssemblage(i) > nSpecies)) then
            lLevelingAssemblageCompoundsOnly = .FALSE.
            exit
        end if
        if (iPhase(iAssemblage(i)) > 0) then
            lLevelingAssemblageCompoundsOnly = .FALSE.
            exit
        end if
    end do

    ! Step 4. Bypass later solution refinement only for a certified true
    ! compound-only set; static-grid discovery never takes this fast exit.
    if(nSolnPhasesSys==0.or. &
        ((.NOT.lGridFrontEndActive).and.lLevelingAssemblageCompoundsOnly.and.nSolnPhases==0.and.&
        MINVAL(dPhasePotential)>=dToleranceLevel.and.dLevelingFunctionNorm<=dPEATol)) then
        lCompbdOnly = .True.
    end if

    return

end subroutine RunLeveling
