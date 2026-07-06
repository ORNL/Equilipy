!> \brief Initialize and run classical Leveling before active-set refinement.
!!
!! \details Builds the initial Leveling active-set candidate and marks true
!! compound-only systems that do not need later PEA/Lagrangian refinement.
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
!   06/26/2026      S.Y. Kwon           Documented the Leveling wrapper.
!   06/26/2026      S.Y. Kwon           Kept one-component systems eligible for PEA phase selection.
!   06/28/2026      S.Y. Kwon           Documented that unary active-component systems are handled by
!                                       direct Leveling-to-Lagrangian conversion in MultiPhaseMinimizer.
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
    ! lCompbdOnly        True only when no solution phase refinement is needed.
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
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine RunLeveling

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    lCompbdOnly  = .False.

    call InitGEMSolver
    lPhaseChange       = .True.
    iPhaseChangeReason = PHASE_CHANGE_REASON_LEVELING_INITIAL

    ! Leveling solely based on stoichiometric compounds and solution endmembers.
    call LevelingSolver

    if(nSolnPhasesSys==0.or. &
        (nSolnPhases==0.and.MINVAL(dPhasePotential)>-1D-13.and.dGEMFunctionNorm<1D-10)) then
        lCompbdOnly = .True.
    end if

    return

end subroutine RunLeveling
