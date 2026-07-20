!> \brief Enforce coherent raw status facts at Lagrangian and GEM boundaries.
!!
!! \details A nonconverged solver state cannot carry GEM_EXIT_STATUS_OK.
!! Valid compounds-only Leveling exits are explicitly marked converged before
!! the terminal status is exposed to Python.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    FinalizeGEMStatus.f90
    !> \brief   Enforce coherent raw solver status at public stage boundaries.
    !> \author  S.Y. Kwon
    !> \date    Jul. 19, 2026
    !> \sa      RunLagrangianGEM.f90
    !> \sa      GEMSolver.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Made terminal GEM convergence and failure status physically consistent on every return path.
    !
    ! Purpose:
    ! ========
    !
    !> \details Keep the raw convergence flag and GEM exit status logically
    !! consistent without changing a solver's convergence decision.
    !
    ! Required input variables:
    ! =========================
    !
    ! lConverged          Raw convergence fact from the completed solver stage.
    ! lCompbdOnly         True when Leveling/PEA returned a compounds-only state.
    ! iPEAExitStatus      Strict PEA certificate status for compounds-only exits.
    ! INFOThermo          Thermodynamic workflow error status.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! iGEMExitStatus      Nonzero whenever the terminal state is unconverged.
    ! lConverged          Set true only for a valid compounds-only terminal state.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! None.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! RunLagrangianGEM   Finalizes one fixed-active-set Lagrangian attempt.
    ! GEMSolver          Finalizes every return path from MultiPhaseMinimizer.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - These routines classify status facts only.  They do not apply a residual
    !   tolerance, relabel a failed iteration as converged, or certify a phase set.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine FinalizeLagrangianExitStatus

    USE ModuleGEMSolver

    implicit none

    if (.NOT.lConverged) then
        iGEMExitStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
    end if

    return

end subroutine FinalizeLagrangianExitStatus


subroutine FinalizeGEMStatus

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    ! A true no-solution system is settled by Leveling itself.  When solution
    ! candidates exist, compounds-only is successful only with the strict PEA
    ! certificate already recorded by the compounds-only certification path.
    if (lCompbdOnly.AND.(INFOThermo == 0).AND.&
        ((nSolnPhasesSys == 0).OR.&
        (iPEAExitStatus == PEA_EXIT_STATUS_CERTIFIED_SETTLED))) then
        lConverged = .TRUE.
        iGEMExitStatus = GEM_EXIT_STATUS_OK
        return
    end if

    if ((.NOT.lConverged).AND.(iGEMExitStatus == GEM_EXIT_STATUS_OK)) then
        if (iPEAExitStatus == PEA_EXIT_STATUS_UNSETTLED) then
            iGEMExitStatus = GEM_EXIT_STATUS_PEA_UNCERTIFIED_HANDOFF
        else
            iGEMExitStatus = GEM_EXIT_STATUS_LAGRANGIAN_UNCONVERGED
        end if
    end if

    return

end subroutine FinalizeGEMStatus
