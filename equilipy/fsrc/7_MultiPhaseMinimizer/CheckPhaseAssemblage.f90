!
!
subroutine CheckPhaseAssemblage
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPhaseAssemblage.f90
    !> \brief   Gibbs Energy Minimization solver.
    !> \author  S.Y. Kwon
    !> \date    May 2, 2022
    !
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
!   05/02/2022      S.Y. Kwon            Original code
!   10/21/2023      S.Y. Kwon            Cleaning the code
!   06/26/2026      S.Y. Kwon            Added PEA-internal Lagrangian polish before refreshing
!                                         solution driving-force candidates.
!   06/26/2026      S.Y. Kwon            Allowed Lagrangian-certified PEA convergence once the polished
!                                         plane passes the refreshed driving-force check.
!   06/26/2026      S.Y. Kwon            Recorded per-iteration PEA diagnostics for active-set audits.
!   06/27/2026      S.Y. Kwon            Scaled the accepted-polish driving-force check by the polished
!                                         Lagrangian residual so approximate PEA witnesses do not force
!                                         repeated identical PEA iterations.
!   06/28/2026      S.Y. Kwon            Delayed PEA-Lagrangian polish until elemental potentials reach a
!                                         loose plateau based on dMaxElementPotential.
!   06/30/2026      S.Y. Kwon            Added debug-only max-iteration warnings for CheckPhaseAssemblage/PEA.
!   07/01/2026      S.Y. Kwon            Allowed converged PEA-Lagrangian handoffs to exit when
!                                         the refreshed driving-force check adds no phase.
!   07/01/2026      S.Y. Kwon            Added direct inactive-witness handoff after accepted PEA-Lagrangian
!                                         polish to avoid duplicate pseudo-compound Leveling oscillations.
!   07/01/2026      S.Y. Kwon            Broke two-state PEA cycles once a direct-witness Lagrangian handoff
!                                         repeats after both active-set alternatives were polished.
!   07/01/2026      S.Y. Kwon            Stopped repeated accepted handoff cycles after multiple repeats when
!                                         intervening PEA candidates do not yield an accepted Lagrangian state.
!   07/01/2026      S.Y. Kwon            Prevented solution endmember rows from being classified as
!                                         compound-only before PostProcessPEA converts them to solution phases.
!   07/02/2026      S.Y. Kwon            Added driving-force backoff after rejected embedded PEA-Lagrangian
!                                         polish attempts.
!   07/04/2026      S.Y. Kwon            Added passive C4-a counters for polish gating, direct-witness
!                                         handoff, and repeated-handoff exits.
!   07/05/2026      S.Y. Kwon            Added passive timing/counting for inactive-candidate
!                                         certification sweeps used by WS-a2.
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the quantities of species and phases at thermodynamic
    !! equilibrium using the Gibbs Energy Minimization (GEM) method.  This subroutine uses the theory developed by
    !! Capitani 1983
    !! Issue: Immiscibility distinguising problem, Subminimization issue
    !
    ! Pertinent variables:
    ! ====================
    !
    ! nConPhases            The number of pure condensed phases in the assemblage
    ! nSolnPhases           The number of solution phases in the assemblage
    ! nSolnPhasesSys        The number of solution phases in the system
    ! iAssemblage           Integer vector containing the indices of phases in the assemblage
    !                        (1:nConphases represent pure condensed phases and (nElements-nSolnPhases:nSolnPhases)
    !                        represent solution phases.
    ! INFOThermo            An integer scalar identifying whether the program exits successfully or if
    !                        it encounters an error.
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver
!
    implicit none
!
    integer                             :: i, iterLevel, nPolishAttemptBefore, iPolishAttempted
    integer                             :: iMinPotentialRow, iDirectWitnessPhaseNext
    real(8)                             :: dGibbsGlobal, dPolishDrivingForceTol
    real(8)                             :: dCertStart, dCertStop
    real(8)                             :: dPreviousMinPhasePotential
    logical                             :: lPolishEnabledSave, lHitPEAMax
    logical                             :: lPolishedRepeatedHandoff, lDirectWitnessHandoffNext
    logical                             :: lDirectWitnessHandoff
    logical                             :: lRequireDrivingForceReadyForPolish, lDrivingForceReadyForPolish
    logical                             :: lMinPotentialActive, lRepeatedHandoffExit
    logical                             :: MinPotentialPhaseIsActive
    integer                             :: PhaseIndexForLevelingRow

    ! Initialize variables
    iterPEA         = 0
    dGibbsGlobal    = 0
    dMaxElementPotential = 10
    lDirectWitnessHandoffNext = .FALSE.
    iDirectWitnessPhaseNext = 0
    dPreviousMinPhasePotential = -HUGE(1D0)
    lRequireDrivingForceReadyForPolish = .FALSE.

    ! Initializatrion: Allocate variables and calculate minimum points of each solution
    call InitCheckPhaseAssemblage
!
!
!
    LOOP_GEM: do iterPEA = 1,iterPEAMax
        ! Initialize variables related to history tracking
        iterHistoryLevel = 0
        dMolesPhaseHistory = 0d0
        dElementPotentialLast=dElementPotential
        dToleranceLevel = -1D-8
        iPEAGatePass = 0
        iPEAGateBlock = 0
        iPEAGatePotBlock = 0
        iPEAGateDFBlock = 0
        iPEADirectHandoff = 0
        iPEADirectPhase = 0
        iPEARepeatExit = 0
        iPEARepeatExitReason = 0
        
!
        ! Part 1: Leveling.  If the previous accepted PEA-Lagrangian polish
        ! exposed a favorable inactive phase after solution-minimum refresh,
        ! try that witness directly before the duplicate pseudo-compound
        ! Leveling plane can cycle back to the old collapsed active set.
        lDirectWitnessHandoff = lDirectWitnessHandoffNext
        lDirectWitnessHandoffNext = .FALSE.
        if (lDirectWitnessHandoff) then
            iPEADirectHandoff = 1
            iPEADirectPhase = iDirectWitnessPhaseNext
            nPEADirectHandoff = nPEADirectHandoff + 1
            iterLevel = 1
            call GetNewAssemblage(iterLevel)
            dPhasePotential = dLevelingChemicalPotential - MATMUL(dLevelingCompositionSpecies,dElementPotential)
            dMinPhasePotential = MINVAL(dPhasePotential)
            if (INFOThermo /= 0) then
                call cpu_time(dCertStart)
                call CompMinSolnPoint
                dPhasePotential = dLevelingChemicalPotential - MATMUL(dLevelingCompositionSpecies,dElementPotential)
                dMinPhasePotential = MINVAL(dPhasePotential)
                call cpu_time(dCertStop)
                dGEMTimingCertification = dGEMTimingCertification + MAX(0D0, dCertStop - dCertStart)
                nGEMCertificationSweep = nGEMCertificationSweep + 1
                exit LOOP_GEM
            end if
        else
            LOOP_Leveling: do iterLevel = 1, 500

                ! Calculate phase potential of each species
                dPhasePotential = dLevelingChemicalPotential - MATMUL(dLevelingCompositionSpecies,dElementPotential)
                dMinPhasePotential = MINVAL(dPhasePotential)
                ! Check global minimum: if all elements in phase potential are positive, the system is in global minimum
                if ((dMinPhasePotential> dToleranceLevel)) exit LOOP_Leveling
!
                ! Determine the next phase assemblage to be tested:
                call GetNewAssemblage(iterLevel)
!
                ! Exit if an error has been encountered:
                if (INFOThermo /= 0) then

                    call cpu_time(dCertStart)
                    call CompMinSolnPoint
                    dPhasePotential = dLevelingChemicalPotential - MATMUL(dLevelingCompositionSpecies,dElementPotential)
                    dMinPhasePotential = MINVAL(dPhasePotential)
                    call cpu_time(dCertStop)
                    dGEMTimingCertification = dGEMTimingCertification + MAX(0D0, dCertStop - dCertStart)
                    nGEMCertificationSweep = nGEMCertificationSweep + 1
                    exit LOOP_GEM
                end if
!
            end do LOOP_Leveling
        end if
       
        
!
        ! Part2: Polish the PEA Leveling plane with the current fixed active
        ! set only after the PEA elemental potentials have reached a loose
        ! plateau.  If an embedded polish rejects the current active set, PEA
        ! waits until the refreshed driving-force pass no longer shows a
        ! strongly favorable inactive phase before trying another polish.
        ! Early PEA iterations are still discovering the coarse active set;
        ! repeated rejected polishes can lock the loop into expensive
        ! handoffs.  Then calculate the minimum point of each solution phase
        ! based on the current elemental potentials.
        dMaxElementPotential = MAXVAL(DABS(dElementPotential-dElementPotentialLast))
        nPolishAttemptBefore = nPEALagrangianPolishAttempt
        lDrivingForceReadyForPolish = (.NOT.lRequireDrivingForceReadyForPolish).OR.&
            (dPreviousMinPhasePotential >= -dPEALagrangianPolishMaxDrivingForce)
        if ((dMaxElementPotential < dPEALagrangianPolishMaxElementPotential).AND.&
            lDrivingForceReadyForPolish) then
            iPEAGatePass = 1
            nPEAGatePassed = nPEAGatePassed + 1
            call PEALagrangianPolish
        else
            iPEAGateBlock = 1
            nPEAGateBlocked = nPEAGateBlocked + 1
            if (dMaxElementPotential >= dPEALagrangianPolishMaxElementPotential) then
                iPEAGatePotBlock = 1
                nPEAGatePotBlocked = nPEAGatePotBlocked + 1
            end if
            if (.NOT.lDrivingForceReadyForPolish) then
                iPEAGateDFBlock = 1
                nPEAGateDFBlocked = nPEAGateDFBlocked + 1
            end if
            lPolishEnabledSave = lPEALagrangianPolishEnabled
            lPEALagrangianPolishEnabled = .FALSE.
            call PEALagrangianPolish
            lPEALagrangianPolishEnabled = lPolishEnabledSave
        end if
        iPolishAttempted = 0
        if (nPEALagrangianPolishAttempt > nPolishAttemptBefore) iPolishAttempted = 1
        call cpu_time(dCertStart)
        call CompMinSolnPoint
        dPhasePotential = dLevelingChemicalPotential - MATMUL(dLevelingCompositionSpecies,dElementPotential)
        dMinPhasePotential = MINVAL(dPhasePotential)
        call cpu_time(dCertStop)
        dGEMTimingCertification = dGEMTimingCertification + MAX(0D0, dCertStop - dCertStart)
        nGEMCertificationSweep = nGEMCertificationSweep + 1
!
        
        ! Part3: Check convergence of Capitani algorithm
        dGibbsGlobal = dot_product(dMolesElement,dElementPotential)
        dPreviousMinPhasePotential = dMinPhasePotential
        dMaxElementPotential = MAXVAL(DABS(dElementPotential-dElementPotentialLast))
        dPolishDrivingForceTol = DMAX1(dPEATol, 10D0*dPEALagrangianPolishNormAfter)
        if (iPolishAttempted == 1) then
            lRequireDrivingForceReadyForPolish = .NOT.lPEALagrangianPolishAccepted
        end if
        lPolishedRepeatedHandoff = lPEALagrangianPolishAccepted.AND.&
            (iPEALagrangianHandoffRepeat == 1).AND.&
            (dPEALagrangianPolishNormAfter < 1D-5)
        lMinPotentialActive = MinPotentialPhaseIsActive()
        lDirectWitnessHandoffNext = lPEALagrangianPolishAccepted.AND.&
            (dPEALagrangianPolishNormAfter < 1D-5).AND.&
            (dMinPhasePotential < -dPolishDrivingForceTol).AND.&
            (.NOT.lMinPotentialActive)
        if (lDirectWitnessHandoffNext) then
            iMinPotentialRow = MAXVAL(MINLOC(dPhasePotential))
            iDirectWitnessPhaseNext = PhaseIndexForLevelingRow(iMinPotentialRow)
        else
            iDirectWitnessPhaseNext = 0
        end if
        lRepeatedHandoffExit = lPolishedRepeatedHandoff.AND.&
            ((nPEALagrangianHandoffRepeated >= 2).OR.&
            lDirectWitnessHandoff.OR.&
            (.NOT.lDirectWitnessHandoffNext).OR.&
            lMinPotentialActive)
        if (lRepeatedHandoffExit) then
            iPEARepeatExit = 1
            nPEARepeatExit = nPEARepeatExit + 1
            if (nPEALagrangianHandoffRepeated >= 2) then
                iPEARepeatExitReason = 1
            else if (lDirectWitnessHandoff) then
                iPEARepeatExitReason = 2
            else if (.NOT.lDirectWitnessHandoffNext) then
                iPEARepeatExitReason = 3
            else if (lMinPotentialActive) then
                iPEARepeatExitReason = 4
            else
                iPEARepeatExitReason = 9
            end if
        end if
        call RecordPEAIterationDiagnostics(iterLevel, iPolishAttempted)

        ! When stable phases are only composed of stoichiomatric compounds
        CheckCompdOnly: if ((dMinPhasePotential>=-1D-7).AND.(dGEMFunctionNorm<1D-7)) then 
            do i = 1, nElements
                if(iAssemblage(i)<0) exit CheckCompdOnly
                if(iAssemblage(i)>nSpecies) exit CheckCompdOnly
                if(iPhase(iAssemblage(i))>0) exit CheckCompdOnly
            end do 
            lCompbdOnly = .True.
            exit LOOP_GEM
        end if CheckCompdOnly

        ! When a single phase appears and the amount converges to unity, exit the loop
        if((MINVAL(iAssemblage)==MAXVAL(iAssemblage)).AND.(MAXVAL(dMolesPhase)>=1000)) exit LOOP_GEM

        ! Global minimum from the original PEA residual.
        if ((dMinPhasePotential >=-dPEATol).AND.&
        (dMaxElementPotential<dPEATol).AND.(dGEMFunctionNorm<dPEATol)) exit LOOP_GEM
!
        ! Lagrangian-polished PEA convergence.  If the embedded Lagrangian
        ! solve converged and the refreshed solution-phase minima do not sit
        ! below the polished plane, PEA has no more active-set evidence to add.
        if (lPEALagrangianPolishAccepted.AND.&
            (dMinPhasePotential >= -dPolishDrivingForceTol).AND.&
            (dPEALagrangianPolishNormAfter < 1D-5)) exit LOOP_GEM
!
        ! If the embedded Lagrangian solve converges and Leveling hands the
        ! same active set back to Lagrangian, PEA has no new assemblage to add.
        ! Requiring elemental-potential plateauing in this case can force the
        ! loop to repeat the same certified active set until iterPEAMax.
        if (lRepeatedHandoffExit) exit LOOP_GEM
!
    end do LOOP_GEM
    lHitPEAMax = (INFOThermo == 0).AND.(iterPEA > iterPEAMax)
    if (lDebugMode.AND.lHitPEAMax) then
        print *, 'WARNING: CheckPhaseAssemblage/PEA hit max iterations ', &
            ' iterPEAMax=', iterPEAMax, &
            ' T=', dTemperature, &
            ' P=', dPressure, &
            ' dMinPhasePotential=', dMinPhasePotential, &
            ' dMaxElementPotential=', dMaxElementPotential, &
            ' dGEMFunctionNorm=', dGEMFunctionNorm, &
            ' dPEATol=', dPEATol, &
            ' nPEARecorded=', nPEARecorded
        print *, 'WARNING: CheckPhaseAssemblage/PEA polish summary ', &
            ' attempts=', nPEALagrangianPolishAttempt, &
            ' accepted=', nPEALagrangianPolishAccepted, &
            ' rejected=', nPEALagrangianPolishRejected, &
            ' lastReason=', iPEALagrangianPolishReason, &
            ' lastIterGlobal=', iPEALagrangianPolishIterGlobal, &
            ' lastNormAfter=', dPEALagrangianPolishNormAfter
    end if
    call PostProcessPEA
    
!
!
    return
!
end subroutine CheckPhaseAssemblage
!
!

!> \brief Return whether the most favorable Leveling row already belongs to the active Lagrangian set.
logical function MinPotentialPhaseIsActive()
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPhaseAssemblage.f90
    !> \brief   Check whether the current minimum driving-force witness is already active.
    !> \author  S.Y. Kwon
    !> \date    Jul. 01, 2026
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/01/2026      S.Y. Kwon            Original helper for direct inactive-witness PEA handoff logic.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details PEA sometimes receives a negative driving-force witness after an accepted embedded Lagrangian
    !! polish.  This helper distinguishes a genuinely inactive phase from another pseudo-compound row that
    !! already maps to a phase in the current Lagrangian active set.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! dPhasePotential                   Current Leveling-row driving-force vector.
    ! iLevel2LagrangeOutputAssemblage   Current Leveling-to-Lagrangian row mapping.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! MinPotentialPhaseIsActive         TRUE when the lowest-potential row maps to an active phase.
    !
    !
    ! Important called routines:
    ! ==========================
    !
    ! PhaseIndexForLevelingRow          Maps Leveling rows/pseudo-compound ids to real phase ids.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! CheckPhaseAssemblage              Uses this helper for PEA convergence and direct-witness handoff checks.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - A row id can be a stoichiometric compound, pseudo-compound candidate, or negative phase marker.
    ! - Phase identity, not candidate identity, determines whether the driving-force witness is already active.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: i, iMinPotentialRow, iWitnessPhase, iActivePhase
    integer :: PhaseIndexForLevelingRow

    MinPotentialPhaseIsActive = .FALSE.
    if (.NOT.allocated(dPhasePotential)) return
    if (.NOT.allocated(iLevel2LagrangeOutputAssemblage)) return

    iMinPotentialRow = MAXVAL(MINLOC(dPhasePotential))
    iWitnessPhase = PhaseIndexForLevelingRow(iMinPotentialRow)
    if (iWitnessPhase <= 0) return

    do i = 1, SIZE(iLevel2LagrangeOutputAssemblage)
        iActivePhase = PhaseIndexForLevelingRow(iLevel2LagrangeOutputAssemblage(i))
        if (iActivePhase == iWitnessPhase) then
            MinPotentialPhaseIsActive = .TRUE.
            return
        end if
    end do

end function MinPotentialPhaseIsActive


!> \brief Map a Leveling row or phase marker to the corresponding real phase id.
integer function PhaseIndexForLevelingRow(iRow)
!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CheckPhaseAssemblage.f90
    !> \brief   Map Leveling candidate rows to real thermodynamic phase ids.
    !> \author  S.Y. Kwon
    !> \date    Jul. 01, 2026
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/01/2026      S.Y. Kwon            Original helper for comparing PEA pseudo-candidates by phase id.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details Leveling rows can refer to real stoichiometric compounds, pseudo-compound solution minima, or
    !! negative phase markers.  PEA cycle detection must compare these rows by their underlying phase identity
    !! so that duplicate pseudo-compounds from the same solution phase are not treated as different active sets.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! iRow          Leveling row id or negative phase marker.
    ! iPhaseLevel   Phase id for each Leveling pseudo-compound row.
    ! iPhase        Fallback phase id vector.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! PhaseIndexForLevelingRow  Real phase id for the row, or zero if no phase can be identified.
    !
    !
    ! Important called routines:
    ! ==========================
    !
    ! None.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! MinPotentialPhaseIsActive          Checks whether a driving-force witness is active.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Negative row ids already encode real phase ids.
    ! - `iPhaseLevel` is preferred because PEA pseudo-compounds can have row ids unrelated to `iPhase`.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo

    implicit none

    integer, intent(in) :: iRow

    PhaseIndexForLevelingRow = 0
    if (iRow == 0) return

    if (iRow < 0) then
        PhaseIndexForLevelingRow = -iRow
        return
    end if

    if (allocated(iPhaseLevel)) then
        if ((iRow >= 1).AND.(iRow <= SIZE(iPhaseLevel))) then
            PhaseIndexForLevelingRow = iPhaseLevel(iRow)
            return
        end if
    end if

    if (allocated(iPhase)) then
        if ((iRow >= 1).AND.(iRow <= SIZE(iPhase))) then
            PhaseIndexForLevelingRow = iPhase(iRow)
            return
        end if
    end if

end function PhaseIndexForLevelingRow
