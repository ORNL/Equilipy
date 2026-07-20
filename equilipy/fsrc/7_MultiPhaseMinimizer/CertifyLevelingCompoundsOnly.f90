!> \brief Certify a provisional compounds-only result from classical Leveling.
!!
!! \details Builds the standard PEA solution-candidate rows with one
!! CompInitMinSolnPoint sweep and tests every resulting phase potential
!! against the same named Leveling tolerance used by the PEA exit certificate.
!! A below-plane solution minimum clears the provisional fast-exit flag and
!! leaves its candidate row available to the normal PEA pipeline.  Certification
!! also requires typed coverage of every physical solution-candidate identity.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CertifyLevelingCompoundsOnly.f90
    !> \brief   Certify a provisional compounds-only classical Leveling result.
    !> \author  S.Y. Kwon
    !> \date    Jul. 16, 2026
    !> \sa      RunLeveling.f90
    !> \sa      InitCheckPhaseAssemblage.f90
    !> \sa      CompInitMinSolnPoint.f90
    !> \sa      CheckPhaseAssemblage.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Certified compounds-only Leveling exits with complete typed solution-phase coverage.
    !
    ! Purpose:
    ! ========
    !
    !> \details Classical Leveling compares compounds with solution endmembers
    !! but cannot prove that a mixed solution constitution lies above the
    !! resulting elemental-potential plane.  This routine performs exactly one
    !! initial all-solution minimum-point sweep before allowing that provisional
    !! compounds-only result to become a terminal minimizer exit.
    !
    ! Required input variables:
    ! =========================
    !
    ! lCompbdOnly                 Provisional compounds-only flag from RunLeveling.
    ! dElementPotential           Elemental-potential plane from classical Leveling.
    ! dLevelingChemicalPotential  Classical Leveling candidate energies.
    ! dLevelingCompositionSpecies Classical Leveling candidate compositions.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    ! lCompbdOnly                 Remains true only when the solution sweep certifies the plane.
    ! iPEAExitStatus              Certified-settled for a valid fast return.
    ! iPEAExitReason              Compounds-only for a valid fast return.
    ! iPEAExitFreshMinPointSweep  Records that the certificate used a fresh sweep.
    ! dPEAExitMinPhasePotential   Minimum refreshed phase potential.
    ! dPEAExitTolerance           Named Leveling tolerance used by the certificate.
    ! dLeveling* candidate rows   Retain any below-plane solution witness for PEA.
    ! lConverged/lPhaseChange     Coherent terminal state for a certified fast return.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! InitCheckPhaseAssemblage      Allocates the standard PEA candidate pool and
    !                               calls CompInitMinSolnPoint once.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! MultiPhaseMinimizer           Guards the provisional compounds-only fast return.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - RunLeveling has already produced a compounds-only classical Leveling
    !   assemblage while at least one solution phase exists in the system.
    ! - InitCheckPhaseAssemblage sets dPEATol before generating the initial
    !   solution minima.  The strict exit threshold is the corresponding
    !   signed Leveling tolerance, -dPEATol.
    ! - Candidate identity and composition storage are owned by the existing
    !   CompInitMinSolnPoint/SetLevelingSolutionCandidateRow machinery.
    ! - Rejected and max-iteration statuses are unresolved.  They may prevent
    !   certification but can never prove that the Leveling plane is settled.
    ! - Typed helper aliases, duplicate rows, and unary ordered parents may be
    !   covered structurally because they are redundant or inapplicable.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine CertifyLevelingCompoundsOnly

    USE ModuleThermo
    USE ModuleThermoIO
    USE ModuleGEMSolver

    implicit none

    real(8) :: dTimerStart, dTimerStop
    logical :: lCandidateCoverageComplete

    call cpu_time(dTimerStart)
    call InitCheckPhaseAssemblage
    dToleranceLevel = -dPEATol
    dPhasePotential = dLevelingChemicalPotential - &
        MATMUL(dLevelingCompositionSpecies, dElementPotential)
    dMinPhasePotential = MINVAL(dPhasePotential)
    lCandidateCoverageComplete = CandidateCoverageIsComplete()
    call cpu_time(dTimerStop)

    dGEMTimingCertification = dGEMTimingCertification + MAX(0D0, dTimerStop - dTimerStart)
    nGEMCertificationSweep = nGEMCertificationSweep + 1
    iPEAExitFreshMinPointSweep = 1
    dPEAExitMinPhasePotential = dMinPhasePotential
    dPEAExitTolerance = dToleranceLevel

    if ((INFOThermo == 0).AND.lCandidateCoverageComplete.AND.&
        (dMinPhasePotential >= dToleranceLevel)) then
        lCompbdOnly = .TRUE.
        iPEAExitStatus = PEA_EXIT_STATUS_CERTIFIED_SETTLED
        iPEAExitReason = PEA_EXIT_REASON_COMPOUNDS_ONLY
        lConverged = .TRUE.
        lPhaseChange = .FALSE.
        iPhaseChangeReason = PHASE_CHANGE_REASON_NONE
        iGEMExitStatus = GEM_EXIT_STATUS_OK
    else
        lCompbdOnly = .FALSE.
        call ResetRejectedCertificateTrace
    end if

    return

contains

    logical function CandidateCoverageIsComplete()

        integer :: iSolnPhase, iCandidateStatus, iTopologyClass
        logical :: OrderDisorderPhaseIsEligible, IsOrderDisorderHelperAliasPhase

        CandidateCoverageIsComplete = .FALSE.
        if (.NOT.allocated(iSubMinCandidateStatusSoln)) return
        if (SIZE(iSubMinCandidateStatusSoln) < nSolnPhasesSys) return

        do iSolnPhase = 1, nSolnPhasesSys
            ! Typed helper aliases are nonphysical representations.  The parent
            ! classification below proves whether their physical branch was resolved.
            if (IsOrderDisorderHelperAliasPhase(iSolnPhase)) cycle
            if (.NOT.OrderDisorderPhaseIsEligible(iSolnPhase)) cycle

            iCandidateStatus = iSubMinCandidateStatusSoln(iSolnPhase)
            select case (iCandidateStatus)
            case (SUBMIN_CANDIDATE_CONVERGED, SUBMIN_CANDIDATE_NEGATIVE_WITNESS)
                continue
            case (SUBMIN_CANDIDATE_DUPLICATE)
                ! The typed duplicate status is structural coverage by an
                ! already evaluated physical candidate identity.
                cycle
            case default
                return
            end select

            if (TRIM(cSolnPhaseType(iSolnPhase)) /= 'SUBOM') cycle
            if (.NOT.allocated(iDisorderedPhase)) cycle
            if (iSolnPhase > SIZE(iDisorderedPhase)) return
            if (iDisorderedPhase(iSolnPhase) <= 0) cycle
            if (.NOT.allocated(iODTopologyClass)) return
            if (iSolnPhase > SIZE(iODTopologyClass)) return

            iTopologyClass = iODTopologyClass(iSolnPhase)
            select case (iTopologyClass)
            case (OD_TOPOLOGY_HELPER_STANDALONE, OD_TOPOLOGY_DIRECT_TARGET, &
                OD_TOPOLOGY_HELPER_ONLY)
                if (.NOT.CanonicalPartitionIdentityIsResolved(iSolnPhase)) return
            case (OD_TOPOLOGY_LEGACY_NO_METADATA)
                ! Legacy DAT phases retain their accepted per-phase status.
                continue
            case default
                ! Not-evaluated, ambiguous, and unsupported typed topology
                ! classes are unresolved identities and cannot confirm an exit.
                return
            end select
        end do

        CandidateCoverageIsComplete = .TRUE.

        return
    end function CandidateCoverageIsComplete


    logical function CanonicalPartitionIdentityIsResolved(iSolnPhase)

        integer, intent(in) :: iSolnPhase
        integer             :: iCandidateClass

        CanonicalPartitionIdentityIsResolved = .FALSE.
        if (.NOT.allocated(iODCandidateClass)) return
        if (iSolnPhase > SIZE(iODCandidateClass)) return

        iCandidateClass = iODCandidateClass(iSolnPhase)
        CanonicalPartitionIdentityIsResolved = &
            (iCandidateClass == OD_CANDIDATE_DISORDERED).OR.&
            (iCandidateClass == OD_CANDIDATE_ORDERED).OR.&
            (iCandidateClass == OD_CANDIDATE_ORDERED_COMPANION_UNSTABLE).OR.&
            (iCandidateClass == OD_CANDIDATE_DISORDERED_PROJECTED)

        return
    end function CanonicalPartitionIdentityIsResolved


    subroutine ResetRejectedCertificateTrace

        ! The normal PEA entry rebuilds this candidate pool, so a rejected
        ! certificate still pays a second initial sweep.  Discard the first
        ! sweep's passive SUBOM trace so that recomputation records each event once.
        nSUBOMTwoSetTrace = 0
        if (allocated(iSUBOMTwoSetTraceStage)) iSUBOMTwoSetTraceStage = 0
        if (allocated(iSUBOMTwoSetTraceIterPEA)) iSUBOMTwoSetTraceIterPEA = 0
        if (allocated(iSUBOMTwoSetTraceIterGlobal)) iSUBOMTwoSetTraceIterGlobal = 0
        if (allocated(iSUBOMTwoSetTracePhase)) iSUBOMTwoSetTracePhase = 0
        if (allocated(iSUBOMTwoSetTraceOrdinal)) iSUBOMTwoSetTraceOrdinal = 0
        if (allocated(iSUBOMTwoSetTraceSlot)) iSUBOMTwoSetTraceSlot = 0
        if (allocated(dSUBOMTwoSetTraceAmount)) dSUBOMTwoSetTraceAmount = 0D0
        if (allocated(dSUBOMTwoSetTraceMol)) dSUBOMTwoSetTraceMol = 0D0
        if (allocated(dSUBOMTwoSetTraceSite)) dSUBOMTwoSetTraceSite = 0D0

        return
    end subroutine ResetRejectedCertificateTrace

end subroutine CertifyLevelingCompoundsOnly
