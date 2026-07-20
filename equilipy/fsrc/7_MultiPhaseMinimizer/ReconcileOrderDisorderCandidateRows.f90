!> \brief Coalesce a symmetric ordered minimum with its disordered companion.
!!
!! \details Treats one DIS_PART pair as one candidate identity when the
!! ordered parent's selected minimum is no lower than its exact disordered
!! projection and the projected ordering curvature is positive.  Exact
!! ordering nodes and one-representation-step energy overlaps remain typed
!! while using the structural disordered projection.  For an
!! ordered minimum, an independent standalone disordered phase remains a
!! separate candidate.  Ambiguous or failed classifications are recorded and
!! never silently coalesced.
!!
subroutine ReconcileOrderDisorderCandidateRows
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    ReconcileOrderDisorderCandidateRows.f90
    !> \brief   Reconcile ordered/disordered PEA candidate identity.
    !> \author  S.Y. Kwon
    !> \date    Jul. 16, 2026
    !> \sa      CompOrderingModeSUBOM.f90
    !> \sa      SetLevelingSolutionCandidateRow.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Reconciled order/disorder candidate ownership from structural stability and propagated certified owners into selected PEA slots.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The ordered SUBOM model contains the disordered manifold.  Its
    !! pseudo-compound row is therefore redundant when the selected minimum is
    !! the symmetric state already represented by the canonical disordered
    !! companion.  The classification uses same-parent energy and projected
    !! ordering curvature, not a site-fraction distance tolerance.
    !
    ! Required input variables:
    ! =========================
    !
    ! iDisorderedPhase             Thermodynamic DIS_PART helper mapping.
    ! iODCompanionPhase            Parser-declared DIS_PART helper mapping used by the parent.
    ! iODStandalonePhase           Separately declared standalone disordered phase, if present.
    ! iSubMinCandidateStatusSoln    Raw selected subminimization status per solution phase.
    ! dMolFraction                 Selected phase-local constitutions.
    !
    ! Output/updated variables:
    ! =========================
    !
    ! iODCandidateClass            Typed identity verdict per ordered SUBOM phase.
    ! iODCandidateCompanionPhase   Canonical companion used for the verdict.
    ! dLevelingChemicalPotential   Ordered pseudo-row disabled only for a proven symmetric minimum.
    ! dODCompanionEigenMin
    !                              Ordering curvature at the companion's own composition.
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! CompOrderingModeSUBOM        Evaluates current/disordered energy and ordering curvature.
    ! ProjectOrderDisorderCompanionFraction
    !                              Maps the companion state into the ordered parent.
    ! SetLevelingSolutionCandidateRow
    !                              Invalidates a redundant ordered pseudo-row.
    !
    ! Primary callers:
    ! ================
    !
    ! CompInitMinSolnPoint         Reconciles the initial PEA candidate sweep.
    ! CompMinSolnPoint             Reconciles every refreshed PEA candidate sweep.
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Equality is evaluated through same-parent Gibbs arithmetic.  No
    !   ordering-degree or composition-distance tolerance is introduced.
    ! - A negative projected ordering curvature without a lower selected
    !   ordered minimum is an ambiguous failed branch search, not disordered
    !   evidence.
    !
    !-------------------------------------------------------------------------------------------------------------
    USE ModuleThermo
    USE ModuleGEMSolver
    USE GridDiscovery, ONLY: GridTangentRowIndex

    implicit none

    integer :: iOrderedPhase, iCompanionPhase, iRawHelperPhase, iStandalonePhase
    integer :: nModeCapacity, nModeOut, iInfo
    integer :: iCompanionStatus, iCompanionInfo, iProjectionStatus
    real(8) :: dOrderingEigenMin, dCompanionOrderingEigenMin
    real(8), dimension(:), allocatable :: dOrderingGradient
    real(8), dimension(:,:), allocatable :: dOrderingHessian
    logical :: lCommittedProjection, lRecoveredAmbiguousProjection

    if (.NOT.lODPartitionUnifiedActive) return
    if (.NOT.allocated(iDisorderedPhase)) return
    if (.NOT.allocated(iODCompanionPhase)) return
    if (.NOT.allocated(iODTopologyClass)) return
    if (.NOT.allocated(iODCandidateClass)) return
    if (.NOT.allocated(iSubMinCandidateStatusSoln)) return

    nModeCapacity = MAX(1, nMaxSublatticeSys*nMaxConstituentSys)
    allocate(dOrderingGradient(nModeCapacity))
    allocate(dOrderingHessian(nModeCapacity,nModeCapacity))

    do iOrderedPhase = 1, MIN(nSolnPhasesSys, SIZE(iDisorderedPhase))
        lCommittedProjection = .FALSE.
        lRecoveredAmbiguousProjection = .FALSE.
        if (TRIM(cSolnPhaseType(iOrderedPhase)) /= 'SUBOM') cycle
        if ((iODTopologyClass(iOrderedPhase) < OD_TOPOLOGY_HELPER_STANDALONE).OR.&
            (iODTopologyClass(iOrderedPhase) > OD_TOPOLOGY_HELPER_ONLY)) cycle

        iRawHelperPhase = iDisorderedPhase(iOrderedPhase)
        iCompanionPhase = iODCompanionPhase(iOrderedPhase)
        iStandalonePhase = StandaloneDisorderedPhase(iOrderedPhase)
        if ((iRawHelperPhase <= 0).OR.(iCompanionPhase <= 0)) cycle
        if (iCompanionPhase > nSolnPhasesSys) cycle

        iODCandidateCompanionPhase(iOrderedPhase) = iCompanionPhase
        if ((iODTopologyClass(iOrderedPhase) == OD_TOPOLOGY_HELPER_STANDALONE).OR.&
            (iODTopologyClass(iOrderedPhase) == OD_TOPOLOGY_HELPER_ONLY)) then
            call SuppressCandidatePseudoRow(iRawHelperPhase)
        end if

        if (.NOT.CandidateStatusIsUsable(iSubMinCandidateStatusSoln(iOrderedPhase))) then
            iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_EVALUATION_FAILED
            call SuppressDisPartHelperRow(iOrderedPhase, iCompanionPhase)
            cycle
        end if

        lODCandidateClassifyActive = .TRUE.
        call CompOrderingModeSUBOM(iOrderedPhase, nModeCapacity, dOrderingHessian, &
            dOrderingGradient, dOrderingEigenMin, nModeOut, iInfo)
        lODCandidateClassifyActive = .FALSE.
        if ((iInfo /= 0).OR.(nModeOut <= 0)) then
            iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_EVALUATION_FAILED
            call SuppressDisPartHelperRow(iOrderedPhase, iCompanionPhase)
            cycle
        end if

        if (iODCandidateClass(iOrderedPhase) == OD_CANDIDATE_ORDERED) then
            call EvaluateCompanionOrderingCurvature(iOrderedPhase, iCompanionPhase, &
                dCompanionOrderingEigenMin, iCompanionInfo)
            if (iCompanionInfo /= 0) then
                iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_AMBIGUOUS_COMPANION
                call SuppressDisPartHelperRow(iOrderedPhase, iCompanionPhase)
                cycle
            end if
            dODCompanionEigenMin(iOrderedPhase) = &
                dCompanionOrderingEigenMin
            if (dCompanionOrderingEigenMin < 0D0) then
                iODCandidateClass(iOrderedPhase) = &
                    OD_CANDIDATE_ORDERED_COMPANION_UNSTABLE
            else if (dCompanionOrderingEigenMin == 0D0) then
                iODCandidateClass(iOrderedPhase) = &
                    OD_CANDIDATE_AMBIGUOUS_COMPANION
            end if
            call MarkCandidateIdentity(GridTangentRowIndex(iOrderedPhase), iOrderedPhase, &
                iOrderedPhase, 1)
            cycle
        end if

        if (.NOT.CandidateUsesDisorderedProjection(&
            iODCandidateClass(iOrderedPhase))) then
            call MarkCandidateIdentity(GridTangentRowIndex(iOrderedPhase), iOrderedPhase, &
                iOrderedPhase, 1)
            cycle
        end if

        if (iStandalonePhase > 0) then
            iCompanionStatus = iSubMinCandidateStatusSoln(iStandalonePhase)
            if (.NOT.CandidateStatusIsUsable(iCompanionStatus)) then
                iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_AMBIGUOUS_COMPANION
                cycle
            end if

            if (iODCandidateClass(iOrderedPhase) == OD_CANDIDATE_DISORDERED_PROJECTED) then
                lCommittedProjection = CandidateConstitutionIsOnDisorderedManifold(&
                    iOrderedPhase,iStandalonePhase)
                call CommitProjectedParentConstitution(iOrderedPhase, iStandalonePhase, &
                    iCompanionInfo)
                if ((iCompanionInfo /= 0).OR.&
                    (iODCandidateClass(iOrderedPhase) /= OD_CANDIDATE_DISORDERED)) then
                    if (lCommittedProjection.AND.(iCompanionInfo == 0)) then
                        iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_DISORDERED
                        lRecoveredAmbiguousProjection = .TRUE.
                    else
                        iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_AMBIGUOUS_COMPANION
                        cycle
                    end if
                end if
            end if

            ! The parent-disordered state is redundant when a separately
            ! declared standalone phase represents the same physical branch.
            ! Keep the standalone candidate unchanged; its own Gibbs model
            ! decides the tie, including any database-author offset.
            call SuppressCandidatePseudoRow(iOrderedPhase)
            call MarkCandidateIdentity(GridTangentRowIndex(iStandalonePhase), &
                iStandalonePhase, iStandalonePhase, 1)
            if (lGridFrontEndActive.AND.lRecoveredAmbiguousProjection) then
                lODCommittedOwnershipApplied = .TRUE.
                call ReplaceSelectedOrderedCandidatesWithOwner(iOrderedPhase,iStandalonePhase)
            end if
            cycle
        end if

        iCompanionStatus = iSubMinCandidateStatusSoln(iCompanionPhase)
        if (.NOT.CandidateStatusIsUsable(iCompanionStatus)) then
            iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_AMBIGUOUS_COMPANION
            call SuppressCandidatePseudoRow(iCompanionPhase)
            cycle
        end if

        if (iODCandidateClass(iOrderedPhase) == OD_CANDIDATE_DISORDERED_PROJECTED) then
            call CommitProjectedParentConstitution(iOrderedPhase, iCompanionPhase, &
                iCompanionInfo)
            if ((iCompanionInfo /= 0).OR.&
                (iODCandidateClass(iOrderedPhase) /= OD_CANDIDATE_DISORDERED)) then
                iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_AMBIGUOUS_COMPANION
                cycle
            end if
        end if

        call WriteProjectedParentCandidate(iOrderedPhase, iCompanionPhase, &
            GridTangentRowIndex(iOrderedPhase), iCompanionPhase, 1, iCompanionInfo)
        if (iCompanionInfo /= 0) then
            iODCandidateClass(iOrderedPhase) = OD_CANDIDATE_AMBIGUOUS_COMPANION
        end if
    end do

    deallocate(dOrderingGradient, dOrderingHessian)

    return

contains

    logical function CandidateConstitutionIsOnDisorderedManifold(iOrderedPhaseIn, iOwnerPhaseIn)

        integer, intent(in) :: iOrderedPhaseIn, iOwnerPhaseIn
        integer :: iFirst, iLast, nLocal, iProjectionStatusLocal, iPhaseID
        real(8), dimension(:), allocatable :: dProjectedFraction, dCurrentFraction
        real(8), dimension(nMaxSublatticeSys,nMaxConstituentSys) :: dCurrentSite, dProjectedSite
        logical :: lProjected

        CandidateConstitutionIsOnDisorderedManifold = .FALSE.
        iFirst = nSpeciesPhase(iOrderedPhaseIn-1)+1
        iLast = nSpeciesPhase(iOrderedPhaseIn)
        nLocal = iLast-iFirst+1
        if (nLocal <= 0) return
        allocate(dProjectedFraction(nLocal),dCurrentFraction(nLocal))
        call ProjectOrderDisorderCompanionFraction(iOwnerPhaseIn,iOrderedPhaseIn,nLocal, &
            dProjectedFraction,lProjected,iProjectionStatusLocal)
        if (.NOT.lProjected) then
            deallocate(dProjectedFraction,dCurrentFraction)
            return
        end if
        iPhaseID = iPhaseSublattice(iOrderedPhaseIn)
        if (.NOT.allocated(dSiteFraction)) then
            deallocate(dProjectedFraction,dCurrentFraction)
            return
        end if
        if ((iPhaseID < 1).OR.(iPhaseID > SIZE(dSiteFraction,1))) then
            deallocate(dProjectedFraction,dCurrentFraction)
            return
        end if
        dCurrentFraction = dMolFraction(iFirst:iLast)
        call CompExcessGibbsEnergy(iOrderedPhaseIn)
        dCurrentSite = dSiteFraction(iPhaseID,:,:)
        dMolFraction(iFirst:iLast) = dProjectedFraction
        call CompExcessGibbsEnergy(iOrderedPhaseIn)
        dProjectedSite = dSiteFraction(iPhaseID,:,:)
        dMolFraction(iFirst:iLast) = dCurrentFraction
        call CompExcessGibbsEnergy(iOrderedPhaseIn)
        CandidateConstitutionIsOnDisorderedManifold = &
            MAXVAL(DABS(dCurrentSite-dProjectedSite)) <= dTolerance(1)
        deallocate(dProjectedFraction,dCurrentFraction)

        return
    end function CandidateConstitutionIsOnDisorderedManifold

    subroutine ReplaceSelectedOrderedCandidatesWithOwner(iOrderedPhaseIn, iOwnerPhaseIn)

        integer, intent(in) :: iOrderedPhaseIn, iOwnerPhaseIn
        integer :: iSlot, iOwnerRow, iFirst, iLast

        iOwnerRow = GridTangentRowIndex(iOwnerPhaseIn)
        if ((iOwnerRow < 1).OR.(iOwnerRow > nSpeciesLevel)) return
        iFirst = nSpeciesPhase(iOwnerPhaseIn-1)+1
        iLast = nSpeciesPhase(iOwnerPhaseIn)
        do iSlot = 1, nElements
            if (iPhaseGEM(iSlot) /= iOrderedPhaseIn) cycle
            iAssemblage(iSlot) = iOwnerRow
            iPhaseGEM(iSlot) = iOwnerPhaseIn
            dMolFractionGEM(iSlot,:) = 0D0
            dMolFractionGEM(iSlot,iFirst:iLast) = dMolFraction(iFirst:iLast)
            dStoichSpeciesGEM(iSlot,:) = dStoichSpeciesLevel(iOwnerRow,:)
            dAtomFractionSpeciesGEM(iSlot,:) = dLevelingCompositionSpecies(iOwnerRow,:)
            dChemicalPotentialGEM(iSlot) = dLevelingChemicalPotential(iOwnerRow)
        end do

        return
    end subroutine ReplaceSelectedOrderedCandidatesWithOwner

    logical function CandidateStatusIsUsable(iStatus)

        integer, intent(in) :: iStatus

        CandidateStatusIsUsable = (iStatus == SUBMIN_CANDIDATE_CONVERGED).OR.&
            (iStatus == SUBMIN_CANDIDATE_NEGATIVE_WITNESS)

        return
    end function CandidateStatusIsUsable


    subroutine SuppressDisPartHelperRow(iOrderedPhaseIn, iCompanionPhaseIn)

        integer, intent(in) :: iOrderedPhaseIn, iCompanionPhaseIn

        if ((iODTopologyClass(iOrderedPhaseIn) /= OD_TOPOLOGY_HELPER_STANDALONE).AND.&
            (iODTopologyClass(iOrderedPhaseIn) /= OD_TOPOLOGY_HELPER_ONLY)) return
        call SuppressCandidatePseudoRow(iCompanionPhaseIn)

        return
    end subroutine SuppressDisPartHelperRow


    logical function CandidateUsesDisorderedProjection(iClass)

        integer, intent(in) :: iClass

        CandidateUsesDisorderedProjection = &
            (iClass == OD_CANDIDATE_DISORDERED).OR.&
            (iClass == OD_CANDIDATE_DISORDERED_PROJECTED).OR.&
            (iClass == OD_CANDIDATE_AMBIGUOUS_ROUNDOFF).OR.&
            (iClass == OD_CANDIDATE_AMBIGUOUS_NODE)

        return
    end function CandidateUsesDisorderedProjection


    subroutine EvaluateCompanionOrderingCurvature(iOrderedPhaseIn, iCompanionPhaseIn, &
        dEigenMinOut, iInfoOut)

        integer, intent(in) :: iOrderedPhaseIn, iCompanionPhaseIn
        real(8), intent(out) :: dEigenMinOut
        integer, intent(out) :: iInfoOut

        integer :: iFirst, iLast, nLocal, nModeLocal, iModeInfo
        integer :: iClassSave
        real(8) :: dCurrentSave, dDisorderedSave, dEigenSave, dCompanionEigenSave
        real(8) :: dGibbsSolnPhaseSave
        real(8), dimension(:), allocatable :: dProjectedFraction
        real(8), dimension(:), allocatable :: dMolFractionSave, dMolesSpeciesSave
        logical :: lProjected

        dEigenMinOut = 0D0
        iInfoOut = 1
        iFirst = nSpeciesPhase(iOrderedPhaseIn-1) + 1
        iLast = nSpeciesPhase(iOrderedPhaseIn)
        nLocal = iLast - iFirst + 1
        if (nLocal <= 0) return

        allocate(dProjectedFraction(nLocal), dMolFractionSave(nLocal), &
            dMolesSpeciesSave(nLocal))
        call ProjectOrderDisorderCompanionFraction(iCompanionPhaseIn, iOrderedPhaseIn, &
            nLocal, dProjectedFraction, lProjected, iProjectionStatus)
        if (.NOT.lProjected) then
            deallocate(dProjectedFraction, dMolFractionSave, dMolesSpeciesSave)
            return
        end if

        dMolFractionSave = dMolFraction(iFirst:iLast)
        dMolesSpeciesSave = dMolesSpecies(iFirst:iLast)
        dGibbsSolnPhaseSave = dGibbsSolnPhase(iOrderedPhaseIn)
        iClassSave = iODCandidateClass(iOrderedPhaseIn)
        dCurrentSave = dODCandidateCurrentGibbs(iOrderedPhaseIn)
        dDisorderedSave = dODCandidateDisorderedGibbs(iOrderedPhaseIn)
        dEigenSave = dODCandidateOrderingEigenMin(iOrderedPhaseIn)
        dCompanionEigenSave = dODCompanionEigenMin(iOrderedPhaseIn)

        dMolFraction(iFirst:iLast) = dProjectedFraction
        dMolesSpecies(iFirst:iLast) = dProjectedFraction
        dGibbsSolnPhase(iOrderedPhaseIn) = 0D0
        lODCandidateClassifyActive = .FALSE.
        call CompOrderingModeSUBOM(iOrderedPhaseIn, nModeCapacity, dOrderingHessian, &
            dOrderingGradient, dEigenMinOut, nModeLocal, iModeInfo)

        dMolFraction(iFirst:iLast) = dMolFractionSave
        dMolesSpecies(iFirst:iLast) = dMolesSpeciesSave
        dGibbsSolnPhase(iOrderedPhaseIn) = dGibbsSolnPhaseSave
        iODCandidateClass(iOrderedPhaseIn) = iClassSave
        dODCandidateCurrentGibbs(iOrderedPhaseIn) = dCurrentSave
        dODCandidateDisorderedGibbs(iOrderedPhaseIn) = dDisorderedSave
        dODCandidateOrderingEigenMin(iOrderedPhaseIn) = dEigenSave
        dODCompanionEigenMin(iOrderedPhaseIn) = dCompanionEigenSave

        if ((iModeInfo == 0).AND.(nModeLocal > 0)) iInfoOut = 0

        return
    end subroutine EvaluateCompanionOrderingCurvature


    subroutine CommitProjectedParentConstitution(iOrderedPhaseIn, iCompanionPhaseIn, &
        iInfoOut)

        integer, intent(in) :: iOrderedPhaseIn, iCompanionPhaseIn
        integer, intent(out) :: iInfoOut

        integer :: iFirst, iLast, nLocal, nModeLocal, iModeInfo
        real(8), dimension(:), allocatable :: dProjectedFraction
        logical :: lProjected

        iInfoOut = 1
        iFirst = nSpeciesPhase(iOrderedPhaseIn-1) + 1
        iLast = nSpeciesPhase(iOrderedPhaseIn)
        nLocal = iLast - iFirst + 1
        if (nLocal <= 0) return

        allocate(dProjectedFraction(nLocal))
        call ProjectOrderDisorderCompanionFraction(iCompanionPhaseIn, iOrderedPhaseIn, &
            nLocal, dProjectedFraction, lProjected, iProjectionStatus)
        if (.NOT.lProjected) then
            deallocate(dProjectedFraction)
            return
        end if

        ! Class 7 means the carried asymmetric state is not the candidate
        ! minimum: its exact stable disordered projection is lower.  Commit
        ! that structural state before classifying again so the final identity
        ! is a property of the selected constitution, not of the start that
        ! happened to reach the rejected branch.
        dMolFraction(iFirst:iLast) = dProjectedFraction
        dMolesSpecies(iFirst:iLast) = dProjectedFraction
        dGibbsSolnPhase(iOrderedPhaseIn) = 0D0
        call CompExcessGibbsEnergy(iOrderedPhaseIn)

        lODCandidateClassifyActive = .TRUE.
        call CompOrderingModeSUBOM(iOrderedPhaseIn, nModeCapacity, dOrderingHessian, &
            dOrderingGradient, dOrderingEigenMin, nModeLocal, iModeInfo)
        lODCandidateClassifyActive = .FALSE.
        if ((iModeInfo == 0).AND.(nModeLocal > 0)) iInfoOut = 0

        deallocate(dProjectedFraction)

        return
    end subroutine CommitProjectedParentConstitution


    subroutine WriteProjectedParentCandidate(iOrderedPhaseIn, iCompanionPhaseIn, &
        iLevelRowIn, iDisplayPhaseIn, iOrdinalIn, iInfoOut)

        integer, intent(in) :: iOrderedPhaseIn, iCompanionPhaseIn
        integer, intent(in) :: iLevelRowIn, iDisplayPhaseIn, iOrdinalIn
        integer, intent(out) :: iInfoOut

        integer :: iFirst, iLast, nLocal
        real(8), dimension(:), allocatable :: dProjectedFraction, dFractionSave
        logical :: lProjected

        iInfoOut = 1
        iFirst = nSpeciesPhase(iOrderedPhaseIn-1) + 1
        iLast = nSpeciesPhase(iOrderedPhaseIn)
        nLocal = iLast - iFirst + 1
        if (nLocal <= 0) return

        allocate(dProjectedFraction(nLocal), dFractionSave(nLocal))
        call ProjectOrderDisorderCompanionFraction(iCompanionPhaseIn, iOrderedPhaseIn, &
            nLocal, dProjectedFraction, lProjected, iProjectionStatus)
        if (.NOT.lProjected) then
            deallocate(dProjectedFraction, dFractionSave)
            return
        end if

        dFractionSave = dMolFraction(iFirst:iLast)
        dMolFraction(iFirst:iLast) = dProjectedFraction
        dGibbsSolnPhase(iOrderedPhaseIn) = 0D0
        call CompExcessGibbsEnergy(iOrderedPhaseIn)
        call SetLevelingSolutionCandidateRow(iLevelRowIn, iOrderedPhaseIn, &
            dProjectedFraction, .TRUE.)
        call MarkCandidateIdentity(iLevelRowIn, iOrderedPhaseIn, &
            iDisplayPhaseIn, iOrdinalIn)

        dMolFraction(iFirst:iLast) = dFractionSave
        dGibbsSolnPhase(iOrderedPhaseIn) = 0D0
        call CompExcessGibbsEnergy(iOrderedPhaseIn)
        iInfoOut = 0
        deallocate(dProjectedFraction, dFractionSave)

        return
    end subroutine WriteProjectedParentCandidate


    subroutine MarkCandidateIdentity(iLevelRowIn, iParentPhaseIn, iDisplayPhaseIn, iOrdinalIn)

        integer, intent(in) :: iLevelRowIn, iParentPhaseIn, iDisplayPhaseIn, iOrdinalIn
        integer :: iCandidate

        if (.NOT.allocated(iLevelCandidateFromLevel)) return
        if ((iLevelRowIn < 1).OR.(iLevelRowIn > SIZE(iLevelCandidateFromLevel))) return
        iCandidate = iLevelCandidateFromLevel(iLevelRowIn)
        if ((iCandidate < 1).OR.(iCandidate > nLevelCandidate)) return

        iLevelCandidateParentPhase(iCandidate) = iParentPhaseIn
        iLevelCandidateDisplayPhase(iCandidate) = iDisplayPhaseIn
        iLevelCandidateIdentityOrdinal(iCandidate) = iOrdinalIn

        return
    end subroutine MarkCandidateIdentity


    subroutine SuppressCandidatePseudoRow(iPhase)

        integer, intent(in) :: iPhase
        integer :: iFirst, iLast

        if ((iPhase < 1).OR.(iPhase > nSolnPhasesSys)) return
        iFirst = nSpeciesPhase(iPhase-1) + 1
        iLast = nSpeciesPhase(iPhase)
        call SetLevelingSolutionCandidateRow(GridTangentRowIndex(iPhase), iPhase, &
            dMolFraction(iFirst:iLast), .FALSE.)
        dPhasePotential(GridTangentRowIndex(iPhase)) = 5D9

        return
    end subroutine SuppressCandidatePseudoRow


    integer function StandaloneDisorderedPhase(iOrderedPhaseIn)

        integer, intent(in) :: iOrderedPhaseIn

        StandaloneDisorderedPhase = 0
        if (.NOT.allocated(iODStandalonePhase)) return
        if ((iOrderedPhaseIn < 1).OR.(iOrderedPhaseIn > SIZE(iODStandalonePhase))) return
        if (iODStandalonePhase(iOrderedPhaseIn) == iOrderedPhaseIn) return
        StandaloneDisorderedPhase = iODStandalonePhase(iOrderedPhaseIn)

        return
    end function StandaloneDisorderedPhase

end subroutine ReconcileOrderDisorderCandidateRows
