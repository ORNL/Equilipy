!> \brief Seed active SUBOM phases away from unstable random stationary points.
!!
!! \details A fixed-active-set CEF Lagrangian solve can satisfy first-order
!! stationarity at the random ordered-phase state even when the ordered
!! branch is lower in Gibbs energy.  This routine audits already-active SUBOM
!! phases against the current elemental-potential plane by reusing the
!! endpoint-start phase-local minimization strategy used by PEA.  A lower
!! point is accepted only as a new seed; RunLagrangianGEM then reruns the full
!! mass-conserving Lagrangian solve.
!!
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    StabilizeActiveOrderDisorderPhases.f90
    !> \brief   Seed active SUBOM phases away from unstable random stationary points.
    !> \author  S.Y. Kwon
    !> \date    Jul. 01, 2026
    !> \sa      RunLagrangianGEM.f90
    !> \sa      Subminimization.f90
    !> \sa      CompMinSolnPoint.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   07/20/2026      S.Y. Kwon           Coalesced degenerate order/disorder phases and synchronized accepted constitutions into active slots.
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to prevent a converged CEF
    !! Lagrangian solve from accepting a SUBOM random-state saddle as a true
    !! local minimum.  For each active SUBOM phase, it tries endpoint-like
    !! phase-local starts under the current elemental-potential plane.  If a
    !! start finds a negative driving-force point, the best constitution is
    !! copied back into the active phase and the caller reruns Lagrangian GEM.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! iAssemblage              Current fixed active phase assemblage.
    ! dMolesPhase              Current phase amounts.
    ! dElementPotential        Current Lagrangian elemental-potential plane.
    ! dMolFraction             Current active solution constitutions.
    ! dLevelingChemicalPotential
    !                          Leveling/PEA endmember Gibbs values used only to rank endpoint starts.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] lChanged    True when one active SUBOM phase was reseeded.
    !
    ! dMolFraction             Updated for the reseeded active SUBOM phase.
    ! dMolesSpecies            Updated consistently with the current phase amount.
    ! dActiveSlotMolFraction   Updated with the accepted active-slot product fractions when allocated.
    ! dActiveSlotSiteFraction  Updated with the accepted active-slot site fractions when allocated.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! Subminimization          Finds a lower phase-local point against the current Gibbs plane.
    ! Qsort                    Ranks endpoint starts by endmember grand potential.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! RunLagrangianGEM         Calls after CEF Lagrangian convergence and before accepting the active set.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - The routine does not add or remove phases.  It only changes the
    !   constitution of an already-active SUBOM phase.
    ! - A negative phase-local driving force is interpreted as proof that the
    !   current stationary point is not the minimum for that active phase under
    !   the current elemental-potential plane.
    ! - Endpoint starts use the same low-grand-potential subset strategy as
    !   PEA, including exact ties at the cutoff.  This keeps the audit from
    !   becoming a hidden grid search over ordered endmembers.
    !
    !-------------------------------------------------------------------------------------------------------------



subroutine StabilizeActiveOrderDisorderPhases(lChanged)

    USE ModuleThermo
    USE ModuleSubMin
    USE ModuleGEMSolver

    implicit none

    logical, intent(out) :: lChanged

    real(8), parameter :: dEndpointTieTolerance = 1D-12
    real(8), parameter :: dStartMinMoleFraction = 1D-5

    integer :: iSoln, iSlot, iSolnPhase, iFirst, iLast, nConstituents
    integer :: iStart, nStart, nEndpointStartCount, iEndmember, iLocal
    real(8) :: dStartMaxMoleFraction, dBestDrivingForce, dEndpointCutoffPotential
    logical :: lAddPhase, lAcceptedCandidate

    integer, allocatable :: iEndmemberPotential(:)
    real(8), allocatable :: dEndmemberPotential(:), dBestMolFraction(:)
    real(8), allocatable :: dMolFractionSave(:), dMolesSpeciesSave(:), dMolesPhaseSave(:)
    real(8), allocatable :: dChemicalPotentialSave(:), dPhasePotentialSave(:)
    real(8), allocatable :: dDrivingForceSolnSave(:), dGibbsSolnPhaseSave(:)

    lChanged = .FALSE.
    if (nSolnPhases <= 0) return
    if (.NOT.allocated(dDrivingForceSoln)) return

    allocate(dMolFractionSave(nSpecies), dMolesSpeciesSave(nSpecies), dMolesPhaseSave(nElements))
    allocate(dChemicalPotentialSave(nSpecies), dPhasePotentialSave(nSpecies))
    allocate(dDrivingForceSolnSave(MAX(1,nSolnPhasesSys)), dGibbsSolnPhaseSave(MAX(1,nSolnPhasesSys)))

    dMolFractionSave = dMolFraction
    dMolesSpeciesSave = dMolesSpecies
    dMolesPhaseSave = dMolesPhase
    dChemicalPotentialSave = dChemicalPotential
    dPhasePotentialSave = dPhasePotential
    dDrivingForceSolnSave = dDrivingForceSoln
    dGibbsSolnPhaseSave = dGibbsSolnPhase

    LOOP_ActiveSoln: do iSoln = 1, nSolnPhases
        iSlot = nElements - iSoln + 1
        iSolnPhase = -iAssemblage(iSlot)
        if ((iSolnPhase <= 0).OR.(iSolnPhase > nSolnPhasesSys)) cycle LOOP_ActiveSoln
        if (TRIM(cSolnPhaseType(iSolnPhase)) /= 'SUBOM') cycle LOOP_ActiveSoln
        if (dMolesPhase(iSlot) <= dTolerance(8)) cycle LOOP_ActiveSoln

        iFirst = nSpeciesPhase(iSolnPhase-1) + 1
        iLast = nSpeciesPhase(iSolnPhase)
        nConstituents = iLast - iFirst + 1
        if (nConstituents <= 1) cycle LOOP_ActiveSoln

        if (allocated(iEndmemberPotential)) deallocate(iEndmemberPotential)
        if (allocated(dEndmemberPotential)) deallocate(dEndmemberPotential)
        if (allocated(dBestMolFraction)) deallocate(dBestMolFraction)
        allocate(iEndmemberPotential(nConstituents), dEndmemberPotential(nConstituents))
        allocate(dBestMolFraction(nConstituents))

        call BuildEndpointRanking(iFirst, iLast, nConstituents, dEndmemberPotential, iEndmemberPotential)

        nEndpointStartCount = MIN(nElements, nConstituents)
        dEndpointCutoffPotential = dEndmemberPotential(nEndpointStartCount)
        do while (nEndpointStartCount < nConstituents)
            if (DABS(dEndmemberPotential(nEndpointStartCount+1) - &
                dEndpointCutoffPotential) > dEndpointTieTolerance) exit
            nEndpointStartCount = nEndpointStartCount + 1
        end do
        nStart = nEndpointStartCount

        dBestDrivingForce = 9D5
        dBestMolFraction = dMolFractionSave(iFirst:iLast)

        LOOP_EndpointStarts: do iStart = 1, nStart
            dMolFraction = dMolFractionSave
            dMolesSpecies = dMolesSpeciesSave
            dMolesPhase = dMolesPhaseSave
            dChemicalPotential = dChemicalPotentialSave
            dPhasePotential = dPhasePotentialSave
            dDrivingForceSoln = dDrivingForceSolnSave
            dGibbsSolnPhase = dGibbsSolnPhaseSave

            dStartMaxMoleFraction = 1D0 - dStartMinMoleFraction * DFLOAT(nConstituents-1)
            dMolFraction(iFirst:iLast) = dStartMinMoleFraction
            iEndmember = iFirst + iEndmemberPotential(iStart) - 1
            dMolFraction(iEndmember) = DMAX1(dStartMaxMoleFraction, 0.99D0)
            dMolFraction(iFirst:iLast) = dMolFraction(iFirst:iLast) / SUM(dMolFraction(iFirst:iLast))
            dMolesSpecies(iFirst:iLast) = dMolesPhase(iSlot) * dMolFraction(iFirst:iLast)

            lAddPhase = .FALSE.
            call Subminimization(iSolnPhase, lAddPhase)
            lAcceptedCandidate = (lSubMinConverged.OR.(dDrivingForce < dTolerance(4))).AND.&
                (dDrivingForce < dBestDrivingForce)
            if (lAcceptedCandidate) then
                dBestDrivingForce = dDrivingForce
                dBestMolFraction = dMolFraction(iFirst:iLast)
            end if
        end do LOOP_EndpointStarts

        dMolFraction = dMolFractionSave
        dMolesSpecies = dMolesSpeciesSave
        dMolesPhase = dMolesPhaseSave
        dChemicalPotential = dChemicalPotentialSave
        dPhasePotential = dPhasePotentialSave
        dDrivingForceSoln = dDrivingForceSolnSave
        dGibbsSolnPhase = dGibbsSolnPhaseSave

        if (dBestDrivingForce < dTolerance(4)) then
            dMolFraction(iFirst:iLast) = dBestMolFraction
            do iLocal = iFirst, iLast
                dMolesSpecies(iLocal) = dMolesPhase(iSlot) * dMolFraction(iLocal)
            end do
            call SyncAcceptedSubomRestartToActiveSlot(iSlot, iSolnPhase, iFirst, iLast)
            lChanged = .TRUE.
            exit LOOP_ActiveSoln
        end if
    end do LOOP_ActiveSoln

    if (allocated(iEndmemberPotential)) deallocate(iEndmemberPotential)
    if (allocated(dEndmemberPotential)) deallocate(dEndmemberPotential)
    if (allocated(dBestMolFraction)) deallocate(dBestMolFraction)
    if (allocated(dMolFractionSave)) deallocate(dMolFractionSave)
    if (allocated(dMolesSpeciesSave)) deallocate(dMolesSpeciesSave)
    if (allocated(dMolesPhaseSave)) deallocate(dMolesPhaseSave)
    if (allocated(dChemicalPotentialSave)) deallocate(dChemicalPotentialSave)
    if (allocated(dPhasePotentialSave)) deallocate(dPhasePotentialSave)
    if (allocated(dDrivingForceSolnSave)) deallocate(dDrivingForceSolnSave)
    if (allocated(dGibbsSolnPhaseSave)) deallocate(dGibbsSolnPhaseSave)

    return

contains

    subroutine BuildEndpointRanking(iFirstIn, iLastIn, nConstituentsIn, dPotentialOut, iIndexOut)

        implicit none

        integer, intent(in) :: iFirstIn, iLastIn, nConstituentsIn
        real(8), intent(out) :: dPotentialOut(nConstituentsIn)
        integer, intent(out) :: iIndexOut(nConstituentsIn)

        integer :: iSpecies, iLocal, iElement
        logical :: lUseLevelingBasis

        lUseLevelingBasis = allocated(dLevelingChemicalPotential).AND.&
            allocated(dLevelingCompositionSpecies)
        if (lUseLevelingBasis) then
            lUseLevelingBasis = (SIZE(dLevelingChemicalPotential) >= iLastIn).AND.&
                (SIZE(dLevelingCompositionSpecies,1) >= iLastIn)
        end if

        if (lUseLevelingBasis) then
            dPotentialOut = dLevelingChemicalPotential(iFirstIn:iLastIn) - &
                MATMUL(dLevelingCompositionSpecies(iFirstIn:iLastIn,1:nElements), &
                dElementPotential(1:nElements))
        else
            do iSpecies = iFirstIn, iLastIn
                iLocal = iSpecies - iFirstIn + 1
                dPotentialOut(iLocal) = dChemicalPotential(iSpecies)
                do iElement = 1, nElements
                    dPotentialOut(iLocal) = dPotentialOut(iLocal) - &
                        dElementPotential(iElement) * dStoichSpecies(iSpecies,iElement) / &
                        DFLOAT(iParticlesPerMole(iSpecies))
                end do
            end do
        end if

        call Qsort(dPotentialOut, iIndexOut, nConstituentsIn)

        return

    end subroutine BuildEndpointRanking

    subroutine SyncAcceptedSubomRestartToActiveSlot(iSlotIn, iSolnPhaseIn, iFirstIn, iLastIn)

        implicit none

        integer, intent(in) :: iSlotIn, iSolnPhaseIn, iFirstIn, iLastIn

        integer :: iPhaseIDLocal
        real(8) :: dSumLocal
        real(8) :: dAcceptedSite(nMaxSublatticeSys,nMaxConstituentSys)

        if ((iSlotIn < 1).OR.(iSlotIn > nElements)) return
        if ((iSolnPhaseIn < 1).OR.(iSolnPhaseIn > nSolnPhasesSys)) return
        if ((iFirstIn < 1).OR.(iLastIn > nSpecies).OR.(iFirstIn > iLastIn)) return

        if (allocated(dActiveSlotMolFraction)) then
            if ((iSlotIn <= SIZE(dActiveSlotMolFraction,1)).AND.&
                (iLastIn <= SIZE(dActiveSlotMolFraction,2))) then
                dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) = DMAX1(dMolFraction(iFirstIn:iLastIn), 0D0)
                dSumLocal = SUM(dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn))
                if (dSumLocal > 1D-300) then
                    dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) = &
                        dActiveSlotMolFraction(iSlotIn,iFirstIn:iLastIn) / dSumLocal
                end if
            end if
        end if

        if (allocated(dMolFractionGEM)) then
            if ((iSlotIn <= SIZE(dMolFractionGEM,1)).AND.(iLastIn <= SIZE(dMolFractionGEM,2))) then
                dMolFractionGEM(iSlotIn,iFirstIn:iLastIn) = DMAX1(dMolFraction(iFirstIn:iLastIn), 0D0)
                dSumLocal = SUM(dMolFractionGEM(iSlotIn,iFirstIn:iLastIn))
                if (dSumLocal > 1D-300) then
                    dMolFractionGEM(iSlotIn,iFirstIn:iLastIn) = &
                        dMolFractionGEM(iSlotIn,iFirstIn:iLastIn) / dSumLocal
                end if
            end if
        end if

        if (allocated(iActiveSlotThermoPhase)) then
            if (iSlotIn <= SIZE(iActiveSlotThermoPhase)) iActiveSlotThermoPhase(iSlotIn) = iSolnPhaseIn
        end if
        if (allocated(iActiveSlotDisplayPhase)) then
            if (iSlotIn <= SIZE(iActiveSlotDisplayPhase)) iActiveSlotDisplayPhase(iSlotIn) = iSolnPhaseIn
        end if

        iPhaseIDLocal = iPhaseSublattice(iSolnPhaseIn)
        if (iPhaseIDLocal <= 0) return
        if (.NOT.allocated(dActiveSlotSiteFraction)) return
        if (iSlotIn > SIZE(dActiveSlotSiteFraction,1)) return

        call BuildAcceptedSiteFraction(iSolnPhaseIn, iPhaseIDLocal, dAcceptedSite)
        dActiveSlotSiteFraction(iSlotIn,:,:) = dAcceptedSite

        return

    end subroutine SyncAcceptedSubomRestartToActiveSlot

    subroutine BuildAcceptedSiteFraction(iSolnPhaseIn, iPhaseIDIn, dSiteOut)

        implicit none

        integer, intent(in)  :: iSolnPhaseIn, iPhaseIDIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iFirstLocal, iLastLocal, iSpeciesLocal, iLocal
        integer :: iSub, iCon

        dSiteOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIn-1) + 1
        iLastLocal  = nSpeciesPhase(iSolnPhaseIn)
        do iSpeciesLocal = iFirstLocal, iLastLocal
            iLocal = iSpeciesLocal - iFirstLocal + 1
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                    DMAX1(dMolFraction(iSpeciesLocal), 0D0)
            end do
        end do
        call NormalizeAcceptedSiteFraction(iPhaseIDIn, dSiteOut)

        return

    end subroutine BuildAcceptedSiteFraction

    subroutine NormalizeAcceptedSiteFraction(iPhaseIDIn, dSiteInOut)

        implicit none

        integer, intent(in)    :: iPhaseIDIn
        real(8), intent(inout) :: dSiteInOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iSub, iCon
        real(8) :: dSumLocal

        do iSub = 1, nSublatticePhase(iPhaseIDIn)
            dSumLocal = 0D0
            do iCon = 1, nConstituentSublattice(iPhaseIDIn,iSub)
                dSumLocal = dSumLocal + dSiteInOut(iSub,iCon)
            end do
            if (dSumLocal <= 0D0) cycle
            do iCon = 1, nConstituentSublattice(iPhaseIDIn,iSub)
                dSiteInOut(iSub,iCon) = dSiteInOut(iSub,iCon) / dSumLocal
            end do
        end do

        return

    end subroutine NormalizeAcceptedSiteFraction

end subroutine StabilizeActiveOrderDisorderPhases
