!-------------------------------------------------------------------------------------------------------------
!
!> \file    GEMLineSearchCEF.f90
!> \brief   Apply a CEF site-fraction Lagrangian Newton direction.
!> \author  S.Y. Kwon
!> \date    Jun. 25, 2026
!> \sa      GEMNewtonCEF.f90
!> \sa      CompFunctionNorm.f90
!
! Revisions:
! ==========
!
!   Date            Programmer          Description of change
!   ----            ----------          ---------------------
!   06/25/2026      S.Y. Kwon           Added positivity-preserving CEF site-fraction GEM line search.
!   06/25/2026      S.Y. Kwon           Rebuilt active CEF species moles before saving line-search state.
!   06/25/2026      S.Y. Kwon           Limited the initial CEF line-search step by analytical phase and
!                                        site-fraction positivity bounds.
!   06/28/2026      S.Y. Kwon           Extended the guarded line search to mixed CEF/non-CEF active sets.
!   07/01/2026      S.Y. Kwon           Clarified that CEF trial directions are applied in site-fraction
!                                        variables before product endmember fractions are rebuilt.
!   07/02/2026      S.Y. Kwon           Added passive Gibbs/merit line-search diagnostics and
!                                        residual-only no-descent classification.
!   07/02/2026      S.Y. Kwon           Synchronized slot-local CEF site fractions during each feasible
!                                        trial before residual evaluation.
!   07/03/2026      S.Y. Kwon           Preserved bound-active phase-amount retry directions instead of
!                                        overwriting them with mass-balance projection.
!   07/03/2026      S.Y. Kwon           Preserved the lowest Gibbs+mass merit trial for pre-residual-LM
!                                        class-1 diagnostics without accepting it as a state.
!   07/04/2026      S.Y. Kwon           Removed rejected C3-c2 funnel acceptance while retaining
!                                        passive merit diagnostics.
!
!
! Purpose:
! ========
!
!> \details This routine applies the mixed-coordinate Newton direction from
!! `GEMNewtonCEF`.  Unlike the legacy pseudo-endmember lambda correction, this
!! line search updates CEF site fractions directly and supported non-CEF
!! solution phases by independent phase-local fractions.  It rejects a trial
!! step if a phase amount, site fraction, or non-CEF local fraction would
!! become non-positive.  Product endmember fractions are reconstructed only
!! after a feasible CEF site-fraction trial is accepted for evaluation.
!
!
! Required input variables:
! =========================
!
! dGEMCEFPhaseDirection      Phase amount direction from GEMNewtonCEF.
! dGEMCEFSiteDirection       Independent site-fraction direction from GEMNewtonCEF.
! dGEMCEFElementDirection    Element-potential direction from GEMNewtonCEF.
! dGEMCEFSiteLast            Site fractions at the start of the line search.
!
!
! Output/updated variables:
! =========================
!
! dMolesPhase                Accepted phase amounts.
! dMolesSpecies              Accepted species amounts.
! dMolFraction               Accepted solution endmember fractions.
! dElementPotential          Accepted elemental potentials.
! dGEMFunctionNorm           Residual norm at the accepted state.
!
!
! Called subroutines/functions:
! =============================
!
! CompChemicalPotential      Recomputes phase/species chemical potentials.
! CompFunctionNorm           Recomputes the mixed CEF residual norm.
!
!
! Primary callers:
! ================
!
! RunLagrangianGEM          Calls this after GEMNewtonCEF.
!
!
! Numerical assumptions:
! ======================
!
! - Site normalization is enforced by using one reference constituent per
!   sublattice; the reference update is the negative sum of non-reference
!   updates on that sublattice.
! - dGEMCEFSiteDirection is the CEF search coordinate.  dMolFraction is rebuilt
!   from accepted site fractions only for shared thermodynamic model calls and
!   later PEA/Level2Lagrange handoff.
! - Failure to find descent leaves the previous state intact and lets
!   RunLagrangianGEM classify stagnation or nonconvergence.
!
!-------------------------------------------------------------------------------------------------------------



subroutine GEMLineSearchCEF

    USE ModuleThermo
    USE ModuleGEMSolver

    implicit none

    integer :: iterWolfe
    real(8) :: dStepLength, dBestFunctionNorm, dBestStepLength
    real(8) :: dLineSearchNormTolerance
    real(8) :: dTrialGibbs, dTrialMerit, dBestGibbs, dBestMerit
    real(8) :: dMeritFunctionNorm, dMeritStepLength, dMeritGibbs, dMeritMerit
    real(8) :: dInitialMassNorm, dMeritMassNorm
    logical :: lSlopeRecorded
    logical :: lCompEverything, lFeasible
    real(8), allocatable :: dBestMolesSpecies(:), dBestMolesPhase(:), dBestElementPotential(:)
    real(8), allocatable :: dBestPhaseSite(:,:,:), dBasePhaseSite(:,:,:)

    if (.NOT.lGEMCEFSiteDirectionActive) return

    lCompEverything = .FALSE.
    call RebuildActiveCEFMolesFromStoredSites
    call CompChemicalPotential(lCompEverything)
    call CompFunctionNorm
    dGEMFunctionNormLast = dGEMFunctionNorm
    dInitialMassNorm = dGEMMassBalanceNorm
    dGEMLineSearchInitialGibbs = CurrentActiveGibbs()
    dGEMLineSearchInitialMerit = ActiveMassMerit(dGEMLineSearchInitialGibbs)
    dLineSearchNormTolerance = 1D-5
    if (lPostProcess) dLineSearchNormTolerance = 1D-10

    if (allocated(dMolesSpeciesLast)) deallocate(dMolesSpeciesLast)
    if (allocated(dMolesPhaseLast)) deallocate(dMolesPhaseLast)
    if (allocated(dElementPotentialLast)) deallocate(dElementPotentialLast)
    allocate(dMolesSpeciesLast(nSpecies), dMolesPhaseLast(nElements), dElementPotentialLast(nElements))
    dMolesSpeciesLast = dMolesSpecies
    dMolesPhaseLast = dMolesPhase
    dElementPotentialLast = dElementPotential

    allocate(dBestMolesSpecies(nSpecies), dBestMolesPhase(nElements), dBestElementPotential(nElements))
    allocate(dBestPhaseSite(MAX(1,nGEMCEFPhaseVariables),nMaxSublatticeSys,nMaxConstituentSys), &
        dBasePhaseSite(MAX(1,nGEMCEFPhaseVariables),nMaxSublatticeSys,nMaxConstituentSys))
    dBestMolesSpecies = dMolesSpeciesLast
    dBestMolesPhase = dMolesPhaseLast
    dBestElementPotential = dElementPotentialLast
    dBestPhaseSite = dGEMCEFPhaseSiteLast
    dBasePhaseSite = dGEMCEFPhaseSiteLast
    dBestFunctionNorm = dGEMFunctionNormLast
    dBestGibbs = dGEMLineSearchInitialGibbs
    dBestMerit = dGEMLineSearchInitialMerit
    dBestStepLength = 0D0
    dMeritFunctionNorm = dGEMFunctionNormLast
    dMeritStepLength = 0D0
    dMeritGibbs = dGEMLineSearchInitialGibbs
    dMeritMerit = dGEMLineSearchInitialMerit
    dMeritMassNorm = dInitialMassNorm
    lSlopeRecorded = .FALSE.

    iGEMLineSearchIterationCount = 0
    iGEMLineSearchNegativeFactorCount = 0
    iGEMLineSearchFloorCount = 0
    iGEMLineSearchNoDescent = 0
    iGEMLineSearchNoDescentClass = 0
    dGEMLineSearchInitialNorm = dGEMFunctionNormLast
    dGEMLineSearchBestNorm = dGEMFunctionNormLast
    dGEMLineSearchFinalNorm = dGEMFunctionNormLast
    dGEMLineSearchInitialStep = 1D0
    dGEMLineSearchBestStep = 0D0
    dGEMLineSearchFinalStep = 0D0
    dGEMLineSearchMinRawPhaseMoles = MINVAL(dMolesPhaseLast)
    dGEMLineSearchMinFinalPhaseMoles = MINVAL(dMolesPhaseLast)
    dGEMLineSearchBestGibbs = dGEMLineSearchInitialGibbs
    dGEMLineSearchFinalGibbs = dGEMLineSearchInitialGibbs
    dGEMLineSearchBestMerit = dGEMLineSearchInitialMerit
    dGEMLineSearchFinalMerit = dGEMLineSearchInitialMerit
    dGEMLineSearchMeritCandNorm = dGEMFunctionNormLast
    dGEMLineSearchMeritCandMass = dInitialMassNorm
    dGEMLineSearchMeritCandStep = 0D0
    dGEMLineSearchMeritCandGibbs = dGEMLineSearchInitialGibbs
    dGEMLineSearchMeritCandMerit = dGEMLineSearchInitialMerit
    dGEMNewtonDirNormSlope = 0D0
    dGEMNewtonDirGibbsSlope = 0D0
    dGEMNewtonDirMeritSlope = 0D0

    if (allocated(dGEMLSRawPhaseMoles)) then
        dGEMLSRawPhaseMoles = dMolesPhaseLast
        call StoreRawCEFPhaseMoles
        dGEMLineSearchMinRawPhaseMoles = MINVAL(dGEMLSRawPhaseMoles)
        iGEMLineSearchNegativePhaseCount = COUNT(dGEMLSRawPhaseMoles < 0D0)
    end if

    dStepLength = 1D0
    call LimitCEFInitialStep(dStepLength)
    dGEMLineSearchInitialStep = dStepLength
    LOOP_WOLFE: do iterWolfe = 1, 12
        iGEMLineSearchIterationCount = iterWolfe
        dGEMLineSearchFinalStep = dStepLength

        call ApplyCEFLineSearchStep(dStepLength, lFeasible)
        if (lFeasible) then
            call CompChemicalPotential(lCompEverything)
            call CompFunctionNorm
            dTrialGibbs = CurrentActiveGibbs()
            dTrialMerit = ActiveMassMerit(dTrialGibbs)

            if ((.NOT.lSlopeRecorded).AND.(dStepLength > 0D0)) then
                call RecordDirectionalSlope(dStepLength, dTrialGibbs, dTrialMerit)
                lSlopeRecorded = .TRUE.
            end if

            if (dTrialGibbs < dBestGibbs) then
                dBestGibbs = dTrialGibbs
            end if

            if (dTrialMerit < dBestMerit) then
                dBestMerit = dTrialMerit
            end if

            if (dTrialMerit < dMeritMerit) then
                dMeritFunctionNorm = dGEMFunctionNorm
                dMeritStepLength = dStepLength
                dMeritGibbs = dTrialGibbs
                dMeritMerit = dTrialMerit
                dMeritMassNorm = dGEMMassBalanceNorm
            end if

            if (dGEMFunctionNorm < dBestFunctionNorm) then
                dBestFunctionNorm = dGEMFunctionNorm
                dBestStepLength = dStepLength
                dBestMolesSpecies = dMolesSpecies
                dBestMolesPhase = dMolesPhase
                dBestElementPotential = dElementPotential
                dBestPhaseSite = dGEMCEFPhaseSiteLast
            end if

            if (dGEMFunctionNorm < dLineSearchNormTolerance) exit LOOP_WOLFE
            if (dGEMFunctionNorm < 0.999D0 * dGEMFunctionNormLast) exit LOOP_WOLFE
        end if

        dStepLength = 0.5D0 * dStepLength
        if (dStepLength < 1D-12) exit LOOP_WOLFE
    end do LOOP_WOLFE

    if (dBestFunctionNorm < dGEMFunctionNormLast) then
        dMolesSpecies = dBestMolesSpecies
        dMolesPhase = dBestMolesPhase
        dElementPotential = dBestElementPotential
        dGEMCEFPhaseSiteLast = dBestPhaseSite
    else
        dMolesSpecies = dMolesSpeciesLast
        dMolesPhase = dMolesPhaseLast
        dElementPotential = dElementPotentialLast
        dGEMCEFPhaseSiteLast = dBasePhaseSite
        iGEMLineSearchNoDescent = 1
    end if

    call RebuildActiveCEFMolesFromStoredSites
    call RebuildActiveCEFMolFractionsFromMoles
    call CompChemicalPotential(lCompEverything)
    call CompFunctionNorm

    dGEMLineSearchBestNorm = dBestFunctionNorm
    dGEMLineSearchBestStep = dBestStepLength
    dGEMLineSearchBestGibbs = dBestGibbs
    dGEMLineSearchBestMerit = dBestMerit
    dGEMLineSearchFinalNorm = dGEMFunctionNorm
    dGEMLineSearchFinalGibbs = CurrentActiveGibbs()
    dGEMLineSearchFinalMerit = ActiveMassMerit(dGEMLineSearchFinalGibbs)
    dGEMLineSearchMinFinalPhaseMoles = MINVAL(dMolesPhase)
    dGEMLineSearchMeritCandNorm = dMeritFunctionNorm
    dGEMLineSearchMeritCandMass = dMeritMassNorm
    dGEMLineSearchMeritCandStep = dMeritStepLength
    dGEMLineSearchMeritCandGibbs = dMeritGibbs
    dGEMLineSearchMeritCandMerit = dMeritMerit
    call ClassifyNoDescent
    if (allocated(dGEMLSFinalPhaseMoles)) dGEMLSFinalPhaseMoles = dMolesPhase

    if (allocated(dBestMolesSpecies)) deallocate(dBestMolesSpecies)
    if (allocated(dBestMolesPhase)) deallocate(dBestMolesPhase)
    if (allocated(dBestElementPotential)) deallocate(dBestElementPotential)
    if (allocated(dBestPhaseSite)) deallocate(dBestPhaseSite)
    if (allocated(dBasePhaseSite)) deallocate(dBasePhaseSite)

    return

contains

    subroutine StoreRawCEFPhaseMoles

        implicit none

        integer :: iPhaseVar, iSlot

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            iSlot = iGEMCEFPhaseSlot(iPhaseVar)
            dGEMLSRawPhaseMoles(iSlot) = dMolesPhaseLast(iSlot) + dGEMCEFPhaseDirection(iPhaseVar)
        end do

        return

    end subroutine StoreRawCEFPhaseMoles

    real(8) function CurrentActiveGibbs()

        implicit none

        integer :: iSlot, iSolnPhase, iSpecies

        CurrentActiveGibbs = 0D0
        do iSlot = 1, nConPhases
            iSpecies = iAssemblage(iSlot)
            if (iSpecies > 0) then
                CurrentActiveGibbs = CurrentActiveGibbs + dMolesPhase(iSlot) * dStdGibbsEnergy(iSpecies)
            end if
        end do

        do iSlot = 1, nSolnPhases
            iSolnPhase = -iAssemblage(nElements - iSlot + 1)
            if ((iSolnPhase > 0).AND.(iSolnPhase <= nSolnPhasesSys)) then
                CurrentActiveGibbs = CurrentActiveGibbs + dGibbsSolnPhase(iSolnPhase)
            end if
        end do

        return

    end function CurrentActiveGibbs

    real(8) function ActiveMassMerit(dGibbsIn)

        implicit none

        real(8), intent(in) :: dGibbsIn

        ActiveMassMerit = dGibbsIn + 0.5D0 * dGEMMassBalanceNorm * dGEMMassBalanceNorm

        return

    end function ActiveMassMerit

    subroutine RecordDirectionalSlope(dStepIn, dTrialGibbsIn, dTrialMeritIn)

        implicit none

        real(8), intent(in) :: dStepIn, dTrialGibbsIn, dTrialMeritIn

        dGEMNewtonDirNormSlope = (dGEMFunctionNorm - dGEMLineSearchInitialNorm) / dStepIn
        dGEMNewtonDirGibbsSlope = (dTrialGibbsIn - dGEMLineSearchInitialGibbs) / dStepIn
        dGEMNewtonDirMeritSlope = (dTrialMeritIn - dGEMLineSearchInitialMerit) / dStepIn

        return

    end subroutine RecordDirectionalSlope

    subroutine ClassifyNoDescent

        implicit none

        real(8) :: dGibbsTol

        if (iGEMLineSearchNoDescent <= 0) return

        dGibbsTol = 1D-12 * DMAX1(1D0,DABS(dGEMLineSearchInitialGibbs))
        if (dGEMLineSearchBestGibbs < dGEMLineSearchInitialGibbs - dGibbsTol) then
            iGEMLineSearchNoDescentClass = 1
        else
            iGEMLineSearchNoDescentClass = 2
        end if

        return

    end subroutine ClassifyNoDescent

    subroutine LimitCEFInitialStep(dStepInOut)

        implicit none

        real(8), intent(inout) :: dStepInOut

        integer :: iPhaseVar, iSiteVar, iSlot, iPhaseID, iSolnPhase, iSub, iCon, iRef
        integer :: iFirstLocal, iSpeciesLocal, iSpeciesRefLocal
        real(8) :: dDirection, dAvailable, dSiteFloorLocal, dCandidate

        dStepInOut = DMIN1(1D0, DMAX1(0D0, dStepInOut))
        dSiteFloorLocal = DMAX1(dTraceSpeciesRemoveFraction, 1D-30)

        if (.NOT.UseMassBalanceProjection()) then
            do iPhaseVar = 1, nGEMCEFPhaseVariables
                iSlot = iGEMCEFPhaseSlot(iPhaseVar)
                if (iSlot <= 0) cycle
                dDirection = dGEMCEFPhaseDirection(iPhaseVar)
                if (dDirection < -1D-300) then
                    dAvailable = dMolesPhaseLast(iSlot) - dTolerance(8)
                    if (dAvailable <= 0D0) then
                        dStepInOut = 0D0
                        return
                    end if
                    dCandidate = 0.99D0 * dAvailable / (-dDirection)
                    dStepInOut = DMIN1(dStepInOut, DMAX1(0D0, dCandidate))
                end if
            end do
        end if

        do iSiteVar = 1, nGEMCEFSiteVariables
            iPhaseID = iGEMCEFVarPhaseID(iSiteVar)
            iSub = iGEMCEFVarSub(iSiteVar)
            iCon = iGEMCEFVarCon(iSiteVar)
            iRef = iGEMCEFVarRef(iSiteVar)
            if ((iCon <= 0).OR.(iRef <= 0)) cycle

            dDirection = dGEMCEFSiteDirection(iSiteVar)
            if (iPhaseID > 0) then
                if (iSub <= 0) cycle
                if (dDirection < -1D-300) then
                    dAvailable = dGEMCEFPhaseSiteLast(iGEMCEFVarPhaseVar(iSiteVar),iSub,iCon) - &
                        dSiteFloorLocal
                    if (dAvailable <= 0D0) then
                        dStepInOut = 0D0
                        return
                    end if
                    dCandidate = 0.99D0 * dAvailable / (-dDirection)
                    dStepInOut = DMIN1(dStepInOut, DMAX1(0D0, dCandidate))
                else if (dDirection > 1D-300) then
                    dAvailable = dGEMCEFPhaseSiteLast(iGEMCEFVarPhaseVar(iSiteVar),iSub,iRef) - &
                        dSiteFloorLocal
                    if (dAvailable <= 0D0) then
                        dStepInOut = 0D0
                        return
                    end if
                    dCandidate = 0.99D0 * dAvailable / dDirection
                    dStepInOut = DMIN1(dStepInOut, DMAX1(0D0, dCandidate))
                end if
            else
                iSolnPhase = iGEMCEFVarSolnPhase(iSiteVar)
                iFirstLocal = nSpeciesPhase(iSolnPhase-1) + 1
                iSpeciesLocal = iFirstLocal + iCon - 1
                iSpeciesRefLocal = iFirstLocal + iRef - 1
                if (dDirection < -1D-300) then
                    dAvailable = dMolFraction(iSpeciesLocal) - dSiteFloorLocal
                    if (dAvailable <= 0D0) then
                        dStepInOut = 0D0
                        return
                    end if
                    dCandidate = 0.99D0 * dAvailable / (-dDirection)
                    dStepInOut = DMIN1(dStepInOut, DMAX1(0D0, dCandidate))
                else if (dDirection > 1D-300) then
                    dAvailable = dMolFraction(iSpeciesRefLocal) - dSiteFloorLocal
                    if (dAvailable <= 0D0) then
                        dStepInOut = 0D0
                        return
                    end if
                    dCandidate = 0.99D0 * dAvailable / dDirection
                    dStepInOut = DMIN1(dStepInOut, DMAX1(0D0, dCandidate))
                end if
            end if
        end do

        return

    end subroutine LimitCEFInitialStep

    subroutine ApplyCEFLineSearchStep(dStepIn, lFeasibleOut)

        implicit none

        real(8), intent(in) :: dStepIn
        logical, intent(out) :: lFeasibleOut

        integer :: iPhaseVar, iSiteVar, iSlot, iSolnPhase, iPhaseID
        integer :: iSub, iCon, iRef, iSpecies, iFirstLocal, iLastLocal
        integer :: iSpeciesLocal, iSpeciesRefLocal, iLocal
        real(8) :: dPhaseAmount, dSumLocal, dFractionFloorLocal
        real(8), allocatable :: dSiteTrial(:,:,:)
        real(8), allocatable :: dFractionTrial(:)

        lFeasibleOut = .TRUE.
        allocate(dSiteTrial(MAX(1,nGEMCEFPhaseVariables),nMaxSublatticeSys,nMaxConstituentSys))
        allocate(dFractionTrial(nSpecies))

        dMolesSpecies = dMolesSpeciesLast
        dMolesPhase = dMolesPhaseLast
        dElementPotential = dElementPotentialLast + dStepIn * dGEMCEFElementDirection
        dSiteTrial = dGEMCEFPhaseSiteLast
        dFractionTrial = dMolFraction
        dFractionFloorLocal = DMAX1(dTraceSpeciesRemoveFraction, 1D-30)

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            iSolnPhase = iGEMCEFPhaseSoln(iPhaseVar)
            if (iSolnPhase <= 0) cycle
            iFirstLocal = nSpeciesPhase(iSolnPhase-1) + 1
            iLastLocal = nSpeciesPhase(iSolnPhase)
            dPhaseAmount = dMolesPhaseLast(iGEMCEFPhaseSlot(iPhaseVar))
            if (dPhaseAmount > dTolerance(8)) then
                dFractionTrial(iFirstLocal:iLastLocal) = &
                    dMolesSpeciesLast(iFirstLocal:iLastLocal) / dPhaseAmount
            end if
        end do

        if (.NOT.UseMassBalanceProjection()) then
            do iPhaseVar = 1, nGEMCEFPhaseVariables
                iSlot = iGEMCEFPhaseSlot(iPhaseVar)
                dPhaseAmount = dMolesPhaseLast(iSlot) + dStepIn * dGEMCEFPhaseDirection(iPhaseVar)
                if (dPhaseAmount <= dTolerance(8)) then
                    lFeasibleOut = .FALSE.
                    exit
                end if
                dMolesPhase(iSlot) = dPhaseAmount
            end do
        end if

        if (lFeasibleOut) then
            do iSiteVar = 1, nGEMCEFSiteVariables
                iPhaseID = iGEMCEFVarPhaseID(iSiteVar)
                iCon = iGEMCEFVarCon(iSiteVar)
                iRef = iGEMCEFVarRef(iSiteVar)
                if (iPhaseID > 0) then
                    iSub = iGEMCEFVarSub(iSiteVar)
                    iPhaseVar = iGEMCEFVarPhaseVar(iSiteVar)
                    dSiteTrial(iPhaseVar,iSub,iCon) = dSiteTrial(iPhaseVar,iSub,iCon) + &
                        dStepIn * dGEMCEFSiteDirection(iSiteVar)
                    dSiteTrial(iPhaseVar,iSub,iRef) = dSiteTrial(iPhaseVar,iSub,iRef) - &
                        dStepIn * dGEMCEFSiteDirection(iSiteVar)
                else
                    iSolnPhase = iGEMCEFVarSolnPhase(iSiteVar)
                    iFirstLocal = nSpeciesPhase(iSolnPhase-1) + 1
                    iSpeciesLocal = iFirstLocal + iCon - 1
                    iSpeciesRefLocal = iFirstLocal + iRef - 1
                    dFractionTrial(iSpeciesLocal) = dFractionTrial(iSpeciesLocal) + &
                        dStepIn * dGEMCEFSiteDirection(iSiteVar)
                    dFractionTrial(iSpeciesRefLocal) = dFractionTrial(iSpeciesRefLocal) - &
                        dStepIn * dGEMCEFSiteDirection(iSiteVar)
                end if
            end do

            call CheckCEFSiteFeasibility(dSiteTrial, lFeasibleOut)
            if (lFeasibleOut) call CheckNonCEFFractionFeasibility(dFractionTrial, lFeasibleOut)
            if (lFeasibleOut .AND. UseMassBalanceProjection()) then
                call ProjectMixedPhaseAmounts(dSiteTrial, dFractionTrial, lFeasibleOut)
            end if
        end if

        if (lFeasibleOut) then
            do iPhaseVar = 1, nGEMCEFPhaseVariables
                iSolnPhase = iGEMCEFPhaseSoln(iPhaseVar)
                if (iSolnPhase > 0) then
                    if (IsCEFLineSearchPhase(iSolnPhase)) then
                        iPhaseID = iPhaseSublattice(iSolnPhase)
                        ! Convert the accepted site trial back to product
                        ! fractions only after the site-coordinate step passes
                        ! positivity and mass-balance feasibility.
                        call SetCEFMolesFromSite(iSolnPhase, iPhaseID, &
                            dSiteTrial(iPhaseVar,:,:), dMolesPhase(iGEMCEFPhaseSlot(iPhaseVar)))
                        dGEMCEFPhaseSiteLast(iPhaseVar,:,:) = dSiteTrial(iPhaseVar,:,:)
                        if (allocated(dActiveSlotSiteFraction)) then
                            iSlot = iGEMCEFPhaseSlot(iPhaseVar)
                            if ((iSlot > 0).AND.(iSlot <= SIZE(dActiveSlotSiteFraction,1))) then
                                dActiveSlotSiteFraction(iSlot,:,:) = dSiteTrial(iPhaseVar,:,:)
                            end if
                        end if
                    else
                        iFirstLocal = nSpeciesPhase(iSolnPhase-1) + 1
                        iLastLocal = nSpeciesPhase(iSolnPhase)
                        dSumLocal = SUM(dFractionTrial(iFirstLocal:iLastLocal))
                        if (dSumLocal <= 0D0) then
                            lFeasibleOut = .FALSE.
                            exit
                        end if
                        do iLocal = iFirstLocal, iLastLocal
                            dMolFraction(iLocal) = DMAX1(dFractionTrial(iLocal), dFractionFloorLocal) / dSumLocal
                            dMolesSpecies(iLocal) = dMolesPhase(iGEMCEFPhaseSlot(iPhaseVar)) * dMolFraction(iLocal)
                        end do
                    end if
                else
                    iSpecies = iGEMCEFPhaseSpecies(iPhaseVar)
                    dMolesSpecies(iSpecies) = dMolesPhase(iGEMCEFPhaseSlot(iPhaseVar))
                    dMolFraction(iSpecies) = 1D0
                end if
            end do
        end if

        if (.NOT.lFeasibleOut) then
            dMolesSpecies = dMolesSpeciesLast
            dMolesPhase = dMolesPhaseLast
            dElementPotential = dElementPotentialLast
        end if

        if (allocated(dSiteTrial)) deallocate(dSiteTrial)
        if (allocated(dFractionTrial)) deallocate(dFractionTrial)

        return

    end subroutine ApplyCEFLineSearchStep

    subroutine CheckCEFSiteFeasibility(dSiteIn, lFeasibleOut)

        implicit none

        real(8), intent(in) :: dSiteIn(MAX(1,nGEMCEFPhaseVariables),nMaxSublatticeSys,nMaxConstituentSys)
        logical, intent(out) :: lFeasibleOut

        integer :: iPhaseVar, iSolnPhase, iPhaseID, iSub, iCon
        real(8) :: dSiteFloorLocal

        lFeasibleOut = .TRUE.
        dSiteFloorLocal = DMAX1(dTraceSpeciesRemoveFraction, 1D-30)

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            iSolnPhase = iGEMCEFPhaseSoln(iPhaseVar)
            if (iSolnPhase <= 0) cycle
            if (.NOT.IsCEFLineSearchPhase(iSolnPhase)) cycle
            iPhaseID = iPhaseSublattice(iSolnPhase)
            do iSub = 1, nSublatticePhase(iPhaseID)
                do iCon = 1, nConstituentSublattice(iPhaseID,iSub)
                    if (dSiteIn(iPhaseVar,iSub,iCon) <= dSiteFloorLocal) then
                        lFeasibleOut = .FALSE.
                        return
                    end if
                end do
            end do
        end do

        return

    end subroutine CheckCEFSiteFeasibility

    subroutine CheckNonCEFFractionFeasibility(dFractionIn, lFeasibleOut)

        implicit none

        real(8), intent(in) :: dFractionIn(nSpecies)
        logical, intent(out) :: lFeasibleOut

        integer :: iPhaseVar, iSolnPhase, iFirstLocal, iLastLocal, iSpeciesLocal
        real(8) :: dFractionFloorLocal

        lFeasibleOut = .TRUE.
        dFractionFloorLocal = DMAX1(dTraceSpeciesRemoveFraction, 1D-30)

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            iSolnPhase = iGEMCEFPhaseSoln(iPhaseVar)
            if (iSolnPhase <= 0) cycle
            if (IsCEFLineSearchPhase(iSolnPhase)) cycle
            iFirstLocal = nSpeciesPhase(iSolnPhase-1) + 1
            iLastLocal = nSpeciesPhase(iSolnPhase)
            do iSpeciesLocal = iFirstLocal, iLastLocal
                if (dFractionIn(iSpeciesLocal) <= dFractionFloorLocal) then
                    lFeasibleOut = .FALSE.
                    return
                end if
            end do
        end do

        return

    end subroutine CheckNonCEFFractionFeasibility

    logical function UseMassBalanceProjection()

        implicit none

        UseMassBalanceProjection = (nGEMCEFPhaseVariables == nElements).AND.(.NOT.lGEMCEFBndPhaseActive)

        return

    end function UseMassBalanceProjection

    subroutine ProjectMixedPhaseAmounts(dSiteIn, dFractionIn, lFeasibleOut)

        implicit none

        real(8), intent(in) :: dSiteIn(MAX(1,nGEMCEFPhaseVariables),nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(in) :: dFractionIn(nSpecies)
        logical, intent(out) :: lFeasibleOut

        integer :: iPhaseVar, iSlot, iSolnPhase, iPhaseID, iSpecies, iElem, INFO_LU
        integer, allocatable :: IPIVLocal(:)
        real(8), allocatable :: ALocal(:,:), BLocal(:), dPhaseComp(:)

        lFeasibleOut = .TRUE.
        if (.NOT.UseMassBalanceProjection()) return

        allocate(ALocal(nElements,nElements), BLocal(nElements), IPIVLocal(nElements))
        allocate(dPhaseComp(nElements))
        ALocal = 0D0
        BLocal = dMolesElement(1:nElements)
        IPIVLocal = 0

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            dPhaseComp = 0D0
            iSolnPhase = iGEMCEFPhaseSoln(iPhaseVar)
            if (iSolnPhase > 0) then
                if (IsCEFLineSearchPhase(iSolnPhase)) then
                    iPhaseID = iPhaseSublattice(iSolnPhase)
                    call ComputeCEFLineSearchPhaseComposition(iSolnPhase, iPhaseID, &
                        dSiteIn(iPhaseVar,:,:), dPhaseComp)
                else
                    call ComputeNonCEFLineSearchPhaseComposition(iSolnPhase, dFractionIn, dPhaseComp)
                end if
            else
                iSpecies = iGEMCEFPhaseSpecies(iPhaseVar)
                dPhaseComp = dStoichSpecies(iSpecies,1:nElements)
            end if

            do iElem = 1, nElements
                ALocal(iElem,iPhaseVar) = dPhaseComp(iElem)
            end do
        end do

        call DGESV(nElements, 1, ALocal, nElements, IPIVLocal, BLocal, nElements, INFO_LU)
        if (INFO_LU /= 0) lFeasibleOut = .FALSE.

        if (lFeasibleOut) then
            do iPhaseVar = 1, nGEMCEFPhaseVariables
                if (BLocal(iPhaseVar) <= dTolerance(8)) then
                    lFeasibleOut = .FALSE.
                    exit
                end if
            end do
        end if

        if (lFeasibleOut) then
            do iPhaseVar = 1, nGEMCEFPhaseVariables
                iSlot = iGEMCEFPhaseSlot(iPhaseVar)
                dMolesPhase(iSlot) = BLocal(iPhaseVar)
            end do
        end if

        if (allocated(ALocal)) deallocate(ALocal)
        if (allocated(BLocal)) deallocate(BLocal)
        if (allocated(IPIVLocal)) deallocate(IPIVLocal)
        if (allocated(dPhaseComp)) deallocate(dPhaseComp)

        return

    end subroutine ProjectMixedPhaseAmounts

    subroutine ComputeCEFLineSearchPhaseComposition(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, dCompositionOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dCompositionOut(nElements)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dProductSum
        real(8), allocatable :: dProductSpecies(:)

        dCompositionOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        allocate(dProductSpecies(iLastLocal-iFirstLocal+1))
        dProductSpecies = 0D0
        dProductSum = 0D0

        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 1D-75)
            end do
            dProductSpecies(iLocal) = dProduct
            dProductSum = dProductSum + dProduct
        end do

        if (dProductSum > 0D0) then
            do iLocalSpecies = iFirstLocal, iLastLocal
                iLocal = iLocalSpecies - iFirstLocal + 1
                dCompositionOut = dCompositionOut + &
                    dProductSpecies(iLocal) * dStoichSpecies(iLocalSpecies,1:nElements) / dProductSum
            end do
        end if

        if (allocated(dProductSpecies)) deallocate(dProductSpecies)

        return

    end subroutine ComputeCEFLineSearchPhaseComposition

    subroutine ComputeNonCEFLineSearchPhaseComposition(iSolnPhaseIndexIn, dFractionIn, dCompositionOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn
        real(8), intent(in) :: dFractionIn(nSpecies)
        real(8), intent(out) :: dCompositionOut(nElements)

        integer :: iFirstLocal, iLastLocal, iSpeciesLocal

        dCompositionOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)

        do iSpeciesLocal = iFirstLocal, iLastLocal
            dCompositionOut = dCompositionOut + dFractionIn(iSpeciesLocal) * &
                dStoichSpecies(iSpeciesLocal,1:nElements) / DBLE(iParticlesPerMole(iSpeciesLocal))
        end do

        return

    end subroutine ComputeNonCEFLineSearchPhaseComposition

    subroutine SetCEFMolesFromSite(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, dPhaseAmountIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(in) :: dPhaseAmountIn

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dSumLocal

        ! The CEF line search state is dSiteIn.  Product fractions are rebuilt
        ! here so existing Gibbs-model and reporting arrays stay synchronized.
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        dSumLocal = 0D0
        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 1D-75)
            end do
            dMolFraction(iLocalSpecies) = dProduct
            dSumLocal = dSumLocal + dProduct
        end do

        if (dSumLocal > 0D0) then
            do iLocalSpecies = iFirstLocal, iLastLocal
                dMolFraction(iLocalSpecies) = dMolFraction(iLocalSpecies) / dSumLocal
                dMolesSpecies(iLocalSpecies) = dPhaseAmountIn * dMolFraction(iLocalSpecies)
            end do
        end if

        return

    end subroutine SetCEFMolesFromSite

    subroutine RebuildActiveCEFMolFractionsFromMoles

        implicit none

        integer :: iPhaseVar, iSolnPhase, iFirstLocal, iLastLocal, iSpeciesLocal
        real(8) :: dSumLocal

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            iSolnPhase = iGEMCEFPhaseSoln(iPhaseVar)
            if (iSolnPhase <= 0) cycle
            if (IsCEFLineSearchPhase(iSolnPhase)) cycle
            iFirstLocal = nSpeciesPhase(iSolnPhase-1) + 1
            iLastLocal = nSpeciesPhase(iSolnPhase)
            dSumLocal = SUM(dMolesSpecies(iFirstLocal:iLastLocal))
            if (dSumLocal > 0D0) then
                do iSpeciesLocal = iFirstLocal, iLastLocal
                    dMolFraction(iSpeciesLocal) = dMolesSpecies(iSpeciesLocal) / dSumLocal
                end do
            end if
        end do

        return

    end subroutine RebuildActiveCEFMolFractionsFromMoles

    subroutine RebuildActiveCEFMolesFromStoredSites

        implicit none

        integer :: iPhaseVar, iSolnPhase, iPhaseID
        integer :: iSlot

        do iPhaseVar = 1, nGEMCEFPhaseVariables
            iSolnPhase = iGEMCEFPhaseSoln(iPhaseVar)
            if (iSolnPhase <= 0) cycle
            if (.NOT.IsCEFLineSearchPhase(iSolnPhase)) cycle
            iPhaseID = iPhaseSublattice(iSolnPhase)
            iSlot = iGEMCEFPhaseSlot(iPhaseVar)
            call SetCEFMolesFromSite(iSolnPhase, iPhaseID, dGEMCEFPhaseSiteLast(iPhaseVar,:,:), &
                dMolesPhase(iSlot))
            if (allocated(dActiveSlotSiteFraction)) then
                if ((iSlot > 0).AND.(iSlot <= SIZE(dActiveSlotSiteFraction,1))) then
                    dActiveSlotSiteFraction(iSlot,:,:) = dGEMCEFPhaseSiteLast(iPhaseVar,:,:)
                end if
            end if
        end do

        return

    end subroutine RebuildActiveCEFMolesFromStoredSites

    logical function IsCEFLineSearchPhase(iSolnPhaseIndexIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn
        character(len=8) :: cTypeLocal

        IsCEFLineSearchPhase = .FALSE.
        if ((iSolnPhaseIndexIn <= 0).OR.(iSolnPhaseIndexIn > nSolnPhasesSys)) return
        if (.NOT.allocated(iPhaseSublattice)) return
        if (iSolnPhaseIndexIn > SIZE(iPhaseSublattice)) return
        if (iPhaseSublattice(iSolnPhaseIndexIn) <= 0) return

        cTypeLocal = TRIM(cSolnPhaseType(iSolnPhaseIndexIn))
        IsCEFLineSearchPhase = (TRIM(cTypeLocal) == 'SUBL').OR.&
            (TRIM(cTypeLocal) == 'SUBLM').OR.&
            (TRIM(cTypeLocal) == 'SUBOM')

        return

    end function IsCEFLineSearchPhase

end subroutine GEMLineSearchCEF
