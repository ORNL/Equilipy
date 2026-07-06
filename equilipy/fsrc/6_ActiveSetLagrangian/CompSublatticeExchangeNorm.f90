!> \brief Compute active CEF site-gradient exchange residuals for GEM convergence.
!!
!! \details Active CEF phases are minimized in site-fraction coordinates.
!! Their stationarity residual is therefore the difference
!! between site-gradient grand potentials on the same sublattice, not the
!! legacy pseudo-endmember chemical-potential residual.
!!
subroutine CompSublatticeExchangeNorm(dExchangeNorm)

    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    CompSublatticeExchangeNorm.f90
    !> \brief   Compute site-gradient exchange residuals for active CEF phases.
    !> \author  S.Y. Kwon
    !> \date    Jun. 25, 2026
    !> \sa      CompFunctionNorm.f90
    !> \sa      CompGradientSUBL.f90
    !
    ! Revisions:
    ! ==========
    !
    !   Date            Programmer          Description of change
    !   ----            ----------          ---------------------
    !   06/25/2026      S.Y. Kwon           Replaced the pseudo-endmember residual for SUBOM phases with
    !                                        analytical site-gradient exchange residuals.
    !   06/25/2026      S.Y. Kwon           Extended the same residual to active SUBL/SUBLM CEF phases
    !                                        for site-fraction Lagrangian GEM.
!   06/25/2026      S.Y. Kwon           Added the CEF phase grand-potential residual so the GEM norm
!                                        matches the mixed-coordinate KKT equations.
!   06/25/2026      S.Y. Kwon           Computed CEF composition projections from numeric species
!                                        stoichiometry to match GEM mass balance.
!   06/25/2026      S.Y. Kwon           Return zero unless the current fixed assemblage is using the
!                                        CEF site-fraction KKT path.
!   07/01/2026      S.Y. Kwon           Applied site-fraction complementarity for trace constituents so
!                                        lower-bound species with positive residuals do not block convergence.
!   07/05/2026      S.Y. Kwon           Excluded zero-amount solution remnants from fixed-assemblage
!                                        postprocess exchange diagnostics.
!   07/05/2026      S.Y. Kwon           Scaled trace-element-edge site-exchange residuals by active
!                                        phase fraction to avoid unit-dependent trace-phase norm failures.
!   07/05/2026      S.Y. Kwon           Applied the same trace-edge phase-fraction scaling to the CEF
!                                        scalar phase-plane residual in the exchange norm.
!   07/05/2026      S.Y. Kwon           Attenuated lower-bound trace-site residuals only in trace-edge
!                                        systems when the carrying active phase is itself trace-scale.
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to compute the Lagrangian
    !! stationarity residual for active CEF phases in the same variables used
    !! by CEF site-fraction subminimization.  For each sublattice, the residual is
    !! phase grand-potential residual
    !! \f$g_p-\lambda\cdot a_p\f$ and the site-exchange residual
    !! \f$(\partial G/\partial y_{s,c} - \lambda\cdot\partial M/\partial y_{s,c})
    !! - (\partial G/\partial y_{s,r} - \lambda\cdot\partial M/\partial y_{s,r})\f$,
    !! where r is the largest site-fraction constituent on that sublattice.
    !
    !
    ! Required input variables:
    ! =========================
    !
    ! iAssemblage           Current active phase assemblage.
    ! dMolFraction          Current active solution endmember fractions.
    ! dElementPotential     Current Lagrangian elemental potentials.
    ! cSolnPhaseType        Solution model type; SUBL, SUBLM, and SUBOM are handled here.
    !
    !
    ! Output/updated variables:
    ! =========================
    !
    !> \param[out] dExchangeNorm  Square-rooted site-gradient exchange norm.
    !
    ! dSiteFraction         Rebuilt by CompGradientSUBL from dMolFraction.
    !
    !
    ! Called subroutines/functions:
    ! =============================
    !
    ! CompGradientSUBL      Computes fixed-state CEF scalar gradients.
    !
    !
    ! Primary callers:
    ! ================
    !
    ! CompFunctionNorm      Combines this residual with mass balance and non-CEF
    !                       chemical-potential residuals.
    !
    !
    ! Numerical assumptions:
    ! ======================
    !
    ! - Lower-bound site constituents satisfy complementarity.  A lower-bound
    !   constituent only contributes when the residual indicates it should grow.
    ! - In trace-element-edge active sets, CEF scalar phase-plane and
    !   site-fraction stationarity residuals are scaled by active phase
    !   fraction so trace phases do not dominate the extensive-system
    !   convergence norm.  Order/disorder mapped active sets remain on the
    !   unscaled norm because their near-degenerate transfer directions are
    !   separately classified and pinned.
    ! - SUBOM uses the fixed ordered-state gradient from CompGradientSUBL, which
    !   includes the disorder correction.  SUBL/SUBLM use the same CEF scalar
    !   gradient without an order/disorder correction.
    !
    !-------------------------------------------------------------------------------------------------------------

    USE ModuleThermo
    USE ModuleGEMSolver
    USE ModuleSubMin, ONLY: dSubMinHenrianTraceThreshold

    implicit none

    real(8), intent(out) :: dExchangeNorm

    integer :: k, l, s, c, cRef, iSlot, iThermoPhase
    integer :: iSublPhase, nSublattice, nConstituent, nSiteCapacity, nSiteOut, iInfo
    integer :: iSiteIndex(nMaxSublatticeSys,nMaxConstituentSys)
    integer :: iSiteCursor
    real(8) :: dSiteFloor, dSiteTraceThreshold, dY, dYRef
    real(8) :: dGrad, dGradRef, dResidual, dWeightedResidual, dPhaseFraction
    real(8) :: dTracePhaseExchangeWeight
    real(8) :: dTotalActivePhaseAmount
    real(8) :: dExchangeNormSquared, dPhaseResidual
    real(8) :: dScalarGibbs, dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity
    real(8), allocatable :: dSiteGradient(:), dSiteGradientH(:)
    real(8), allocatable :: dSiteGradientS(:), dSiteGradientCp(:)
    real(8), allocatable :: dPhaseComp(:), dSiteCompDeriv(:)
    real(8), allocatable :: dSiteLocal(:,:)
    logical :: lIsCEFPhase, lScaleExchangeByPhaseFraction
    logical :: lActiveSetContainsOrderDisorderMapped, lTraceElementEdgeSystem
    real(8), parameter :: dTracePhaseExchangeFraction = 1D-3

    dExchangeNorm = 0D0
    dExchangeNormSquared = 0D0
    iGEMMaxSublExchangeSlot = 0
    iGEMMaxSublExchangePhase = 0
    iGEMMaxSublExchangeSite = 0
    iGEMMaxSublExchangeConstituent = 0
    dGEMMaxSublExchangeResidual = 0D0
    dGEMMaxSublExchangeWeightedResidual = 0D0

    if (.NOT.lGEMCEFSiteLagrangianActive) return
    if (.NOT.allocated(iPhaseSublattice)) return
    if (.NOT.allocated(dSiteFraction)) return
    if (.NOT.allocated(iConstituentSublattice)) return

    dSiteFloor = DMAX1(dTraceSpeciesRemoveFraction, 1D-30)
    dSiteTraceThreshold = DMAX1(dSiteFloor, dSubMinHenrianTraceThreshold)
    nSiteCapacity = nMaxSublatticeSys * nMaxConstituentSys
    dTotalActivePhaseAmount = DMAX1(SUM(DABS(dMolesPhase(1:nElements))), 1D-300)
    lActiveSetContainsOrderDisorderMapped = ActiveSetContainsOrderDisorderMappedForExchange()
    lTraceElementEdgeSystem = TraceElementEdgeSystemForExchange()

    allocate(dSiteGradient(nSiteCapacity), dSiteGradientH(nSiteCapacity), &
        dSiteGradientS(nSiteCapacity), dSiteGradientCp(nSiteCapacity))
    allocate(dPhaseComp(nElements), dSiteCompDeriv(nElements))
    allocate(dSiteLocal(nMaxSublatticeSys,nMaxConstituentSys))

    do l = 1, nSolnPhases
        iSlot = nElements - l + 1
        if (lPostProcess.AND.(dMolesPhase(iSlot) <= DMAX1(dTolerance(8), 1D-300))) cycle
        k = -iAssemblage(iSlot)
        if (k <= 0) cycle
        iThermoPhase = ActiveSlotThermoPhaseForExchange(iSlot, k)
        dPhaseFraction = DMAX1(dMolesPhase(iSlot), 0D0) / dTotalActivePhaseAmount
        if (iThermoPhase > SIZE(iPhaseSublattice)) cycle
        lIsCEFPhase = (TRIM(cSolnPhaseType(iThermoPhase)) == 'SUBL').OR.&
            (TRIM(cSolnPhaseType(iThermoPhase)) == 'SUBLM').OR.&
            (TRIM(cSolnPhaseType(iThermoPhase)) == 'SUBOM')
        if (.NOT.lIsCEFPhase) cycle
        lScaleExchangeByPhaseFraction = lTraceElementEdgeSystem.AND.&
            (.NOT.lActiveSetContainsOrderDisorderMapped)

        iSublPhase = iPhaseSublattice(iThermoPhase)
        if (iSublPhase <= 0) cycle
        call ActiveSlotSiteForExchange(iThermoPhase, iSublPhase, iSlot, dSiteLocal)
        call SetExchangeMolFractionFromSite(iThermoPhase, iSublPhase, dSiteLocal)

        nSublattice = nSublatticePhase(iSublPhase)
        iSiteIndex = 0
        iSiteCursor = 0
        do s = 1, nSublattice
            do c = 1, nConstituentSublattice(iSublPhase,s)
                iSiteCursor = iSiteCursor + 1
                iSiteIndex(s,c) = iSiteCursor
            end do
        end do

        call CompGradientSUBL(iThermoPhase, nSiteCapacity, dSiteGradient, dScalarGibbs, &
            dSiteGradientH, dSiteGradientS, dSiteGradientCp, &
            dScalarEnthalpy, dScalarEntropy, dScalarHeatCapacity, nSiteOut, iInfo)
        if (iInfo /= 0) cycle

        call ComputeCEFPhaseComposition(iThermoPhase, iSublPhase, dSiteLocal, dPhaseComp)
        dPhaseResidual = dScalarGibbs - SUM(dElementPotential(1:nElements) * dPhaseComp(1:nElements))
        if (lScaleExchangeByPhaseFraction) then
            dWeightedResidual = dPhaseFraction * dPhaseResidual
        else
            dWeightedResidual = dPhaseResidual
        end if
        dExchangeNormSquared = dExchangeNormSquared + dWeightedResidual**2
        if (DABS(dWeightedResidual) > DABS(dGEMMaxSublExchangeWeightedResidual)) then
            iGEMMaxSublExchangeSlot = iSlot
            iGEMMaxSublExchangePhase = iThermoPhase
            iGEMMaxSublExchangeSite = 0
            iGEMMaxSublExchangeConstituent = 0
            dGEMMaxSublExchangeResidual = dPhaseResidual
            dGEMMaxSublExchangeWeightedResidual = dWeightedResidual
        end if

        do s = 1, nSublattice
            nConstituent = nConstituentSublattice(iSublPhase,s)
            if (nConstituent <= 1) cycle

            cRef = 1
            do c = 2, nConstituent
                if (dSiteLocal(s,c) > dSiteLocal(s,cRef)) cRef = c
            end do

            dYRef = DMAX1(dSiteLocal(s,cRef), dSiteFloor)
            if (dYRef <= 0D0) cycle
            dGradRef = dSiteGradient(iSiteIndex(s,cRef))

            do c = 1, nConstituent
                if (c == cRef) cycle

                dY = DMAX1(dSiteLocal(s,c), dSiteFloor)
                if (dY <= 0D0) cycle
                dGrad = dSiteGradient(iSiteIndex(s,c))

                call ComputeCEFCompositionDerivative(iThermoPhase, iSublPhase, s, c, cRef, &
                    dSiteLocal, dSiteCompDeriv)
                dResidual = dGrad - dGradRef - &
                    SUM(dElementPotential(1:nElements) * dSiteCompDeriv(1:nElements))

                if ((dSiteLocal(s,c) > dSiteTraceThreshold).OR.&
                    (dResidual < -dTraceSpeciesResidualTolerance)) then
                    if (lScaleExchangeByPhaseFraction) then
                        dWeightedResidual = dPhaseFraction * dResidual
                    else
                        dWeightedResidual = dResidual
                    end if
                    if (lScaleExchangeByPhaseFraction.AND.&
                        (dSiteLocal(s,c) <= dSiteTraceThreshold).AND.&
                        (dResidual < -dTraceSpeciesResidualTolerance).AND.&
                        (dPhaseFraction < dTracePhaseExchangeFraction)) then
                        dTracePhaseExchangeWeight = &
                            DMAX1(dPhaseFraction / dTracePhaseExchangeFraction, 0D0)
                        dWeightedResidual = dWeightedResidual * dTracePhaseExchangeWeight
                    end if
                    dExchangeNormSquared = dExchangeNormSquared + dWeightedResidual**2
                    if (DABS(dWeightedResidual) > DABS(dGEMMaxSublExchangeWeightedResidual)) then
                        iGEMMaxSublExchangeSlot = iSlot
                        iGEMMaxSublExchangePhase = iThermoPhase
                        iGEMMaxSublExchangeSite = s
                        iGEMMaxSublExchangeConstituent = c
                        dGEMMaxSublExchangeResidual = dResidual
                        dGEMMaxSublExchangeWeightedResidual = dWeightedResidual
                    end if
                end if
            end do
        end do
    end do

    dExchangeNorm = DSQRT(dExchangeNormSquared)

    if (allocated(dSiteGradient)) deallocate(dSiteGradient)
    if (allocated(dSiteGradientH)) deallocate(dSiteGradientH)
    if (allocated(dSiteGradientS)) deallocate(dSiteGradientS)
    if (allocated(dSiteGradientCp)) deallocate(dSiteGradientCp)
    if (allocated(dPhaseComp)) deallocate(dPhaseComp)
    if (allocated(dSiteCompDeriv)) deallocate(dSiteCompDeriv)
    if (allocated(dSiteLocal)) deallocate(dSiteLocal)

    return

contains

    integer function ActiveSlotThermoPhaseForExchange(iSlotIn, iDisplayPhaseIn)

        implicit none

        integer, intent(in) :: iSlotIn, iDisplayPhaseIn

        ActiveSlotThermoPhaseForExchange = iDisplayPhaseIn
        if (.NOT.allocated(iActiveSlotThermoPhase)) return
        if ((iSlotIn <= 0).OR.(iSlotIn > SIZE(iActiveSlotThermoPhase))) return
        if (iActiveSlotThermoPhase(iSlotIn) > 0) ActiveSlotThermoPhaseForExchange = iActiveSlotThermoPhase(iSlotIn)

        return

    end function ActiveSlotThermoPhaseForExchange

    logical function PhaseIsOrderDisorderMappedForExchange(iPhaseIn)

        implicit none

        integer, intent(in) :: iPhaseIn
        integer :: iOrderedPhase

        PhaseIsOrderDisorderMappedForExchange = .FALSE.
        if (.NOT.allocated(iDisorderedPhase)) return
        if (iPhaseIn <= 0) return

        if (iPhaseIn <= SIZE(iDisorderedPhase)) then
            if (iDisorderedPhase(iPhaseIn) > 0) then
                PhaseIsOrderDisorderMappedForExchange = .TRUE.
                return
            end if
        end if

        do iOrderedPhase = 1, MIN(nSolnPhasesSys, SIZE(iDisorderedPhase))
            if (iDisorderedPhase(iOrderedPhase) == iPhaseIn) then
                PhaseIsOrderDisorderMappedForExchange = .TRUE.
                return
            end if
        end do

        return

    end function PhaseIsOrderDisorderMappedForExchange

    logical function ActiveSetContainsOrderDisorderMappedForExchange()

        implicit none

        integer :: iSlotLocal, iPhaseLocal

        ActiveSetContainsOrderDisorderMappedForExchange = .FALSE.
        do iSlotLocal = 1, nElements
            if (iAssemblage(iSlotLocal) >= 0) cycle
            iPhaseLocal = -iAssemblage(iSlotLocal)
            iPhaseLocal = ActiveSlotThermoPhaseForExchange(iSlotLocal, iPhaseLocal)
            if (PhaseIsOrderDisorderMappedForExchange(iPhaseLocal)) then
                ActiveSetContainsOrderDisorderMappedForExchange = .TRUE.
                return
            end if
        end do

        return

    end function ActiveSetContainsOrderDisorderMappedForExchange

    logical function TraceElementEdgeSystemForExchange()

        implicit none

        integer :: iElementLocal
        real(8) :: dTotalElementMoles, dMinPositiveElementFraction

        TraceElementEdgeSystemForExchange = .FALSE.
        dTotalElementMoles = SUM(DMAX1(dMolesElement(1:nElements), 0D0))
        if (dTotalElementMoles <= 0D0) return

        dMinPositiveElementFraction = HUGE(1D0)
        do iElementLocal = 1, nElements
            if (dMolesElement(iElementLocal) <= 0D0) cycle
            dMinPositiveElementFraction = DMIN1(dMinPositiveElementFraction, &
                dMolesElement(iElementLocal) / dTotalElementMoles)
        end do

        if (dMinPositiveElementFraction < 1D-4) TraceElementEdgeSystemForExchange = .TRUE.

        return

    end function TraceElementEdgeSystemForExchange

    subroutine ActiveSlotSiteForExchange(iSolnPhaseIndexIn, iPhaseIDIn, iSlotIn, dSiteOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn, iSlotIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        if (allocated(dActiveSlotSiteFraction)) then
            if ((iSlotIn > 0).AND.(iSlotIn <= SIZE(dActiveSlotSiteFraction,1)).AND.&
                (SUM(dActiveSlotSiteFraction(iSlotIn,:,:)) > 0D0)) then
                dSiteOut = dActiveSlotSiteFraction(iSlotIn,:,:)
                return
            end if
        end if

        call BuildExchangeSiteFromMolFraction(iSolnPhaseIndexIn, iPhaseIDIn, dSiteOut)

        return

    end subroutine ActiveSlotSiteForExchange

    subroutine BuildExchangeSiteFromMolFraction(iSolnPhaseIndexIn, iPhaseIDIn, dSiteOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(out) :: dSiteOut(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dSumLocal

        dSiteOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dSiteOut(iSub,iCon) = dSiteOut(iSub,iCon) + &
                    DMAX1(dMolFraction(iLocalSpecies), 0D0)
            end do
        end do

        do iSub = 1, nSublatticePhase(iPhaseIDIn)
            dSumLocal = SUM(dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)))
            if (dSumLocal > 0D0) then
                dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)) = &
                    dSiteOut(iSub,1:nConstituentSublattice(iPhaseIDIn,iSub)) / dSumLocal
            end if
        end do

        return

    end subroutine BuildExchangeSiteFromMolFraction

    subroutine SetExchangeMolFractionFromSite(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dSumLocal

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
            end do
        end if

        return

    end subroutine SetExchangeMolFractionFromSite

    subroutine ComputeCEFPhaseComposition(iSolnPhaseIndexIn, iPhaseIDIn, dSiteIn, dCompositionOut)

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

    end subroutine ComputeCEFPhaseComposition

    subroutine ComputeCEFCompositionDerivative(iSolnPhaseIndexIn, iPhaseIDIn, iSubIn, iConIn, &
        iRefIn, dSiteIn, dDerivativeOut)

        implicit none

        integer, intent(in) :: iSolnPhaseIndexIn, iPhaseIDIn, iSubIn, iConIn, iRefIn
        real(8), intent(in) :: dSiteIn(nMaxSublatticeSys,nMaxConstituentSys)
        real(8), intent(out) :: dDerivativeOut(nElements)

        integer :: iFirstLocal, iLastLocal, iLocalSpecies, iLocal, iSub, iCon
        real(8) :: dProduct, dProductSum, dProductDerivative, dProductDerivativeSum
        real(8) :: dYCon, dYRef, dDerivativeFraction
        real(8), allocatable :: dProductSpecies(:), dProductDerivativeSpecies(:)

        dDerivativeOut = 0D0
        iFirstLocal = nSpeciesPhase(iSolnPhaseIndexIn-1) + 1
        iLastLocal = nSpeciesPhase(iSolnPhaseIndexIn)
        allocate(dProductSpecies(iLastLocal-iFirstLocal+1), &
            dProductDerivativeSpecies(iLastLocal-iFirstLocal+1))
        dProductSpecies = 0D0
        dProductDerivativeSpecies = 0D0
        dProductSum = 0D0
        dProductDerivativeSum = 0D0
        dYCon = DMAX1(dSiteIn(iSubIn,iConIn), 1D-75)
        dYRef = DMAX1(dSiteIn(iSubIn,iRefIn), 1D-75)

        do iLocalSpecies = iFirstLocal, iLastLocal
            iLocal = iLocalSpecies - iFirstLocal + 1
            dProduct = 1D0
            do iSub = 1, nSublatticePhase(iPhaseIDIn)
                iCon = iConstituentSublattice(iPhaseIDIn,iSub,iLocal)
                if (iCon > 0) dProduct = dProduct * DMAX1(dSiteIn(iSub,iCon), 1D-75)
            end do

            dProductDerivative = 0D0
            iCon = iConstituentSublattice(iPhaseIDIn,iSubIn,iLocal)
            if (iCon == iConIn) dProductDerivative = dProductDerivative + dProduct / dYCon
            if (iCon == iRefIn) dProductDerivative = dProductDerivative - dProduct / dYRef

            dProductSpecies(iLocal) = dProduct
            dProductDerivativeSpecies(iLocal) = dProductDerivative
            dProductSum = dProductSum + dProduct
            dProductDerivativeSum = dProductDerivativeSum + dProductDerivative
        end do

        if (dProductSum > 0D0) then
            do iLocalSpecies = iFirstLocal, iLastLocal
                iLocal = iLocalSpecies - iFirstLocal + 1
                dDerivativeFraction = (dProductDerivativeSpecies(iLocal) * dProductSum - &
                    dProductSpecies(iLocal) * dProductDerivativeSum) / dProductSum**2
                dDerivativeOut = dDerivativeOut + &
                    dDerivativeFraction * dStoichSpecies(iLocalSpecies,1:nElements)
            end do
        end if

        if (allocated(dProductSpecies)) deallocate(dProductSpecies)
        if (allocated(dProductDerivativeSpecies)) deallocate(dProductDerivativeSpecies)

        return

    end subroutine ComputeCEFCompositionDerivative

end subroutine CompSublatticeExchangeNorm
